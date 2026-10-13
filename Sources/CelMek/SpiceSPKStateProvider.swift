// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum SpiceSPKAdapterError: Error, Equatable {
    case invalidFrameBinding
    case duplicateFrameBinding(Int32)
    case unsupportedBodyIdentifier
    case originFrameRequired
    case unknownNativeFrame(Int32)
    case frameTreeRequired
    case frameTreeUnavailable
}

/// Explicit interpretation of a kernel axes ID using an exact library frame identity.
/// The frame's origin is ignored for displacement conversion; its orientation/rates are used.
public struct SpiceFrameBinding {
    public let frameID: Int32
    public let axesFrame: Frame
    public let provenance: PhysicalDataProvenance
    public init(frameID: Int32, axesFrame: Frame, provenance: PhysicalDataProvenance) throws {
        guard frameID != 0 else { throw SpiceSPKAdapterError.invalidFrameBinding }
        self.frameID = frameID; self.axesFrame = axesFrame; self.provenance = provenance
    }
    /// Caller confirms the dataset really uses ICRF axes under its J2000 label.
    public static func icrfJ2000(provenance: PhysicalDataProvenance) throws -> SpiceFrameBinding {
        try SpiceFrameBinding(frameID: 1,axesFrame: CelestialFrames.icrf,provenance: provenance)
    }
}
public struct SpiceFrameBindings {
    public let bindings: [SpiceFrameBinding]
    private let lookup: [Int32:SpiceFrameBinding]
    public init(_ bindings: [SpiceFrameBinding]) throws {
        var lookup: [Int32:SpiceFrameBinding] = [:]
        for binding in bindings {
            guard lookup.updateValue(binding,forKey: binding.frameID) == nil else {
                throw SpiceSPKAdapterError.duplicateFrameBinding(binding.frameID)
            }
        }
        self.bindings = bindings; self.lookup = lookup
    }
    public func binding(for frameID: Int32) throws -> SpiceFrameBinding {
        guard let binding = lookup[frameID] else { throw SpiceSPKAdapterError.unknownNativeFrame(frameID) }
        return binding
    }
}
public struct SpiceSPKEphemerisEvaluation {
    public let state: CartesianState
    /// Native contributing states and selection provenance, before unit/time/frame conversion.
    public let chain: SpiceSPKChain
    public let tdbCorrection: TDBCorrection
}

/// Geometric SI ephemeris over a retained kernel snapshot. No implicit fallback.
/// Output frame must physically represent the declared origin. Caller owns that assertion.
/// Retain the supplied FrameTree externally; the weak reference avoids frame-provider cycles.
public final class SpiceSPKStateProvider: AccuracyAwareEphemerisProvider {
    public let metadata: EphemerisMetadata
    public let accuracy: AccuracyLevel
    public let effectiveAccuracy: AccuracyLevel = .full
    public var availableAccuracyLevels: Set<AccuracyLevel> { [.full] }
    public let snapshot: SpiceKernelSnapshot
    public let frameBindings: SpiceFrameBindings
    public let resolutionOptions: SpiceSPKResolutionOptions
    private let targetID: Int32
    private let centreID: Int32
    private let tdbModel: any TDBTimeModel
    private weak var tree: FrameTree?

    /// A non-SSB origin requires an explicit origin-centred frame. For SSB, nil uses ICRF.
    /// Bindings must describe every native axes ID needed by the selected chain.
    /// Coverage is caller-selected; endpoints are validated, while interior gaps fail on query.
    public init(snapshot: SpiceKernelSnapshot, body: BodyIdentifier,
                origin: BodyIdentifier = .solarSystemBarycentre, frame: Frame? = nil,
                frameBindings: SpiceFrameBindings, frameTree: FrameTree? = nil,
                coverage: ClosedRange<Epoch>, tdbModel: any TDBTimeModel,
                dataProvenance: PhysicalDataProvenance, accuracy: AccuracyLevel = .full,
                resolutionOptions: SpiceSPKResolutionOptions? = nil) throws {
        guard case .naif(let target) = body, case .naif(let centre) = origin else {
            throw SpiceSPKAdapterError.unsupportedBodyIdentifier
        }
        guard frame != nil || centre == 0 else { throw SpiceSPKAdapterError.originFrameRequired }
        let output = frame ?? CelestialFrames.icrf
        for binding in frameBindings.bindings where binding.axesFrame !== output {
            guard let frameTree else { throw SpiceSPKAdapterError.frameTreeRequired }
            guard frameTree.isRegistered(binding.axesFrame), frameTree.isRegistered(output) else {
                throw FrameTransformError.disconnectedFrames
            }
        }
        guard coverage.lowerBound >= tdbModel.coverage.lowerBound,
              coverage.upperBound <= tdbModel.coverage.upperBound else { throw EphemerisError.outsideCoverage }
        self.snapshot = snapshot; self.targetID = target; self.centreID = centre
        self.frameBindings = frameBindings; self.tree = frameTree; self.tdbModel = tdbModel
        self.resolutionOptions = try resolutionOptions ?? SpiceSPKResolutionOptions(); self.accuracy = accuracy
        let assets = snapshot.binaryKernels.map { $0.source.path }.joined(separator: "; ")
        let axes = frameBindings.bindings.map { "\($0.frameID): \($0.provenance.source) [\($0.provenance.version)]" }.joined(separator: "; ")
        metadata = try EphemerisMetadata(body: body,origin: origin,frame: output,coverage: coverage,
            coverageBasis: .callerSelected,modelTimeScale: .tdb,
            velocityMethod: .supplied(description: "SPK type 2 position derivative / type 3 supplied velocity; TDB rate and explicit frame kinematics"),
            accelerationMethod: .unavailable,
            provenance: PhysicalDataProvenance(source: "\(dataProvenance.source); SPK snapshot assets: \(assets); axes: \(axes); TDB: \(tdbModel.provenance)",
                version: dataProvenance.version,reference: dataProvenance.reference),
            accuracy: accuracy,effectiveAccuracy: .full)
        _ = try evaluation(at: coverage.lowerBound)
        if coverage.upperBound != coverage.lowerBound { _ = try evaluation(at: coverage.upperBound) }
    }
    public func state(at epoch: Epoch) throws -> CartesianState { try evaluation(at: epoch).state }

    public func evaluation(at epoch: Epoch) throws -> SpiceSPKEphemerisEvaluation {
        guard metadata.coverage.contains(epoch) else { throw EphemerisError.outsideCoverage }
        let correction = try epoch.tdbCorrection(using: tdbModel)
        // Avoid subtracting large Julian dates. SPK ultimately consumes one Double ET.
        let seconds = try epoch.offsetFromJ2000.adding(correction.tdbMinusTT).seconds
        let chain = try snapshot.spkChain(targetID: targetID,centreID: centreID,atSecondsTDB: seconds,options: resolutionOptions)
        var p = SPKCompensatedVector(), v = SPKCompensatedVector()
        for term in chain.contributions {
            let binding = try frameBindings.binding(for: term.state.descriptor.frameID)
            let native = try CartesianState(epoch: epoch,frame: binding.axesFrame,
                position: term.state.positionKilometres*1000,
                velocity: term.state.velocityKilometresPerTDBSecond*1000*correction.rate)
            let converted: CartesianState
            if binding.axesFrame === frame { converted = native }
            else {
                guard let tree else { throw SpiceSPKAdapterError.frameTreeUnavailable }
                let transform = try tree.transform(from: binding.axesFrame,to: frame,at: epoch)
                // Centre subtraction is already in the chain. These are relative vectors:
                // use only R and its derivatives, never the tree's origin translation.
                let orientation = try FrameTransform(epoch: epoch,rotation: transform.rotation,
                    angularVelocity: transform.angularVelocity,angularAcceleration: transform.angularAcceleration)
                converted = try orientation.applying(to: native,destination: frame)
            }
            let sign = term.side == .target ? 1.0 : -1.0
            try p.add(sign*converted.position); try v.add(sign*converted.velocity)
        }
        let state = try CartesianState(epoch: epoch,frame: frame,position: p.value(),velocity: v.value())
        return SpiceSPKEphemerisEvaluation(state: state,chain: chain,tdbCorrection: correction)
    }
}
