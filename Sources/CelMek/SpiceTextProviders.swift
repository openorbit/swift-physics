// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum SpiceTextProviderError: Error, Equatable {
    case outsideCoverage
    case invalidReferenceFrame
    case conflictingReferenceParameters
    case invalidPhaseParameters
}

extension SpiceKernelSnapshot {
    internal func physicalProvenance(names: [String], supplied: PhysicalDataProvenance) throws -> PhysicalDataProvenance {
        try PhysicalDataProvenance(source: supplied.source + "; " + sources(for: names).map(\.path).joined(separator: "; "),
                                   version: supplied.version, reference: supplied.reference)
    }
    public func bodyReferenceRadii(bodyID: Int32, axesDescription: String,
                                   provenance: PhysicalDataProvenance) throws -> BodyReferenceRadii {
        let radii = try bodyRadiiMetres(bodyID: bodyID)
        return try BodyReferenceRadii(metres: SIMD3(radii[0], radii[1], radii[2]), axesDescription: axesDescription,
            provenance: physicalProvenance(names: ["BODY\(bodyID)_RADII"], supplied: provenance))
    }
    /// Caller confirms conventional km³/s² input and its coordinate-time convention.
    public func bodyGravitationalParameter(bodyID: Int32, convention: GravitationalParameterConvention,
                                            provenance: PhysicalDataProvenance) throws -> GravitationalParameter {
        let name = "BODY\(bodyID)_GM"
        return try GravitationalParameter(cubicMetresPerSecondSquared: numeric(name, count: 1)[0] * 1e9,
            convention: convention, provenance: physicalProvenance(names: [name], supplied: provenance))
    }
}

/// Full text PCK polynomial/periodic model using a caller-declared system ID.
/// Coverage is caller-selected, not a certified physical accuracy interval.
public struct SpicePCKOrientationProvider: OrientationProvider {
    public let bodyID: Int32
    public let systemID: Int32
    public let referenceFrame: Frame
    public let referenceFrameID: Int32
    public let referenceJulianDate: Double
    public let coverage: ClosedRange<Epoch>
    public let provenance: PhysicalDataProvenance
    public let cacheScope: EphemerisCacheScope
    private let tdbModel: any TDBTimeModel
    private let ra: IAUAngleSeries
    private let dec: IAUAngleSeries
    private let pm: IAUAngleSeries

    /// systemID is the planetary barycentre for planets/satellites, otherwise bodyID.
    /// This explicit argument avoids guessing relationships for unusual NAIF IDs.
    public init(snapshot: SpiceKernelSnapshot, bodyID: Int32, systemID: Int32,
                frameBindings: SpiceFrameBindings, coverage: ClosedRange<Epoch>,
                tdbModel: any TDBTimeModel, provenance: PhysicalDataProvenance) throws {
        guard tdbModel.coverage.contains(coverage.lowerBound), tdbModel.coverage.contains(coverage.upperBound) else {
            throw TDBModelError.outsideCoverage
        }
        let parameters = try snapshot.bodyOrientationParameters(bodyID: bodyID)
        let prefix = "BODY\(systemID)_"
        var names = ["BODY\(bodyID)_POLE_RA", "BODY\(bodyID)_POLE_DEC", "BODY\(bodyID)_PM"]
        func reference(_ suffix: String, fallback: Double) throws -> Double {
            let long = prefix + "CONSTANTS_" + suffix, short = prefix + "CONSTS_" + suffix
            let a = try snapshot[long].map { _ in try snapshot.numeric(long, count: 1)[0] }
            let b = try snapshot[short].map { _ in try snapshot.numeric(short, count: 1)[0] }
            if a != nil { names.append(long) }; if b != nil { names.append(short) }
            if let a, let b, a != b { throw SpiceTextProviderError.conflictingReferenceParameters }
            return a ?? b ?? fallback
        }
        let frameNumber = try reference("REF_FRAME", fallback: 1)
        guard let frameID = Int32(exactly: frameNumber), (1...21).contains(frameID) else {
            throw SpiceTextProviderError.invalidReferenceFrame
        }
        let epoch = try reference("JED_EPOCH", fallback: 2451545)
        let binding = try frameBindings.binding(for: frameID)
        guard binding.axesFrame.pseudoInertial else { throw SpiceTextProviderError.invalidReferenceFrame }
        // Ensure representability in the library's supported continuous date range.
        _ = try Epoch(julianDate: epoch, scale: .tt)
        let degreeName = prefix + "MAX_PHASE_DEGREE", phasesName = prefix + "NUT_PREC_ANGLES"
        var phaseDegree = 1
        if snapshot[degreeName] != nil {
            let degree = try snapshot.numeric(degreeName, count: 1)[0]
            guard let value = Int(exactly: degree), (0...3).contains(value) else {
                throw SpiceTextProviderError.invalidPhaseParameters
            }
            phaseDegree = value; names.append(degreeName)
        }
        var phases: [[Double]] = []
        if snapshot[phasesName] != nil {
            let values = try snapshot.numeric(phasesName)
            guard values.count % (phaseDegree + 1) == 0 else { throw SpiceTextProviderError.invalidPhaseParameters }
            phases = stride(from: 0, to: values.count, by: phaseDegree + 1).map { Array(values[$0...($0 + phaseDegree)]) }
            names.append(phasesName)
        }
        func series(_ polynomial: [Double], _ suffix: String, _ variable: IAUAngleSeries.Variable,
                    _ function: IAUHarmonic.Function) throws -> IAUAngleSeries {
            let key = "BODY\(bodyID)_NUT_PREC_" + suffix
            let amplitudes = parameters.additionalParameters["NUT_PREC_" + suffix] ?? []
            guard amplitudes.count <= phases.count else { throw SpiceTextProviderError.invalidPhaseParameters }
            if !amplitudes.isEmpty { names.append(key) }
            return IAUAngleSeries(polynomialDegrees: polynomial, variable: variable,
                harmonics: amplitudes.enumerated().map { IAUHarmonic(amplitudeDegrees: $0.element, phaseDegrees: phases[$0.offset], function: function) })
        }
        ra = try series(parameters.poleRightAscension, "RA", .centuries, .sine)
        dec = try series(parameters.poleDeclination, "DEC", .centuries, .cosine)
        pm = try series(parameters.primeMeridian, "PM", .days, .sine)
        self.bodyID = bodyID; self.systemID = systemID; self.coverage = coverage; self.tdbModel = tdbModel
        referenceFrame = binding.axesFrame; referenceFrameID = frameID; referenceJulianDate = epoch
        cacheScope = snapshot.cacheScope
        let data = try snapshot.physicalProvenance(names: names, supplied: provenance)
        self.provenance = try PhysicalDataProvenance(source: data.source + "; axes: " + binding.provenance.source + " (" + binding.provenance.version + "); time: " + tdbModel.provenance,
            version: data.version, reference: data.reference)
        _ = try orientation(at: coverage.lowerBound)
        _ = try orientation(at: coverage.upperBound)
    }
    public func orientation(at epoch: Epoch) throws -> BodyOrientation {
        guard coverage.contains(epoch) else { throw SpiceTextProviderError.outsideCoverage }
        let clock = try epoch.tdbCorrection(using: tdbModel)
        let seconds = try epoch.offsetFromJ2000.adding(clock.tdbMinusTT).seconds
        let days = seconds / 86400 + (2451545 - referenceJulianDate)
        let r = try ra.evaluate(days: days), d = try dec.evaluate(days: days), w = try pm.evaluate(days: days)
        func rotation(_ angle: Double, _ first: Double, _ second: Double, _ axis: SIMD3<Double>) throws -> FrameTransform {
            guard angle.isFinite else { throw FrameTransformError.nonFinite }
            return try FrameTransform(epoch: epoch,
                rotation: Quaternion(angle: angle.truncatingRemainder(dividingBy: 2 * .pi), axis: axis),
                angularVelocity: axis * (first * clock.rate),
                angularAcceleration: axis * (second * clock.rate * clock.rate + first * clock.acceleration))
        }
        let node = try rotation(-(r.radians + .pi/2), -r.rate, -r.acceleration, SIMD3(0,0,1))
        let tilt = try rotation(d.radians - .pi/2, d.rate, d.acceleration, SIMD3(1,0,0))
        let spin = try rotation(-w.radians, -w.rate, -w.acceleration, SIMD3(0,0,1))
        let transform = try node.followed(by: tilt).followed(by: spin)
        return try BodyOrientation(referenceFrame: referenceFrame, epoch: epoch, rotation: transform.rotation,
                                   angularVelocity: transform.angularVelocity, angularAcceleration: transform.angularAcceleration)
    }
}

/// Constant FK axes relative to the exact caller-bound parent axes.
public struct SpiceFixedOrientationProvider: OrientationProvider {
    public let referenceFrame: Frame
    public let definition: SpiceFixedFrameDefinition
    public let provenance: PhysicalDataProvenance
    private let rotation: Quaternion
    /// Caller asserts relativeAxes physically represents definition.relativeFrameName.
    public init(definition: SpiceFixedFrameDefinition, relativeAxes: Frame, provenance: PhysicalDataProvenance) throws {
        self.definition = definition; referenceFrame = relativeAxes
        self.provenance = try PhysicalDataProvenance(source: provenance.source + "; " + definition.sources.map(\.path).joined(separator: "; "),
            version: provenance.version, reference: provenance.reference)
        // Decode a validated proper column-major matrix using its largest quaternion component.
        let m = definition.matrix
        let diagonal = [1+m[0]+m[4]+m[8], 1+m[0]-m[4]-m[8], 1-m[0]+m[4]-m[8], 1-m[0]-m[4]+m[8]]
        let i = diagonal.indices.max(by: { diagonal[$0] < diagonal[$1] })!
        let s = sqrt(max(0, diagonal[i])) * 2
        let q: Quaternion
        switch i {
        case 0: q = Quaternion(x: (m[5]-m[7])/s, y: (m[6]-m[2])/s, z: (m[1]-m[3])/s, w: s/4)
        case 1: q = Quaternion(x: s/4, y: (m[3]+m[1])/s, z: (m[6]+m[2])/s, w: (m[5]-m[7])/s)
        case 2: q = Quaternion(x: (m[3]+m[1])/s, y: s/4, z: (m[7]+m[5])/s, w: (m[6]-m[2])/s)
        default: q = Quaternion(x: (m[6]+m[2])/s, y: (m[7]+m[5])/s, z: s/4, w: (m[1]-m[3])/s)
        }
        // FK matrix is local-to-relative; OrientationProvider is reference-to-local.
        rotation = try FrameTransform(epoch: .j2000, rotation: q.conj).rotation
    }
    public func orientation(at epoch: Epoch) throws -> BodyOrientation {
        try BodyOrientation(referenceFrame: referenceFrame, epoch: epoch, rotation: rotation)
    }
}

extension SpiceFixedFrameDefinition {
    /// Caller supplies the physical centreID's state in relativeAxes coordinates.
    /// Translation and its derivatives are handled by the existing body-fixed provider.
    public func makeFrame(in tree: FrameTree, relativeAxes: Frame, centre: any StateProvider,
                          provenance: PhysicalDataProvenance,
                          accelerationPolicy: OriginAccelerationPolicy = .requireProvided) throws -> Frame {
        let orientation = try SpiceFixedOrientationProvider(definition: self, relativeAxes: relativeAxes, provenance: provenance)
        return try tree.makeBodyFixedFrame(name: name, centre: centre, orientation: orientation, accelerationPolicy: accelerationPolicy)
    }
}
