// SPDX-License-Identifier: Apache-2.0
import Math

public enum EphemerisOriginError: Error, Equatable {
  case incompatibleOrigins, incompatibleFrames, disjointCoverage, invalidComponents, invalidMassRatio
}

/// Target relative to a new origin, expressed in a registered translated child frame.
/// Construct through FrameTree.makeOriginConvertedEphemeris so the frame's physical
/// translation and identity are established together. Native queries subtract SI states.
public struct OriginConvertedEphemeris: EphemerisProvider {
  public let metadata: EphemerisMetadata
  private let target: any EphemerisProvider
  private let newOrigin: any EphemerisProvider
  internal init(target: any EphemerisProvider, newOrigin: any EphemerisProvider,
                frame: Frame, coverage: ClosedRange<Epoch>) throws {
    self.target = target; self.newOrigin = newOrigin
    metadata = try EphemerisMetadata(body: target.metadata.body,origin: newOrigin.metadata.body,frame: frame,
      coverage: coverage,coverageBasis: .callerSelected,modelTimeScale: target.metadata.modelTimeScale,
      velocityMethod: combinedDerivative([target.metadata.velocityMethod,newOrigin.metadata.velocityMethod]),
      accelerationMethod: combinedDerivative([target.metadata.accelerationMethod,newOrigin.metadata.accelerationMethod]),
      provenance: PhysicalDataProvenance(source: "Explicit origin subtraction; target: \(target.metadata.provenance.source) [\(target.metadata.provenance.version)]; origin: \(newOrigin.metadata.provenance.source) [\(newOrigin.metadata.provenance.version)]",
        version: "SI state composition 1"),
      modelTimeScales: target.metadata.modelTimeScales.union(newOrigin.metadata.modelTimeScales))
  }
  public func state(at epoch: Epoch) throws -> CartesianState {
    guard metadata.coverage.contains(epoch) else { throw EphemerisError.outsideCoverage }
    let a = try target.checkedState(at: epoch), b = try newOrigin.checkedState(at: epoch)
    let acceleration: SIMD3<Double>?
    if let aa = a.acceleration, let ba = b.acceleration { acceleration = aa-ba } else { acceleration = nil }
    return try CartesianState(epoch: epoch,frame: frame,position: a.position-b.position,
      velocity: a.velocity-b.velocity,acceleration: acceleration)
  }
}

extension FrameTree {
  /// Both ephemerides must use the exact same native frame and declared origin.
  /// For differently oriented frames, convert them explicitly before composition.
  /// The edge's acceleration policy governs later frame transforms, not subtraction:
  /// recentered acceleration is available only when both input states supply it.
  public func makeOriginConvertedEphemeris(target: any EphemerisProvider,
    newOrigin: any EphemerisProvider, frameName: String,
    accelerationPolicy: OriginAccelerationPolicy = .requireProvided) throws -> OriginConvertedEphemeris {
    guard target.metadata.origin == newOrigin.metadata.origin else { throw EphemerisOriginError.incompatibleOrigins }
    guard target.frame === target.metadata.frame, newOrigin.frame === newOrigin.metadata.frame,
          target.frame === newOrigin.frame else { throw EphemerisOriginError.incompatibleFrames }
    let coverage = try intersectCoverage([target.metadata.coverage,newOrigin.metadata.coverage])
    let axes = try FixedAxisOrientationProvider(referenceFrame: target.frame,axis: SIMD3(1,0,0),rate: 0)
    // Use a checked wrapper even for custom implementations retained by the frame tree.
    let checkedOrigin = CheckedEphemerisState(provider: newOrigin)
    let child = try makeBodyFixedFrame(name: frameName,centre: checkedOrigin,orientation: axes,
      accelerationPolicy: accelerationPolicy)
    return try OriginConvertedEphemeris(target: target,newOrigin: newOrigin,frame: child,coverage: coverage)
  }
}

/// A constant positive mass ratio and a constituent ephemeris. Ratios are relative
/// weights, not independently assumed GM values or automatically inferred masses.
public struct BarycentreComponent {
  public let ephemeris: any EphemerisProvider
  public let relativeMass: Double
  public let provenance: PhysicalDataProvenance
  public init(ephemeris: any EphemerisProvider, relativeMass: Double,
              provenance: PhysicalDataProvenance) throws {
    guard relativeMass.isFinite, relativeMass > 0 else { throw EphemerisOriginError.invalidMassRatio }
    self.ephemeris = ephemeris; self.relativeMass = relativeMass; self.provenance = provenance
  }
}

/// Constant-mass weighted barycentre. Constituents must share an origin and frame;
/// every constituent must be present. No satellite reconstruction is implied.
public struct BarycentreEphemeris: EphemerisProvider {
  public let metadata: EphemerisMetadata
  public let components: [BarycentreComponent]
  private let weights: [Double]
  public init(bodyIdentifier: BodyIdentifier, components: [BarycentreComponent],
              provenance: PhysicalDataProvenance) throws {
    try bodyIdentifier.validate()
    guard components.count >= 2, let first = components.first else { throw EphemerisOriginError.invalidComponents }
    let ids = components.map { $0.ephemeris.metadata.body }
    guard Set(ids).count == ids.count, !ids.contains(bodyIdentifier) else { throw EphemerisOriginError.invalidComponents }
    let metadata = first.ephemeris.metadata
    for component in components {
      guard component.ephemeris.metadata.origin == metadata.origin else { throw EphemerisOriginError.incompatibleOrigins }
      guard component.ephemeris.frame === metadata.frame,
            component.ephemeris.metadata.frame === metadata.frame else { throw EphemerisOriginError.incompatibleFrames }
    }
    let coverage = try intersectCoverage(components.map { $0.ephemeris.metadata.coverage })
    let maximum = components.map(\.relativeMass).max()!
    let scaled = components.map { $0.relativeMass/maximum }, total = scaled.reduce(0,+)
    weights = scaled.map { $0/total }; self.components = components
    let details = components.map { "\($0.ephemeris.metadata.provenance.source) [\($0.ephemeris.metadata.provenance.version)]; mass ratio: \($0.provenance.source) [\($0.provenance.version)]" }.joined(separator: "; ")
    let source = try PhysicalDataProvenance(source: "\(provenance.source); components: \(details)",
      version: provenance.version,reference: provenance.reference)
    self.metadata = try EphemerisMetadata(body: bodyIdentifier,origin: metadata.origin,frame: metadata.frame,
      coverage: coverage,coverageBasis: .callerSelected,modelTimeScale: metadata.modelTimeScale,
      velocityMethod: combinedDerivative(components.map { $0.ephemeris.metadata.velocityMethod }),
      accelerationMethod: combinedDerivative(components.map { $0.ephemeris.metadata.accelerationMethod }),
      provenance: source,modelTimeScales: components.reduce(into: Set<TimeScale>()) { $0.formUnion($1.ephemeris.metadata.modelTimeScales) })
  }
  public func state(at epoch: Epoch) throws -> CartesianState {
    guard metadata.coverage.contains(epoch) else { throw EphemerisError.outsideCoverage }
    let samples = try components.map { try $0.ephemeris.checkedState(at: epoch) }
    var p = SIMD3<Double>.zero, v = p, a = p
    for (i,sample) in samples.enumerated() {
      p += weights[i]*sample.position; v += weights[i]*sample.velocity
      if let supplied = sample.acceleration { a += weights[i]*supplied }
    }
    return try CartesianState(epoch: epoch,frame: frame,position: p,velocity: v,
      acceleration: metadata.accelerationMethod == .unavailable ? nil : a)
  }
}

private struct CheckedEphemerisState: StateProvider {
  let provider: any EphemerisProvider
  var frame: Frame { provider.frame }
  func state(at epoch: Epoch) throws -> CartesianState { try provider.checkedState(at: epoch) }
}

private func intersectCoverage(_ ranges: [ClosedRange<Epoch>]) throws -> ClosedRange<Epoch> {
  guard let lower = ranges.map(\.lowerBound).max(), let upper = ranges.map(\.upperBound).min(),
        lower <= upper else { throw EphemerisOriginError.disjointCoverage }
  return lower...upper
}

private func combinedDerivative(_ methods: [EphemerisDerivativeMethod]) -> EphemerisDerivativeMethod {
  if methods.contains(.unavailable) { return .unavailable }
  if methods.allSatisfy({ $0 == .analytic }) { return .analytic }
  return .supplied(description: "Constant linear combination of component derivatives: \(methods)")
}
