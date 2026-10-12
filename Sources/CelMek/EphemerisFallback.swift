// SPDX-License-Identifier: Apache-2.0
import Math

public enum EphemerisFallbackError: Error, Equatable {
  case noProviders, incompatibleBody, incompatibleOrigin, coverageGap
}

/// Ordered tiers, independent of requested numerical accuracy.
public enum EphemerisFallbackTier: Sendable { case de441, model, kepler }

public struct EphemerisFallbackResult {
  public let state: CartesianState
  public let tier: EphemerisFallbackTier
  /// Metadata of the selected source, before the chain's derivative normalization.
  public let source: EphemerisMetadata
}

/// DE441 -> analytic/fitted model -> Kepler. Missing tiers are skipped. Selection
/// uses declared coverage only; errors from an eligible source are never swallowed.
/// Inputs must already describe the same target, origin and frame identity.
public struct FallbackEphemeris: EphemerisProvider {
  public let metadata: EphemerisMetadata
  private let sources: [(EphemerisFallbackTier, any EphemerisProvider)]

  public init(de441: (any EphemerisProvider)? = nil,
              model: (any EphemerisProvider)? = nil,
              kepler: (any EphemerisProvider)? = nil) throws {
    var sources: [(EphemerisFallbackTier, any EphemerisProvider)] = []
    if let de441 { sources.append((.de441, de441)) }
    if let model { sources.append((.model, model)) }
    if let kepler { sources.append((.kepler, kepler)) }
    guard let first = sources.first?.1.metadata else { throw EphemerisFallbackError.noProviders }
    for (_, p) in sources {
      guard p.metadata.body == first.body else { throw EphemerisFallbackError.incompatibleBody }
      guard p.metadata.origin == first.origin else { throw EphemerisFallbackError.incompatibleOrigin }
      guard p.frame === first.frame, p.metadata.frame === first.frame else { throw EphemerisError.frameMismatch }
    }
    // A single metadata interval must not claim coverage across an uncovered gap.
    let intervals = sources.map { $0.1.metadata.coverage }.sorted { $0.lowerBound < $1.lowerBound }
    var end = intervals[0].upperBound
    for interval in intervals.dropFirst() {
      guard interval.lowerBound <= end else { throw EphemerisFallbackError.coverageGap }
      end = max(end, interval.upperBound)
    }
    let hasAcceleration = sources.allSatisfy { $0.1.metadata.accelerationMethod != .unavailable }
    metadata = try EphemerisMetadata(body: first.body, origin: first.origin, frame: first.frame,
      coverage: intervals[0].lowerBound...end, coverageBasis: .callerSelected,
      modelTimeScale: first.modelTimeScale,
      velocityMethod: .supplied(description: "Selected fallback source; see result source metadata"),
      accelerationMethod: hasAcceleration ? .supplied(description: "Selected fallback source") : .unavailable,
      provenance: PhysicalDataProvenance(source: sources.map { $0.1.metadata.provenance.source }.joined(separator: " -> "),
        version: "DE441 -> model -> Kepler coverage policy"),
      modelTimeScales: sources.reduce(into: Set<TimeScale>()) { $0.formUnion($1.1.metadata.modelTimeScales) })
    self.sources = sources
  }

  public func selection(at epoch: Epoch) throws -> EphemerisFallbackResult {
    guard let (tier, provider) = sources.first(where: { $0.1.metadata.coverage.contains(epoch) }) else {
      throw EphemerisError.outsideCoverage
    }
    let native = try provider.checkedState(at: epoch)
    let state = try CartesianState(epoch: epoch, frame: frame, position: native.position,
      velocity: native.velocity, acceleration: metadata.accelerationMethod == .unavailable ? nil : native.acceleration)
    return EphemerisFallbackResult(state: state, tier: tier, source: provider.metadata)
  }
  public func state(at epoch: Epoch) throws -> CartesianState { try selection(at: epoch).state }
}

/// Explicit fixed-seed Kepler fallback. Coverage is a caller's usable interval,
/// not a claim about the physical accuracy of unperturbed two-body motion.
public struct KeplerEphemeris: EphemerisProvider {
  public let metadata: EphemerisMetadata
  public let propagator: TwoBodyPropagator
  public init(body: BodyIdentifier, propagator: TwoBodyPropagator,
              coverage: ClosedRange<Epoch>, provenance: PhysicalDataProvenance) throws {
    self.propagator = propagator
    metadata = try EphemerisMetadata(body: body, origin: propagator.centralBody, frame: propagator.frame,
      coverage: coverage, coverageBasis: .callerSelected, modelTimeScale: .tt,
      velocityMethod: .analytic, accelerationMethod: .analytic, provenance: provenance)
  }
  public func state(at epoch: Epoch) throws -> CartesianState {
    guard metadata.coverage.contains(epoch) else { throw EphemerisError.outsideCoverage }
    return try propagator.state(at: epoch)
  }
}

/// Built-in defaults use exact body identities: a planet centre is not silently
/// replaced with its system barycentre. Caller-supplied fitted models can override these.
public enum DefaultEphemerisModel: Equatable, Sendable {
  case vsop87, meeusMoon, goffinPluto
  public static func model(for body: BodyIdentifier) -> Self? {
    switch body {
    case .sun, .mercury, .venus, .earth, .marsBarycentre, .jupiterBarycentre,
         .saturnBarycentre, .uranusBarycentre, .neptuneBarycentre: return .vsop87
    case .moon: return .meeusMoon
    case .plutoBarycentre: return .goffinPluto
    default: return nil
    }
  }
  /// Native-frame model. Convert frames and origins explicitly before composing
  /// with DE441 or a central-body Kepler ephemeris. Lunar/Pluto frames are required
  /// only for those bodies. Their query coverage is shortened by derivative margins.
  public static func makeProvider(for body: BodyIdentifier, coverage: ClosedRange<Epoch>,
    tdbModel: any TDBTimeModel, frames: AnalyticEphemerisFrames? = nil,
    velocityStep: Duration = try! Duration(seconds: 60), accuracy: AccuracyLevel = .full) throws -> (any EphemerisProvider)? {
    let coefficients: vsop87_body?
    switch body {
    case .sun: coefficients = sun
    case .mercury: coefficients = mercury
    case .venus: coefficients = venus
    case .earth: coefficients = earth
    case .marsBarycentre: coefficients = mars
    case .jupiterBarycentre: coefficients = jupiter
    case .saturnBarycentre: coefficients = saturn
    case .uranusBarycentre: coefficients = uranus
    case .neptuneBarycentre: coefficients = neptune
    default: coefficients = nil
    }
    if let coefficients { return try VSOP87Ephemeris(body: coefficients, coverage: coverage, tdbModel: tdbModel, accuracy: accuracy) }
    guard body == .moon || body == .plutoBarycentre else { return nil }
    guard let frames else { throw EphemerisError.frameMismatch }
    return try SphericalAnalyticEphemeris(model: body == .moon ? .meeusMoon : .goffinPluto,
      frames: frames, positionCoverage: coverage, velocityStep: velocityStep, accuracy: accuracy)
  }
}
