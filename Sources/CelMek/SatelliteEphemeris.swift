import Foundation
import Math

/// Native, planet-centred identities. These are deliberately distinct from ICRF
/// and body-fixed surface frames. Retain this bundle to share identities between
/// moons. A FrameTree may adopt either identity as its root.
public struct SatelliteEphemerisFrames {
  public let jupiterEquatorialOfDate: Frame
  public let uranusGUST86: Frame
  public init() {
    jupiterEquatorialOfDate = Frame(name: "Jupiter equator, ascending node on ecliptic of date (E2x3)", pseudoInertial: false)
    uranusGUST86 = Frame(name: "Uranus equator, GUST86 B1950 node", pseudoInertial: true)
  }
  internal init(jupiter: Frame, uranus: Frame) {
    jupiterEquatorialOfDate = jupiter; uranusGUST86 = uranus
  }
}

extension FrameTree {
  /// Register native satellite frames under caller-owned planet-centred parents.
  /// Providers must rotate the parent axes into the documented native axes,
  /// including their time derivatives. No pole/precession convention is guessed.
  public func makeSatelliteEphemerisFrames(jupiterParent: Frame, jupiterTransform: any TransformProvider,
    uranusParent: Frame, uranusTransform: any TransformProvider) throws -> SatelliteEphemerisFrames {
    let j = try makeFrame(parent: jupiterParent, name: "Jupiter E2x3 equatorial of date", provider: jupiterTransform)
    let u = try makeFrame(parent: uranusParent, name: "Uranus GUST86 equatorial", pseudoInertial: true, provider: uranusTransform)
    return SatelliteEphemerisFrames(jupiter: j, uranus: u)
  }
}

/// Adapter for the supplied historical satellite series. Coverage is caller
/// selected, not an accuracy claim. Uses TT as the legacy ephemeris-time argument.
/// Immutable local evaluation state permits independent/concurrent moon queries.
public struct SatelliteEphemeris: EphemerisProvider {
  public enum Model: Sendable {
    case lieskeE2x3(LieskeE2x3.Moon)
    case gust86(GUST86.Moon)
  }
  public let model: Model
  public let metadata: EphemerisMetadata
  public let positionCoverage: ClosedRange<Epoch>
  public let velocityStep: Duration

  public init(model: Model, frames: SatelliteEphemerisFrames, positionCoverage: ClosedRange<Epoch>,
              velocityStep: Duration = try! Duration(seconds: 30)) throws {
    guard velocityStep > .zero else { throw EphemerisError.invalidDerivativeMethod }
    let margin = try velocityStep.adding(velocityStep)
    let lower = try positionCoverage.lowerBound.advanced(by: margin)
    let upper = try positionCoverage.upperBound.advanced(by: margin.negated())
    guard lower <= upper else { throw EphemerisError.outsideCoverage }
    self.model = model; self.positionCoverage = positionCoverage; self.velocityStep = velocityStep
    let body: BodyIdentifier, origin: BodyIdentifier, frame: Frame, source: String
    switch model {
    case .lieskeE2x3(let moon):
      body = .naif(Int32(501+moon.rawValue)); origin = .jupiter
      frame = frames.jupiterEquatorialOfDate
      source = "Lieske E2x3 supplied legacy coefficients; corrected degree trig, Io evaluation, Europa latitude arguments, fractional radius and atan latitude; Jupiter radius 71492 km"
    case .gust86(let moon):
      let codes: [Int32] = [705,701,702,703,704]
      body = .naif(codes[moon.rawValue]); origin = .uranus; frame = frames.uranusGUST86
      source = "GUST86 supplied unfinished coefficient set; corrected eccentricity phase lookup; Cartesian conversion uses table 5 a0 and n0; not precision validated"
    }
    metadata = try EphemerisMetadata(body: body, origin: origin, frame: frame,
      coverage: lower...upper, coverageBasis: .callerSelected, modelTimeScale: .tt,
      velocityMethod: .finiteDifference(step: velocityStep), accelerationMethod: .unavailable,
      provenance: PhysicalDataProvenance(source: source, version: "2012/2013 legacy series, Swift port"))
  }
  private func position(at epoch: Epoch) throws -> SIMD3<Double> {
    guard positionCoverage.contains(epoch) else { throw EphemerisError.outsideCoverage }
    switch model {
    case .lieskeE2x3(let moon): return try LieskeE2x3.position(of: moon, at: epoch)
    case .gust86(let moon): return try GUST86.position(of: moon, at: epoch)
    }
  }
  public func state(at epoch: Epoch) throws -> CartesianState {
    guard metadata.coverage.contains(epoch) else { throw EphemerisError.outsideCoverage }
    let p = try position(at: epoch), h = velocityStep.seconds
    let samples = try [-2.0,-1,1,2].map { try position(at: epoch.advanced(by: Duration(seconds: $0*h))) }
    let v = (samples[0]-8*samples[1]+8*samples[2]-samples[3])/(12*h)
    return try CartesianState(epoch: epoch, frame: frame, position: p, velocity: v)
  }
}

extension BodyIdentifier {
  public static let io = Self.naif(501)
  public static let europa = Self.naif(502)
  public static let ganymede = Self.naif(503)
  public static let callisto = Self.naif(504)
  public static let miranda = Self.naif(705)
  public static let ariel = Self.naif(701)
  public static let umbriel = Self.naif(702)
  public static let titania = Self.naif(703)
  public static let oberon = Self.naif(704)
}
