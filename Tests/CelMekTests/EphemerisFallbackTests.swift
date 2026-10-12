import Testing
import Math
@testable import CelMek

@Suite struct EphemerisFallbackTests {
  private enum Failure: Error { case corrupt }
  private struct Source: EphemerisProvider {
    let metadata: EphemerisMetadata
    let value: Double
    var fails = false
    var wrongEpoch = false
    func state(at epoch: Epoch) throws -> CartesianState {
      if fails { throw Failure.corrupt }
      return try CartesianState(epoch: wrongEpoch ? .j2000 : epoch, frame: frame,
        position: SIMD3(value,0,0), velocity: .zero,
        acceleration: metadata.accelerationMethod == .unavailable ? nil : .zero)
    }
  }
  private func epoch(_ t: Double) throws -> Epoch { try Epoch.j2000.advanced(by: CelMek.Duration(seconds: t)) }
  private func source(_ frame: Frame, _ lower: Double, _ upper: Double, _ value: Double,
    body: BodyIdentifier = .earth, origin: BodyIdentifier = .sun, acceleration: Bool = false) throws -> Source {
    Source(metadata: try EphemerisMetadata(body: body, origin: origin, frame: frame,
      coverage: epoch(lower)...epoch(upper), coverageBasis: .callerSelected, modelTimeScale: .tt,
      velocityMethod: .analytic, accelerationMethod: acceleration ? .analytic : .unavailable,
      provenance: PhysicalDataProvenance(source: "source \(value)", version: "1")), value: value)
  }
  @Test func priorityBoundariesAndOrderIndependence() throws {
    let frame = Frame(name: "test", pseudoInertial: true)
    let chain = try FallbackEphemeris(de441: source(frame,-1,1,1), model: source(frame,-2,2,2),
      kepler: source(frame,-3,3,3,acceleration: true))
    for (t, expected, tier) in [(3.0,3.0,EphemerisFallbackTier.kepler), (1,1,.de441),
      (-2,2,.model), (0,1,.de441), (-3,3,.kepler), (2,2,.model)] {
      let result = try chain.selection(at: epoch(t))
      #expect(result.state.position.x == expected)
      #expect(result.tier == tier)
      #expect(result.source.provenance.source == "source \(expected)")
      #expect(result.state.acceleration == nil)
      #expect(try chain.checkedState(at: epoch(t)) == result.state)
    }
    #expect(throws: EphemerisError.outsideCoverage) { try chain.state(at: epoch(4)) }
  }
  @Test func missingTiersAndErrors() throws {
    let frame = Frame(name: "test", pseudoInertial: true)
    var primary = try source(frame,-1,1,1); primary.fails = true
    let backup = try source(frame,-3,3,3)
    let chain = try FallbackEphemeris(de441: primary, kepler: backup)
    #expect(throws: Failure.corrupt) { try chain.state(at: .j2000) }
    #expect(try chain.selection(at: epoch(2)).tier == .kepler)
    primary.fails = false; primary.wrongEpoch = true
    let invalid = try FallbackEphemeris(de441: primary, model: backup)
    #expect(throws: FrameTransformError.epochMismatch) { try invalid.state(at: epoch(1)) }
    #expect(try FallbackEphemeris(kepler: backup).selection(at: .j2000).tier == .kepler)
    #expect(throws: EphemerisFallbackError.noProviders) { try FallbackEphemeris() }
  }
  @Test func rejectsIncompatibleSourcesAndGaps() throws {
    let frame = Frame(name: "test", pseudoInertial: true)
    let p = try source(frame,-1,1,1)
    #expect(throws: EphemerisFallbackError.incompatibleBody) {
      try FallbackEphemeris(de441: p, model: source(frame,-2,2,2,body: .moon))
    }
    #expect(throws: EphemerisFallbackError.incompatibleOrigin) {
      try FallbackEphemeris(de441: p, model: source(frame,-2,2,2,origin: .earth))
    }
    #expect(throws: EphemerisError.frameMismatch) {
      try FallbackEphemeris(de441: p, model: source(Frame(name: "other", pseudoInertial: true),-2,2,2))
    }
    #expect(throws: EphemerisFallbackError.coverageGap) {
      try FallbackEphemeris(de441: p, model: source(frame,2,3,2))
    }
    let touching = try FallbackEphemeris(de441: p, model: source(frame,1,3,2))
    #expect(try touching.selection(at: epoch(1)).tier == .de441)
    #expect(try touching.selection(at: epoch(3)).tier == .model)
  }
  @Test func defaultBodyIdentityMapping() {
    for body: BodyIdentifier in [.sun,.mercury,.venus,.earth,.marsBarycentre,.jupiterBarycentre,
      .saturnBarycentre,.uranusBarycentre,.neptuneBarycentre] {
      #expect(DefaultEphemerisModel.model(for: body) == .vsop87)
    }
    #expect(DefaultEphemerisModel.model(for: .moon) == .meeusMoon)
    #expect(DefaultEphemerisModel.model(for: .plutoBarycentre) == .goffinPluto)
    #expect(DefaultEphemerisModel.model(for: .mars) == nil)
    #expect(DefaultEphemerisModel.model(for: .naif(123456)) == nil)
  }
  @Test func realKeplerFallbackAndBuiltInFactory() throws {
    let frame = Frame(name: "Earth-centred inertial", pseudoInertial: true)
    let provenance = try PhysicalDataProvenance(source: "explicit orbit seed", version: "1")
    let gm = try GravitationalParameter(cubicMetresPerSecondSquared: 3.986004418e14,
      convention: .newtonian, provenance: provenance)
    let initial = try CartesianState(epoch: .j2000, frame: frame,
      position: SIMD3(7e6,0,0), velocity: SIMD3(0,7500,0))
    let propagator = try TwoBodyPropagator(initialState: initial, centralBody: .earth, gravitationalParameter: gm)
    let body = BodyIdentifier.naif(123456)
    let kepler = try KeplerEphemeris(body: body, propagator: propagator,
      coverage: epoch(-100)...epoch(100), provenance: provenance)
    let chain = try FallbackEphemeris(kepler: kepler)
    for t in [100.0,-100,0,50] {
      #expect(try chain.checkedState(at: epoch(t)) == propagator.state(at: epoch(t)))
      #expect(try chain.selection(at: epoch(t)).tier == .kepler)
    }
    #expect(chain.metadata.origin == .earth)
    let tdb = try FairheadTDBModel()
    #expect(try DefaultEphemerisModel.makeProvider(for: body, coverage: epoch(-100)...epoch(100), tdbModel: tdb) == nil)
    let optionalBuiltin = try DefaultEphemerisModel.makeProvider(for: .earth,
      coverage: epoch(-100)...epoch(100), tdbModel: tdb, accuracy: .low)
    let builtin = try #require(optionalBuiltin)
    #expect(builtin.metadata.body == .earth)
    #expect(builtin.metadata.accuracy == .low)
    let expected = try VSOP87Ephemeris(body: earth, coverage: epoch(-100)...epoch(100), tdbModel: tdb, accuracy: .low)
    #expect(try builtin.checkedState(at: .j2000) == expected.checkedState(at: .j2000))
    #expect(throws: EphemerisError.frameMismatch) {
      try DefaultEphemerisModel.makeProvider(for: .moon, coverage: epoch(-100)...epoch(100), tdbModel: tdb)
    }
  }

}
