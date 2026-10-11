import Testing
import Math
@testable import CelMek

@Suite struct EphemerisProviderTests {
  private func metadata(frame: Frame, body: BodyIdentifier = .earth,
    acceleration: EphemerisDerivativeMethod = .unavailable) throws -> EphemerisMetadata {
    try EphemerisMetadata(body: body,origin: .solarSystemBarycentre,frame: frame,
      coverage: Epoch.j2000.advanced(by: CelMek.Duration(seconds: -10))...Epoch.j2000.advanced(by: CelMek.Duration(seconds: 10)),
      coverageBasis: .callerSelected,modelTimeScale: .tdb,velocityMethod: .analytic,
      accelerationMethod: acceleration,provenance: PhysicalDataProvenance(source: "test trajectory",version: "1"))
  }
  private struct Linear: StateProvider {
    let frame: Frame
    var wrongEpoch = false
    var outputFrame: Frame? = nil
    var acceleration: SIMD3<Double>? = nil
    func state(at epoch: Epoch) throws -> CartesianState {
      let t = try epoch.duration(since: .j2000).seconds
      return try CartesianState(epoch: wrongEpoch ? epoch.advanced(by: CelMek.Duration(seconds: 1)) : epoch,
        frame: outputFrame ?? frame,position: SIMD3(1+t,2+2*t,3+3*t),velocity: SIMD3(1,2,3),acceleration: acceleration)
    }
  }
  @Test func queryAndCoverageThroughBothInterfaces() throws {
    let frame = Frame(name: "barycentric",pseudoInertial: true)
    let p = try StateProviderEphemeris(provider: Linear(frame: frame),metadata: metadata(frame: frame))
    let collection = try EphemerisCollection(providers: [p])
    for t in [-10.0,0,10] {
      let epoch = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: t))
      let state = try collection.state(of: .earth,at: epoch)
      #expect(state.position == SIMD3(1+t,2+2*t,3+3*t))
      let direct = try p.state(at: epoch)
      #expect(state == direct)
      #expect(state.frame === frame)
    }
    #expect(collection.bodies == [.earth])
    #expect(try collection.metadata(for: .earth).origin == .solarSystemBarycentre)
    #expect(throws: EphemerisError.outsideCoverage) {
      try p.state(at: Epoch.j2000.advanced(by: CelMek.Duration(seconds: 11)))
    }
    #expect(throws: EphemerisError.unsupportedBody) { try collection.state(of: .moon,at: .j2000) }
    #expect(throws: EphemerisError.duplicateBody) { try EphemerisCollection(providers: [p,p]) }
  }
  @Test func providerContractFailures() throws {
    let frame = Frame(name: "one",pseudoInertial: true), other = Frame(name: "two",pseudoInertial: true)
    let m = try metadata(frame: frame)
    #expect(throws: EphemerisError.frameMismatch) { try StateProviderEphemeris(provider: Linear(frame: other),metadata: m) }
    let wrongTime = try StateProviderEphemeris(provider: Linear(frame: frame,wrongEpoch: true),metadata: m)
    #expect(throws: FrameTransformError.epochMismatch) { try wrongTime.state(at: .j2000) }
    let wrongFrame = try StateProviderEphemeris(provider: Linear(frame: frame,outputFrame: other),metadata: m)
    #expect(throws: EphemerisError.frameMismatch) { try wrongFrame.state(at: .j2000) }
    let unexpected = try StateProviderEphemeris(provider: Linear(frame: frame,acceleration: .zero),metadata: m)
    #expect(throws: EphemerisError.accelerationMismatch) { try unexpected.state(at: .j2000) }
    let missing = try StateProviderEphemeris(provider: Linear(frame: frame),metadata: metadata(frame: frame,acceleration: .analytic))
    #expect(throws: EphemerisError.accelerationMismatch) { try missing.state(at: .j2000) }
  }
  @Test func existingVSOPStateIntegration() throws {
    let m = try metadata(frame: earth.frame)
    let wrapped = try StateProviderEphemeris(provider: earth,metadata: m)
    let checked = try wrapped.checkedState(at: .j2000)
    let native = try earth.state(at: .j2000)
    #expect(checked == native)
    let tree = try CelestialFrames.makeTree()
    let converted = try wrapped.state(at: .j2000,in: CelestialFrames.icrf,using: tree)
    #expect(converted.frame === CelestialFrames.icrf)
    #expect(converted.acceleration == nil)
  }
  @Test func metadataRejectsInvalidDerivativeClaims() throws {
    let frame = Frame(name: "reference",pseudoInertial: true)
    #expect(throws: EphemerisError.invalidDerivativeMethod) {
      try metadata(frame: frame,acceleration: .finiteDifference(step: .zero))
    }
    #expect(throws: EphemerisError.invalidDerivativeMethod) {
      try metadata(frame: frame,acceleration: .supplied(description: " "))
    }
  }

  private final class CustomProvider: EphemerisProvider {
    let metadata: EphemerisMetadata
    var evaluations = 0
    init(metadata: EphemerisMetadata) { self.metadata = metadata }
    func state(at epoch: Epoch) throws -> CartesianState {
      evaluations += 1
      return try CartesianState(epoch: epoch,frame: metadata.frame,position: .zero,velocity: .zero)
    }
  }
  @Test func collectionValidatesCustomProviderBeforeEvaluation() throws {
    let frame = Frame(name: "reference",pseudoInertial: true)
    let p = CustomProvider(metadata: try metadata(frame: frame,acceleration: .analytic))
    let collection = try EphemerisCollection(providers: [p])
    #expect(throws: EphemerisError.outsideCoverage) {
      try collection.state(of: .earth,at: Epoch.j2000.advanced(by: CelMek.Duration(seconds: 11)))
    }
    #expect(p.evaluations == 0)
    #expect(throws: EphemerisError.accelerationMismatch) { try collection.state(of: .earth,at: .j2000) }
    #expect(p.evaluations == 1)
  }
}
