import Math
import Testing
@testable import CelMek

@Suite struct EphemerisOriginTests {
  private func source() throws -> PhysicalDataProvenance {
    try PhysicalDataProvenance(source: "synthetic trajectories and mass ratios",version: "1")
  }
  private struct Motion: EphemerisProvider {
    let metadata: EphemerisMetadata
    let base: Double
    let hasAcceleration: Bool
    func state(at epoch: Epoch) throws -> CartesianState {
      let t = try epoch.duration(since: .j2000).seconds
      return try CartesianState(epoch: epoch,frame: frame,position: SIMD3(base+base*t+0.5*base*t*t,0,0),
        velocity: SIMD3(base+base*t,0,0),acceleration: hasAcceleration ? SIMD3(base,0,0) : nil)
    }
  }
  private func motion(body: BodyIdentifier, frame: Frame, base: Double, origin: BodyIdentifier = .solarSystemBarycentre,
    acceleration: Bool = true, lower: Double = -10, upper: Double = 10, scale: TimeScale = .tdb) throws -> Motion {
    let m = try EphemerisMetadata(body: body,origin: origin,frame: frame,
      coverage: Epoch.j2000.advanced(by: CelMek.Duration(seconds: lower))...Epoch.j2000.advanced(by: CelMek.Duration(seconds: upper)),
      coverageBasis: .callerSelected,modelTimeScale: scale,velocityMethod: .analytic,
      accelerationMethod: acceleration ? .analytic : .unavailable,provenance: source())
    return Motion(metadata: m,base: base,hasAcceleration: acceleration)
  }
  @Test func recenteringAndFrameRoundTrip() throws {
    let tree = FrameTree()
    let a = try motion(body: .moon,frame: tree.root,base: 10)
    let b = try motion(body: .earth,frame: tree.root,base: 2,lower: -5,scale: .tt)
    let relative = try tree.makeOriginConvertedEphemeris(target: a,newOrigin: b,frameName: "Earth-centred")
    #expect(relative.metadata.origin == .earth)
    #expect(relative.metadata.body == .moon)
    #expect(relative.metadata.modelTimeScales == [.tt,.tdb])
    #expect(relative.metadata.coverage == b.metadata.coverage)
    let epoch = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: 2))
    let local = try relative.checkedState(at: epoch)
    #expect(local.position == SIMD3(40,0,0))
    #expect(local.velocity == SIMD3(24,0,0))
    #expect(local.acceleration == SIMD3(8,0,0))
    let throughTree = try a.state(at: epoch,in: relative.frame,using: tree)
    #expect(local == throughTree)
    let parent = try relative.state(at: epoch,in: tree.root,using: tree)
    let original = try a.state(at: epoch)
    #expect(parent == original)
    #expect(throws: EphemerisError.outsideCoverage) {
      try relative.state(at: Epoch.j2000.advanced(by: CelMek.Duration(seconds: -6)))
    }
  }
  @Test func barycentreIndependentMassWeightedValues() throws {
    let tree = FrameTree()
    let earth = try motion(body: .earth,frame: tree.root,base: 10)
    let moon = try motion(body: .moon,frame: tree.root,base: 2)
    let barycentre = try BarycentreEphemeris(bodyIdentifier: .earthMoonBarycentre,components: [
      BarycentreComponent(ephemeris: earth,relativeMass: 3,provenance: source()),
      BarycentreComponent(ephemeris: moon,relativeMass: 1,provenance: source())],provenance: source())
    let s = try barycentre.checkedState(at: .j2000)
    #expect(s.position == SIMD3(8,0,0))
    #expect(s.velocity == SIMD3(8,0,0))
    #expect(s.acceleration == SIMD3(8,0,0))
    let localEarth = try tree.makeOriginConvertedEphemeris(target: earth,newOrigin: barycentre,frameName: "EMB")
    #expect(try localEarth.state(at: .j2000).position == SIMD3(2,0,0))
    let scaled = try BarycentreEphemeris(bodyIdentifier: .earthMoonBarycentre,components: [
      BarycentreComponent(ephemeris: earth,relativeMass: 3e300,provenance: source()),
      BarycentreComponent(ephemeris: moon,relativeMass: 1e300,provenance: source())],provenance: source())
    #expect(try scaled.state(at: .j2000).position == s.position)
  }
  @Test func missingAccelerationRemainsUnknown() throws {
    let tree = FrameTree()
    let a = try motion(body: .moon,frame: tree.root,base: 10)
    let b = try motion(body: .earth,frame: tree.root,base: 2,acceleration: false)
    let relative = try tree.makeOriginConvertedEphemeris(target: a,newOrigin: b,frameName: "Earth")
    #expect(try relative.checkedState(at: .j2000).acceleration == nil)
    #expect(relative.metadata.accelerationMethod == .unavailable)
    #expect(throws: BodyFixedFrameError.accelerationUnavailable) {
      try tree.transform(from: tree.root,to: relative.frame,at: .j2000)
    }
    let barycentre = try BarycentreEphemeris(bodyIdentifier: .earthMoonBarycentre,components: [
      BarycentreComponent(ephemeris: a,relativeMass: 1,provenance: source()),
      BarycentreComponent(ephemeris: b,relativeMass: 1,provenance: source())],provenance: source())
    #expect(try barycentre.checkedState(at: .j2000).acceleration == nil)
  }
  @Test func incompatibleInputsAndDisjointCoverageFail() throws {
    let tree = FrameTree()
    let a = try motion(body: .earth,frame: tree.root,base: 1)
    let wrongOrigin = try motion(body: .moon,frame: tree.root,base: 1,origin: .sun)
    let wrongFrame = try motion(body: .moon,frame: Frame(name: "other",pseudoInertial: true),base: 1)
    let noOverlap = try motion(body: .moon,frame: tree.root,base: 1,lower: 11,upper: 20)
    for (b,error) in [(wrongOrigin,EphemerisOriginError.incompatibleOrigins),(wrongFrame,.incompatibleFrames),(noOverlap,.disjointCoverage)] {
      #expect(throws: error) { try tree.makeOriginConvertedEphemeris(target: a,newOrigin: b,frameName: "invalid") }
      #expect(throws: error) {
        try BarycentreEphemeris(bodyIdentifier: .earthMoonBarycentre,components: [
          BarycentreComponent(ephemeris: a,relativeMass: 1,provenance: source()),
          BarycentreComponent(ephemeris: b,relativeMass: 1,provenance: source())],provenance: source())
      }
    }
    for weight in [0.0,-1,.nan,.infinity] {
      #expect(throws: EphemerisOriginError.invalidMassRatio) {
        try BarycentreComponent(ephemeris: a,relativeMass: weight,provenance: source())
      }
    }
    let component = try BarycentreComponent(ephemeris: a,relativeMass: 1,provenance: source())
    #expect(throws: EphemerisOriginError.invalidComponents) {
      try BarycentreEphemeris(bodyIdentifier: .earthMoonBarycentre,components: [component,component],provenance: source())
    }
  }
}
