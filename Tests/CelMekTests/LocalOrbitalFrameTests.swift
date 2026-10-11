import Foundation
import Testing
import Math
@testable import CelMek

@Suite struct LocalOrbitalFrameTests {
  private struct Orbit: StateProvider {
    let frame: Frame
    var offset = SIMD3<Double>(1000,2000,3000)
    var drift = SIMD3<Double>(10,20,30)
    var radialRate = 0.0
    func state(at epoch: Epoch) throws -> CartesianState {
      let t = try epoch.duration(since: .j2000).seconds
      let n = 0.001, a = n*t, r = 7e6+radialRate*t
      let u = SIMD3(cos(a),sin(a),0), tangent = SIMD3(-sin(a),cos(a),0)
      return try CartesianState(epoch: epoch,frame: frame,position: offset+drift*t+r*u,
        velocity: drift+radialRate*u+r*n*tangent,acceleration: -r*n*n*u+2*radialRate*n*tangent)
    }
  }
  private struct Centre: StateProvider {
    let frame: Frame
    func state(at epoch: Epoch) throws -> CartesianState {
      let t = try epoch.duration(since: .j2000).seconds
      return try CartesianState(epoch: epoch,frame: frame,position: SIMD3(1000,2000,3000)+SIMD3(10,20,30)*t,
        velocity: SIMD3(10,20,30),acceleration: .zero)
    }
  }
  @Test func circularAxesRatesAndSpacecraftOrigin() throws {
    let tree = FrameTree(), orbit = Orbit(frame: tree.root), centre = Centre(frame: tree.root)
    for axes in [LocalOrbitalAxes.rtn,.tnw,.lvlhCCSDS] {
      let p = try LocalOrbitalOrientationProvider(spacecraft: orbit,centralBody: centre,axes: axes,
        derivativeStep: CelMek.Duration(seconds: 0.5))
      let attitude = try p.orientation(at: .j2000)
      let expected: [[Double]]
      let omega: SIMD3<Double>
      switch axes {
      case .rtn: expected = [[1,0,0],[0,1,0],[0,0,1]]; omega = SIMD3(0,0,-0.001)
      case .tnw: expected = [[0,1,0],[-1,0,0],[0,0,1]]; omega = SIMD3(0,0,-0.001)
      case .lvlhCCSDS: expected = [[0,1,0],[0,0,-1],[-1,0,0]]; omega = SIMD3(0,0.001,0)
      }
      for i in 0..<3 { for j in 0..<3 { #expect(abs(attitude.rotation.matrix[i,j]-expected[i][j]) < 1e-14) } }
      for i in 0..<3 {
        #expect(abs(attitude.angularVelocity[i]-omega[i]) < 1e-15)
        #expect(abs(attitude.angularAcceleration[i]) < 1e-14)
      }
      let frame = try tree.makeLocalOrbitalFrame(name: "local",spacecraft: orbit,centralBody: centre,
        axes: axes,derivativeStep: CelMek.Duration(seconds: 0.5))
      let state = try orbit.state(at: .j2000)
      let local = try tree.transform(from: tree.root,to: frame,at: .j2000).applying(to: state,destination: frame)
      #expect(local.position == .zero)
      #expect(local.velocity == .zero)
      #expect(local.acceleration == .zero)
    }
  }

  @Test func nonCircularVelocityAxesAndMountedAttitude() throws {
    let tree = FrameTree(), orbit = Orbit(frame: tree.root,radialRate: 1000), centre = Centre(frame: tree.root)
    let p = try LocalOrbitalOrientationProvider(spacecraft: orbit,centralBody: centre,axes: .tnw,
      derivativeStep: CelMek.Duration(seconds: 0.5))
    let a = try p.orientation(at: .j2000)
    #expect(abs(a.rotation.matrix[0,0]-1/sqrt(50)) < 1e-14)
    #expect(abs(a.rotation.matrix[0,1]-7/sqrt(50)) < 1e-14)
    // atan2(r*n, radialRate): frame angle and its exact first/second derivatives.
    let n = 0.001, r = 7e6, k = 1000.0, denominator = k*k+r*r*n*n
    let rate = n+n*k*k/denominator
    let acceleration = -2*n*n*n*k*k*k*r/(denominator*denominator)
    #expect(abs(a.angularVelocity.z+rate) < 2e-14)
    #expect(abs(a.angularAcceleration.z+acceleration) < 2e-14)
    let mount = Quaternion(angle: .pi/2,axis: SIMD3(1,0,0))
    let mounted = try LocalOrbitalOrientationProvider(spacecraft: orbit,centralBody: centre,axes: .tnw,
      derivativeStep: CelMek.Duration(seconds: 0.5),mountingRotation: mount)
    let b = try mounted.orientation(at: .j2000)
    for i in 0..<3 { #expect(abs(b.angularVelocity[i]-(mount.matrix*a.angularVelocity)[i]) < 1e-15) }
    let frame = try tree.makeSpacecraftFrame(name: "mounted",spacecraft: orbit,attitude: mounted)
    let transform = try tree.transform(from: tree.root,to: frame,at: .j2000)
    let local = try transform.applying(to: orbit.state(at: .j2000),destination: frame)
    #expect(local.position == .zero)
    #expect(local.velocity == .zero)
  }
  private struct Static: StateProvider {
    let frame: Frame
    let r: SIMD3<Double>, v: SIMD3<Double>
    var wrongEpoch = false
    func state(at epoch: Epoch) throws -> CartesianState {
      try CartesianState(epoch: wrongEpoch ? epoch.advanced(by: CelMek.Duration(seconds: 1)) : epoch,
        frame: frame,position: r,velocity: v)
    }
  }
  @Test func degenerateInputsAndConfiguration() throws {
    let frame = Frame(name: "reference",pseudoInertial: true)
    for (r,v,error) in [
      (SIMD3<Double>.zero,SIMD3(0,1,0),LocalOrbitalFrameError.zeroPosition),
      (SIMD3(1,0,0),SIMD3<Double>.zero,.zeroVelocity),
      (SIMD3(1,0,0),SIMD3(1,0,0),.undefinedOrbitalPlane),
      (SIMD3(1,0,0),SIMD3(1,1e-14,0),.undefinedOrbitalPlane)] {
      let p = try LocalOrbitalOrientationProvider(spacecraft: Static(frame: frame,r: r,v: v),axes: .rtn,
        derivativeStep: CelMek.Duration(seconds: 1))
      #expect(throws: error) { try p.orientation(at: .j2000) }
    }
    let state = Static(frame: frame,r: SIMD3(1,0,0),v: SIMD3(0,1,0),wrongEpoch: true)
    let p = try LocalOrbitalOrientationProvider(spacecraft: state,axes: .rtn,derivativeStep: CelMek.Duration(seconds: 1))
    #expect(throws: FrameTransformError.epochMismatch) { try p.orientation(at: .j2000) }
    #expect(throws: LocalOrbitalFrameError.invalidConfiguration) {
      try LocalOrbitalOrientationProvider(spacecraft: state,axes: .rtn,derivativeStep: .zero)
    }
    let rotating = Static(frame: Frame(name: "rotating",pseudoInertial: false),r: SIMD3(1,0,0),v: SIMD3(0,1,0))
    #expect(throws: LocalOrbitalFrameError.nonInertialReference) {
      try LocalOrbitalOrientationProvider(spacecraft: rotating,axes: .rtn,derivativeStep: CelMek.Duration(seconds: 1))
    }
    #expect(throws: LocalOrbitalFrameError.frameMismatch) {
      try LocalOrbitalOrientationProvider(spacecraft: state,centralBody: rotating,axes: .rtn,
        derivativeStep: CelMek.Duration(seconds: 1))
    }
  }

  private struct RotatedOrbit: StateProvider {
    let frame: Frame
    let rotation: Quaternion
    func state(at epoch: Epoch) throws -> CartesianState {
      let s = try Orbit(frame: frame,offset: .zero,drift: .zero).state(at: epoch)
      return try CartesianState(epoch: epoch,frame: frame,position: rotation.matrix*s.position,
        velocity: rotation.matrix*s.velocity,acceleration: s.acceleration.map { rotation.matrix*$0 })
    }
  }
  @Test func inclinedHalfTurnAndTransformedStateDerivatives() throws {
    let tree = FrameTree()
    let orbit = RotatedOrbit(frame: tree.root,rotation: Quaternion(angle: .pi,axis: SIMD3(1,2,3)))
    let frame = try tree.makeLocalOrbitalFrame(name: "inclined",spacecraft: orbit,axes: .rtn,
      derivativeStep: CelMek.Duration(seconds: 0.5))
    let t = try tree.transform(from: tree.root,to: frame,at: .j2000)
    let expected = orbit.rotation.matrix
    for i in 0..<3 { for j in 0..<3 { #expect(abs(t.rotation.matrix[i,j]-expected[j,i]) < 1e-14) } }
    // A stationary parent point moves in the rotating, translating local frame.
    func local(_ seconds: Double) throws -> CartesianState {
      let epoch = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: seconds))
      let point = try CartesianState(epoch: epoch,frame: tree.root,position: SIMD3(2e6,-1e6,3e6),
        velocity: .zero,acceleration: .zero)
      return try tree.transform(from: tree.root,to: frame,at: epoch).applying(to: point,destination: frame)
    }
    let centre = try local(0), before = try local(-0.1), after = try local(0.1)
    let v = (after.position-before.position)/0.2
    let a = (after.velocity-before.velocity)/0.2
    for i in 0..<3 {
      #expect(abs(v[i]-centre.velocity[i]) < 1e-5)
      #expect(abs(a[i]-centre.acceleration![i]) < 1e-7)
    }
  }
  private enum CoverageError: Error { case outside }
  private struct LimitedOrbit: StateProvider {
    let frame: Frame
    func state(at epoch: Epoch) throws -> CartesianState {
      guard abs(try epoch.duration(since: .j2000).seconds) < 1 else { throw CoverageError.outside }
      return try Orbit(frame: frame).state(at: epoch)
    }
  }
  @Test func derivativeStencilRequiresNeighbourCoverage() throws {
    let frame = Frame(name: "reference",pseudoInertial: true)
    let p = try LocalOrbitalOrientationProvider(spacecraft: LimitedOrbit(frame: frame),axes: .rtn,
      derivativeStep: CelMek.Duration(seconds: 1))
    #expect(throws: CoverageError.outside) { try p.orientation(at: .j2000) }
  }
}
