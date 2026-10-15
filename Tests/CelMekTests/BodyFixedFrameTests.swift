// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Math
import CelMek

@Suite
struct BodyFixedFrameTests {
  struct Centre: StateProvider {
    let frame: Frame
    var suppliesAcceleration = true
    func state(at epoch: Epoch) throws -> CartesianState {
      let t = epoch.offsetFromJ2000.seconds
      return try CartesianState(epoch: epoch, frame: frame,
        position: SIMD3(10 + 2*t + t*t, -3 + t, 4), velocity: SIMD3(2 + 2*t, 1, 0),
        acceleration: suppliesAcceleration ? SIMD3(2, 0, 0) : nil)
    }
  }
  func close(_ actual: SIMD3<Double>, _ expected: SIMD3<Double>, tolerance: Double = 1e-10) {
    for axis in 0..<3 { #expect(abs(actual[axis] - expected[axis]) < tolerance) }
  }

  @Test
  func centreIsStationaryAtOriginAndConversionsRoundTrip() throws {
    let tree = FrameTree()
    let centre = Centre(frame: tree.root)
    let spin = try FixedAxisOrientationProvider(referenceFrame: tree.root,
      axis: SIMD3(0, 0, 1), rate: 0.3, initialRotation: Quaternion(angle: 0.4, axis: SIMD3(1, 0, 0)))
    let child = try tree.makeBodyFixedFrame(name: "Body", centre: centre, orientation: spin)
    let epoch = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: 3))
    let origin = try centre.state(at: epoch)
    let mapped = try tree.convert(origin, to: child)
    close(mapped.position, .zero)
    close(mapped.velocity, .zero)
    close(mapped.acceleration!, .zero)
    let state = try CartesianState(epoch: epoch, frame: tree.root,
      position: SIMD3(3, 8, -2), velocity: SIMD3(2, -1, 4), acceleration: SIMD3(1, 3, -2))
    let restored = try tree.convert(tree.convert(state, to: child), to: tree.root)
    close(restored.position, state.position)
    close(restored.velocity, state.velocity)
    close(restored.acceleration!, state.acceleration!)
    #expect(child.parent === tree.root && !child.pseudoInertial)
  }

  @Test
  func movingOriginAndRotatingAxesMatchIndependentDifferences() throws {
    let tree = FrameTree()
    let child = try tree.makeBodyFixedFrame(name: "Body", centre: Centre(frame: tree.root),
      orientation: FixedAxisOrientationProvider(referenceFrame: tree.root, axis: SIMD3(0, 0, 1), rate: 0.3))
    func point(_ t: Double) -> SIMD3<Double> {
      let x = 5 + t - (10 + 2*t + t*t)
      let y = 2 - 2*t - (-3 + t)
      return SIMD3(cos(0.3*t)*x + sin(0.3*t)*y,
                   -sin(0.3*t)*x + cos(0.3*t)*y, -3)
    }
    let t = 0.7
    let epoch = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: t))
    let state = try CartesianState(epoch: epoch, frame: tree.root,
      position: SIMD3(5+t, 2-2*t, 1), velocity: SIMD3(1, -2, 0), acceleration: .zero)
    let result = try tree.convert(state, to: child)
    let h = 0.0001
    close(result.position, point(t))
    close(result.velocity, (point(t+h)-point(t-h))/(2*h), tolerance: 1e-7)
    close(result.acceleration!, (point(t+h)-2*point(t)+point(t-h))/(h*h), tolerance: 1e-6)
  }

  @Test
  func missingOriginAccelerationRequiresExplicitPolicy() throws {
    let tree = FrameTree()
    let centre = Centre(frame: tree.root, suppliesAcceleration: false)
    let spin = try FixedAxisOrientationProvider(referenceFrame: tree.root, axis: SIMD3(0, 0, 1), rate: 0)
    let strict = try BodyFixedTransformProvider(centre: centre, orientation: spin)
    #expect(throws: BodyFixedFrameError.accelerationUnavailable) { try strict.transform(at: .j2000) }
    let estimated = try BodyFixedTransformProvider(centre: centre, orientation: spin,
      accelerationPolicy: .estimate(step: CelMek.Duration(seconds: 1)))
    close(try estimated.transform(at: .j2000).acceleration, SIMD3(-2, 0, 0))
    let zero = try BodyFixedTransformProvider(centre: centre, orientation: spin, accelerationPolicy: .assumeZero)
    close(try zero.transform(at: .j2000).acceleration, .zero)
    #expect(throws: BodyFixedFrameError.invalidDifferenceStep) {
      try BodyFixedTransformProvider(centre: centre, orientation: spin, accelerationPolicy: .estimate(step: .zero))
    }
    let child = try tree.makeBodyFixedFrame(name: "Estimated", centre: centre, orientation: spin,
      accelerationPolicy: .estimate(step: CelMek.Duration(seconds: 1)))
    #expect(try tree.convert(centre.state(at: .j2000), to: child).acceleration == nil)
  }

  @Test
  func validatesProviderConventionsAndSupportsVSOP() throws {
    let tree = try FrameTree(root: CelestialFrames.vsop87E)
    let spin = try FixedAxisOrientationProvider(referenceFrame: tree.root, axis: SIMD3(0, 0, 1), rate: 0.001)
    let child = try tree.makeBodyFixedFrame(name: "Prescribed Earth axes", centre: earth, orientation: spin,
      accelerationPolicy: .estimate(step: CelMek.Duration(seconds: 60)))
    let centre = try earth.state(at: .j2000)
    let result = try tree.convert(centre, to: child)
    close(result.position, .zero, tolerance: 1e-4)
    close(result.velocity, .zero, tolerance: 1e-7)
    let other = Frame(name: tree.root.name, pseudoInertial: true)
    #expect(throws: BodyFixedFrameError.frameMismatch) {
      try BodyFixedTransformProvider(centre: Centre(frame: other), orientation: spin)
    }
    #expect(throws: FrameTransformError.disconnectedFrames) {
      try FrameTree().makeBodyFixedFrame(name: "Disconnected", centre: earth, orientation: spin)
    }
  }
  struct StaticCentre: StateProvider {
    let frame: Frame
    let result: CartesianState
    func state(at epoch: Epoch) throws -> CartesianState { result }
  }
  struct StaticOrientation: OrientationProvider {
    let referenceFrame: Frame
    let result: BodyOrientation
    func orientation(at epoch: Epoch) throws -> BodyOrientation { result }
  }

  @Test
  func rejectsIncorrectReturnedEpochsAndFrames() throws {
    let frame = Frame(name: "Reference", pseudoInertial: true)
    let other = Frame(name: "Reference", pseudoInertial: true)
    let spin = try FixedAxisOrientationProvider(referenceFrame: frame, axis: SIMD3(0, 0, 1), rate: 0)
    let wrongState = try CartesianState(epoch: .j2000, frame: other, position: .zero, velocity: .zero,
      acceleration: .zero)
    let wrongCentre = StaticCentre(frame: frame, result: wrongState)
    #expect(throws: BodyFixedFrameError.frameMismatch) {
      try BodyFixedTransformProvider(centre: wrongCentre, orientation: spin).transform(at: .j2000)
    }
    let centre = Centre(frame: frame)
    let later = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: 1))
    let stale = StaticCentre(frame: frame, result: try centre.state(at: .j2000))
    #expect(throws: FrameTransformError.epochMismatch) {
      try BodyFixedTransformProvider(centre: stale, orientation: spin).transform(at: later)
    }
    let wrongOrientation = StaticOrientation(referenceFrame: frame,
      result: try BodyOrientation(referenceFrame: other, epoch: .j2000, rotation: .identity))
    #expect(throws: BodyFixedFrameError.frameMismatch) {
      try BodyFixedTransformProvider(centre: centre, orientation: wrongOrientation).transform(at: .j2000)
    }
    let staleOrientation = StaticOrientation(referenceFrame: frame,
      result: try BodyOrientation(referenceFrame: frame, epoch: .j2000, rotation: .identity))
    #expect(throws: FrameTransformError.epochMismatch) {
      try BodyFixedTransformProvider(centre: centre, orientation: staleOrientation).transform(at: later)
    }
    // A provider supplying only the central sample cannot satisfy estimation.
    let fixedSample = StaticCentre(frame: frame,
      result: try Centre(frame: frame, suppliesAcceleration: false).state(at: .j2000))
    #expect(throws: FrameTransformError.epochMismatch) {
      try BodyFixedTransformProvider(centre: fixedSample, orientation: spin,
        accelerationPolicy: .estimate(step: CelMek.Duration(seconds: 1))).transform(at: .j2000)
    }
  }

}
