// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Math
import CelMek

@Suite
struct FrameTransformTests {
  func close(_ actual: SIMD3<Double>, _ expected: SIMD3<Double>, tolerance: Double = 1e-10) {
    for axis in 0..<3 { #expect(abs(actual[axis] - expected[axis]) < tolerance) }
  }

  @Test
  func rotatingCoordinatesIncludeCentripetalAndCoriolisTerms() throws {
    let tree = FrameTree()
    let child = try tree.makeFrame(parent: tree.root, name: "Rotating",
      provider: UniformRotationProvider(axis: SIMD3(0, 0, 2), rate: 0.5))
    let state = try CartesianState(epoch: .j2000, frame: tree.root,
      position: SIMD3(2, 0, 0), velocity: SIMD3(0, 3, 0), acceleration: .zero)
    let converted = try tree.convert(state, to: child)
    close(converted.position, SIMD3(2, 0, 0))
    close(converted.velocity, SIMD3(0, 2, 0))
    close(converted.acceleration!, SIMD3(2.5, 0, 0))
    let restored = try tree.convert(converted, to: tree.root)
    close(restored.position, state.position)
    close(restored.velocity, state.velocity)
    close(restored.acceleration!, state.acceleration!)
    let quarterTurn = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: .pi))
    let transform = try tree.transform(from: tree.root, to: child, at: quarterTurn)
    close(transform.rotation.matrix * SIMD3(2, 0, 0), SIMD3(0, -2, 0))
  }

  @Test
  func compositionAndInverseRecoverAllDerivatives() throws {
    let frame = Frame(name: "Source", pseudoInertial: true)
    let destination = Frame(name: "Destination", pseudoInertial: false)
    let a = try FrameTransform(epoch: .j2000, rotation: Quaternion(angle: 0.7, axis: SIMD3(1, 2, 3)),
      translation: SIMD3(5, -2, 7), velocity: SIMD3(2, 3, -1), acceleration: SIMD3(-1, 2, 4),
      angularVelocity: SIMD3(0.1, -0.2, 0.3), angularAcceleration: SIMD3(-0.2, 0.05, 0.1))
    let b = try FrameTransform(epoch: .j2000, rotation: Quaternion(angle: -0.4, axis: SIMD3(2, -1, 1)),
      translation: SIMD3(-3, 8, 1), velocity: SIMD3(4, -1, 2), acceleration: SIMD3(3, 1, -2),
      angularVelocity: SIMD3(-0.3, 0.2, 0.1), angularAcceleration: SIMD3(0.1, 0.2, -0.1))
    let state = try CartesianState(epoch: .j2000, frame: frame,
      position: SIMD3(3, 7, -4), velocity: SIMD3(1, -2, 5), acceleration: SIMD3(-2, 1, 3))
    let composed = try a.followed(by: b)
    let sequential = try b.applying(to: a.applying(to: state, destination: destination), destination: destination)
    let direct = try composed.applying(to: state, destination: destination)
    close(direct.position, sequential.position)
    close(direct.velocity, sequential.velocity)
    close(direct.acceleration!, sequential.acceleration!)
    let restored = try composed.inverted().applying(to: direct, destination: frame)
    close(restored.position, state.position)
    close(restored.velocity, state.velocity)
    close(restored.acceleration!, state.acceleration!)
  }

  @Test
  func derivativesMatchIndependentMovingTrajectory() throws {
    // Independent explicit trajectory: rotate a polynomial inertial path into
    // clockwise-turning axes, then add a polynomial translation. Angle has
    // nonzero acceleration to exercise the Euler term.
    let frame = Frame(name: "Test", pseudoInertial: true)
    func point(_ t: Double) -> SIMD3<Double> {
      let theta = -(0.3*t + 0.02*t*t)
      let x = 2 + 3*t + 0.5*t*t
      let y = -1 + 2*t - t*t
      return SIMD3(cos(theta)*x - sin(theta)*y + t*t,
                   sin(theta)*x + cos(theta)*y - 3*t, 4 + t + 2*t*t)
    }
    let t = 0.7
    let epoch = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: t))
    let transform = try FrameTransform(epoch: epoch,
      rotation: Quaternion(angle: -(0.3*t + 0.02*t*t), axis: SIMD3(0, 0, 1)),
      translation: SIMD3(t*t, -3*t, 2*t*t), velocity: SIMD3(2*t, -3, 4*t),
      acceleration: SIMD3(2, 0, 4), angularVelocity: SIMD3(0, 0, -(0.3 + 0.04*t)),
      angularAcceleration: SIMD3(0, 0, -0.04))
    let state = try CartesianState(epoch: epoch, frame: frame,
      position: SIMD3(2 + 3*t + 0.5*t*t, -1 + 2*t - t*t, 4 + t),
      velocity: SIMD3(3 + t, 2 - 2*t, 1), acceleration: SIMD3(1, -2, 0))
    let result = try transform.applying(to: state, destination: frame)
    let h = 0.0001
    close(result.position, point(t))
    close(result.velocity, (point(t+h)-point(t-h))/(2*h), tolerance: 1e-7)
    close(result.acceleration!, (point(t+h)-2*point(t)+point(t-h))/(h*h), tolerance: 1e-6)
  }

  @Test
  func treeHandlesSiblingsIdentityOwnershipAndUnknownAcceleration() throws {
    let tree = FrameTree()
    let parent = try tree.makeFrame(parent: tree.root, name: "Offset",
      provider: FixedTransformProvider(translation: SIMD3(5, 0, 0)))
    let a = try tree.makeFrame(parent: parent, name: "Same",
      provider: FixedTransformProvider(translation: SIMD3(0, 2, 0)))
    let b = try tree.makeFrame(parent: parent, name: "Same",
      provider: FixedTransformProvider(translation: SIMD3(0, 0, 3)))
    #expect(a.parent === parent && a.depth == 2)
    let state = try CartesianState(epoch: .j2000, frame: a, position: SIMD3(1, 2, 3), velocity: .zero)
    let converted = try tree.convert(state, to: b)
    close(converted.position, SIMD3(1, 0, 6))
    #expect(converted.acceleration == nil)
    #expect(try tree.convert(state, to: a) == state)
    let outsider = Frame(name: "Same", pseudoInertial: true)
    #expect(throws: FrameTransformError.disconnectedFrames) { try tree.convert(state, to: outsider) }
    #expect(throws: FrameTransformError.disconnectedFrames) {
      try tree.makeFrame(parent: outsider, name: "Bad", provider: FixedTransformProvider())
    }
    // A child retains its parent even without a tree or another owner.
    var temporary: Frame? = Frame(name: "Retained", pseudoInertial: true)
    weak var weakParent = temporary
    let child = Frame(parent: temporary, name: "Child", pseudoInertial: false)
    temporary = nil
    #expect(weakParent != nil && child.parent === weakParent)
  }

  @Test
  func releasingTreeDoesNotCycleThroughFramesOrBuiltInProviders() throws {
    var tree: FrameTree? = FrameTree()
    weak var observedTree = tree
    weak var observedRoot = tree!.root
    var child: Frame? = try tree!.makeFrame(parent: tree!.root, name: "Child",
      provider: UniformRotationProvider(axis: SIMD3(0, 0, 1), rate: 1))
    tree = nil
    #expect(observedTree == nil)
    #expect(observedRoot != nil && child!.parent === observedRoot)
    child = nil
    #expect(observedRoot == nil)
  }

  struct FailingProvider: TransformProvider {
    func transform(at epoch: Epoch) throws -> FrameTransform { throw FrameTransformError.nonFinite }
  }
  struct WrongEpochProvider: TransformProvider {
    func transform(at epoch: Epoch) throws -> FrameTransform { .identity(at: .j2000) }
  }

  @Test
  func sharedEdgesAreSkippedAndProviderErrorsPropagate() throws {
    let tree = FrameTree()
    let parent = try tree.makeFrame(parent: tree.root, name: "Unavailable", provider: FailingProvider())
    let child = try tree.makeFrame(parent: parent, name: "Available", provider: FixedTransformProvider())
    _ = try tree.transform(from: parent, to: child, at: .j2000)
    #expect(throws: FrameTransformError.nonFinite) { try tree.transform(from: tree.root, to: child, at: .j2000) }
    let wrong = try tree.makeFrame(parent: tree.root, name: "Wrong epoch", provider: WrongEpochProvider())
    let later = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: 1))
    #expect(throws: FrameTransformError.epochMismatch) { try tree.transform(from: tree.root, to: wrong, at: later) }
  }

  @Test
  func validatesInputsAndAdoptsCelestialRoot() throws {
    #expect(throws: FrameTransformError.invalidRotation) { try FrameTransform(epoch: .j2000, rotation: Quaternion()) }
    #expect(throws: FrameTransformError.nonFinite) {
      try FrameTransform(epoch: .j2000, angularVelocity: SIMD3(.infinity, 0, 0))
    }
    #expect(throws: FrameTransformError.nonFinite) { try UniformRotationProvider(axis: .zero, rate: 1) }
    let tree = try FrameTree(root: CelestialFrames.vsop87E)
    let child = try tree.makeFrame(parent: tree.root, name: "Fixed",
      provider: FixedTransformProvider(translation: SIMD3(100, 0, 0)))
    let state = try earth.state(at: .j2000)
    let converted = try tree.convert(state, to: child)
    close(converted.position, state.position + SIMD3(100, 0, 0), tolerance: 1e-5)
    close(converted.velocity, state.velocity)
    let later = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: 1))
    #expect(throws: FrameTransformError.epochMismatch) {
      try FrameTransform.identity(at: later).applying(to: state, destination: child)
    }
  }
}
