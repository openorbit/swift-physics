// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct RotatingDynamicsTests {
    private func epoch(_ t: Double) throws -> Epoch { try Epoch.j2000.advanced(by: Duration(seconds: t)) }
    private func options() throws -> NumericalPropagationOptions {
        try NumericalPropagationOptions(absolutePositionToleranceMetres: 1e-10,
            absoluteVelocityToleranceMetresPerSecond: 1e-11, relativeTolerance: 1e-12,
            initialStepSeconds: 0.1, maximumStepSeconds: 0.5)
    }
    private func distance(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> Double { sqrt(((a-b)*(a-b)).sum()) }
    private struct AcceleratingOrigin: TransformProvider {
        func transform(at epoch: Epoch) throws -> FrameTransform {
            let t = try epoch.duration(since: .j2000).seconds
            return try FrameTransform(epoch: epoch, translation: SIMD3(-0.15*t*t,0,0),
                velocity: SIMD3(-0.3*t,0,0), acceleration: SIMD3(-0.3,0,0))
        }
    }
    private struct AcceleratingRotation: TransformProvider {
        func transform(at epoch: Epoch) throws -> FrameTransform {
            let t = try epoch.duration(since: .j2000).seconds
            return try FrameTransform(epoch: epoch, rotation: Quaternion(angle: -0.2*t-0.01*t*t, axis: SIMD3(0,0,1)),
                angularVelocity: SIMD3(0,0,-0.2-0.02*t), angularAcceleration: SIMD3(0,0,-0.02))
        }
    }
    @Test func translatingAndAcceleratingRotationMatchesAnalyticInertialMotion() throws {
        let tree = FrameTree()
        let translated = try tree.makeFrame(parent: tree.root, name: "Accelerating origin", provider: AcceleratingOrigin())
        let rotating = try tree.makeFrame(parent: translated, name: "Accelerating rotation", provider: AcceleratingRotation())
        let inertial = try SpacecraftDynamics(models: [], context: ForceEvaluationContext(frame: tree.root, origin: .earth))
        let dynamics = try SpacecraftDynamics(rotatingFrame: rotating, origin: .custom(namespace: "test", name: "moving-origin"), inertialDynamics: inertial, using: tree)
        let original = try CartesianState(epoch: .j2000, frame: tree.root, position: SIMD3(3,1,0), velocity: SIMD3(0.2,0.5,0), acceleration: .zero)
        let initial = SpacecraftState(cartesian: try tree.convert(original, to: rotating), mass: try SpacecraftMassProperties(dryMassKilograms: 1))
        let provider = try NumericalSpacecraftPropagator(initialState: initial, dynamics: dynamics, options: options())
        for t in [-10.0,10.0] {
            let end = try epoch(t)
            let reference = try tree.convert(CartesianState(epoch: end, frame: tree.root,
                position: original.position+original.velocity*t, velocity: original.velocity, acceleration: .zero), to: rotating)
            let actual = try provider.spacecraftState(at: end)
            #expect(distance(actual.cartesian.position,reference.position) < 1e-8)
            #expect(distance(actual.cartesian.velocity,reference.velocity) < 1e-8)
            #expect(distance(try #require(actual.cartesian.acceleration),try #require(reference.acceleration)) < 1e-8)
            #expect(actual.mass == initial.mass)
        }
    }
    @Test func centrifugalAndCoriolisSigns() throws {
        let tree = FrameTree()
        let frame = try tree.makeFrame(parent: tree.root, name: "Spinning", provider: UniformRotationProvider(axis: SIMD3(0,0,1), rate: 2))
        let inertial = try SpacecraftDynamics(models: [], context: ForceEvaluationContext(frame: tree.root, origin: .earth))
        let dynamics = try SpacecraftDynamics(rotatingFrame: frame, origin: .earth, inertialDynamics: inertial, using: tree)
        let stationary = SpacecraftState(cartesian: try CartesianState(epoch: .j2000, frame: frame,
            position: SIMD3(3,0,0), velocity: .zero), mass: try SpacecraftMassProperties(dryMassKilograms: 1))
        #expect(distance(try dynamics.derivative(at: stationary).velocityRateMetresPerSecondSquared,SIMD3(12,0,0)) < 1e-14)
        let moving = SpacecraftState(cartesian: try CartesianState(epoch: .j2000, frame: frame,
            position: SIMD3(3,0,0), velocity: SIMD3(0,1,0)), mass: stationary.mass)
        #expect(distance(try dynamics.derivative(at: moving).velocityRateMetresPerSecondSquared,SIMD3(16,0,0)) < 1e-14)
    }
    @Test func centralGravityOrbitMatchesTwoBodyInRotatingAxes() throws {
        let tree = FrameTree()
        let frame = try tree.makeFrame(parent: tree.root, name: "Orbit view", provider: UniformRotationProvider(axis: SIMD3(0,0,1), rate: 0.3))
        let provenance = try PhysicalDataProvenance(source: "Synthetic GM", version: "1")
        let gm = try GravitationalParameter(cubicMetresPerSecondSquared: 10, convention: .newtonian, provenance: provenance)
        let inertial = try SpacecraftDynamics(models: [CentralGravity(body: .earth, gravitationalParameter: gm)],
            context: ForceEvaluationContext(frame: tree.root, origin: .earth))
        let dynamics = try SpacecraftDynamics(rotatingFrame: frame, origin: .earth, inertialDynamics: inertial, using: tree)
        let seed = try CartesianState(epoch: .j2000, frame: tree.root, position: SIMD3(2,0,0), velocity: SIMD3(0,sqrt(5),0))
        let initial = SpacecraftState(cartesian: try tree.convert(seed, to: frame), mass: try SpacecraftMassProperties(dryMassKilograms: 1))
        let provider = try NumericalSpacecraftPropagator(initialState: initial, dynamics: dynamics, options: options())
        let t = 20.0, n = sqrt(10.0/8), angle = (n-0.3)*t
        let actual = try provider.state(at: epoch(t))
        #expect(distance(actual.position,SIMD3(2*cos(angle),2*sin(angle),0)) < 1e-7)
        #expect(distance(actual.velocity,SIMD3(-2*(n-0.3)*sin(angle),2*(n-0.3)*cos(angle),0)) < 1e-7)
        #expect(throws: FrameTransformError.disconnectedFrames) {
            try SpacecraftDynamics(rotatingFrame: Frame(name: "Unregistered", pseudoInertial: false), origin: .earth, inertialDynamics: inertial, using: tree)
        }
    }
    @Test func fuelCutoffPreservesFrameDynamics() throws {
        let tree = FrameTree()
        let frame = try tree.makeFrame(parent: tree.root, name: "Powered view",
            provider: UniformRotationProvider(axis: SIMD3(0,0,1), rate: 0.01))
        let engine = try FiniteThrust(frame: tree.root, thrustNewtons: 2,
            effectiveExhaustVelocityMetresPerSecond: 100, direction: SIMD3(1,0,0))
        let inertial = try SpacecraftDynamics(models: [engine],
            context: ForceEvaluationContext(frame: tree.root, origin: .earth))
        let dynamics = try SpacecraftDynamics(rotatingFrame: frame, origin: .earth, inertialDynamics: inertial, using: tree)
        let initial = SpacecraftState(cartesian: try CartesianState(epoch: .j2000, frame: frame, position: .zero, velocity: .zero),
            mass: try SpacecraftMassProperties(dryMassKilograms: 9, propellantMassKilograms: 1))
        let result = try FuelLimitedNumericalPropagator(initialState: initial, dynamics: dynamics,
            numericalOptions: options(), eventOptions: PropagationEventOptions(timeToleranceSeconds: 1e-8))
            .propagation(at: epoch(60))
        let recovered = try tree.convert(result.state.cartesian, to: tree.root)
        let speed = 100*log(10.0/9)
        let position = 100*(50-450*log(10.0/9))+10*speed
        #expect(abs(recovered.velocity.x-speed) < 1e-7)
        #expect(abs(recovered.position.x-position) < 1e-6)
        #expect(distance(try #require(recovered.acceleration),.zero) < 1e-9)
        #expect(result.state.mass.propellantMassKilograms == 0)
    }
}
