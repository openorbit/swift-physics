// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct FiniteThrustTests {
    private let frame = Frame(name: "Thrust inertial", pseudoInertial: true)
    private func epoch(_ seconds: Double) throws -> Epoch { try Epoch.j2000.advanced(by: Duration(seconds: seconds)) }
    private func thrust(_ force: Double = 20, direction: SIMD3<Double> = SIMD3(1,0,0)) throws -> FiniteThrust {
        try FiniteThrust(frame: frame, thrustNewtons: force, effectiveExhaustVelocityMetresPerSecond: 100, direction: direction)
    }
    private func seed(propellant: Double = 10) throws -> SpacecraftState {
        SpacecraftState(cartesian: try CartesianState(epoch: .j2000, frame: frame, position: .zero, velocity: .zero),
            mass: try SpacecraftMassProperties(dryMassKilograms: 90, propellantMassKilograms: propellant))
    }
    private func context() throws -> ForceEvaluationContext { try ForceEvaluationContext(frame: frame, origin: .earth) }
    private func options() throws -> NumericalPropagationOptions {
        try NumericalPropagationOptions(absolutePositionToleranceMetres: 1e-9,
            absoluteVelocityToleranceMetresPerSecond: 1e-10, relativeTolerance: 1e-12,
            initialStepSeconds: 10, maximumStepSeconds: 300)
    }
    private func propagator(_ burns: [FiniteBurn], initial: SpacecraftState? = nil) throws -> ScheduledThrustPropagator {
        try ScheduledThrustPropagator(initialState: initial ?? seed(),
            schedule: ScheduledThrust(frame: frame, burns: burns), context: context(), options: options())
    }
    @Test func rocketEquationAndPosition() throws {
        let engine = try thrust()
        let initial = try seed()
        let dynamics = try SpacecraftDynamics(models: [engine], context: context())
        let provider = try NumericalSpacecraftPropagator(initialState: initial, dynamics: dynamics, options: options())
        let result = try provider.spacecraftState(at: epoch(20))
        let finalMass = 96.0, exhaust = 100.0, flow = 0.2
        let velocity = exhaust * log(100/finalMass)
        let position = exhaust * (20-finalMass/flow*log(100/finalMass))
        #expect(abs(result.cartesian.velocity.x-velocity) < 1e-8)
        #expect(abs(result.cartesian.position.x-position) < 1e-7)
        #expect(abs(result.mass.propellantMassKilograms-6) < 1e-12)
        #expect(result.mass.dryMassKilograms == 90)
        #expect(result.cartesian.velocity.y == 0 && result.cartesian.position.z == 0)
        let derivative = try engine.contribution(for: result, in: context())
        #expect(abs(derivative.accelerationMetresPerSecondSquared.x-20/96.0) < 1e-14)
        #expect(derivative.propellantMassRateKilogramsPerSecond == -0.2)
    }
    @Test func shortBurnBoundariesCoastAndReplay() throws {
        let burn = try FiniteBurn(start: epoch(10), end: epoch(10.1), thrust: thrust())
        let provider = try propagator([burn])
        let ignition = try provider.spacecraftState(at: burn.start)
        #expect(ignition.cartesian.position == .zero)
        #expect(ignition.cartesian.acceleration == SIMD3(0.2,0,0))
        let cutoff = try provider.spacecraftState(at: burn.end)
        #expect(cutoff.cartesian.acceleration == .zero)
        let late = try provider.spacecraftState(at: epoch(100))
        let expectedVelocity = 100*log(100/99.98)
        #expect(abs(late.cartesian.velocity.x-expectedVelocity) < 1e-9)
        #expect(abs(late.mass.propellantMassKilograms-9.98) < 1e-12)
        #expect(abs(late.cartesian.position.x-cutoff.cartesian.position.x-expectedVelocity*89.9) < 1e-7)
        #expect(try provider.spacecraftState(at: epoch(100)) == late)
        let reverse = try propagator([burn], initial: late).spacecraftState(at: .j2000)
        #expect(abs(reverse.cartesian.position.x) < 1e-7)
        #expect(abs(reverse.cartesian.velocity.x) < 1e-9)
        #expect(abs(reverse.mass.propellantMassKilograms-10) < 1e-12)
    }
    @Test func overlappingAndAdjacentEngines() throws {
        let a = try FiniteBurn(start: epoch(1), end: epoch(3), thrust: thrust())
        let b = try FiniteBurn(start: epoch(2), end: epoch(3), thrust: thrust(10, direction: SIMD3(0,1,0)))
        let c = try FiniteBurn(start: epoch(3), end: epoch(4), thrust: thrust(10, direction: SIMD3(0,0,1)))
        let provider = try propagator([c,b,a])
        let atThree = try provider.spacecraftState(at: epoch(3))
        #expect(abs(atThree.mass.propellantMassKilograms-9.5) < 1e-12)
        #expect(atThree.cartesian.acceleration?.x == 0 && atThree.cartesian.acceleration?.y == 0)
        let end = try provider.spacecraftState(at: epoch(4))
        #expect(abs(end.mass.propellantMassKilograms-9.4) < 1e-12)
        #expect(end.cartesian.velocity.x > 0 && end.cartesian.velocity.y > 0 && end.cartesian.velocity.z > 0)
        #expect(end.cartesian.acceleration == .zero)
    }
    @Test func invalidInputsAndDepletion() throws {
        #expect(throws: FiniteThrustError.invalidThrust) { try thrust(-1) }
        #expect(throws: FiniteThrustError.invalidDirection) { try thrust(direction: .zero) }
        #expect(throws: FiniteThrustError.invalidExhaustVelocity) {
            try FiniteThrust(frame: frame, thrustNewtons: 1, effectiveExhaustVelocityMetresPerSecond: 0, direction: SIMD3(1,0,0))
        }
        #expect(throws: FiniteThrustError.invalidInterval) { try FiniteBurn(start: epoch(2), end: epoch(1), thrust: thrust()) }
        let burn = try FiniteBurn(start: epoch(0), end: epoch(100), thrust: thrust())
        #expect(throws: NumericalPropagationError.propellantDepleted) { try propagator([burn]).state(at: epoch(80)) }
        let empty = try seed(propellant: 0)
        #expect(throws: NumericalPropagationError.propellantDepleted) { try thrust().contribution(for: empty, in: context()) }
        #expect(try thrust(0).contribution(for: empty, in: context()).accelerationMetresPerSecondSquared == .zero)
        let normalized = try thrust(direction: SIMD3(1e308,1e308,0))
        #expect(abs(normalized.direction.x-1/sqrt(2)) < 1e-15)
        #expect(throws: SpacecraftDynamicsError.frameMismatch) {
            try thrust().contribution(for: seed(), in: ForceEvaluationContext(frame: Frame(name: "Other", pseudoInertial: true), origin: .earth))
        }
    }
}
