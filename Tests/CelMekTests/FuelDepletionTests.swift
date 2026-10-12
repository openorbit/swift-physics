// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct FuelDepletionTests {
    private let frame = Frame(name: "Fuel inertial", pseudoInertial: true)
    private func epoch(_ t: Double) throws -> Epoch { try Epoch.j2000.advanced(by: Duration(seconds: t)) }
    private func make(burnEnd: Double = 100, background: [any ForceModel] = [], fuel: Double = 10) throws -> ScheduledThrustPropagator {
        let initial = SpacecraftState(cartesian: try CartesianState(epoch: .j2000, frame: frame, position: .zero, velocity: .zero),
            mass: try SpacecraftMassProperties(dryMassKilograms: 90, propellantMassKilograms: fuel))
        let engine = try FiniteThrust(frame: frame, thrustNewtons: 20, effectiveExhaustVelocityMetresPerSecond: 100,
            direction: SIMD3(1,0,0))
        let schedule = try ScheduledThrust(frame: frame, burns: [FiniteBurn(start: .j2000, end: epoch(burnEnd), thrust: engine)])
        return try ScheduledThrustPropagator(initialState: initial, schedule: schedule, backgroundModels: background,
            context: ForceEvaluationContext(frame: frame, origin: .earth),
            options: NumericalPropagationOptions(absolutePositionToleranceMetres: 1e-9,
                absoluteVelocityToleranceMetresPerSecond: 1e-10, relativeTolerance: 1e-12), fuelDepletionPolicy: .cutoff)
    }
    @Test func exhaustionAndCoastAgainstRocketEquation() throws {
        let provider = try make()
        let result = try provider.propagation(at: epoch(80))
        let dv = 100 * log(100.0/90)
        let position = 100 * (50-450*log(100.0/90)) + 30*dv
        #expect(try result.depletionEpoch == epoch(50))
        #expect(result.state.mass.propellantMassKilograms == 0)
        #expect(result.state.mass.dryMassKilograms == 90)
        #expect(abs(result.state.cartesian.velocity.x-dv) < 1e-8)
        #expect(abs(result.state.cartesian.position.x-position) < 1e-7)
        #expect(result.state.cartesian.acceleration == .zero)
        #expect(try provider.spacecraftState(at: epoch(80)) == result.state)
        let boundary = try provider.propagation(at: epoch(50))
        #expect(boundary.state.mass.propellantMassKilograms == 0)
        #expect(boundary.state.cartesian.acceleration == .zero)
        #expect(try boundary.depletionEpoch == epoch(50))
        let before = try provider.propagation(at: epoch(20))
        #expect(before.depletionEpoch == nil && abs(before.state.mass.propellantMassKilograms-6) < 1e-12)
        let reverse = try ScheduledThrustPropagator(initialState: before.state, schedule: provider.schedule,
            context: provider.context, options: provider.options, fuelDepletionPolicy: .cutoff)
            .spacecraftState(at: .j2000)
        #expect(abs(reverse.mass.propellantMassKilograms-10) < 1e-12)
        #expect(abs(reverse.cartesian.velocity.x) < 1e-8)
        #expect(abs(reverse.cartesian.position.x) < 1e-7)
    }
    @Test func exactScheduledCutoffAndEmptySeed() throws {
        let result = try make(burnEnd: 50).propagation(at: epoch(70))
        #expect(try result.depletionEpoch == epoch(50))
        #expect(result.state.mass.propellantMassKilograms == 0)
        let empty = try make(fuel: 0).propagation(at: epoch(70))
        #expect(empty.depletionEpoch == .j2000)
        #expect(empty.state.cartesian.position == .zero && empty.state.cartesian.velocity == .zero)
    }
    private struct SideForce: ForceModel {
        let massRate: Double
        func contribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
            try ForceContribution(epoch: state.epoch, frame: state.frame,
                accelerationMetresPerSecondSquared: SIMD3(0,10/state.mass.totalMassKilograms,0),
                propellantMassRateKilogramsPerSecond: massRate)
        }
    }
    @Test func backgroundSeesEvolvingMassAndContinuesAfterCutoff() throws {
        let result = try make(background: [SideForce(massRate: 0)]).spacecraftState(at: epoch(80))
        let expected = 50*log(100.0/90) + 30*10/90.0
        #expect(abs(result.cartesian.velocity.y-expected) < 1e-8)
        #expect(abs((result.cartesian.acceleration?.y ?? 0)-10/90.0) < 1e-14)
        #expect(throws: ScheduledThrustPropagationError.backgroundMassFlowUnsupported) {
            try make(background: [SideForce(massRate: -1)]).spacecraftState(at: epoch(1))
        }
    }
    @Test func overlappingEnginesAndLaterIgnitionsStayOff() throws {
        let base = try make()
        let engine = base.schedule.burns[0].thrust
        let other = try FiniteThrust(frame: frame, thrustNewtons: 20, effectiveExhaustVelocityMetresPerSecond: 100,
            direction: SIMD3(0,1,0))
        let schedule = try ScheduledThrust(frame: frame, burns: [
            FiniteBurn(start: .j2000, end: epoch(40), thrust: engine),
            FiniteBurn(start: .j2000, end: epoch(40), thrust: other),
            FiniteBurn(start: epoch(60), end: epoch(80), thrust: engine)])
        let provider = try ScheduledThrustPropagator(initialState: base.initialState, schedule: schedule,
            context: base.context, options: base.options, fuelDepletionPolicy: .cutoff)
        let result = try provider.propagation(at: epoch(90))
        #expect(try result.depletionEpoch == epoch(25))
        #expect(abs(result.state.cartesian.velocity.x-50*log(100.0/90)) < 1e-8)
        #expect(abs(result.state.cartesian.velocity.y-result.state.cartesian.velocity.x) < 1e-12)
        #expect(result.state.mass.propellantMassKilograms == 0)
        let dryReplay = try ScheduledThrustPropagator(initialState: result.state, schedule: schedule,
            context: base.context, fuelDepletionPolicy: .cutoff)
        #expect(throws: NumericalPropagationError.propellantDepleted) { try dryReplay.spacecraftState(at: .j2000) }
    }
}
