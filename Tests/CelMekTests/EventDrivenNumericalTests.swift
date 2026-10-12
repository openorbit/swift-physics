// SPDX-License-Identifier: Apache-2.0
import Math
import Testing
@testable import CelMek

@Suite struct EventDrivenNumericalTests {
    private let frame = Frame(name: "Event reset axes", pseudoInertial: true)
    private func epoch(_ t: Double) throws -> Epoch { try Epoch.j2000.advanced(by: Duration(seconds: t)) }
    private func seed() throws -> SpacecraftState {
        SpacecraftState(cartesian: try CartesianState(epoch: .j2000, frame: frame, position: .zero, velocity: SIMD3(1,0,0)),
            mass: try SpacecraftMassProperties(dryMassKilograms: 1))
    }
    private func dynamics(_ acceleration: Double = 0) throws -> SpacecraftDynamics {
        try SpacecraftDynamics(models: [ConstantAcceleration(value: acceleration)],
            context: ForceEvaluationContext(frame: frame, origin: .earth))
    }
    private struct ConstantAcceleration: ForceModel {
        let value: Double
        func contribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
            try ForceContribution(epoch: state.epoch, frame: state.frame, accelerationMetresPerSecondSquared: SIMD3(value,0,0))
        }
    }
    private struct Plane: PropagationEventDetector {
        let identifier = "plane"
        let action: PropagationEventAction = .stop
        let direction: EventCrossingDirection = .increasing
        func value(at state: SpacecraftState) throws -> Double { state.cartesian.position.x-2.25 }
    }
    @Test func acceptedStepTerminalCrossing() throws {
        let result = try EventDrivenNumericalPropagator(initialState: seed(), dynamics: dynamics(), detectors: [Plane()],
            eventOptions: PropagationEventOptions(maximumCheckIntervalSeconds: 1, timeToleranceSeconds: 1e-9)).propagation(at: epoch(10))
        #expect(!result.reachedTarget && result.events.count == 1)
        #expect(abs(result.state.cartesian.position.x-2.25) < 1e-8)
        #expect(try abs(result.state.epoch.duration(since: epoch(2.25)).seconds) < 1e-8)
    }
    @Test func impulseAndInverseReplay() throws {
        let event = try EpochEvent(identifier: "impulse", epoch: epoch(2), action: .record)
        let reset: (PropagationEventOccurrence,Bool,SpacecraftDynamics) throws -> EventStateTransition = { occurrence, forward, force in
            let s = occurrence.state
            return EventStateTransition(state: SpacecraftState(cartesian: try CartesianState(epoch: s.epoch, frame: s.frame,
                position: s.cartesian.position, velocity: s.cartesian.velocity+SIMD3(forward ? 2 : -2,0,0)), mass: s.mass), dynamics: force)
        }
        let provider = try EventDrivenNumericalPropagator(initialState: seed(), dynamics: dynamics(), detectors: [event], transition: reset)
        let result = try provider.propagation(at: epoch(5))
        #expect(result.reachedTarget && result.events.count == 1)
        #expect(abs(result.state.cartesian.position.x-11) < 1e-9)
        #expect(abs(result.state.cartesian.velocity.x-3) < 1e-9)
        #expect(try provider.propagation(at: epoch(5)).state == result.state)
        let reverse = try EventDrivenNumericalPropagator(initialState: result.state, dynamics: dynamics(), detectors: [event], transition: reset)
            .propagation(at: .j2000)
        #expect(abs(reverse.state.cartesian.position.x) < 1e-8)
        #expect(abs(reverse.state.cartesian.velocity.x-1) < 1e-9)
    }
    @Test func forceChangesAtTwoEvents() throws {
        let on = try EpochEvent(identifier: "on", epoch: epoch(2), action: .record)
        let off = try EpochEvent(identifier: "off", epoch: epoch(4), action: .record)
        let provider = try EventDrivenNumericalPropagator(initialState: seed(), dynamics: dynamics(), detectors: [off,on]) { event, _, _ in
            EventStateTransition(state: event.state, dynamics: try dynamics(event.identifier == "on" ? 2 : 0))
        }
        let result = try provider.propagation(at: epoch(6))
        #expect(result.events.map(\.identifier) == ["on","off"])
        #expect(abs(result.state.cartesian.position.x-18) < 1e-8)
        #expect(abs(result.state.cartesian.velocity.x-5) < 1e-9)
        #expect(result.state.cartesian.acceleration == .zero)
    }
    @Test func resetValidationAndOneShotInitialRoot() throws {
        let event = try EpochEvent(identifier: "initial", epoch: .j2000, action: .record)
        let result = try EventDrivenNumericalPropagator(initialState: seed(), dynamics: dynamics(), detectors: [event]).propagation(at: epoch(3))
        #expect(result.events.count == 1 && result.reachedTarget)
        let bad = try EventDrivenNumericalPropagator(initialState: seed(), dynamics: dynamics(), detectors: [event]) { occurrence, _, force in
            let wrong = try CartesianState(epoch: epoch(1), frame: frame, position: .zero, velocity: .zero)
            return EventStateTransition(state: SpacecraftState(cartesian: wrong, mass: occurrence.state.mass), dynamics: force)
        }
        #expect(throws: SpacecraftDynamicsError.epochMismatch) { try bad.propagation(at: epoch(3)) }
    }
}
