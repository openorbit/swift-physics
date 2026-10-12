// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct PropagationEventTests {
    private let frame = Frame(name: "Event inertial", pseudoInertial: true)
    private func epoch(_ t: Double) throws -> Epoch { try Epoch.j2000.advanced(by: Duration(seconds: t)) }
    private func trajectory() throws -> NumericalSpacecraftPropagator {
        let initial = SpacecraftState(cartesian: try CartesianState(epoch: .j2000, frame: frame,
            position: SIMD3(10,0,0), velocity: SIMD3(-1,0,0)), mass: try SpacecraftMassProperties(dryMassKilograms: 1))
        return try NumericalSpacecraftPropagator(initialState: initial,
            dynamics: SpacecraftDynamics(models: [], context: ForceEvaluationContext(frame: frame, origin: .earth)))
    }
    @Test func impactLocalizationAndBackwardDirection() throws {
        let impact = try SphericalImpactEvent(identifier: "surface", frame: frame, radiusMetres: 2)
        let search = try EventDetectedPropagation(trajectory: trajectory(), detectors: [impact],
            options: PropagationEventOptions(maximumCheckIntervalSeconds: 3, timeToleranceSeconds: 1e-8))
        let result = try search.propagation(from: .j2000, to: epoch(9))
        #expect(!result.reachedTarget)
        #expect(result.events.count == 1)
        #expect(try abs(result.state.epoch.duration(since: epoch(8)).seconds) <= 1e-8)
        #expect(abs(result.state.cartesian.position.x-2) <= 1e-8)
        #expect(result.events[0].bracketDurationSeconds <= 1e-8)
        let reverse = try search.propagation(from: epoch(9), to: .j2000)
        #expect(try abs(reverse.state.epoch.duration(since: epoch(8)).seconds) <= 1e-8)
    }
    @Test func orderedRecordsTerminalAndInitialRoots() throws {
        let a = try EpochEvent(identifier: "first", epoch: epoch(2), action: .record)
        let b = try EpochEvent(identifier: "stop", epoch: epoch(5))
        let c = try EpochEvent(identifier: "later", epoch: epoch(7), action: .record)
        let search = try EventDetectedPropagation(trajectory: trajectory(), detectors: [c,b,a],
            options: PropagationEventOptions(maximumCheckIntervalSeconds: 10, timeToleranceSeconds: 1e-8))
        let result = try search.propagation(from: .j2000, to: epoch(9))
        #expect(result.events.map(\.identifier) == ["first","stop"])
        #expect(try result.state.epoch == epoch(5))
        #expect(result.events[1].bracketDurationSeconds == 0)
        let initial = try search.propagation(from: epoch(5), to: epoch(9))
        #expect(try initial.state.epoch == epoch(5))
        #expect(initial.events.count == 1)
        let ignored = try search.propagation(from: epoch(5), to: epoch(9), includeInitialEvents: false)
        #expect(ignored.reachedTarget && ignored.events.map(\.identifier) == ["later"])
        let records = try EventDetectedPropagation(trajectory: trajectory(), detectors: [a,c])
        #expect(try records.propagation(from: .j2000, to: epoch(9)).events.map(\.identifier) == ["first","later"])
    }
    private struct BadDetector: PropagationEventDetector {
        let identifier = "bad"
        let direction: EventCrossingDirection = .any
        let action: PropagationEventAction = .record
        func value(at state: SpacecraftState) throws -> Double { .nan }
    }
    @Test func limitsAndInvalidContracts() throws {
        #expect(throws: PropagationEventError.invalidOptions) { try PropagationEventOptions(timeToleranceSeconds: 0) }
        let event = try EpochEvent(identifier: "time", epoch: epoch(2))
        #expect(throws: PropagationEventError.invalidDetector) { try EventDetectedPropagation(trajectory: trajectory(), detectors: [event,event]) }
        let budget = try EventDetectedPropagation(trajectory: trajectory(), detectors: [event],
            options: PropagationEventOptions(maximumStateEvaluations: 1))
        #expect(throws: PropagationEventError.evaluationLimitExceeded) { try budget.propagation(from: .j2000, to: epoch(3)) }
        let bad = try EventDetectedPropagation(trajectory: trajectory(), detectors: [BadDetector()])
        #expect(throws: PropagationEventError.nonFiniteValue) { try bad.propagation(from: .j2000, to: epoch(3)) }
        let wrong = try SphericalImpactEvent(identifier: "wrong", frame: Frame(name: "Other", pseudoInertial: true), radiusMetres: 2)
        #expect(throws: SpacecraftDynamicsError.frameMismatch) {
            try EventDetectedPropagation(trajectory: trajectory(), detectors: [wrong]).propagation(from: .j2000, to: epoch(3))
        }
    }
    @Test func scheduledBurnSearch() throws {
        let base = try trajectory()
        let engine = try FiniteThrust(frame: frame, thrustNewtons: 1, effectiveExhaustVelocityMetresPerSecond: 100, direction: SIMD3(1,0,0))
        let initial = SpacecraftState(cartesian: base.initialState.cartesian,
            mass: try SpacecraftMassProperties(dryMassKilograms: 1, propellantMassKilograms: 1))
        let burn = try FiniteBurn(start: epoch(1), end: epoch(2), thrust: engine)
        let scheduled = try ScheduledThrustPropagator(initialState: initial,
            schedule: ScheduledThrust(frame: frame, burns: [burn]), context: base.dynamics.context)
        let end = try EpochEvent(identifier: "cutoff", epoch: burn.end)
        let result = try EventDetectedPropagation(trajectory: scheduled, detectors: [end])
            .propagation(from: .j2000, to: epoch(4))
        #expect(result.events.count == 1)
        #expect(abs(result.state.mass.propellantMassKilograms-0.99) < 1e-8)
        #expect(result.state.epoch == burn.end)
        #expect(result.state.cartesian.acceleration == .zero)
    }
}
