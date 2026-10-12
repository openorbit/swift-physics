// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct VariableFuelAndRearmingTests {
    private let frame = Frame(name: "Variable fuel axes", pseudoInertial: true)
    private func epoch(_ t: Double) throws -> Epoch { try Epoch.j2000.advanced(by: Duration(seconds: t)) }
    private func seed(fuel: Double = 10) throws -> SpacecraftState {
        SpacecraftState(cartesian: try CartesianState(epoch: .j2000, frame: frame, position: .zero, velocity: .zero),
            mass: try SpacecraftMassProperties(dryMassKilograms: 90, propellantMassKilograms: fuel))
    }
    private func options() throws -> NumericalPropagationOptions {
        try NumericalPropagationOptions(absolutePositionToleranceMetres: 1e-9,
            absoluteVelocityToleranceMetresPerSecond: 1e-10, absolutePropellantToleranceKilograms: 1e-11,
            relativeTolerance: 1e-12, initialStepSeconds: 0.1, minimumStepSeconds: 1e-9, maximumStepSeconds: 1)
    }
    @Test func variableFlowAnalyticDepletionAndCoast() throws {
        let engine = try VariableThrust(frame: frame) { state, _ in
            let t = try state.epoch.duration(since: .j2000).seconds
            return try FiniteThrust(frame: frame, thrustNewtons: 10*(1+t),
                effectiveExhaustVelocityMetresPerSecond: 100, direction: SIMD3(1,0,0))
        }
        let dynamics = try SpacecraftDynamics(models: [engine], context: ForceEvaluationContext(frame: frame, origin: .earth))
        let provider = try FuelLimitedNumericalPropagator(initialState: seed(), dynamics: dynamics,
            numericalOptions: options(), eventOptions: PropagationEventOptions(timeToleranceSeconds: 1e-8))
        let result = try provider.propagation(at: epoch(20))
        let exhaustion = sqrt(201.0)-1
        let measured = try #require(result.depletionEpoch).duration(since: .j2000).seconds
        #expect(abs(measured-exhaustion) < 2e-8)
        #expect(try #require(result.depletionBracketSeconds) <= 1e-8)
        #expect(result.state.mass.propellantMassKilograms == 0)
        #expect(abs(result.state.cartesian.velocity.x-100*log(100.0/90)) < 1e-7)
        // Independent Simpson quadrature of analytic v(t)=ve*ln(m0/m(t)).
        let n = 2000, h = exhaustion/Double(n)
        var sum = 0.0
        for i in 0...n {
            let t = Double(i)*h, mass = 100-0.1*(t+t*t/2)
            sum += (i == 0 || i == n ? 1 : i % 2 == 0 ? 2 : 4)*100*log(100/mass)
        }
        let expectedPosition = sum*h/3 + (20-exhaustion)*100*log(100.0/90)
        #expect(abs(result.state.cartesian.position.x-expectedPosition) < 2e-6)
        #expect(result.state.cartesian.acceleration == .zero)
        #expect(try provider.propagation(at: epoch(20)).state == result.state)
        let before = try provider.propagation(at: epoch(5))
        #expect(before.depletionEpoch == nil)
        #expect(abs(before.state.mass.propellantMassKilograms-8.25) < 1e-9)
        let limited = try FuelLimitedNumericalPropagator(initialState: seed(), dynamics: dynamics,
            numericalOptions: options(), eventOptions: PropagationEventOptions(maximumStateEvaluations: 1))
        #expect(throws: PropagationEventError.evaluationLimitExceeded) { try limited.propagation(at: epoch(20)) }
    }
    private struct Recurring: PropagationEventDetector {
        let identifier = "periodic"
        let direction: EventCrossingDirection = .any
        let action: PropagationEventAction = .record
        func value(at state: SpacecraftState) throws -> Double {
            let t = try state.epoch.duration(since: .j2000).seconds
            return sin(Double.pi*(t-0.25))
        }
    }
    @Test func recurringRootsRearmAfterLeavingBand() throws {
        let dynamics = try SpacecraftDynamics(models: [], context: ForceEvaluationContext(frame: frame, origin: .earth))
        let provider = try EventDrivenNumericalPropagator(initialState: seed(), dynamics: dynamics, detectors: [Recurring()],
            eventOptions: PropagationEventOptions(maximumCheckIntervalSeconds: 0.1, timeToleranceSeconds: 1e-9),
            rearmingThresholds: ["periodic": 0.05])
        let result = try provider.propagation(at: epoch(3.5))
        #expect(result.events.count == 4)
        for (i,event) in result.events.enumerated() {
            #expect(try abs(event.state.epoch.duration(since: epoch(Double(i)+0.25)).seconds) < 1e-8)
        }
        let reverse = try EventDrivenNumericalPropagator(initialState: result.state, dynamics: dynamics, detectors: [Recurring()],
            eventOptions: PropagationEventOptions(maximumCheckIntervalSeconds: 0.1), rearmingThresholds: ["periodic": 0.05])
            .propagation(at: .j2000)
        #expect(reverse.events.count == 4)
        let once = try EventDrivenNumericalPropagator(initialState: seed(), dynamics: dynamics, detectors: [Recurring()])
            .propagation(at: epoch(0.5))
        #expect(once.events.count == 1)
        let bounded = try EventDrivenNumericalPropagator(initialState: seed(), dynamics: dynamics, detectors: [Recurring()],
            eventOptions: PropagationEventOptions(maximumCheckIntervalSeconds: 0.1, maximumEvents: 2),
            rearmingThresholds: ["periodic": 0.05])
        #expect(throws: PropagationEventError.eventLimitExceeded) { try bounded.propagation(at: epoch(3.5)) }
    }
    private struct UnclassifiedConsumption: ForceModel {
        func contribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
            try ForceContribution(epoch: state.epoch, frame: state.frame, accelerationMetresPerSecondSquared: .zero,
                propellantMassRateKilogramsPerSecond: -1)
        }
    }
    @Test func unsupportedFlowAndInvalidRearming() throws {
        let dynamics = try SpacecraftDynamics(models: [UnclassifiedConsumption()], context: ForceEvaluationContext(frame: frame, origin: .earth))
        #expect(throws: FuelBoundaryError.unsupportedMassFlow) {
            try FuelLimitedNumericalPropagator(initialState: seed(), dynamics: dynamics).propagation(at: epoch(20))
        }
        #expect(throws: PropagationEventError.invalidOptions) {
            try EventDrivenNumericalPropagator(initialState: seed(), dynamics: dynamics, detectors: [Recurring()],
                rearmingThresholds: ["periodic": 0])
        }
        let empty = try FuelLimitedNumericalPropagator(initialState: seed(fuel: 0),
            dynamics: SpacecraftDynamics(models: [], context: ForceEvaluationContext(frame: frame, origin: .earth)))
            .propagation(at: epoch(10))
        #expect(empty.state.mass.propellantMassKilograms == 0 && empty.state.cartesian.velocity == .zero)
    }
}
