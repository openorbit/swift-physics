// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct WaypointTimingOptimizationTests {
    let frame = Frame(name: "Optimized timings", pseudoInertial: true)
    func epoch(_ seconds: Double) throws -> Epoch { try Epoch.j2000.advanced(by: Duration(seconds: seconds)) }
    func optimize(seed: Double = 10, budget: Int = 1000, iterations: Int = 256) throws -> WaypointTimingOptimizationResult {
        let limits = try WaypointRouteEngineLimits(maximumAccelerationMetresPerSecondSquared: 1,
            maximumThrustNewtons: 1000, effectiveExhaustVelocityMetresPerSecond: 100)
        return try WaypointRouteTimingPlanner.optimizeTiming(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: [BrachistochroneWaypoint(positionMetres: .zero), BrachistochroneWaypoint(positionMetres: SIMD3(4,0,0))],
            legAccelerationsMetresPerSecondSquared: [1], arrivalWindows: [epoch(4)...epoch(12)], seedArrivalEpochs: [epoch(seed)],
            optimizationOptions: WaypointTimingOptimizationOptions(initialStepSeconds: 2, minimumStepSeconds: 0.001,
                maximumIterations: iterations, maximumCandidateEvaluations: budget),
            initialMass: SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50),
            backgroundDynamics: SpacecraftDynamics(models: [], context: ForceEvaluationContext(frame: frame, origin: .sun)),
            legEngineLimits: [limits], numericalOptions: .preset(for: .full))
    }
    @Test func convergesToIndependentSingleCubicAccelerationThreshold() throws {
        let result = try optimize()
        #expect(result.status == .stepLimitReached)
        let route = try #require(result.route)
        // For rest-to-rest cubic motion max |a| = 6*d/T², so T >= sqrt(24).
        // This is a threshold for this cubic family, not the ideal bang-bang optimum.
        let threshold = sqrt(24.0)
        #expect(route.duration.seconds >= threshold)
        #expect(route.duration.seconds-threshold < 0.002)
        #expect(result.iterations.contains { $0.improved })
        let epochs = result.iterations.map(\.bestArrivalEpoch)
        #expect(zip(epochs, epochs.dropFirst()).allSatisfy { $0.0 >= $0.1 })
        #expect(try optimize().propagation?.state == result.propagation?.state)
    }
    @Test func budgetsPreserveFeasibleIncumbentAndBadSeedIsExplicit() throws {
        let limited = try optimize(budget: 1)
        #expect(limited.status == .evaluationLimitReached)
        #expect(limited.evaluatedCandidates == 1)
        #expect(limited.route?.duration.seconds == 10)
        #expect(limited.propagation != nil)
        let iteration = try optimize(iterations: 1)
        #expect(iteration.status == .iterationLimitReached)
        #expect(iteration.route?.duration.seconds == 8)
        let rejected = try optimize(seed: 4)
        #expect(rejected.status == .seedRejected)
        #expect(rejected.route == nil)
        #expect(rejected.evaluatedCandidates == 1)
        #expect(rejected.failedAttempts[0].failure == .geometry(.timingAccelerationLimitExceeded(leg: 0)))
    }
    @Test func preservesExactPassageAndTimedEngineConstraints() throws {
        let passage = try epoch(20)
        let limits = try WaypointRouteEngineLimits(maximumAccelerationMetresPerSecondSquared: 0.1,
            maximumThrustNewtons: 1000, effectiveExhaustVelocityMetresPerSecond: 100)
        let phase = try WaypointRouteEnginePhase(startEpoch: .j2000, endEpoch: epoch(15), limits: limits)
        let result = try WaypointRouteTimingPlanner.optimizeTiming(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: [BrachistochroneWaypoint(positionMetres: .zero),
                BrachistochroneWaypoint(positionMetres: SIMD3(4,0,0), arrivalEpoch: passage),
                BrachistochroneWaypoint(positionMetres: SIMD3(8,0,0))],
            legAccelerationsMetresPerSecondSquared: [1,1],
            arrivalWindows: [epoch(10)...epoch(30), epoch(30)...epoch(50)], seedArrivalEpochs: [passage, epoch(50)],
            optimizationOptions: WaypointTimingOptimizationOptions(initialStepSeconds: 2, minimumStepSeconds: 0.01),
            initialMass: SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50),
            backgroundDynamics: SpacecraftDynamics(models: [], context: ForceEvaluationContext(frame: frame, origin: .sun)),
            legEngineLimits: [limits, limits], timedEnginePhases: [phase], numericalOptions: .preset(for: .full))
        #expect(result.status == .stepLimitReached)
        let route = try #require(result.route)
        #expect(route.waypointArrivalEpochs[1] == passage)
        #expect(try route.arrivalEpoch < epoch(50))
        let execution = try #require(result.propagation)
        #expect(execution.checkpoints[0].state.cartesian.velocity.x > 0)
        #expect(execution.checkpoints[0].state.mass == execution.segments[2].initialState.mass)
        #expect(execution.checkpoints.allSatisfy { $0.positionErrorMetres < 1e-8 })
    }
}
