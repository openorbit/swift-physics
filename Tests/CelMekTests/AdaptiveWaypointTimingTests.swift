// SPDX-License-Identifier: Apache-2.0
import Math
import Testing
@testable import CelMek

@Suite struct AdaptiveWaypointTimingTests {
    let frame = Frame(name: "Adaptive timing", pseudoInertial: true)
    func epoch(_ seconds: Double) throws -> Epoch { try Epoch.j2000.advanced(by: Duration(seconds: seconds)) }
    func search(minimumSpacing: Double = 15, budget: Int = 10, rounds: Int = 8, generated: Int = 100,
                window: ClosedRange<Epoch>? = nil) throws -> WaypointTimingRefinementResult {
        let passage = try window ?? epoch(10)...epoch(40), arrival = try epoch(50)
        let limits = try WaypointRouteEngineLimits(maximumAccelerationMetresPerSecondSquared: 0.1,
            maximumThrustNewtons: 1000, effectiveExhaustVelocityMetresPerSecond: 100)
        return try WaypointRouteTimingPlanner.adaptiveSearch(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: [BrachistochroneWaypoint(positionMetres: .zero),
                BrachistochroneWaypoint(positionMetres: SIMD3(4,0,0), arrivalWindow: passage),
                BrachistochroneWaypoint(positionMetres: SIMD3(8,0,0), arrivalEpoch: arrival)],
            legAccelerationsMetresPerSecondSquared: [1,1], arrivalWindows: [passage, arrival...arrival],
            refinementOptions: WaypointTimingRefinementOptions(initialSpacingSeconds: 30,
                minimumSpacingSeconds: minimumSpacing, maximumRounds: rounds,
                maximumCandidateEvaluations: budget, maximumGeneratedCandidates: generated),
            initialMass: SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50),
            backgroundDynamics: SpacecraftDynamics(models: [], context: ForceEvaluationContext(frame: frame, origin: .sun)),
            legEngineLimits: [limits, limits], numericalOptions: .preset(for: .full))
    }
    @Test func refinementFindsInteriorFeasibilityWithoutRepeatingFailures() throws {
        let result = try search()
        #expect(result.status == .feasible)
        #expect(result.evaluatedCandidates == 3)
        #expect(result.rounds.map(\.spacingSeconds) == [30,15])
        #expect(result.rounds.map(\.evaluatedCandidates) == [2,1])
        #expect(result.failedAttempts.count == 2)
        #expect(result.failedAttempts[0].failure == .engine(.accelerationLimitExceeded(leg: 0)))
        #expect(result.failedAttempts[1].failure == .engine(.accelerationLimitExceeded(leg: 1)))
        let route = try #require(result.route)
        #expect(try route.waypointArrivalEpochs[1] == epoch(25))
        let propagation = try #require(result.propagation)
        #expect(propagation.checkpoints[0].state.cartesian.velocity.x > 0)
        #expect(propagation.checkpoints[0].state.mass == propagation.segments[1].initialState.mass)
        #expect(try search().propagation?.state == propagation.state)
    }
    @Test func terminationReasonsRemainDistinct() throws {
        #expect(try search(budget: 2).status == .evaluationLimitReached)
        #expect(try search(rounds: 1).status == .refinementLimitReached)
        #expect(try search(generated: 3).status == .candidateGenerationLimitReached)
        let coarse = try search(minimumSpacing: 30)
        #expect(coarse.status == .resolutionLimitReached)
        #expect(coarse.route == nil)
        #expect(coarse.evaluatedCandidates == 2)
        // An interior feasible route exists even though the requested resolution
        // is exhausted: this result must not be called global infeasibility.
        #expect(try search().status == .feasible)
    }
    @Test func exactFailedWindowsTerminateWithoutUnnecessaryRefinement() throws {
        let instant = try epoch(10)
        let result = try search(minimumSpacing: 0.01, window: instant...instant)
        #expect(result.status == .resolutionLimitReached)
        #expect(result.rounds.count == 1)
        #expect(result.evaluatedCandidates == 1)
        #expect(throws: WaypointTimingSearchError.invalidCandidateSpacing) {
            try WaypointTimingRefinementOptions(initialSpacingSeconds: 1, minimumSpacingSeconds: 2)
        }
        #expect(throws: WaypointTimingSearchError.invalidBudget) {
            try WaypointTimingRefinementOptions(initialSpacingSeconds: 1, minimumSpacingSeconds: 0.1, maximumRounds: 33)
        }
    }
}
