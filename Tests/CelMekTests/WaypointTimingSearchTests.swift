// SPDX-License-Identifier: Apache-2.0
import Math
import Testing
@testable import CelMek

@Suite struct WaypointTimingSearchTests {
    let frame = Frame(name: "Timing search", pseudoInertial: true)
    func epoch(_ seconds: Double) throws -> Epoch { try Epoch.j2000.advanced(by: Duration(seconds: seconds)) }
    func points() throws -> [BrachistochroneWaypoint] {
        try [BrachistochroneWaypoint(positionMetres: .zero),
            BrachistochroneWaypoint(positionMetres: SIMD3(4,0,0), arrivalWindow: epoch(10)...epoch(20)),
            BrachistochroneWaypoint(positionMetres: SIMD3(8,0,0), arrivalEpoch: epoch(50))]
    }
    func search(_ midpoints: [Double], budget: Int = 10) throws -> WaypointTimingSearchResult {
        let normal = try WaypointRouteEngineLimits(maximumAccelerationMetresPerSecondSquared: 1,
            maximumThrustNewtons: 1000, effectiveExhaustVelocityMetresPerSecond: 100)
        let reduced = try WaypointRouteEngineLimits(maximumAccelerationMetresPerSecondSquared: 0.1,
            maximumThrustNewtons: 1000, effectiveExhaustVelocityMetresPerSecond: 100)
        return try WaypointRouteTimingPlanner.search(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: points(), legAccelerationsMetresPerSecondSquared: [1,1],
            candidateArrivalEpochs: [midpoints.map { try epoch($0) }, [epoch(50)]],
            initialMass: SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50),
            backgroundDynamics: SpacecraftDynamics(models: [], context: ForceEvaluationContext(frame: frame, origin: .sun)),
            legEngineLimits: [normal, normal],
            timedEnginePhases: [WaypointRouteEnginePhase(startEpoch: .j2000, endEpoch: epoch(15), limits: reduced)],
            maximumCandidateEvaluations: budget, numericalOptions: .preset(for: .full))
    }
    @Test func nonuniformTimingMovesPassageAcrossTemporalBoundary() throws {
        #expect(throws: BrachistochroneRouteError.uniformRetimingUnavailable) {
            try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
                waypoints: points(), legAccelerationsMetresPerSecondSquared: [1,1])
        }
        let result = try search([20,10,10])
        #expect(result.status == .feasible)
        #expect(result.evaluatedCandidates == 2)
        #expect(result.failedAttempts.count == 1)
        #expect(result.failedAttempts[0].failure == .engine(.accelerationLimitExceeded(leg: 0)))
        let route = try #require(result.route)
        #expect(route.legs.map { $0.duration.seconds } == [20,30])
        #expect(try route.waypointArrivalEpochs[1] == epoch(20))
        let propagation = try #require(result.propagation)
        #expect(propagation.segments.count == 3)
        #expect(propagation.checkpoints[0].state.cartesian.velocity.x > 0)
        #expect(propagation.checkpoints[0].state.mass == propagation.segments[2].initialState.mass)
        #expect(propagation.state.mass.propellantMassKilograms < 50)
    }
    @Test func budgetAndExhaustionHaveDistinctMeanings() throws {
        let bounded = try search([10,20], budget: 1)
        #expect(bounded.status == .evaluationLimitReached)
        #expect(bounded.evaluatedCandidates == 1)
        #expect(bounded.route == nil)
        let exhausted = try search([10], budget: 1)
        #expect(exhausted.status == .candidateSetExhausted)
        #expect(exhausted.evaluatedCandidates == 1)
        #expect(exhausted.propagation == nil)
        #expect(throws: WaypointTimingSearchError.invalidCandidates) { try search([]) }
        #expect(throws: WaypointTimingSearchError.invalidBudget) { try search([20], budget: 0) }
    }
    @Test func explicitDurationsPreserveEpochsAndValidateConstraints() throws {
        let route = try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: points(), legAccelerationsMetresPerSecondSquared: [1,1],
            legDurations: [Duration(seconds: 20), Duration(seconds: 30)])
        #expect(try route.arrivalEpoch == epoch(50))
        #expect(throws: BrachistochroneRouteError.arrivalWindowNotSatisfied(index: 1)) {
            try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
                waypoints: points(), legAccelerationsMetresPerSecondSquared: [1,1],
                legDurations: [Duration(seconds: 25), Duration(seconds: 25)])
        }
        #expect(throws: BrachistochroneRouteError.invalidLegDuration) {
            try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
                waypoints: points(), legAccelerationsMetresPerSecondSquared: [1,1],
                legDurations: [.zero, Duration(seconds: 50)])
        }
    }
}
