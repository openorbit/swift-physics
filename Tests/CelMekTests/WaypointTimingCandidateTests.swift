// SPDX-License-Identifier: Apache-2.0
import Math
import Testing
@testable import CelMek

@Suite struct WaypointTimingCandidateTests {
    let frame = Frame(name: "Generated timings", pseudoInertial: true)
    func epoch(_ seconds: Double) throws -> Epoch { try Epoch.j2000.advanced(by: Duration(seconds: seconds)) }
    func limits(_ acceleration: Double) throws -> WaypointRouteEngineLimits {
        try WaypointRouteEngineLimits(maximumAccelerationMetresPerSecondSquared: acceleration,
            maximumThrustNewtons: 1000, effectiveExhaustVelocityMetresPerSecond: 100)
    }
    @Test func fractionalGridIncludesOffGridPhaseBoundariesAndExactWindows() throws {
        let phase = try WaypointRouteEnginePhase(startEpoch: epoch(10.375), endEpoch: epoch(11.5), limits: limits(1))
        let windows = try [epoch(10.125)...epoch(12), epoch(50.25)...epoch(50.25)]
        let generated = try WaypointRouteTimingPlanner.arrivalCandidates(in: windows, spacing: Duration(seconds: 1),
            including: [phase, phase])
        #expect(try generated[0] == [epoch(10.125), epoch(10.375), epoch(11.125), epoch(11.5), epoch(12)])
        #expect(try generated[1] == [epoch(50.25)])
    }
    @Test func generatedCandidatesFindNonuniformFeasibleRoute() throws {
        let window = try epoch(10)...epoch(20), arrival = try epoch(50)
        let phase = try WaypointRouteEnginePhase(startEpoch: .j2000, endEpoch: epoch(15), limits: limits(0.1))
        let candidates = try WaypointRouteTimingPlanner.arrivalCandidates(in: [window, arrival...arrival],
            spacing: Duration(seconds: 10), including: [phase])
        #expect(try candidates[0] == [epoch(10), epoch(15), epoch(20)])
        let result = try WaypointRouteTimingPlanner.search(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: [BrachistochroneWaypoint(positionMetres: .zero),
                BrachistochroneWaypoint(positionMetres: SIMD3(4,0,0), arrivalWindow: window),
                BrachistochroneWaypoint(positionMetres: SIMD3(8,0,0), arrivalEpoch: arrival)],
            legAccelerationsMetresPerSecondSquared: [1,1], candidateArrivalEpochs: candidates,
            initialMass: SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50),
            backgroundDynamics: SpacecraftDynamics(models: [], context: ForceEvaluationContext(frame: frame, origin: .sun)),
            legEngineLimits: [limits(1), limits(1)], timedEnginePhases: [phase], numericalOptions: .preset(for: .full))
        #expect(result.status == .feasible)
        #expect(result.evaluatedCandidates == 2)
        let route = try #require(result.route)
        #expect(try route.waypointArrivalEpochs[1] == epoch(15))
        #expect(route.legs.map { $0.duration.seconds } == [15,35])
        let propagation = try #require(result.propagation)
        #expect(propagation.checkpoints[0].state.cartesian.velocity.x > 0)
    }
    @Test func spacingAndGenerationBudgetsAreExplicit() throws {
        let window = try epoch(0)...epoch(1_000_000)
        #expect(throws: WaypointTimingSearchError.candidateGenerationLimitExceeded) {
            try WaypointRouteTimingPlanner.arrivalCandidates(in: [window], spacing: Duration(seconds: 1),
                maximumTotalCandidates: 3)
        }
        #expect(throws: WaypointTimingSearchError.invalidCandidateSpacing) {
            try WaypointRouteTimingPlanner.arrivalCandidates(in: [window], spacing: .zero)
        }
        #expect(throws: WaypointTimingSearchError.invalidCandidates) {
            try WaypointRouteTimingPlanner.arrivalCandidates(in: [], spacing: Duration(seconds: 1))
        }
        let single = try epoch(2)...epoch(2)
        #expect(throws: WaypointTimingSearchError.candidateGenerationLimitExceeded) {
            try WaypointRouteTimingPlanner.arrivalCandidates(in: [single, single], spacing: Duration(seconds: 1),
                maximumTotalCandidates: 1)
        }
    }
}
