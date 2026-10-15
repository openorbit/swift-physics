// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct MovingWaypointPassageTests {
    let frame = Frame(name: "Moving passage", pseudoInertial: true)
    func epoch(_ seconds: Double) throws -> Epoch { try Epoch.j2000.advanced(by: Duration(seconds: seconds)) }
    func passage(speed: Double = 0.1, matching: Bool = false, origin: BodyIdentifier = .sun) throws -> MovingWaypointPassage {
        MovingWaypointPassage(provider: try MovingWaypointRouteTarget(provider: LinearPoint(frame: frame, speed: speed),
            origin: origin, coverage: epoch(8)...epoch(12)), velocityMode: matching ? .matchProviderVelocity : .preserveRouteVelocity)
    }
    func points(velocity: SIMD3<Double>? = nil) throws -> [BrachistochroneWaypoint] {
        try [BrachistochroneWaypoint(positionMetres: .zero),
            BrachistochroneWaypoint(positionMetres: SIMD3(-999,0,0), velocityMetresPerSecond: velocity),
            BrachistochroneWaypoint(positionMetres: SIMD3(6,0,0))]
    }
    func limits() throws -> WaypointRouteEngineLimits {
        try WaypointRouteEngineLimits(maximumAccelerationMetresPerSecondSquared: 1,
            maximumThrustNewtons: 1000, effectiveExhaustVelocityMetresPerSecond: 100)
    }
    func background() throws -> SpacecraftDynamics {
        try SpacecraftDynamics(models: [], context: ForceEvaluationContext(frame: frame, origin: .sun))
    }
    func search(_ passage: MovingWaypointPassage, index: Int = 1, velocity: SIMD3<Double>? = nil) throws -> WaypointTimingSearchResult {
        try WaypointRouteTimingPlanner.search(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: points(velocity: velocity), legAccelerationsMetresPerSecondSquared: [1,1],
            candidateArrivalEpochs: [[epoch(5),epoch(10)],[epoch(20)]],
            initialMass: SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50),
            backgroundDynamics: background(), legEngineLimits: [limits(),limits()], numericalOptions: .preset(for: .full),
            movingWaypoints: [index: passage])
    }
    @Test func movingAndStationaryPositionsRemainFlyThroughByDefault() throws {
        let found = try search(passage())
        #expect(found.status == .feasible)
        #expect(found.failedAttempts[0].failure == .movingWaypointOutsideCoverage(index: 1))
        let route = try #require(found.route)
        #expect(route.waypoints[1].positionMetres == SIMD3(2,0,0))
        #expect(abs(route.waypointVelocitiesMetresPerSecond[1].x-0.3) < 1e-12)
        let result = try #require(found.propagation)
        #expect(abs(result.checkpoints[0].state.cartesian.position.x-2) < 1e-10)
        #expect(result.checkpoints[0].state.mass == result.segments[1].initialState.mass)
        let stationary = try search(passage(speed: 0))
        #expect(stationary.status == .feasible)
        #expect(try #require(stationary.route).waypointVelocitiesMetresPerSecond[1].x > 0)
        let matched = try search(passage(matching: true))
        #expect(try #require(matched.route).waypointVelocitiesMetresPerSecond[1] == SIMD3(0.1,0,0))
        let stopped = try search(passage(speed: 0, matching: true))
        #expect(stopped.status == .candidateSetExhausted)
        #expect(stopped.failedAttempts.last?.failure == .geometry(.degeneratePassage(index: 1)))
    }
    @Test func adaptiveAndOptimizedPassagesRespectProviderCoverage() throws {
        let moving = try passage()
        let mass = try SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50)
        let found = try WaypointRouteTimingPlanner.adaptiveSearch(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: points(), legAccelerationsMetresPerSecondSquared: [1,1],
            arrivalWindows: [epoch(0)...epoch(20), epoch(24)...epoch(24)],
            refinementOptions: WaypointTimingRefinementOptions(initialSpacingSeconds: 10, minimumSpacingSeconds: 1),
            initialMass: mass, backgroundDynamics: background(), legEngineLimits: [limits(),limits()], movingWaypoints: [1: moving])
        #expect(found.status == .feasible)
        let seed = try #require(found.route)
        #expect(try seed.waypointArrivalEpochs[1] == epoch(8))
        let optimized = try WaypointRouteTimingPlanner.optimizeTiming(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: points(), legAccelerationsMetresPerSecondSquared: [1,1],
            arrivalWindows: [epoch(0)...epoch(20), epoch(18)...epoch(24)],
            seedArrivalEpochs: Array(seed.waypointArrivalEpochs.dropFirst()),
            optimizationOptions: WaypointTimingOptimizationOptions(initialStepSeconds: 2, minimumStepSeconds: 0.01),
            initialMass: mass, backgroundDynamics: background(), legEngineLimits: [limits(),limits()], movingWaypoints: [1: moving])
        let route = try #require(optimized.route)
        #expect(moving.provider.coverage.contains(route.waypointArrivalEpochs[1]))
        #expect(try route.waypoints[1].positionMetres == moving.provider.state(at: route.waypointArrivalEpochs[1]).position)
        #expect(route.arrivalEpoch < seed.arrivalEpoch)
        #expect(route.waypointVelocitiesMetresPerSecond[1].x > 0)
    }
    @Test func rejectsIndexOriginAndVelocityConflicts() throws {
        #expect(throws: MovingWaypointTargetError.invalidWaypointIndex) { try search(passage(), index: 0) }
        #expect(throws: MovingWaypointTargetError.invalidWaypointIndex) { try search(passage(), index: 2) }
        #expect(throws: MovingWaypointTargetError.originMismatch) { try search(passage(origin: .earth)) }
        #expect(throws: MovingWaypointTargetError.conflictingWaypointVelocity(index: 1)) {
            try search(passage(matching: true), velocity: SIMD3(1,0,0))
        }
        let explicit = try search(passage(), velocity: SIMD3(0.7,0,0))
        #expect(try #require(explicit.route).waypointVelocitiesMetresPerSecond[1] == SIMD3(0.7,0,0))
    }
    private struct LinearPoint: StateProvider {
        let frame: Frame
        let speed: Double
        func state(at epoch: Epoch) throws -> CartesianState {
            try CartesianState(epoch: epoch, frame: frame,
                position: SIMD3(1+speed*epoch.duration(since: .j2000).seconds,0,0), velocity: SIMD3(speed,0,0))
        }
    }
}
