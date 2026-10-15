// SPDX-License-Identifier: Apache-2.0
import Math
import Testing
@testable import CelMek

@Suite struct MovingObstacleTests {
    let frame = Frame(name: "Moving obstacle", pseudoInertial: true)
    func epoch(_ t: Double) throws -> Epoch { try Epoch.j2000.advanced(by: Duration(seconds: t)) }
    struct Centre: StateProvider {
        let frame: Frame
        let offset: Double
        func state(at epoch: Epoch) throws -> CartesianState {
            let t = try epoch.duration(since: .j2000).seconds
            return try CartesianState(epoch: epoch, frame: frame,
                position: SIMD3(5, t-5+offset, 0), velocity: SIMD3(0,1,0))
        }
    }
    func obstacle(offset: Double = 0, coverageEnd: Double = 30, bound: Double = 1) throws -> MovingSphericalExclusionRegion {
        try MovingSphericalExclusionRegion(identifier: "crossing", centreProvider:
            MovingWaypointRouteTarget(provider: Centre(frame: frame, offset: offset), origin: .sun,
                coverage: epoch(0)...epoch(coverageEnd)), radiusMetres: 0.2,
            maximumCentreSpeedMetresPerSecond: bound)
    }
    func points() throws -> [BrachistochroneWaypoint] {
        try [BrachistochroneWaypoint(positionMetres: .zero),
             BrachistochroneWaypoint(positionMetres: SIMD3(10,0,0))]
    }
    func route(_ obstacle: MovingSphericalExclusionRegion, time: Double = 10) throws -> WaypointBrachistochroneRoute {
        try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: points(), legAccelerationsMetresPerSecondSquared: [1],
            legDurations: [Duration(seconds: time)], movingObstacles: [obstacle])
    }
    @Test func crossingCollisionAndConservativeClearance() throws {
        #expect(throws: BrachistochroneRouteError.obstructedLeg(index: 0, obstacle: "crossing")) {
            try route(obstacle())
        }
        let safe = try route(obstacle(offset: 10))
        let report = try #require(safe.clearances.first)
        #expect(report.minimumClearanceMetres > 0)
        for i in 0...1000 {
            let t = Double(i)/100
            let p = try safe.state(at: epoch(t)).position
            let c = try Centre(frame: frame, offset: 10).state(at: epoch(t)).position
            let delta = p-c
            let distance = (delta.x*delta.x+delta.y*delta.y).squareRoot()-0.2
            #expect(distance >= report.minimumClearanceMetres)
        }
        #expect(throws: BrachistochroneRouteError.movingObstacleOutsideCoverage(index: 0, obstacle: "crossing")) {
            try route(obstacle(coverageEnd: 9))
        }
        #expect(throws: BrachistochroneRouteError.invalidObstacle) { try obstacle(bound: -1) }
        #expect(throws: BrachistochroneRouteError.unresolvedClearance(index: 0, obstacle: "crossing")) {
            try route(obstacle(offset: 10, bound: 1e15))
        }
    }
    @Test func metadataAndDwellValidation() throws {
        let other = Frame(name: "Other", pseudoInertial: true)
        let mismatch = try MovingSphericalExclusionRegion(identifier: "other",
            centreProvider: MovingWaypointRouteTarget(provider: Centre(frame: other, offset: 10),
                origin: .sun, coverage: epoch(0)...epoch(30)),
            radiusMetres: 1, maximumCentreSpeedMetresPerSecond: 1)
        #expect(throws: SpacecraftDynamicsError.frameMismatch) { try route(mismatch) }
        #expect(throws: BrachistochroneRouteError.invalidObstacle) {
            try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
                waypoints: points(), legAccelerationsMetresPerSecondSquared: [1],
                movingObstacles: [obstacle(offset: 10), obstacle(offset: 10)])
        }
        #expect(throws: BrachistochroneRouteError.unsupportedMovingObstacleDwell) {
            try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
                waypoints: [BrachistochroneWaypoint(positionMetres: .zero, dwellDuration: Duration(seconds: 1)),
                    BrachistochroneWaypoint(positionMetres: SIMD3(10,0,0))],
                legAccelerationsMetresPerSecondSquared: [1], movingObstacles: [obstacle(offset: 10)])
        }
    }
    @Test func candidateTimingChangesClearance() throws {
        let engine = try WaypointRouteEngineLimits(maximumAccelerationMetresPerSecondSquared: 1,
            maximumThrustNewtons: 1000, effectiveExhaustVelocityMetresPerSecond: 100)
        let result = try WaypointRouteTimingPlanner.search(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: points(), legAccelerationsMetresPerSecondSquared: [1],
            candidateArrivalEpochs: [[epoch(10),epoch(20)]],
            initialMass: SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50),
            backgroundDynamics: SpacecraftDynamics(models: [], context: ForceEvaluationContext(frame: frame, origin: .sun)),
            legEngineLimits: [engine], movingObstacles: [obstacle()])
        #expect(result.status == .feasible)
        #expect(result.failedAttempts.first?.failure == .geometry(.obstructedLeg(index: 0, obstacle: "crossing")))
        #expect(try #require(result.route).arrivalEpoch == epoch(20))
    }
}
