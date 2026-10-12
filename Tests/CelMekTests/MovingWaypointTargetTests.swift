// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct MovingWaypointTargetTests {
    let frame = Frame(name: "Target axes", pseudoInertial: true)
    func epoch(_ seconds: Double) throws -> Epoch { try Epoch.j2000.advanced(by: Duration(seconds: seconds)) }
    func target(origin: BodyIdentifier = .sun) throws -> MovingWaypointRouteTarget {
        try MovingWaypointRouteTarget(provider: LinearTarget(frame: frame), origin: origin, coverage: epoch(3)...epoch(20))
    }
    func points() throws -> [BrachistochroneWaypoint] {
        try [BrachistochroneWaypoint(positionMetres: .zero), BrachistochroneWaypoint(positionMetres: SIMD3(-999,0,0))]
    }
    func limits() throws -> WaypointRouteEngineLimits {
        try WaypointRouteEngineLimits(maximumAccelerationMetresPerSecondSquared: 1,
            maximumThrustNewtons: 1000, effectiveExhaustVelocityMetresPerSecond: 100)
    }
    func background() throws -> SpacecraftDynamics {
        try SpacecraftDynamics(models: [], context: ForceEvaluationContext(frame: frame, origin: .sun))
    }
    func search(_ target: MovingWaypointRouteTarget) throws -> WaypointTimingSearchResult {
        try WaypointRouteTimingPlanner.search(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: points(), legAccelerationsMetresPerSecondSquared: [1],
            candidateArrivalEpochs: [[epoch(2),epoch(4),epoch(8)]],
            initialMass: SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50),
            backgroundDynamics: background(), legEngineLimits: [limits()], numericalOptions: .preset(for: .full), movingTarget: target)
    }
    @Test func movingRendezvousAndIndependentTimingThreshold() throws {
        let moving = try target()
        let found = try search(moving)
        #expect(found.status == .feasible)
        #expect(found.evaluatedCandidates == 3)
        #expect(found.failedAttempts[0].failure == .targetOutsideCoverage)
        let route = try #require(found.route)
        #expect(try route.arrivalEpoch == epoch(8))
        #expect(route.waypoints.last?.positionMetres == SIMD3(8,0,0))
        let final = try #require(found.propagation).state.cartesian
        #expect(abs(final.position.x-8) < 1e-10)
        #expect(abs(final.velocity.x-0.5) < 1e-10)
        let optimized = try WaypointRouteTimingPlanner.optimizeTiming(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: points(), legAccelerationsMetresPerSecondSquared: [1], arrivalWindows: [epoch(4)...epoch(12)],
            seedArrivalEpochs: [epoch(10)],
            optimizationOptions: WaypointTimingOptimizationOptions(initialStepSeconds: 2, minimumStepSeconds: 0.001),
            initialMass: SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50),
            backgroundDynamics: background(), legEngineLimits: [limits()], numericalOptions: .preset(for: .full), movingTarget: moving)
        let improved = try #require(optimized.route)
        // p_target(T)=4+0.5*T, v_target=0.5. Cubic a_depart=24/T²+2/T;
        // acceleration <=1 requires T>=6. This differs from a frozen target.
        #expect(improved.duration.seconds >= 6)
        #expect(improved.duration.seconds-6 < 0.002)
        let expected = try moving.state(at: improved.arrivalEpoch)
        let actual = try #require(optimized.propagation).state.cartesian
        #expect(abs(actual.position.x-expected.position.x) < 1e-10)
        #expect(abs(actual.velocity.x-expected.velocity.x) < 1e-10)
    }
    @Test func ephemerisAdapterAndAdaptiveCoverageIntersection() throws {
        let provider = try TargetEphemeris(metadata: EphemerisMetadata(body: .moon, origin: .sun, frame: frame,
            coverage: epoch(3)...epoch(8), coverageBasis: .callerSelected, modelTimeScale: .tt,
            velocityMethod: .analytic, accelerationMethod: .analytic,
            provenance: PhysicalDataProvenance(source: "Synthetic moving target", version: "1")))
        let moving = try MovingWaypointRouteTarget(ephemeris: provider)
        let waypoint = try moving.waypoint(at: epoch(8))
        #expect(waypoint.positionMetres == SIMD3(8,0,0))
        #expect(waypoint.velocityMetresPerSecond == SIMD3(0.5,0,0))
        let result = try WaypointRouteTimingPlanner.adaptiveSearch(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: points(), legAccelerationsMetresPerSecondSquared: [1], arrivalWindows: [epoch(2)...epoch(20)],
            refinementOptions: WaypointTimingRefinementOptions(initialSpacingSeconds: 10, minimumSpacingSeconds: 1),
            initialMass: SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50),
            backgroundDynamics: background(), legEngineLimits: [limits()], movingTarget: moving)
        #expect(result.status == .feasible)
        #expect(result.evaluatedCandidates == 2)
        #expect(try result.route?.arrivalEpoch == epoch(8))
        #expect(throws: MovingWaypointTargetError.outsideCoverage) { try moving.state(at: epoch(9)) }
    }
    @Test func targetContractsAndPhysicalOriginAreExplicit() throws {
        #expect(throws: MovingWaypointTargetError.originMismatch) { try search(target(origin: .earth)) }
        let other = Frame(name: "Other target axes", pseudoInertial: true)
        let mismatch = try MovingWaypointRouteTarget(provider: LinearTarget(frame: other), origin: .sun, coverage: epoch(3)...epoch(20))
        #expect(throws: MovingWaypointTargetError.frameMismatch) { try search(mismatch) }
        let wrongEpoch = try MovingWaypointRouteTarget(provider: LinearTarget(frame: frame, wrongEpoch: true),
            origin: .sun, coverage: epoch(3)...epoch(20))
        #expect(throws: MovingWaypointTargetError.epochMismatch) { try search(wrongEpoch) }
        let wrongFrame = try MovingWaypointRouteTarget(provider: LinearTarget(frame: frame, returnedFrame: other),
            origin: .sun, coverage: epoch(3)...epoch(20))
        #expect(throws: MovingWaypointTargetError.frameMismatch) { try search(wrongFrame) }
    }
    private struct LinearTarget: StateProvider {
        let frame: Frame
        var wrongEpoch = false
        var returnedFrame: Frame? = nil
        func state(at epoch: Epoch) throws -> CartesianState {
            let elapsed = try epoch.duration(since: .j2000).seconds
            return try CartesianState(epoch: wrongEpoch ? .j2000 : epoch, frame: returnedFrame ?? frame,
                position: SIMD3(4+0.5*elapsed,0,0), velocity: SIMD3(0.5,0,0), acceleration: .zero)
        }
    }
    private struct TargetEphemeris: EphemerisProvider {
        let metadata: EphemerisMetadata
        func state(at epoch: Epoch) throws -> CartesianState { try LinearTarget(frame: frame).state(at: epoch) }
    }
}
