// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct WaypointVelocityConstraintTests {
    let frame = Frame(name: "Velocity constraints", pseudoInertial: true)
    func epoch(_ seconds: Double) throws -> Epoch { try Epoch.j2000.advanced(by: Duration(seconds: seconds)) }
    func limits() throws -> WaypointRouteEngineLimits {
        try WaypointRouteEngineLimits(maximumAccelerationMetresPerSecondSquared: 1,
            maximumThrustNewtons: 1000, effectiveExhaustVelocityMetresPerSecond: 100)
    }
    func background() throws -> SpacecraftDynamics {
        try SpacecraftDynamics(models: [], context: ForceEvaluationContext(frame: frame, origin: .sun))
    }
    func route(_ points: [BrachistochroneWaypoint], times: [Duration]? = nil) throws -> WaypointBrachistochroneRoute {
        try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: points, legAccelerationsMetresPerSecondSquared: Array(repeating: 1, count: points.count-1), legDurations: times)
    }
    @Test func independentCoastingAndConstantDecelerationReferences() throws {
        let coast = try route([
            BrachistochroneWaypoint(positionMetres: .zero, velocityMetresPerSecond: SIMD3(1,0,0)),
            BrachistochroneWaypoint(positionMetres: SIMD3(3,0,0), velocityMetresPerSecond: SIMD3(1,0,0))
        ], times: [Duration(seconds: 3)])
        let coasting = try WaypointRoutePropagator(route: coast,
            initialMass: SpacecraftMassProperties(dryMassKilograms: 50), backgroundDynamics: background(), legEngineLimits: [limits()])
        #expect(try coasting.spacecraftState(at: coast.arrivalEpoch).mass.propellantMassKilograms == 0)
        #expect(try coast.state(at: coast.arrivalEpoch).velocity == SIMD3(1,0,0))
        #expect(try abs(coasting.state(at: epoch(1)).position.x-1) < 1e-12)

        let braking = try route([
            BrachistochroneWaypoint(positionMetres: .zero, velocityMetresPerSecond: SIMD3(2,0,0)),
            BrachistochroneWaypoint(positionMetres: SIMD3(4,0,0), velocityMetresPerSecond: .zero)
        ], times: [Duration(seconds: 4)])
        let middle = try braking.legs[0].state(at: epoch(2))
        #expect(abs(middle.position.x-3) < 1e-12)
        #expect(abs(middle.velocity.x-1) < 1e-12)
        let acceleration = try #require(middle.acceleration)
        #expect(abs(acceleration.x+0.5) < 1e-12)
        let execution = try WaypointRoutePropagator(route: braking,
            initialMass: SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50),
            backgroundDynamics: background(), legEngineLimits: [limits()], numericalOptions: .preset(for: .full))
        let result = try execution.propagation(at: braking.arrivalEpoch)
        #expect(abs(result.state.mass.totalMassKilograms-100*exp(-2/100)) < 1e-8)
        #expect(abs(result.state.cartesian.velocity.x) < 1e-12)
    }
    @Test func timingSearchAndOptimizationRetainEndpointConstraints() throws {
        let points = try [BrachistochroneWaypoint(positionMetres: .zero, velocityMetresPerSecond: SIMD3(0.5,0,0)),
            BrachistochroneWaypoint(positionMetres: SIMD3(4,0,0), velocityMetresPerSecond: SIMD3(0.25,0,0))]
        let mass = try SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50)
        let found = try WaypointRouteTimingPlanner.search(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: points, legAccelerationsMetresPerSecondSquared: [1], candidateArrivalEpochs: [[epoch(10)]],
            initialMass: mass, backgroundDynamics: background(), legEngineLimits: [limits()])
        #expect(found.status == .feasible)
        let optimized = try WaypointRouteTimingPlanner.optimizeTiming(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: points, legAccelerationsMetresPerSecondSquared: [1],
            arrivalWindows: [epoch(3)...epoch(12)], seedArrivalEpochs: [epoch(10)],
            optimizationOptions: WaypointTimingOptimizationOptions(initialStepSeconds: 1, minimumStepSeconds: 0.01),
            initialMass: mass, backgroundDynamics: background(), legEngineLimits: [limits()], numericalOptions: .preset(for: .full))
        let plan = try #require(optimized.route)
        #expect(try plan.state(at: .j2000).velocity == SIMD3(0.5,0,0))
        #expect(try plan.state(at: plan.arrivalEpoch).velocity == SIMD3(0.25,0,0))
        let execution = try #require(optimized.propagation)
        #expect(abs(execution.state.cartesian.velocity.x-0.25) < 1e-10)
        #expect(plan.duration.seconds < 10)
    }
    @Test func automaticTimingAndExplicitPassagesValidateVelocityAndDwell() throws {
        let points = try [BrachistochroneWaypoint(positionMetres: .zero, velocityMetresPerSecond: SIMD3(10,0,0)),
            BrachistochroneWaypoint(positionMetres: SIMD3(4,1,0), velocityMetresPerSecond: SIMD3(1,0,0)),
            BrachistochroneWaypoint(positionMetres: SIMD3(8,0,0), velocityMetresPerSecond: SIMD3(2,0,0))]
        let plan = try route(points)
        #expect(plan.waypointVelocitiesMetresPerSecond == [SIMD3(10,0,0), SIMD3(1,0,0), SIMD3(2,0,0)])
        let boundary = plan.legs[0].arrivalEpoch
        #expect(try plan.legs[0].state(at: boundary).velocity == plan.legs[1].state(at: boundary).velocity)
        for leg in plan.legs {
            for end in [leg.departureEpoch, leg.arrivalEpoch] {
                let a = try #require(leg.state(at: end).acceleration)
                #expect(hypot(hypot(a.x,a.y),a.z) <= 1)
            }
        }
        #expect(throws: BrachistochroneRouteError.velocityDwellConflict(index: 0)) {
            try route([BrachistochroneWaypoint(positionMetres: .zero, dwellDuration: Duration(seconds: 1),
                velocityMetresPerSecond: SIMD3(1,0,0)), points[2]])
        }
        #expect(throws: BrachistochroneRouteError.degeneratePassage(index: 1)) {
            try route([points[0], BrachistochroneWaypoint(positionMetres: points[1].positionMetres, velocityMetresPerSecond: .zero), points[2]])
        }
        #expect(throws: BrachistochroneError.nonFiniteEndpoint) {
            try BrachistochroneWaypoint(positionMetres: .zero, velocityMetresPerSecond: SIMD3(.nan,0,0))
        }
    }
}
