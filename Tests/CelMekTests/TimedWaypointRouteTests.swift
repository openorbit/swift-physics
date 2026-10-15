// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct TimedWaypointRouteTests {
    let frame = Frame(name: "Timed route", pseudoInertial: true)
    func plan(_ positions: [SIMD3<Double>]) throws -> WaypointBrachistochroneRoute {
        try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: positions.map { try BrachistochroneWaypoint(positionMetres: $0) },
            legAccelerationsMetresPerSecondSquared: Array(repeating: 1, count: positions.count-1))
    }
    func limits(_ g: Double = 1, exhaust: Double = 100) throws -> WaypointRouteEngineLimits {
        try WaypointRouteEngineLimits(maximumAccelerationMetresPerSecondSquared:
            g*WaypointRouteEngineLimits.standardGravityMetresPerSecondSquared,
            maximumThrustNewtons: 10_000, effectiveExhaustVelocityMetresPerSecond: exhaust)
    }
    func propagator(_ route: WaypointBrachistochroneRoute, phases: [WaypointRouteEnginePhase]) throws -> WaypointRoutePropagator {
        try WaypointRoutePropagator(route: route,
            initialMass: SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50),
            backgroundDynamics: SpacecraftDynamics(models: [], context: ForceEvaluationContext(frame: frame, origin: .sun)),
            legEngineLimits: Array(repeating: limits(), count: route.legs.count),
            timedEnginePhases: phases, numericalOptions: .preset(for: .full))
    }
    @Test func independentTemporalBoundariesAndAnalyticFuel() throws {
        let route = try plan([.zero, SIMD3(4,0,0)])
        let t = route.duration.seconds
        let first = try WaypointRouteEnginePhase(departureEpoch: .j2000,
            startOffset: .zero, duration: Duration(seconds: t/2), limits: limits(1, exhaust: 50))
        let second = try WaypointRouteEnginePhase(startEpoch: first.endEpoch, endEpoch: route.arrivalEpoch,
            limits: limits(0.3, exhaust: 200))
        let execution = try propagator(route, phases: [second, first])
        #expect(execution.timedEnginePhases == [first, second])
        let result = try execution.propagation(at: route.arrivalEpoch)
        #expect(result.segments.count == 2)
        #expect(result.segments.map(\.timedEnginePhaseIndex) == [0, 1])
        let before = result.segments[0].finalState, after = result.segments[1].initialState
        #expect(before.cartesian.position == after.cartesian.position)
        #expect(before.cartesian.velocity == after.cartesian.velocity)
        #expect(before.mass == after.mass)
        #expect(try execution.spacecraftState(at: first.endEpoch).mass == before.mass)
        // Each half contributes peak speed 1.5*d/T to integral |a| dt.
        let expected = 100*exp(-(6/t)/50-(6/t)/200)
        #expect(abs(result.state.mass.totalMassKilograms-expected) < 1e-7)
        #expect(result.checkpoints[0].diagnostics.acceptedSteps == result.segments.reduce(0) { $0+$1.diagnostics.acceptedSteps })
    }
    @Test func spatialPassageAndTemporalPhasesShareOneTimeline() throws {
        let route = try plan([.zero, SIMD3(4,0,0), SIMD3(8,2,0)])
        let passage = route.legs[0].arrivalEpoch
        let early = try route.departureEpoch.advanced(by: Duration(seconds: route.legs[0].duration.seconds/2))
        let late = try passage.advanced(by: Duration(seconds: route.legs[1].duration.seconds/2))
        let phase = try WaypointRouteEnginePhase(startEpoch: early, endEpoch: late, limits: limits(0.3))
        let coincident = try WaypointRouteEnginePhase(startEpoch: passage, endEpoch: late, limits: limits(0.3))
        let execution = try propagator(route, phases: [phase])
        let result = try execution.propagation(at: route.arrivalEpoch)
        #expect(result.segments.count == 4)
        #expect(result.segments.map(\.timedEnginePhaseIndex) == [nil, 0, 0, nil])
        #expect(result.segments[1].finalState.epoch == passage)
        #expect(result.segments[1].finalState.cartesian.velocity.x > 0)
        #expect(result.segments[1].finalState.mass == result.segments[2].initialState.mass)
        let shared = try propagator(route, phases: [coincident]).propagation(at: route.arrivalEpoch)
        #expect(shared.segments.count == 3)
        #expect(shared.checkpoints.count == 2)
        #expect(shared.segments.allSatisfy { $0.finalState.epoch > $0.initialState.epoch })
    }
    @Test func rejectsConflictsAndActivatesLimitsAtExactBoundary() throws {
        let route = try plan([.zero, SIMD3(4,0,0)])
        let quarter = try route.departureEpoch.advanced(by: Duration(seconds: route.duration.seconds/4))
        let phase = try WaypointRouteEnginePhase(startEpoch: quarter, endEpoch: route.arrivalEpoch, limits: limits(0.00001))
        #expect(throws: WaypointRoutePropagationError.accelerationLimitExceeded(leg: 0)) {
            try propagator(route, phases: [phase]).spacecraftState(at: quarter)
        }
        #expect(throws: WaypointRoutePropagationError.overlappingTimedPhases) {
            try propagator(route, phases: [phase, phase])
        }
        #expect(throws: WaypointRoutePropagationError.invalidTimedPhase) {
            try WaypointRouteEnginePhase(startEpoch: quarter, endEpoch: quarter, limits: limits())
        }
        #expect(throws: WaypointRoutePropagationError.timedPhaseOutsideCoverage) {
            let outside = try WaypointRouteEnginePhase(startEpoch: route.departureEpoch,
                endEpoch: route.arrivalEpoch.advanced(by: Duration(seconds: 1)), limits: limits())
            _ = try propagator(route, phases: [outside])
        }
        #expect(throws: WaypointRoutePropagationError.invalidTimedPhase) {
            try WaypointRouteEnginePhase(departureEpoch: .j2000, startOffset: Duration(seconds: -1),
                duration: Duration(seconds: 1), limits: limits())
        }
    }
}
