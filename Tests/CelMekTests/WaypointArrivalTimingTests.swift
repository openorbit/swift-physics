// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct WaypointArrivalTimingTests {
    let frame = Frame(name: "Arrival timing", pseudoInertial: true)
    func epoch(_ seconds: Double) throws -> Epoch { try Epoch.j2000.advanced(by: Duration(seconds: seconds)) }
    func route(_ waypoints: [BrachistochroneWaypoint]) throws -> WaypointBrachistochroneRoute {
        try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: waypoints, legAccelerationsMetresPerSecondSquared: Array(repeating: 1, count: waypoints.count-1))
    }
    @Test func exactArrivalAndEndpointDwellUsePhysicalEpochs() throws {
        let target = try epoch(20.125)
        let plan = try route([
            BrachistochroneWaypoint(positionMetres: .zero, dwellDuration: Duration(seconds: 2)),
            BrachistochroneWaypoint(positionMetres: SIMD3(4,0,0), arrivalEpoch: target, dwellDuration: Duration(seconds: 3))
        ])
        #expect(plan.arrivalEpoch == target)
        #expect(try plan.completionEpoch == target.advanced(by: Duration(seconds: 3)))
        #expect(plan.legs[0].duration.seconds == 18.125)
        let middle = try plan.legs[0].departureEpoch.advanced(by: Duration(seconds: 18.125/2))
        let state = try plan.state(at: middle)
        #expect(abs(state.position.x-2) < 1e-12)
        #expect(abs(state.velocity.x-6/18.125) < 1e-12)
        let initialAcceleration = try #require(plan.legs[0].state(at: plan.legs[0].departureEpoch).acceleration)
        #expect(abs(initialAcceleration.x-24/(18.125*18.125)) < 1e-12)
    }
    @Test func windowsRetimingAndIndependentEngineTimeline() throws {
        let passageWindow = try epoch(15)...epoch(20)
        let arrivalWindow = try epoch(30)...epoch(40)
        let plan = try route([
            BrachistochroneWaypoint(positionMetres: .zero),
            BrachistochroneWaypoint(positionMetres: SIMD3(4,0,0), arrivalWindow: passageWindow),
            BrachistochroneWaypoint(positionMetres: SIMD3(8,0,0), arrivalWindow: arrivalWindow)
        ])
        #expect(passageWindow.contains(plan.waypointArrivalEpochs[1]))
        #expect(arrivalWindow.contains(plan.arrivalEpoch))
        #expect(try abs(plan.waypointArrivalEpochs[1].duration(since: .j2000).seconds-15) < 1e-12)
        #expect(plan.waypointVelocitiesMetresPerSecond[1].x > 0)
        #expect(try plan.legs[0].state(at: plan.waypointArrivalEpochs[1]).velocity == plan.legs[1].state(at: plan.waypointArrivalEpochs[1]).velocity)
        let limits = try WaypointRouteEngineLimits(maximumAccelerationMetresPerSecondSquared: 0.3*9.80665,
            maximumThrustNewtons: 10_000, effectiveExhaustVelocityMetresPerSecond: 100)
        let phase = try WaypointRouteEnginePhase(startEpoch: epoch(5), endEpoch: epoch(25), limits: limits)
        let execution = try WaypointRoutePropagator(route: plan,
            initialMass: SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50),
            backgroundDynamics: SpacecraftDynamics(models: [], context: ForceEvaluationContext(frame: frame, origin: .sun)),
            legEngineLimits: [limits, limits], timedEnginePhases: [phase], numericalOptions: .preset(for: .full))
        let result = try execution.propagation(at: plan.arrivalEpoch)
        #expect(result.segments.count == 4)
        #expect(result.checkpoints[0].state.cartesian.velocity.x > 0)
        #expect(result.checkpoints[0].state.mass == result.segments[2].initialState.mass)
        #expect(result.checkpoints.allSatisfy { $0.positionErrorMetres < 1e-8 && $0.velocityErrorMetresPerSecond < 1e-8 })
    }
    @Test func incompatibleWindowsAreFailureOfUniformRetiming() throws {
        #expect(throws: BrachistochroneRouteError.uniformRetimingUnavailable) {
            try route([
                BrachistochroneWaypoint(positionMetres: .zero),
                BrachistochroneWaypoint(positionMetres: SIMD3(4,0,0), arrivalEpoch: epoch(20)),
                BrachistochroneWaypoint(positionMetres: SIMD3(8,0,0), arrivalEpoch: epoch(21))
            ])
        }
        #expect(throws: BrachistochroneRouteError.uniformRetimingUnavailable) {
            try route([BrachistochroneWaypoint(positionMetres: .zero),
                BrachistochroneWaypoint(positionMetres: SIMD3(4,0,0), arrivalEpoch: epoch(0.01))])
        }
        #expect(throws: BrachistochroneRouteError.departureOutsideArrivalWindow) {
            try route([BrachistochroneWaypoint(positionMetres: .zero, arrivalEpoch: epoch(1)),
                BrachistochroneWaypoint(positionMetres: SIMD3(4,0,0))])
        }
    }
}
