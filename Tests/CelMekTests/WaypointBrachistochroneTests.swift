// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct WaypointBrachistochroneTests {
    private let frame = Frame(name: "Route axes", pseudoInertial: true)
    private func waypoint(_ p: SIMD3<Double>, dwell: Double = 0) throws -> BrachistochroneWaypoint {
        try BrachistochroneWaypoint(positionMetres: p, dwellDuration: Duration(seconds: dwell))
    }
    private func norm(_ p: SIMD3<Double>) -> Double { hypot(hypot(p.x,p.y),p.z) }
    private func obstacle(radius: Double = 2, margin: Double = 1) throws -> SphericalExclusionRegion {
        try SphericalExclusionRegion(identifier: "Sun", frame: frame, origin: .sun, centreMetres: .zero,
            radiusMetres: radius, safetyMarginMetres: margin)
    }
    @Test func flyThroughVelocityContinuityAndIndependentLegLimits() throws {
        let route = try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: [waypoint(.zero),waypoint(SIMD3(100,0,0)),waypoint(SIMD3(100,25,0))],
            legAccelerationsMetresPerSecondSquared: [4,0.25])
        let boundary = route.waypointArrivalEpochs[1]
        let incoming = try route.legs[0].state(at: boundary), outgoing = try route.legs[1].state(at: boundary)
        #expect(incoming.position == SIMD3(100,0,0) && outgoing.position == incoming.position)
        #expect(incoming.velocity == outgoing.velocity && norm(incoming.velocity) > 0)
        #expect(try route.state(at: boundary).velocity == incoming.velocity)
        for leg in route.legs {
            for fraction in [0.0,0.2,0.5,0.8,1] {
                let time = try leg.departureEpoch.advanced(by: Duration(seconds: fraction*leg.duration.seconds))
                let acceleration = try #require(leg.state(at: time).acceleration)
                #expect(norm(acceleration) <= leg.maximumAccelerationMetresPerSecondSquared*(1+1e-12))
            }
        }
        #expect(try route.state(at: route.departureEpoch).velocity == .zero)
        #expect(try route.state(at: route.completionEpoch).velocity == .zero)
        #expect(try route.phase(at: boundary) == .leg(index: 1))
        #expect(try route.state(at: boundary) == route.state(at: boundary))
    }
    @Test func solarAvoidanceAndContinuousCurveCollision() throws {
        let points = try [waypoint(SIMD3(-10,0,0)),waypoint(SIMD3(0,10,0)),waypoint(SIMD3(10,0,0))]
        let route = try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: points, legAccelerationsMetresPerSecondSquared: [1,0.25], obstacles: [obstacle()])
        #expect(route.clearances.count == 2)
        #expect(route.clearances.allSatisfy { $0.minimumClearanceMetres > 0 })
        #expect(norm(route.waypointVelocitiesMetresPerSecond[1]) > 0)
        #expect(throws: BrachistochroneRouteError.obstructedLeg(index: 0, obstacle: "Sun")) {
            try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
                waypoints: [points[0],points[2]], legAccelerationsMetresPerSecondSquared: [1], obstacles: [obstacle()])
        }
        // The Hermite curve bends away from the straight chord: place an obstacle
        // on its interior, while both endpoint positions remain outside.
        let unblocked = try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: points, legAccelerationsMetresPerSecondSquared: [1,0.25])
        let leg = unblocked.legs[0]
        let middle = try leg.state(at: leg.departureEpoch.advanced(by: Duration(seconds: leg.duration.seconds/2))).position
        let blocker = try SphericalExclusionRegion(identifier: "curve", frame: frame, origin: .sun,
            centreMetres: middle, radiusMetres: 0.01)
        #expect(norm(middle-points[0].positionMetres) > 0.01 && norm(middle-points[1].positionMetres) > 0.01)
        #expect(throws: BrachistochroneRouteError.obstructedLeg(index: 0, obstacle: "curve")) {
            try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
                waypoints: points, legAccelerationsMetresPerSecondSquared: [1,0.25], obstacles: [blocker])
        }
    }
    @Test func independentCubicReferenceAndEndpointDwell() throws {
        let route = try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: [waypoint(.zero,dwell: 2),waypoint(SIMD3(4,0,0),dwell: 3)], legAccelerationsMetresPerSecondSquared: [1])
        let leg = route.legs[0], t = leg.duration.seconds
        // With endpoint velocities zero: x=4*(3u²-2u³), v=24*u*(1-u)/T.
        let u = 0.25
        let sample = try leg.state(at: leg.departureEpoch.advanced(by: Duration(seconds: u*t)))
        #expect(abs(sample.position.x-4*(3*u*u-2*u*u*u)) < 1e-13)
        #expect(abs(sample.velocity.x-24*u*(1-u)/t) < 1e-13)
        #expect(try route.phase(at: .j2000) == .dwelling(waypointIndex: 0))
        #expect(try route.phase(at: route.arrivalEpoch) == .dwelling(waypointIndex: 1))
        #expect(try route.phase(at: route.completionEpoch) == .arrived)
    }
    @Test func invalidInputsAndIntermediateStopsRejected() throws {
        #expect(throws: BrachistochroneRouteError.invalidDwell) { try waypoint(.zero,dwell: -1) }
        let points = try [waypoint(SIMD3(10,0,0)),waypoint(SIMD3(20,0,0))]
        #expect(throws: BrachistochroneRouteError.legCountMismatch) {
            try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
                waypoints: points, legAccelerationsMetresPerSecondSquared: [])
        }
        #expect(throws: BrachistochroneRouteError.originMismatch) {
            try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .earth,
                waypoints: points, legAccelerationsMetresPerSecondSquared: [1], obstacles: [obstacle()])
        }
        #expect(throws: BrachistochroneRouteError.unsupportedWaypointDwell) {
            try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
                waypoints: [points[0],waypoint(SIMD3(15,0,0),dwell: 1),points[1]], legAccelerationsMetresPerSecondSquared: [1,1])
        }
        #expect(throws: BrachistochroneRouteError.degeneratePassage(index: 1)) {
            try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
                waypoints: [points[0],points[1],points[0]], legAccelerationsMetresPerSecondSquared: [1,1])
        }
    }
}
