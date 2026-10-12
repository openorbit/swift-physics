// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum BrachistochroneRouteError: Error, Equatable {
    case invalidWaypointCount, legCountMismatch, invalidDwell, invalidObstacle
    case originMismatch, clearanceLimitExceeded, unsupportedWaypointDwell, degeneratePassage(index: Int)
    case unresolvedClearance(index: Int, obstacle: String), obstructedLeg(index: Int, obstacle: String)
}

/// Fixed inertial fly-through waypoint. Dwell is permitted only at route endpoints.
public struct BrachistochroneWaypoint {
    public let positionMetres: SIMD3<Double>
    public let dwellDuration: Duration
    public init(positionMetres: SIMD3<Double>, dwellDuration: Duration = .zero) throws {
        guard [positionMetres.x,positionMetres.y,positionMetres.z].allSatisfy(\.isFinite) else {
            throw BrachistochroneError.nonFiniteEndpoint
        }
        guard dwellDuration >= .zero else { throw BrachistochroneRouteError.invalidDwell }
        self.positionMetres = positionMetres; self.dwellDuration = dwellDuration
    }
}

/// Fixed spherical exclusion region, including a caller-selected safety margin.
public struct SphericalExclusionRegion {
    public let identifier: String
    public let frame: Frame
    public let origin: BodyIdentifier
    public let centreMetres: SIMD3<Double>
    public let radiusMetres: Double
    public let safetyMarginMetres: Double
    public var exclusionRadiusMetres: Double { radiusMetres+safetyMarginMetres }
    public init(identifier: String, frame: Frame, origin: BodyIdentifier, centreMetres: SIMD3<Double>,
                radiusMetres: Double, safetyMarginMetres: Double = 0) throws {
        try origin.validate()
        guard !identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              [centreMetres.x,centreMetres.y,centreMetres.z].allSatisfy(\.isFinite),
              radiusMetres.isFinite, radiusMetres > 0, safetyMarginMetres.isFinite, safetyMarginMetres >= 0,
              (radiusMetres+safetyMarginMetres).isFinite else { throw BrachistochroneRouteError.invalidObstacle }
        self.identifier = identifier; self.frame = frame; self.origin = origin; self.centreMetres = centreMetres
        self.radiusMetres = radiusMetres; self.safetyMarginMetres = safetyMarginMetres
    }
}
public struct BrachistochroneRouteClearance {
    public let legIndex: Int
    public let obstacleIdentifier: String
    /// Conservative lower bound for the entire curved leg's surface clearance.
    public let minimumClearanceMetres: Double
}
public enum BrachistochroneRoutePhase: Equatable, Sendable {
    case leg(index: Int)
    case dwelling(waypointIndex: Int)
    case arrived
}

/// Cubic Hermite leg with continuous endpoint velocity and a certified acceleration
/// magnitude bound (acceleration is affine in time). No gravity/fuel model.
public struct WaypointBrachistochroneLeg: StateProvider {
    public let frame: Frame
    public let departureEpoch: Epoch
    public let arrivalEpoch: Epoch
    public let duration: Duration
    public let maximumAccelerationMetresPerSecondSquared: Double
    public let departurePositionMetres: SIMD3<Double>
    public let arrivalPositionMetres: SIMD3<Double>
    public let departureVelocityMetresPerSecond: SIMD3<Double>
    public let arrivalVelocityMetresPerSecond: SIMD3<Double>
    let controls: [SIMD3<Double>]
    public var coverage: ClosedRange<Epoch> { departureEpoch...arrivalEpoch }
    public func state(at epoch: Epoch) throws -> CartesianState {
        guard coverage.contains(epoch) else { throw BrachistochroneError.outsideCoverage }
        let u = epoch == arrivalEpoch ? 1 : try epoch.duration(since: departureEpoch).seconds/duration.seconds
        let t = duration.seconds, v = 1-u
        let d0 = controls[1]-controls[0], d1 = controls[2]-controls[1], d2 = controls[3]-controls[2]
        let velocity = epoch == departureEpoch ? departureVelocityMetresPerSecond : epoch == arrivalEpoch ? arrivalVelocityMetresPerSecond :
            (d0*(v*v)+d1*(2*u*v)+d2*(u*u))*(3/t)
        let acceleration = ((d1-d0)*v+(d2-d1)*u)*(6/t/t)
        return try CartesianState(epoch: epoch, frame: frame,
            position: epoch == departureEpoch ? departurePositionMetres : epoch == arrivalEpoch ? arrivalPositionMetres : bezierPoint(controls,u),
            velocity: velocity, acceleration: acceleration)
    }
}

/// Feasible gravity-free fly-through route with independently supplied leg limits.
/// Cubic curves carry continuous velocity through waypoints; endpoints are at rest.
/// Durations are conservatively scaled to satisfy acceleration bounds. This is a
/// smooth approximation, not a minimum-time/global brachistochrone optimizer.
public struct WaypointBrachistochroneRoute: StateProvider {
    public let frame: Frame
    public let origin: BodyIdentifier
    public let departureEpoch: Epoch
    public let arrivalEpoch: Epoch
    public let completionEpoch: Epoch
    public let waypoints: [BrachistochroneWaypoint]
    public let legs: [WaypointBrachistochroneLeg]
    public let waypointArrivalEpochs: [Epoch]
    public let waypointDepartureEpochs: [Epoch]
    public let waypointVelocitiesMetresPerSecond: [SIMD3<Double>]
    public let clearances: [BrachistochroneRouteClearance]
    public let duration: Duration
    public var coverage: ClosedRange<Epoch> { departureEpoch...completionEpoch }

    public init(departureEpoch: Epoch, frame: Frame, origin: BodyIdentifier,
                waypoints: [BrachistochroneWaypoint], legAccelerationsMetresPerSecondSquared: [Double],
                obstacles: [SphericalExclusionRegion] = []) throws {
        try origin.validate()
        guard frame.pseudoInertial else { throw BrachistochroneError.nonInertialFrame }
        guard (2...10_001).contains(waypoints.count), obstacles.count <= 10_000 else {
            throw BrachistochroneRouteError.invalidWaypointCount
        }
        guard legAccelerationsMetresPerSecondSquared.count == waypoints.count-1 else { throw BrachistochroneRouteError.legCountMismatch }
        guard waypoints.dropFirst().dropLast().allSatisfy({ $0.dwellDuration == .zero }) else {
            throw BrachistochroneRouteError.unsupportedWaypointDwell
        }
        guard obstacles.isEmpty || waypoints.count-1 <= 1_000_000/obstacles.count else { throw BrachistochroneRouteError.clearanceLimitExceeded }
        guard obstacles.allSatisfy({ $0.frame === frame }) else { throw SpacecraftDynamicsError.frameMismatch }
        guard obstacles.allSatisfy({ $0.origin == origin }) else { throw BrachistochroneRouteError.originMismatch }
        guard Set(obstacles.map(\.identifier)).count == obstacles.count else { throw BrachistochroneRouteError.invalidObstacle }
        let points = waypoints.map(\.positionMetres)
        var times: [Double] = []
        for i in 0..<points.count-1 {
            let leg = try IdealBrachistochroneTrajectory(departureEpoch: departureEpoch, frame: frame,
                departurePositionMetres: points[i], arrivalPositionMetres: points[i+1],
                maximumAccelerationMetresPerSecondSquared: legAccelerationsMetresPerSecondSquared[i])
            guard leg.distanceMetres > 0 else { throw BrachistochroneRouteError.degeneratePassage(index: i) }
            times.append(leg.duration.seconds)
        }
        func velocities(_ durations: [Double]) throws -> [SIMD3<Double>] {
            var result = [SIMD3<Double>.zero]
            for i in 1..<points.count-1 {
                let value = (points[i+1]-points[i-1])/(durations[i-1]+durations[i])
                guard routeNorm(value).isFinite else { throw BrachistochroneError.numericalRangeExceeded }
                guard routeNorm(value) > 0 else { throw BrachistochroneRouteError.degeneratePassage(index: i) }
                result.append(value)
            }
            result.append(.zero); return result
        }
        let nominalVelocities = try velocities(times)
        var scale = 1.0
        for i in times.indices {
            let c = routeControls(points[i],points[i+1],nominalVelocities[i],nominalVelocities[i+1],times[i])
            guard c.allSatisfy({ $0.x.isFinite && $0.y.isFinite && $0.z.isFinite }) else {
                throw BrachistochroneError.numericalRangeExceeded
            }
            let bound = max(routeNorm((c[2]-c[1]-(c[1]-c[0]))*(6/times[i]/times[i])),
                routeNorm((c[3]-c[2]-(c[2]-c[1]))*(6/times[i]/times[i])))
            guard bound.isFinite else { throw BrachistochroneError.numericalRangeExceeded }
            scale = max(scale,sqrt(bound)/sqrt(legAccelerationsMetresPerSecondSquared[i]))
        }
        guard scale.isFinite else { throw BrachistochroneError.numericalRangeExceeded }
        // Margin protects the analytic bound from final floating-point rounding.
        times = times.map { $0*scale*1.01 }
        let passageVelocities = try velocities(times)
        var arrivals = [departureEpoch], departures: [Epoch] = [], paths: [WaypointBrachistochroneLeg] = []
        var reports: [BrachistochroneRouteClearance] = [], clearanceChecks = 0
        for i in waypoints.indices {
            let leave = try arrivals[i].advanced(by: waypoints[i].dwellDuration)
            departures.append(leave)
            guard i < times.count else { break }
            let interval = try Duration(seconds: times[i]), end = try leave.advanced(by: interval)
            guard end > leave else { throw BrachistochroneError.numericalRangeExceeded }
            let c = routeControls(points[i],points[i+1],passageVelocities[i],passageVelocities[i+1],interval.seconds)
            guard c.allSatisfy({ $0.x.isFinite && $0.y.isFinite && $0.z.isFinite }) else {
                throw BrachistochroneError.numericalRangeExceeded
            }
            let actualBound = max(routeNorm((c[2]-c[1]-(c[1]-c[0]))*(6/interval.seconds/interval.seconds)),
                routeNorm((c[3]-c[2]-(c[2]-c[1]))*(6/interval.seconds/interval.seconds)))
            guard actualBound.isFinite, actualBound <= legAccelerationsMetresPerSecondSquared[i] else {
                throw BrachistochroneError.numericalRangeExceeded
            }
            let leg = WaypointBrachistochroneLeg(frame: frame, departureEpoch: leave, arrivalEpoch: end, duration: interval,
                maximumAccelerationMetresPerSecondSquared: legAccelerationsMetresPerSecondSquared[i],
                departurePositionMetres: points[i], arrivalPositionMetres: points[i+1],
                departureVelocityMetresPerSecond: passageVelocities[i], arrivalVelocityMetresPerSecond: passageVelocities[i+1], controls: c)
            for obstacle in obstacles {
                let coordinateScale = max(1,max(routeNorm(obstacle.centreMetres),c.map(routeNorm).max()!))
                let roundoffMargin = 1024*Double.ulpOfOne*coordinateScale
                guard roundoffMargin.isFinite else { throw BrachistochroneError.numericalRangeExceeded }
                var pending = [(c,0)], lowerBound = Double.infinity
                while let (polygon,depth) = pending.popLast() {
                    clearanceChecks += 1
                    guard clearanceChecks <= 1_000_000 else { throw BrachistochroneRouteError.clearanceLimitExceeded }
                    let middle = bezierPoint(polygon,0.5)
                    for sample in [polygon[0],middle,polygon[3]] {
                        let distance = routeNorm(sample-obstacle.centreMetres)
                        guard distance.isFinite else { throw BrachistochroneError.numericalRangeExceeded }
                        guard distance > obstacle.exclusionRadiusMetres else {
                            throw BrachistochroneRouteError.obstructedLeg(index: i, obstacle: obstacle.identifier)
                        }
                    }
                    var offset = SIMD3<Double>.zero
                    for axis in 0..<3 {
                        let low = polygon.map { $0[axis] }.min()!, high = polygon.map { $0[axis] }.max()!
                        offset[axis] = max(0,max(low-obstacle.centreMetres[axis],obstacle.centreMetres[axis]-high))
                    }
                    let bound = routeNorm(offset)-obstacle.exclusionRadiusMetres
                    guard bound.isFinite else { throw BrachistochroneError.numericalRangeExceeded }
                    if bound > roundoffMargin { lowerBound = min(lowerBound,bound-roundoffMargin); continue }
                    guard depth < 24 else { throw BrachistochroneRouteError.unresolvedClearance(index: i, obstacle: obstacle.identifier) }
                    let split = splitBezier(polygon)
                    pending.append((split.0,depth+1)); pending.append((split.1,depth+1))
                }
                reports.append(BrachistochroneRouteClearance(legIndex: i, obstacleIdentifier: obstacle.identifier,
                    minimumClearanceMetres: lowerBound))
            }
            paths.append(leg); arrivals.append(end)
        }
        self.frame = frame; self.origin = origin; self.departureEpoch = departureEpoch
        arrivalEpoch = arrivals.last!; completionEpoch = departures.last!
        duration = try completionEpoch.duration(since: departureEpoch)
        self.waypoints = waypoints; legs = paths; waypointArrivalEpochs = arrivals; waypointDepartureEpochs = departures
        waypointVelocitiesMetresPerSecond = passageVelocities; clearances = reports
    }
    public func phase(at epoch: Epoch) throws -> BrachistochroneRoutePhase {
        guard coverage.contains(epoch) else { throw BrachistochroneError.outsideCoverage }
        if epoch == completionEpoch { return .arrived }
        for i in waypoints.indices {
            if epoch >= waypointArrivalEpochs[i] && epoch < waypointDepartureEpochs[i] { return .dwelling(waypointIndex: i) }
            if i < legs.count && epoch >= legs[i].departureEpoch && epoch < legs[i].arrivalEpoch { return .leg(index: i) }
        }
        throw BrachistochroneError.outsideCoverage
    }
    public func state(at epoch: Epoch) throws -> CartesianState {
        switch try phase(at: epoch) {
        case .leg(let index): return try legs[index].state(at: epoch)
        case .dwelling(let index): return try CartesianState(epoch: epoch, frame: frame, position: waypoints[index].positionMetres, velocity: .zero, acceleration: .zero)
        case .arrived: return try CartesianState(epoch: epoch, frame: frame, position: waypoints.last!.positionMetres, velocity: .zero, acceleration: .zero)
        }
    }
}
private func routeNorm(_ v: SIMD3<Double>) -> Double { hypot(hypot(v.x,v.y),v.z) }
private func routeControls(_ p0: SIMD3<Double>, _ p1: SIMD3<Double>, _ v0: SIMD3<Double>, _ v1: SIMD3<Double>, _ t: Double) -> [SIMD3<Double>] {
    [p0,p0+v0*(t/3),p1-v1*(t/3),p1]
}
private func bezierPoint(_ c: [SIMD3<Double>], _ u: Double) -> SIMD3<Double> {
    let a = c[0]*(1-u)+c[1]*u, b = c[1]*(1-u)+c[2]*u, d = c[2]*(1-u)+c[3]*u
    return (a*(1-u)+b*u)*(1-u)+(b*(1-u)+d*u)*u
}
private func splitBezier(_ c: [SIMD3<Double>]) -> ([SIMD3<Double>],[SIMD3<Double>]) {
    let a = c[0]*0.5+c[1]*0.5, b = c[1]*0.5+c[2]*0.5, d = c[2]*0.5+c[3]*0.5
    let e = a*0.5+b*0.5, f = b*0.5+d*0.5, m = e*0.5+f*0.5
    return ([c[0],a,e,m],[m,f,d,c[3]])
}
