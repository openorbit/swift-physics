//
// SPDX-License-Identifier: Apache-2.0
//
// Copyright 2026 Mattias Holm
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
import Math
import Foundation

/// Moving sphere in the route's native frame and origin. The speed bound is a
/// caller guarantee over the provider's entire declared coverage, not a sample.
public struct MovingSphericalExclusionRegion {
    public let identifier: String
    public let centreProvider: MovingWaypointRouteTarget
    public let radiusMetres: Double
    public let safetyMarginMetres: Double
    public let maximumCentreSpeedMetresPerSecond: Double
    public var exclusionRadiusMetres: Double { radiusMetres + safetyMarginMetres }

    public init(identifier: String, centreProvider: MovingWaypointRouteTarget,
                radiusMetres: Double, safetyMarginMetres: Double = 0,
                maximumCentreSpeedMetresPerSecond: Double) throws {
        guard !identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              radiusMetres.isFinite, radiusMetres > 0,
              safetyMarginMetres.isFinite, safetyMarginMetres >= 0,
              (radiusMetres + safetyMarginMetres).isFinite,
              maximumCentreSpeedMetresPerSecond.isFinite, maximumCentreSpeedMetresPerSecond >= 0 else {
            throw BrachistochroneRouteError.invalidObstacle
        }
        self.identifier = identifier; self.centreProvider = centreProvider
        self.radiusMetres = radiusMetres; self.safetyMarginMetres = safetyMarginMetres
        self.maximumCentreSpeedMetresPerSecond = maximumCentreSpeedMetresPerSecond
    }
}

func validateMovingObstacles(_ obstacles: [MovingSphericalExclusionRegion],
    staticObstacles: [SphericalExclusionRegion], frame: Frame, origin: BodyIdentifier) throws {
    guard obstacles.count <= 10_000 - staticObstacles.count,
          Set(obstacles.map(\.identifier) + staticObstacles.map(\.identifier)).count ==
            obstacles.count + staticObstacles.count else { throw BrachistochroneRouteError.invalidObstacle }
    guard obstacles.allSatisfy({ $0.centreProvider.frame === frame }) else { throw SpacecraftDynamicsError.frameMismatch }
    guard obstacles.allSatisfy({ $0.centreProvider.origin == origin }) else { throw BrachistochroneRouteError.originMismatch }
}

func movingClearance(leg: WaypointBrachistochroneLeg, index: Int,
                     obstacle: MovingSphericalExclusionRegion, checks: inout Int) throws -> BrachistochroneRouteClearance {
    let provider = obstacle.centreProvider
    guard provider.coverage.contains(leg.departureEpoch), provider.coverage.contains(leg.arrivalEpoch) else {
        throw BrachistochroneRouteError.movingObstacleOutsideCoverage(index: index, obstacle: obstacle.identifier)
    }
    var pending = [(leg.controls, leg.departureEpoch, leg.arrivalEpoch, 0)]
    var lowerBound = Double.infinity
    while let (polygon, start, end, depth) = pending.popLast() {
        checks += 1
        guard checks <= 1_000_000 else { throw BrachistochroneRouteError.clearanceLimitExceeded }
        let half = try end.duration(since: start).seconds / 2
        let middle = try start.advanced(by: Duration(seconds: half))
        let epochs = [start, middle, end]
        let positions = [polygon[0], bezierPoint(polygon, 0.5), polygon[3]]
        var centres: [SIMD3<Double>] = []
        for (epoch, position) in zip(epochs, positions) {
            let centre = try provider.state(at: epoch).position
            centres.append(centre)
            let distance = routeNorm(position-centre)
            guard distance.isFinite else { throw BrachistochroneError.numericalRangeExceeded }
            guard distance > obstacle.exclusionRadiusMetres else {
                throw BrachistochroneRouteError.obstructedLeg(index: index, obstacle: obstacle.identifier)
            }
        }
        let centre = centres[1]
        let reach = obstacle.maximumCentreSpeedMetresPerSecond *
            max(try middle.duration(since: start).seconds, try end.duration(since: middle).seconds)
        let scale = max(1, max(routeNorm(centre), polygon.map(routeNorm).max()!))
        let roundoff = 1024 * Double.ulpOfOne * max(scale, reach)
        var offset = SIMD3<Double>.zero
        for axis in 0..<3 {
            let low = polygon.map { $0[axis] }.min()!, high = polygon.map { $0[axis] }.max()!
            offset[axis] = max(0, max(low-centre[axis], centre[axis]-high))
        }
        let bound = routeNorm(offset) - obstacle.exclusionRadiusMetres - reach
        guard bound.isFinite, roundoff.isFinite else { throw BrachistochroneError.numericalRangeExceeded }
        if bound > roundoff { lowerBound = min(lowerBound, bound-roundoff); continue }
        guard depth < 24, middle > start, middle < end else {
            throw BrachistochroneRouteError.unresolvedClearance(index: index, obstacle: obstacle.identifier)
        }
        let split = splitBezier(polygon)
        pending.append((split.0, start, middle, depth+1))
        pending.append((split.1, middle, end, depth+1))
    }
    return BrachistochroneRouteClearance(legIndex: index, obstacleIdentifier: obstacle.identifier,
        minimumClearanceMetres: lowerBound)
}
