// SPDX-License-Identifier: Apache-2.0
import Math

public enum MovingWaypointTargetError: Error, Equatable {
    case invalidWaypointIndex, conflictingWaypointVelocity(index: Int), waypointOutsideCoverage(index: Int)
    case outsideCoverage, frameMismatch, epochMismatch, originMismatch, conflictingEndpointVelocity
}

public enum MovingWaypointPassageVelocity: Equatable, Sendable {
    /// Keep an explicit waypoint velocity, or generate the usual fly-through value.
    case preserveRouteVelocity
    case matchProviderVelocity
}

/// A moving intermediate position constraint. Matching a body's velocity is
/// optional: passing a stationary body must not implicitly stop the spacecraft.
public struct MovingWaypointPassage {
    public let provider: MovingWaypointRouteTarget
    public let velocityMode: MovingWaypointPassageVelocity
    public init(provider: MovingWaypointRouteTarget,
                velocityMode: MovingWaypointPassageVelocity = .preserveRouteVelocity) {
        self.provider = provider; self.velocityMode = velocityMode
    }
    public func waypoint(at epoch: Epoch, template: BrachistochroneWaypoint) throws -> BrachistochroneWaypoint {
        let state = try provider.state(at: epoch)
        return try BrachistochroneWaypoint(positionMetres: state.position, dwellDuration: template.dwellDuration,
            arrivalWindow: template.arrivalWindow,
            velocityMetresPerSecond: velocityMode == .matchProviderVelocity ? state.velocity : template.velocityMetresPerSecond)
    }
}

extension WaypointRouteTimingPlanner {
    static func intersectPassageCoverage(_ windows: [ClosedRange<Epoch>],
        movingWaypoints: [Int: MovingWaypointPassage]) throws -> [ClosedRange<Epoch>] {
        var result = windows
        for index in movingWaypoints.keys.sorted() {
            let coverage = movingWaypoints[index]!.provider.coverage
            let lower = max(result[index-1].lowerBound, coverage.lowerBound)
            let upper = min(result[index-1].upperBound, coverage.upperBound)
            guard lower <= upper else { throw WaypointTimingSearchError.invalidCandidates }
            result[index-1] = lower...upper
        }
        return result
    }
}

/// A rendezvous target relative to an explicit physical origin. Matching uses
/// geometric position and velocity, with no light-time or apparent-position correction.
public struct MovingWaypointRouteTarget {
    public let frame: Frame
    public let origin: BodyIdentifier
    public let coverage: ClosedRange<Epoch>
    private let evaluate: (Epoch) throws -> CartesianState

    /// For arbitrary spacecraft/state providers, the caller supplies and guarantees
    /// physical origin and coverage. The provider must remain pure during planning.
    public init(provider: any StateProvider, origin: BodyIdentifier, coverage: ClosedRange<Epoch>) throws {
        try origin.validate()
        self.frame = provider.frame; self.origin = origin; self.coverage = coverage
        evaluate = { try provider.state(at: $0) }
    }
    /// Uses the ephemeris' own physical origin and coverage, with contract checking.
    public init(ephemeris: any EphemerisProvider) throws {
        try ephemeris.metadata.origin.validate()
        guard ephemeris.frame === ephemeris.metadata.frame else { throw MovingWaypointTargetError.frameMismatch }
        frame = ephemeris.frame; origin = ephemeris.metadata.origin; coverage = ephemeris.metadata.coverage
        evaluate = { try ephemeris.checkedState(at: $0) }
    }
    public func state(at epoch: Epoch) throws -> CartesianState {
        guard coverage.contains(epoch) else { throw MovingWaypointTargetError.outsideCoverage }
        let state = try evaluate(epoch)
        guard state.epoch == epoch else { throw MovingWaypointTargetError.epochMismatch }
        guard state.frame === frame else { throw MovingWaypointTargetError.frameMismatch }
        return state
    }
    /// Useful for a single fixed-time route. Timing planners resolve this anew
    /// for every candidate instead of freezing a waypoint sampled at departure.
    public func waypoint(at epoch: Epoch) throws -> BrachistochroneWaypoint {
        let state = try state(at: epoch)
        return try BrachistochroneWaypoint(positionMetres: state.position, arrivalEpoch: epoch,
            velocityMetresPerSecond: state.velocity)
    }
}

extension MovingWaypointRouteTarget {
    /// A constant displacement in the target's native axes, retaining its velocity.
    public func offset(positionMetres: SIMD3<Double>) throws -> MovingWaypointRouteTarget {
        _ = try CartesianState(epoch: coverage.lowerBound, frame: frame,
            position: positionMetres, velocity: .zero)
        return MovingWaypointRouteTarget(frame: frame, origin: origin, coverage: coverage) { epoch in
            let base = try self.state(at: epoch)
            return try CartesianState(epoch: epoch, frame: frame, position: base.position + positionMetres,
                velocity: base.velocity, acceleration: base.acceleration)
        }
    }

    /// A constant target-relative displacement expressed in possibly rotating axes.
    /// Only axes and their derivatives are used; the axes frame's translation is
    /// ignored because the physical centre is supplied by this target.
    public func offset(positionMetres: SIMD3<Double>, in axes: Frame,
                       using tree: FrameTree) throws -> MovingWaypointRouteTarget {
        _ = try CartesianState(epoch: coverage.lowerBound, frame: axes,
            position: positionMetres, velocity: .zero)
        _ = try tree.transform(from: axes, to: frame, at: coverage.lowerBound)
        return MovingWaypointRouteTarget(frame: frame, origin: origin, coverage: coverage) { epoch in
            let base = try self.state(at: epoch)
            let transform = try tree.transform(from: axes, to: frame, at: epoch)
            let r = transform.rotation.matrix * positionMetres
            let v = offsetCross(transform.angularVelocity, r)
            let a = offsetCross(transform.angularAcceleration, r) + offsetCross(transform.angularVelocity, v)
            return try CartesianState(epoch: epoch, frame: frame, position: base.position + r,
                velocity: base.velocity + v, acceleration: base.acceleration.map { $0 + a })
        }
    }

    private init(frame: Frame, origin: BodyIdentifier, coverage: ClosedRange<Epoch>,
                 evaluate: @escaping (Epoch) throws -> CartesianState) {
        self.frame = frame; self.origin = origin; self.coverage = coverage; self.evaluate = evaluate
    }
}

private func offsetCross(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> SIMD3<Double> {
    SIMD3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x)
}
