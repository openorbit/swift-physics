// SPDX-License-Identifier: Apache-2.0
import Math

public enum MovingWaypointTargetError: Error, Equatable {
    case outsideCoverage, frameMismatch, epochMismatch, originMismatch, conflictingEndpointVelocity
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
