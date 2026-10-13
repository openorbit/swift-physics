// SPDX-License-Identifier: Apache-2.0
import Math

/// Conditions at the central-body two-body arc endpoint. Flyby preserves the
/// incoming velocity; it does not solve a body's local gravity-assist geometry.
public enum LambertArrivalCondition: Equatable, Sendable {
    case intercept, flyby, velocityMatchedRendezvous
}
public enum LambertTransferPlanningError: Error, Equatable {
    case invalidTolerance, invalidEpochs, originMismatch, numericalRangeExceeded, outsideFlight
}

/// A bounded transfer arc, including the two instantaneous endpoint state changes.
/// No post-encounter orbit or pre-departure ephemeris history is extrapolated.
public struct LambertTransferTrajectory: Propagator {
    public let coast: TwoBodyPropagator
    public let departureReferenceState: CartesianState
    public let postArrivalState: CartesianState
    public var frame: Frame { coast.frame }
    public var centralBody: BodyIdentifier { coast.centralBody }
    public var referenceEpoch: Epoch { coast.referenceEpoch }
    public var coverage: ClosedRange<Epoch> { referenceEpoch...postArrivalState.epoch }
    public func state(at epoch: Epoch) throws -> CartesianState { try state(at: epoch, side: .after) }
    public func state(at epoch: Epoch, side: ManeuverEpochSide) throws -> CartesianState {
        guard coverage.contains(epoch) else { throw LambertTransferPlanningError.outsideFlight }
        if epoch == referenceEpoch && side == .before { return departureReferenceState }
        if epoch == postArrivalState.epoch && side == .after { return postArrivalState }
        return try coast.state(at: epoch)
    }
}

public struct LambertTransferPlan {
    public let solution: LambertSolution
    /// Physical state before departure ignition, or the sampled departure body's
    /// state when planning ideal co-located body-to-body velocity differences.
    public let departureReferenceState: CartesianState
    public let arrivalTargetState: CartesianState
    public let arrivalCondition: LambertArrivalCondition
    public let departureBurn: ImpulsiveManeuver
    public let arrivalBurn: ImpulsiveManeuver?
    public let preArrivalState: CartesianState
    public let postArrivalState: CartesianState
    /// Relative vectors in the common inertial axes. For sampled bodies these are
    /// patched-conic excess-velocity inputs, not parking-orbit injection costs.
    public let departureRelativeVelocityMetresPerSecond: SIMD3<Double>
    public let incomingRelativeVelocityMetresPerSecond: SIMD3<Double>
    public let outgoingRelativeVelocityMetresPerSecond: SIMD3<Double>
    public let totalDeltaVelocityMetresPerSecond: Double
    public let positionErrorMetres: Double
    public let velocityErrorMetresPerSecond: Double
    public let arrivalConditionSatisfied: Bool
    /// Begins immediately after the departure impulse. Arrival queries are
    /// right-continuous when a rendezvous impulse is present. No pre-departure
    /// spacecraft history is inferred from a body's ephemeris.
    public let trajectory: LambertTransferTrajectory
    public var maneuvers: [ImpulsiveManeuver] { [departureBurn] + (arrivalBurn.map { [$0] } ?? []) }
}

/// Fixed-epoch, central-body two-body transfers, with actual velocity differences
/// and independently propagated endpoint diagnostics. Sources must be geometric
/// states relative to the same physical origin in the same inertial frame.
public enum LambertTransferPlanner {
    public static func plan(departureState: CartesianState, arrivalTargetState: CartesianState,
        centralBody: BodyIdentifier, gravitationalParameter: GravitationalParameter,
        referencePole: SIMD3<Double>, direction: LambertDirection,
        arrivalCondition: LambertArrivalCondition, revolutions: Int = 0,
        multiRevolutionBranch: LambertMultiRevolutionBranch? = nil,
        positionToleranceMetres: Double = 0.01, velocityToleranceMetresPerSecond: Double = 1e-5,
        lambertOptions: LambertOptions = .default, propagationOptions: TwoBodyPropagationOptions = .default) throws -> LambertTransferPlan {
        guard departureState.frame === arrivalTargetState.frame else { throw SpacecraftDynamicsError.frameMismatch }
        guard arrivalTargetState.epoch > departureState.epoch else { throw LambertTransferPlanningError.invalidEpochs }
        guard [positionToleranceMetres,velocityToleranceMetresPerSecond].allSatisfy({ $0.isFinite && $0 > 0 }) else {
            throw LambertTransferPlanningError.invalidTolerance
        }
        let solution = try LambertSolver.solve(departurePosition: departureState.position,
            arrivalPosition: arrivalTargetState.position, departureEpoch: departureState.epoch,
            arrivalEpoch: arrivalTargetState.epoch, frame: departureState.frame, centralBody: centralBody,
            gravitationalParameter: gravitationalParameter, referencePole: referencePole, direction: direction,
            revolutions: revolutions, multiRevolutionBranch: multiRevolutionBranch, options: lambertOptions)
        let departureRelative = solution.departureState.velocity-departureState.velocity
        let departureBurn = try ImpulsiveManeuver(epoch: departureState.epoch, frame: departureState.frame,
            deltaVelocityMetresPerSecond: departureRelative)
        let coast = try solution.trajectory(options: propagationOptions)
        let before = try coast.state(at: arrivalTargetState.epoch)
        let incoming = before.velocity-arrivalTargetState.velocity
        let arrivalBurn: ImpulsiveManeuver?
        if arrivalCondition == .velocityMatchedRendezvous {
            // Correct the propagated velocity, preserving the propagated position.
            arrivalBurn = try ImpulsiveManeuver(epoch: arrivalTargetState.epoch, frame: departureState.frame,
                deltaVelocityMetresPerSecond: -incoming)
        } else { arrivalBurn = nil }
        let after = try CartesianState(epoch: before.epoch, frame: before.frame, position: before.position,
            velocity: before.velocity + (arrivalBurn?.deltaVelocityMetresPerSecond ?? .zero), acceleration: before.acceleration)
        let trajectory = LambertTransferTrajectory(coast: coast, departureReferenceState: departureState, postArrivalState: after)
        let positionError = routeNorm(after.position-arrivalTargetState.position)
        let outgoing = after.velocity-arrivalTargetState.velocity
        let velocityError = routeNorm(outgoing)
        let cost = routeNorm(departureRelative) + (arrivalBurn.map { routeNorm($0.deltaVelocityMetresPerSecond) } ?? 0)
        guard [positionError,velocityError,cost].allSatisfy(\.isFinite) else {
            throw LambertTransferPlanningError.numericalRangeExceeded
        }
        return LambertTransferPlan(solution: solution, departureReferenceState: departureState,
            arrivalTargetState: arrivalTargetState, arrivalCondition: arrivalCondition,
            departureBurn: departureBurn, arrivalBurn: arrivalBurn, preArrivalState: before, postArrivalState: after,
            departureRelativeVelocityMetresPerSecond: departureRelative,
            incomingRelativeVelocityMetresPerSecond: incoming, outgoingRelativeVelocityMetresPerSecond: outgoing,
            totalDeltaVelocityMetresPerSecond: cost, positionErrorMetres: positionError,
            velocityErrorMetresPerSecond: velocityError,
            arrivalConditionSatisfied: positionError <= positionToleranceMetres &&
                (arrivalCondition != .velocityMatchedRendezvous || velocityError <= velocityToleranceMetresPerSecond),
            trajectory: trajectory)
    }

    /// Samples departure only at departureEpoch and the moving target at arrivalEpoch.
    /// Provider coverage, exact epochs, frame identities and physical origins are checked.
    public static func plan(departure: MovingWaypointRouteTarget, arrival: MovingWaypointRouteTarget,
        departureEpoch: Epoch, arrivalEpoch: Epoch, centralBody: BodyIdentifier,
        gravitationalParameter: GravitationalParameter, referencePole: SIMD3<Double>, direction: LambertDirection,
        arrivalCondition: LambertArrivalCondition, revolutions: Int = 0,
        multiRevolutionBranch: LambertMultiRevolutionBranch? = nil,
        positionToleranceMetres: Double = 0.01, velocityToleranceMetresPerSecond: Double = 1e-5,
        lambertOptions: LambertOptions = .default, propagationOptions: TwoBodyPropagationOptions = .default) throws -> LambertTransferPlan {
        guard departure.origin == centralBody, arrival.origin == centralBody else { throw LambertTransferPlanningError.originMismatch }
        guard departure.frame === arrival.frame else { throw SpacecraftDynamicsError.frameMismatch }
        guard arrivalEpoch > departureEpoch else { throw LambertTransferPlanningError.invalidEpochs }
        return try plan(departureState: departure.state(at: departureEpoch), arrivalTargetState: arrival.state(at: arrivalEpoch),
            centralBody: centralBody, gravitationalParameter: gravitationalParameter, referencePole: referencePole,
            direction: direction, arrivalCondition: arrivalCondition, revolutions: revolutions,
            multiRevolutionBranch: multiRevolutionBranch, positionToleranceMetres: positionToleranceMetres,
            velocityToleranceMetresPerSecond: velocityToleranceMetresPerSecond,
            lambertOptions: lambertOptions, propagationOptions: propagationOptions)
    }

    public static func plan(departure: any EphemerisProvider, arrival: any EphemerisProvider,
        departureEpoch: Epoch, arrivalEpoch: Epoch, centralBody: BodyIdentifier,
        gravitationalParameter: GravitationalParameter, referencePole: SIMD3<Double>, direction: LambertDirection,
        arrivalCondition: LambertArrivalCondition, revolutions: Int = 0,
        multiRevolutionBranch: LambertMultiRevolutionBranch? = nil,
        positionToleranceMetres: Double = 0.01, velocityToleranceMetresPerSecond: Double = 1e-5,
        lambertOptions: LambertOptions = .default, propagationOptions: TwoBodyPropagationOptions = .default) throws -> LambertTransferPlan {
        try plan(departure: MovingWaypointRouteTarget(ephemeris: departure), arrival: MovingWaypointRouteTarget(ephemeris: arrival),
            departureEpoch: departureEpoch, arrivalEpoch: arrivalEpoch, centralBody: centralBody,
            gravitationalParameter: gravitationalParameter, referencePole: referencePole, direction: direction,
            arrivalCondition: arrivalCondition, revolutions: revolutions, multiRevolutionBranch: multiRevolutionBranch,
            positionToleranceMetres: positionToleranceMetres, velocityToleranceMetresPerSecond: velocityToleranceMetresPerSecond,
            lambertOptions: lambertOptions, propagationOptions: propagationOptions)
    }
}
