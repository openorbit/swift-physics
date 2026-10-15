// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum CircularHohmannStateError: Error, Equatable {
    case invalidOptions, nonCircularSpeed, nonTangentialVelocity, epochMismatch, targetOutsideCoverage
}
public struct CircularHohmannValidationOptions: Equatable, Sendable {
    public let relativeSpeedTolerance, maximumRadialSpeedFraction, maximumPlaneNormalDifference: Double
    public init(relativeSpeedTolerance: Double = 1e-8, maximumRadialSpeedFraction: Double = 1e-8,
                maximumPlaneNormalDifference: Double = 1e-10) throws {
        guard [relativeSpeedTolerance, maximumRadialSpeedFraction, maximumPlaneNormalDifference]
            .allSatisfy({ $0.isFinite && $0 >= 0 && $0 < 1 }) else {
            throw CircularHohmannStateError.invalidOptions
        }
        self.relativeSpeedTolerance = relativeSpeedTolerance
        self.maximumRadialSpeedFraction = maximumRadialSpeedFraction
        self.maximumPlaneNormalDifference = maximumPlaneNormalDifference
    }
}
extension CircularHohmannOrbit {
    /// Validate a circular approximation and project velocity onto its exact
    /// circular tangent. This initializer does not preserve a residual radial speed.
    public init(referenceState: CartesianState, centralBody: BodyIdentifier,
                gravitationalParameter: GravitationalParameter,
                validationOptions: CircularHohmannValidationOptions = try! CircularHohmannValidationOptions()) throws {
        let radius = checkedHohmannNorm(referenceState.position)
        guard radius.isFinite, radius > 0 else { throw HohmannTransferError.invalidRadius }
        let speed = try HohmannTransfer(departureRadiusMetres: radius, arrivalRadiusMetres: radius,
            gravitationalParameter: gravitationalParameter).departureCircularSpeedMetresPerSecond
        let radial = referenceState.position/radius
        let radialSpeed = (referenceState.velocity*radial).sum()
        guard abs(radialSpeed)/speed <= validationOptions.maximumRadialSpeedFraction else {
            throw CircularHohmannStateError.nonTangentialVelocity
        }
        guard abs(checkedHohmannNorm(referenceState.velocity)/speed-1) <= validationOptions.relativeSpeedTolerance else {
            throw CircularHohmannStateError.nonCircularSpeed
        }
        try self.init(referenceEpoch: referenceState.epoch, frame: referenceState.frame, centralBody: centralBody,
            radiusMetres: radius, radialDirection: radial,
            progradeDirection: referenceState.velocity-radial*radialSpeed,
            gravitationalParameter: gravitationalParameter)
    }
    public init(provider: MovingWaypointRouteTarget, at epoch: Epoch, centralBody: BodyIdentifier,
                gravitationalParameter: GravitationalParameter,
                validationOptions: CircularHohmannValidationOptions = try! CircularHohmannValidationOptions()) throws {
        guard provider.origin == centralBody else { throw HohmannRendezvousError.incompatibleCentralBody }
        try self.init(referenceState: provider.state(at: epoch), centralBody: centralBody,
            gravitationalParameter: gravitationalParameter, validationOptions: validationOptions)
    }
}
public struct CheckedHohmannTransferPlan {
    public let plan: HohmannTransferPlan
    public let departureCircularModel: CircularHohmannOrbit
    public let targetCircularModel: CircularHohmannOrbit
    public let departureVelocityProjectionErrorMetresPerSecond: Double
    public let targetVelocityProjectionErrorMetresPerSecond: Double
    public let targetArrivalState: CartesianState
    public let positionErrorMetres, velocityErrorMetresPerSecond: Double
    /// Actual impulse norms, including corrections for accepted state approximations.
    public let totalDeltaVelocityMetresPerSecond: Double
    public let rendezvousSatisfied: Bool
}

/// Plans from sampled states, preserving actual departure position/velocity.
/// Projection tolerances describe accepted input approximations, not rendezvous accuracy.
public enum CheckedHohmannTransferPlanner {
    public static func plan(departureState: CartesianState, targetState: CartesianState,
        centralBody: BodyIdentifier, gravitationalParameter: GravitationalParameter,
        validationOptions: CircularHohmannValidationOptions = try! CircularHohmannValidationOptions(),
        positionToleranceMetres: Double = 0.01, velocityToleranceMetresPerSecond: Double = 1e-5,
        propagationOptions: TwoBodyPropagationOptions = try! TwoBodyPropagationOptions()) throws -> CheckedHohmannTransferPlan {
        try build(departureState: departureState, targetState: targetState, centralBody: centralBody,
            gravitationalParameter: gravitationalParameter, validationOptions: validationOptions,
            positionToleranceMetres: positionToleranceMetres,
            velocityToleranceMetresPerSecond: velocityToleranceMetresPerSecond,
            propagationOptions: propagationOptions, actualTarget: nil)
    }
    static func build(departureState: CartesianState, targetState: CartesianState,
        centralBody: BodyIdentifier, gravitationalParameter: GravitationalParameter,
        validationOptions: CircularHohmannValidationOptions,
        positionToleranceMetres: Double, velocityToleranceMetresPerSecond: Double,
        propagationOptions: TwoBodyPropagationOptions, actualTarget: MovingWaypointRouteTarget?) throws -> CheckedHohmannTransferPlan {
        guard departureState.epoch == targetState.epoch else { throw CircularHohmannStateError.epochMismatch }
        guard departureState.frame === targetState.frame else { throw SpacecraftDynamicsError.frameMismatch }
        guard [positionToleranceMetres, velocityToleranceMetresPerSecond].allSatisfy({ $0.isFinite && $0 > 0 }) else {
            throw HohmannRendezvousError.invalidTolerance
        }
        let departure = try CircularHohmannOrbit(referenceState: departureState, centralBody: centralBody,
            gravitationalParameter: gravitationalParameter, validationOptions: validationOptions)
        let target = try CircularHohmannOrbit(referenceState: targetState, centralBody: centralBody,
            gravitationalParameter: gravitationalParameter, validationOptions: validationOptions)
        try checkPlane(departure, target, options: validationOptions)
        let transfer = try HohmannTransfer(departureRadiusMetres: departure.radiusMetres,
            arrivalRadiusMetres: target.radiusMetres, gravitationalParameter: gravitationalParameter)
        let axes = try departure.directions(at: departureState.epoch)
        let ideal = try transfer.plan(departureEpoch: departureState.epoch, frame: departureState.frame,
            centralBody: centralBody, radialDirection: axes.radial, progradeDirection: axes.tangent,
            options: propagationOptions)
        if let actualTarget, !actualTarget.coverage.contains(ideal.arrivalEpoch) {
            throw CircularHohmannStateError.targetOutsideCoverage
        }
        let arrivalTarget = try actualTarget?.state(at: ideal.arrivalEpoch) ?? target.state(at: ideal.arrivalEpoch)
        let arrivalModel = try CircularHohmannOrbit(referenceState: arrivalTarget, centralBody: centralBody,
            gravitationalParameter: gravitationalParameter, validationOptions: validationOptions)
        try checkPlane(departure, arrivalModel, options: validationOptions)
        let desiredDeparture = try ideal.trajectory.state(at: departureState.epoch).velocity
        let departureBurn = try ImpulsiveManeuver(epoch: departureState.epoch, frame: departureState.frame,
            deltaVelocityMetresPerSecond: desiredDeparture-departureState.velocity)
        let base = try TwoBodyPropagator(initialState: departureState, centralBody: centralBody,
            gravitationalParameter: gravitationalParameter, options: propagationOptions)
        let coast = try ManeuveredTwoBodyPropagator(base: base, maneuvers: [departureBurn])
        let preArrival = try coast.state(at: ideal.arrivalEpoch)
        // Circularize in the transfer plane. Do not "match" a remote target's
        // velocity when the spacecraft has missed its physical position.
        let desiredArrival = -axes.tangent*transfer.arrivalCircularSpeedMetresPerSecond
        let arrivalBurn = try ImpulsiveManeuver(epoch: ideal.arrivalEpoch, frame: departureState.frame,
            deltaVelocityMetresPerSecond: desiredArrival-preArrival.velocity)
        let trajectory = try ManeuveredTwoBodyPropagator(base: base, maneuvers: [departureBurn, arrivalBurn])
        let arrival = try trajectory.state(at: ideal.arrivalEpoch)
        let p = checkedHohmannNorm(arrival.position-arrivalTarget.position)
        let v = checkedHohmannNorm(arrival.velocity-arrivalTarget.velocity)
        let departureProjection = try departure.state(at: departureState.epoch)
        let targetProjection = try target.state(at: targetState.epoch)
        return CheckedHohmannTransferPlan(plan: HohmannTransferPlan(transfer: transfer,
                arrivalEpoch: ideal.arrivalEpoch, departureBurn: departureBurn, arrivalBurn: arrivalBurn, trajectory: trajectory),
            departureCircularModel: departure, targetCircularModel: target,
            departureVelocityProjectionErrorMetresPerSecond: checkedHohmannNorm(departureProjection.velocity-departureState.velocity),
            targetVelocityProjectionErrorMetresPerSecond: checkedHohmannNorm(targetProjection.velocity-targetState.velocity),
            targetArrivalState: arrivalTarget, positionErrorMetres: p, velocityErrorMetresPerSecond: v,
            totalDeltaVelocityMetresPerSecond: checkedHohmannNorm(departureBurn.deltaVelocityMetresPerSecond) +
                checkedHohmannNorm(arrivalBurn.deltaVelocityMetresPerSecond),
            rendezvousSatisfied: p <= positionToleranceMetres && v <= velocityToleranceMetresPerSecond)
    }
    private static func checkPlane(_ departure: CircularHohmannOrbit, _ target: CircularHohmannOrbit,
                                   options: CircularHohmannValidationOptions) throws {
        if (departure.normal*target.normal).sum() < 0 { throw HohmannRendezvousError.opposingOrbitalDirections }
        guard checkedHohmannNorm(departure.normal-target.normal) <= options.maximumPlaneNormalDifference else {
            throw HohmannRendezvousError.nonCoplanarOrbits
        }
    }
}
public enum CheckedHohmannCandidateFailure: Equatable {
    case departureOutsideCoverage, targetOutsideCoverage, rendezvousResidual, sameOrbitPhasingRequired
    case state(CircularHohmannStateError), geometry(HohmannRendezvousError)
}
public struct CheckedHohmannAttempt {
    public let departureEpoch: Epoch
    public let failure: CheckedHohmannCandidateFailure
    public let plan: CheckedHohmannTransferPlan?
}
public enum CheckedHohmannSearchStatus: Equatable, Sendable {
    case feasible, candidateSetExhausted, evaluationLimitReached, alreadyRendezvoused
}
public struct CheckedHohmannSearchResult {
    public let rendezvousState: CartesianState?
    public let status: CheckedHohmannSearchStatus
    public let evaluatedCandidates: Int
    public let plan: CheckedHohmannTransferPlan?
    public let failedAttempts: [CheckedHohmannAttempt]
}
extension CheckedHohmannTransferPlanner {
    /// Bounded search of explicit dates. Samples both providers afresh at each
    /// departure and checks the target's actual state at arrival, not its seed fit.
    public static func search(departureProvider: MovingWaypointRouteTarget, targetProvider: MovingWaypointRouteTarget,
        centralBody: BodyIdentifier, gravitationalParameter: GravitationalParameter, candidateDepartureEpochs: [Epoch],
        maximumCandidateEvaluations: Int = 1000,
        validationOptions: CircularHohmannValidationOptions = try! CircularHohmannValidationOptions(),
        positionToleranceMetres: Double = 0.01, velocityToleranceMetresPerSecond: Double = 1e-5,
        propagationOptions: TwoBodyPropagationOptions = try! TwoBodyPropagationOptions()) throws -> CheckedHohmannSearchResult {
        guard !candidateDepartureEpochs.isEmpty, candidateDepartureEpochs.count <= 100_000 else {
            throw WaypointTimingSearchError.invalidCandidates
        }
        guard (1...10_000).contains(maximumCandidateEvaluations) else { throw HohmannRendezvousError.invalidBudget }
        guard departureProvider.frame === targetProvider.frame else { throw SpacecraftDynamicsError.frameMismatch }
        guard departureProvider.origin == centralBody, targetProvider.origin == centralBody else {
            throw HohmannRendezvousError.incompatibleCentralBody
        }
        guard [positionToleranceMetres, velocityToleranceMetresPerSecond].allSatisfy({ $0.isFinite && $0 > 0 }) else {
            throw HohmannRendezvousError.invalidTolerance
        }
        // Reject unsupported GM/frame configuration even if every date misses coverage.
        guard gravitationalParameter.convention == .newtonian || gravitationalParameter.convention == .ttCompatible else {
            throw PropagationError.incompatibleGMConvention
        }
        guard departureProvider.frame.pseudoInertial else { throw PropagationError.nonInertialFrame }
        var failures: [CheckedHohmannAttempt] = [], evaluated = 0
        func finish(_ status: CheckedHohmannSearchStatus, plan: CheckedHohmannTransferPlan? = nil, state: CartesianState? = nil) -> CheckedHohmannSearchResult {
            CheckedHohmannSearchResult(rendezvousState: state ?? plan?.targetArrivalState, status: status, evaluatedCandidates: evaluated, plan: plan, failedAttempts: failures)
        }
        let epochs = Array(Set(candidateDepartureEpochs)).sorted()
        for epoch in epochs {
            if evaluated == maximumCandidateEvaluations { return finish(.evaluationLimitReached) }
            evaluated += 1
            var failure: CheckedHohmannCandidateFailure?
            var plan: CheckedHohmannTransferPlan?
            if !departureProvider.coverage.contains(epoch) { failure = .departureOutsideCoverage }
            else if !targetProvider.coverage.contains(epoch) { failure = .targetOutsideCoverage }
            else {
                do {
                    let departureState = try departureProvider.state(at: epoch)
                    let targetState = try targetProvider.state(at: epoch)
                    let departureModel = try CircularHohmannOrbit(referenceState: departureState, centralBody: centralBody,
                        gravitationalParameter: gravitationalParameter, validationOptions: validationOptions)
                    let targetModel = try CircularHohmannOrbit(referenceState: targetState, centralBody: centralBody,
                        gravitationalParameter: gravitationalParameter, validationOptions: validationOptions)
                    try checkPlane(departureModel, targetModel, options: validationOptions)
                    if departureModel.radiusMetres == targetModel.radiusMetres {
                        if checkedHohmannNorm(departureState.position-targetState.position) <= positionToleranceMetres &&
                            checkedHohmannNorm(departureState.velocity-targetState.velocity) <= velocityToleranceMetresPerSecond {
                            return finish(.alreadyRendezvoused, state: targetState)
                        }
                        failures.append(CheckedHohmannAttempt(departureEpoch: epoch, failure: .sameOrbitPhasingRequired, plan: nil))
                        continue
                    }
                    plan = try build(departureState: departureState, targetState: targetState, centralBody: centralBody,
                        gravitationalParameter: gravitationalParameter, validationOptions: validationOptions,
                        positionToleranceMetres: positionToleranceMetres,
                        velocityToleranceMetresPerSecond: velocityToleranceMetresPerSecond,
                        propagationOptions: propagationOptions, actualTarget: targetProvider)
                    if plan!.rendezvousSatisfied { return finish(.feasible, plan: plan) }
                    failure = .rendezvousResidual
                } catch CircularHohmannStateError.targetOutsideCoverage { failure = .targetOutsideCoverage }
                catch let error as CircularHohmannStateError { failure = .state(error) }
                catch let error as HohmannRendezvousError { failure = .geometry(error) }
            }
            failures.append(CheckedHohmannAttempt(departureEpoch: epoch, failure: failure!, plan: plan))
        }
        return finish(.candidateSetExhausted)
    }
}
private func checkedHohmannNorm(_ v: SIMD3<Double>) -> Double { hypot(hypot(v.x,v.y),v.z) }
