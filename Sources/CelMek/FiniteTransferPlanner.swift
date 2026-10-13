// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum FiniteTransferPlanningError: Error, Equatable {
    case invalidInput, insufficientPropellant, burnsDoNotFit
}
public enum FiniteTransferRefinementStatus: Equatable {
    case converged, evaluationLimitReached, stepLimitReached
}
public struct FiniteTransferEvaluation {
    public let schedule: ScheduledThrust
    public let arrivalState: SpacecraftState
    public let positionErrorMetres: Double
    public let velocityErrorMetresPerSecond: Double
    public let propellantUsedKilograms: Double
}
public struct FiniteTransferRefinementResult {
    public let status: FiniteTransferRefinementStatus
    public let evaluations: Int
    public let rejectedTrials: Int
    public let initial: FiniteTransferEvaluation
    public let best: FiniteTransferEvaluation
}

/// Constant inertial engine commands, followed by bounded local coordinate shooting.
/// The caller supplies the physical origin and gravity/background models explicitly.
public enum FiniteTransferPlanner {
    /// Rocket-equation durations with sequential mass depletion. Burns are centered
    /// on impulses where possible and shifted inside the flight at its endpoints.
    /// This is a seed, not an equivalent realization or a rendezvous certificate.
    public static func initialGuess(initialState: SpacecraftState, arrivalEpoch: Epoch,
        maneuvers: [ImpulsiveManeuver], thrustNewtons: Double,
        effectiveExhaustVelocityMetresPerSecond: Double) throws -> ScheduledThrust {
        guard arrivalEpoch > initialState.epoch, maneuvers.count <= 100,
              thrustNewtons.isFinite, thrustNewtons > 0,
              effectiveExhaustVelocityMetresPerSecond.isFinite, effectiveExhaustVelocityMetresPerSecond > 0 else {
            throw FiniteTransferPlanningError.invalidInput
        }
        let ordered = maneuvers.sorted { $0.epoch < $1.epoch }
        guard ordered.allSatisfy({ $0.frame === initialState.frame }) else { throw SpacecraftDynamicsError.frameMismatch }
        guard ordered.allSatisfy({ $0.epoch >= initialState.epoch && $0.epoch <= arrivalEpoch }),
              zip(ordered, ordered.dropFirst()).allSatisfy({ $0.epoch < $1.epoch }) else {
            throw FiniteTransferPlanningError.invalidInput
        }
        var mass = initialState.mass.totalMassKilograms
        var burns: [FiniteBurn] = []
        for impulse in ordered {
            let speed = routeNorm(impulse.deltaVelocityMetresPerSecond)
            if speed == 0 { continue }
            let consumed = mass * -expm1(-speed / effectiveExhaustVelocityMetresPerSecond)
            let seconds = consumed * effectiveExhaustVelocityMetresPerSecond / thrustNewtons
            guard seconds.isFinite, seconds > 0 else { throw FiniteTransferPlanningError.invalidInput }
            mass -= consumed
            guard mass > initialState.mass.dryMassKilograms else { throw FiniteTransferPlanningError.insufficientPropellant }
            let duration = try Duration(seconds: seconds)
            var start = try impulse.epoch.advanced(by: Duration(seconds: -seconds / 2))
            start = max(start, initialState.epoch)
            var end = try start.advanced(by: duration)
            if end > arrivalEpoch {
                end = arrivalEpoch
                start = try end.advanced(by: Duration(seconds: -seconds))
            }
            guard start >= initialState.epoch, end > start,
                  burns.last.map({ $0.end <= start }) ?? true else { throw FiniteTransferPlanningError.burnsDoNotFit }
            burns.append(try FiniteBurn(start: start, end: end, thrust: FiniteThrust(frame: initialState.frame,
                thrustNewtons: thrustNewtons, effectiveExhaustVelocityMetresPerSecond: effectiveExhaustVelocityMetresPerSecond,
                direction: impulse.deltaVelocityMetresPerSecond)))
        }
        return try ScheduledThrust(frame: initialState.frame, burns: burns)
    }

    /// Fixed arrival epoch/state; thrust and exhaust velocity remain fixed per burn.
    /// Varies ignition, duration and inertial direction, preserving burn order and
    /// non-overlap. No global optimum, collision clearance or attitude feasibility is implied.
    public static func refine(initialState: SpacecraftState, initialSchedule: ScheduledThrust,
        target: CartesianState, context: ForceEvaluationContext, backgroundModels: [any ForceModel] = [],
        positionToleranceMetres: Double = 1, velocityToleranceMetresPerSecond: Double = 0.001,
        minimumRemainingPropellantKilograms: Double = 0, maximumEvaluations: Int = 500,
        numericalOptions: NumericalPropagationOptions = try! NumericalPropagationOptions()) throws -> FiniteTransferRefinementResult {
        guard target.epoch > initialState.epoch, !initialSchedule.burns.isEmpty, initialSchedule.burns.count <= 100,
              (1...10_000).contains(maximumEvaluations),
              [positionToleranceMetres, velocityToleranceMetresPerSecond].allSatisfy({ $0.isFinite && $0 > 0 }),
              minimumRemainingPropellantKilograms.isFinite, minimumRemainingPropellantKilograms >= 0 else {
            throw FiniteTransferPlanningError.invalidInput
        }
        guard target.frame === initialState.frame, initialSchedule.frame === initialState.frame,
              context.frame === initialState.frame else { throw SpacecraftDynamicsError.frameMismatch }
        func valid(_ burns: [FiniteBurn]) -> Bool {
            burns.allSatisfy { $0.start >= initialState.epoch && $0.end <= target.epoch } &&
            zip(burns, burns.dropFirst()).allSatisfy { $0.end <= $1.start }
        }
        guard valid(initialSchedule.burns) else { throw FiniteTransferPlanningError.burnsDoNotFit }
        func evaluate(_ schedule: ScheduledThrust) throws -> FiniteTransferEvaluation {
            let consumed = try schedule.burns.reduce(0.0) {
                try $0 + $1.end.duration(since: $1.start).seconds * $1.thrust.propellantConsumptionKilogramsPerSecond
            }
            guard consumed.isFinite,
                  initialState.mass.propellantMassKilograms - consumed > 0,
                  initialState.mass.propellantMassKilograms - consumed >= minimumRemainingPropellantKilograms else {
                throw FiniteTransferPlanningError.insufficientPropellant
            }
            let state = try ScheduledThrustPropagator(initialState: initialState, schedule: schedule,
                backgroundModels: backgroundModels, context: context, options: numericalOptions).spacecraftState(at: target.epoch)
            guard state.mass.propellantMassKilograms >= minimumRemainingPropellantKilograms else {
                throw FiniteTransferPlanningError.insufficientPropellant
            }
            return FiniteTransferEvaluation(schedule: schedule, arrivalState: state,
                positionErrorMetres: routeNorm(state.cartesian.position - target.position),
                velocityErrorMetresPerSecond: routeNorm(state.cartesian.velocity - target.velocity),
                propellantUsedKilograms: initialState.mass.propellantMassKilograms - state.mass.propellantMassKilograms)
        }
        func score(_ e: FiniteTransferEvaluation) -> Double {
            hypot(e.positionErrorMetres / positionToleranceMetres, e.velocityErrorMetresPerSecond / velocityToleranceMetresPerSecond)
        }
        func converged(_ e: FiniteTransferEvaluation) -> Bool {
            e.positionErrorMetres <= positionToleranceMetres && e.velocityErrorMetresPerSecond <= velocityToleranceMetresPerSecond
        }
        let initial = try evaluate(initialSchedule)
        var best = initial, evaluations = 1, rejected = 0, step = 0.25
        let flight = try target.epoch.duration(since: initialState.epoch).seconds
        while !converged(best) && evaluations < maximumEvaluations && step >= 1e-6 {
            var improved = false
            for i in best.schedule.burns.indices {
                for coordinate in 0..<5 {
                    for sign in [-1.0, 1.0] {
                        if evaluations >= maximumEvaluations { break }
                        var burns = best.schedule.burns
                        let b = burns[i], duration = try b.end.duration(since: b.start).seconds
                        var start = b.start, seconds = duration, direction = b.thrust.direction
                        switch coordinate {
                        case 0: start = try start.advanced(by: Duration(seconds: sign * step * min(flight * 0.01, duration)))
                        case 1: seconds *= 1 + sign * step
                        default: direction[coordinate - 2] += sign * step
                        }
                        let end = try start.advanced(by: Duration(seconds: seconds))
                        guard end > start else { rejected += 1; continue }
                        burns[i] = try FiniteBurn(start: start, end: end, thrust: FiniteThrust(frame: b.thrust.frame,
                            thrustNewtons: b.thrust.thrustNewtons,
                            effectiveExhaustVelocityMetresPerSecond: b.thrust.effectiveExhaustVelocityMetresPerSecond,
                            direction: direction))
                        guard valid(burns) else { rejected += 1; continue }
                        evaluations += 1
                        do {
                            let trial = try evaluate(ScheduledThrust(frame: initialState.frame, burns: burns))
                            if score(trial) < score(best) { best = trial; improved = true }
                        } catch is NumericalPropagationError { rejected += 1 }
                        catch FiniteTransferPlanningError.insufficientPropellant { rejected += 1 }
                        if converged(best) { break }
                    }
                    if converged(best) { break }
                }
                if converged(best) { break }
            }
            if !improved { step *= 0.5 }
        }
        return FiniteTransferRefinementResult(status: converged(best) ? .converged :
            (evaluations >= maximumEvaluations ? .evaluationLimitReached : .stepLimitReached),
            evaluations: evaluations, rejectedTrials: rejected, initial: initial, best: best)
    }
}
