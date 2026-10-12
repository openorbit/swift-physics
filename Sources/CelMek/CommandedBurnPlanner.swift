// SPDX-License-Identifier: Apache-2.0
import Math

public enum CommandedBurnPlanningError: Error, Equatable {
    case invalidCandidates, invalidPassages, invalidTolerance, invalidBudget
}
public struct CommandedBurnCandidate {
    public let arrivalEpoch: Epoch
    public let schedule: ScheduledThrust
    public init(arrivalEpoch: Epoch, schedule: ScheduledThrust) {
        self.arrivalEpoch = arrivalEpoch; self.schedule = schedule
    }
}
/// A physical position constraint independent of ignition and cutoff epochs.
public struct CommandedBurnPassage {
    public let epoch: Epoch
    public let target: MovingWaypointRouteTarget
    public let matchVelocity: Bool
    public init(epoch: Epoch, target: MovingWaypointRouteTarget, matchVelocity: Bool = false) {
        self.epoch = epoch; self.target = target; self.matchVelocity = matchVelocity
    }
}
public enum CommandedBurnCandidateFailure: Equatable {
    case scheduleOutsideFlight, targetOutsideCoverage, passageOutsideFlight(index: Int)
    case passageOutsideCoverage(index: Int), passageResidual(index: Int)
    case arrivalResidual, propellantReserve
    case integration(NumericalPropagationError)
}
public struct CommandedBurnAttempt {
    public let checkpoints: [CommandedBurnCheckpoint]
    public let candidateIndex: Int
    public let failure: CommandedBurnCandidateFailure
}
public struct CommandedBurnCheckpoint {
    public let state: SpacecraftState
    public let positionErrorMetres: Double
    public let velocityErrorMetresPerSecond: Double
}
public enum CommandedBurnSearchStatus: Equatable, Sendable {
    case feasible, candidateSetExhausted, evaluationLimitReached
}
public struct CommandedBurnSearchResult {
    public let status: CommandedBurnSearchStatus
    public let evaluatedCandidates: Int
    public let candidateIndex: Int?
    public let candidate: CommandedBurnCandidate?
    /// Passage checkpoints followed by the rendezvous checkpoint, with carried mass.
    public let checkpoints: [CommandedBurnCheckpoint]
    public let failedAttempts: [CommandedBurnAttempt]
}

/// Bounded shooting over explicit inertial burn commands. It does not compensate
/// gravity to follow a prescribed curve: the supplied commands drive propagation.
/// The earliest feasible supplied candidate wins; no continuous optimum is claimed.
public enum CommandedBurnPlanner {
    public static func search(initialState: SpacecraftState, context: ForceEvaluationContext,
        backgroundModels: [any ForceModel] = [], candidates: [CommandedBurnCandidate],
        arrivalTarget: MovingWaypointRouteTarget, passages: [CommandedBurnPassage] = [],
        maximumCandidateEvaluations: Int = 1000, minimumRemainingPropellantKilograms: Double = 0,
        positionToleranceMetres: Double = 0.01, velocityToleranceMetresPerSecond: Double = 1e-5,
        numericalOptions: NumericalPropagationOptions = try! NumericalPropagationOptions()) throws -> CommandedBurnSearchResult {
        guard !candidates.isEmpty, candidates.count <= 10_000 else { throw CommandedBurnPlanningError.invalidCandidates }
        guard (1...10_000).contains(maximumCandidateEvaluations) else { throw CommandedBurnPlanningError.invalidBudget }
        guard passages.count <= 10_000,
              zip(passages, passages.dropFirst()).allSatisfy({ $0.epoch < $1.epoch }),
              passages.allSatisfy({ $0.epoch > initialState.epoch }) else {
            throw CommandedBurnPlanningError.invalidPassages
        }
        guard [positionToleranceMetres, velocityToleranceMetresPerSecond].allSatisfy({ $0.isFinite && $0 > 0 }),
              minimumRemainingPropellantKilograms.isFinite, minimumRemainingPropellantKilograms >= 0 else {
            throw CommandedBurnPlanningError.invalidTolerance
        }
        guard context.frame === initialState.frame,
              candidates.allSatisfy({ $0.schedule.frame === initialState.frame }),
              arrivalTarget.frame === initialState.frame,
              passages.allSatisfy({ $0.target.frame === initialState.frame }) else {
            throw SpacecraftDynamicsError.frameMismatch
        }
        guard arrivalTarget.origin == context.origin, passages.allSatisfy({ $0.target.origin == context.origin }) else {
            throw MovingWaypointTargetError.originMismatch
        }
        _ = try SpacecraftDynamics(models: backgroundModels, context: context)
        var failures: [CommandedBurnAttempt] = [], evaluated = 0
        let ordered = candidates.indices.sorted {
            candidates[$0].arrivalEpoch == candidates[$1].arrivalEpoch ? $0 < $1 :
                candidates[$0].arrivalEpoch < candidates[$1].arrivalEpoch
        }
        func finish(_ status: CommandedBurnSearchStatus) -> CommandedBurnSearchResult {
            CommandedBurnSearchResult(status: status, evaluatedCandidates: evaluated, candidateIndex: nil,
                candidate: nil, checkpoints: [], failedAttempts: failures)
        }
        for index in ordered {
            if evaluated == maximumCandidateEvaluations { return finish(.evaluationLimitReached) }
            evaluated += 1
            let candidate = candidates[index]
            var failure: CommandedBurnCandidateFailure?
            if candidate.arrivalEpoch <= initialState.epoch || candidate.schedule.burns.contains(where: {
                $0.start < initialState.epoch || $0.end > candidate.arrivalEpoch
            }) { failure = .scheduleOutsideFlight }
            else if !arrivalTarget.coverage.contains(candidate.arrivalEpoch) { failure = .targetOutsideCoverage }
            for (i, passage) in passages.enumerated() where failure == nil {
                if passage.epoch >= candidate.arrivalEpoch { failure = .passageOutsideFlight(index: i) }
                else if !passage.target.coverage.contains(passage.epoch) { failure = .passageOutsideCoverage(index: i) }
            }
            if let failure {
                failures.append(CommandedBurnAttempt(checkpoints: [], candidateIndex: index, failure: failure)); continue
            }
            var state = initialState, checkpoints: [CommandedBurnCheckpoint] = []
            do {
                let goals = passages + [CommandedBurnPassage(epoch: candidate.arrivalEpoch,
                    target: arrivalTarget, matchVelocity: true)]
                for (i, goal) in goals.enumerated() {
                    let propagator = try ScheduledThrustPropagator(initialState: state, schedule: candidate.schedule,
                        backgroundModels: backgroundModels, context: context, options: numericalOptions)
                    state = try propagator.spacecraftState(at: goal.epoch)
                    let target = try goal.target.state(at: goal.epoch)
                    let positionError = routeNorm(state.cartesian.position-target.position)
                    let velocityError = routeNorm(state.cartesian.velocity-target.velocity)
                    checkpoints.append(CommandedBurnCheckpoint(state: state, positionErrorMetres: positionError,
                        velocityErrorMetresPerSecond: velocityError))
                    if positionError > positionToleranceMetres ||
                        (goal.matchVelocity && velocityError > velocityToleranceMetresPerSecond) {
                        failure = i == passages.count ? .arrivalResidual : .passageResidual(index: i); break
                    }
                    if state.mass.propellantMassKilograms < minimumRemainingPropellantKilograms {
                        failure = .propellantReserve; break
                    }
                }
            } catch let error as NumericalPropagationError { failure = .integration(error) }
            if let failure { failures.append(CommandedBurnAttempt(checkpoints: checkpoints, candidateIndex: index, failure: failure)); continue }
            return CommandedBurnSearchResult(status: .feasible, evaluatedCandidates: evaluated, candidateIndex: index,
                candidate: candidate, checkpoints: checkpoints, failedAttempts: failures)
        }
        return finish(.candidateSetExhausted)
    }
}
