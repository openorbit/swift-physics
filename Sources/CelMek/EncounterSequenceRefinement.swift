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
import Foundation
import Math

public enum EncounterSequenceRefinementStatus: Equatable {
    case converged, unresolvedConstraints, evaluationLimitReached, proposalLimitReached, stepLimitReached, toleranceUnconverged
}
public struct EncounterSequenceToleranceReport {
    public let maximumPositionDifferenceMetres: Double
    public let maximumVelocityDifferenceMetresPerSecond: Double
    public let maximumMassDifferenceKilograms: Double
    public let maximumClosestApproachEpochDifferenceSeconds: Double
    public let maximumClosestApproachDistanceDifferenceMetres: Double
    public let consistent: Bool
}
public struct EncounterSequenceRefinementResult {
    public let status: EncounterSequenceRefinementStatus
    public let evaluations: Int
    public let rejectedTrials: Int
    public let proposalAttempts: Int
    public let initial: NumericalSequenceEvaluation
    public let best: NumericalSequenceEvaluation
    public let verification: NumericalSequenceEvaluation?
    public let toleranceReport: EncounterSequenceToleranceReport?
    public let numericalOptions: NumericalPropagationOptions
    public let verificationOptions: NumericalPropagationOptions
    public var provenance: String { "Bounded local coordinate shooting of encounter epochs, impulse epochs/vectors and finite ignition/duration/direction; no claim of a global optimum." }
}
extension NumericalEncounterSequence {
    /// Converts explicit impulse windows to engine seeds with sequential ideal fuel
    /// consumption. Windows must have positive duration and cover the impulse.
    public static func finiteSeed(_ seed: EncounterSequenceSeed, burnWindows: [ClosedRange<Epoch>],
        thrustNewtons: Double) throws -> [DeepSpaceManeuverNode] {
        let ordered = seed.maneuvers.sorted { sequenceInterval($0).lowerBound < sequenceInterval($1).lowerBound }
        guard burnWindows.count == ordered.count else { throw NumericalSequenceError.invalidInput }
        var mass = seed.initialState.mass, result: [DeepSpaceManeuverNode] = []
        for (node,window) in zip(ordered,burnWindows) {
            guard case .impulse(let impulse,let ve) = node.command else { throw NumericalSequenceError.invalidInput }
            let state = try SpacecraftState(cartesian: CartesianState(epoch: impulse.epoch,frame: impulse.frame,
                position: seed.initialState.cartesian.position,velocity: seed.initialState.cartesian.velocity),mass: mass)
            let bounded = try DeepSpaceManeuverNode(window: window,command: node.command)
            let initial = SpacecraftState(cartesian: seed.initialState.cartesian,mass: mass)
            result.append(try bounded.finiteBurnSeed(initialState: initial,thrustNewtons: thrustNewtons))
            mass = try EncounterImpulse.apply(impulse,to: state,effectiveExhaustVelocityMetresPerSecond: ve).mass
        }
        let intervals = result.map(sequenceInterval)
        guard intervals.allSatisfy({ $0.lowerBound >= seed.initialState.epoch && $0.upperBound <= seed.arrivalEpoch }),
              zip(intervals,intervals.dropFirst()).allSatisfy({ $0.upperBound < $1.lowerBound }) else { throw NumericalSequenceError.invalidInput }
        return result
    }
    /// Fixed flight bounds; encounter epochs and command times are independent.
    /// Multi-body gravity enters through backgroundModels (e.g. point-mass N-body
    /// or third-body models), with no hidden central-gravity default.
    public static func refine(initialState: SpacecraftState, arrivalEpoch: Epoch, nodes: [DeepSpaceManeuverNode],
        requirements: [NumericalEncounterRequirement], context: ForceEvaluationContext, backgroundModels: [any ForceModel],
        constraints: NumericalSequenceConstraints = try! NumericalSequenceConstraints(), maximumEvaluations: Int = 200, maximumProposalAttempts: Int = 10_000,
        initialTimeStepSeconds: Double = 60, initialDeltaVelocityStepMetresPerSecond: Double = 1,
        minimumRelativeStep: Double = 1e-5,
        numericalOptions: NumericalPropagationOptions = .preset(for: .high),
        verificationOptions: NumericalPropagationOptions? = nil,
        convergencePositionMetres: Double = 0.01, convergenceVelocityMetresPerSecond: Double = 1e-5,
        convergenceMassKilograms: Double = 1e-7,
        eventOptions: PropagationEventOptions = try! PropagationEventOptions(), maximumStateQueriesPerEvaluation: Int = 10_000) throws -> EncounterSequenceRefinementResult {
        guard (1...10_000).contains(maximumEvaluations), (1...100_000).contains(maximumProposalAttempts),
              [initialTimeStepSeconds,initialDeltaVelocityStepMetresPerSecond,minimumRelativeStep,
               convergencePositionMetres,convergenceVelocityMetresPerSecond,convergenceMassKilograms].allSatisfy({ $0.isFinite && $0 > 0 }),
              minimumRelativeStep < 1 else { throw NumericalSequenceError.invalidInput }
        let checkedOptions = try verificationOptions ?? NumericalPropagationOptions(
            absolutePositionToleranceMetres: numericalOptions.absolutePositionToleranceMetres*0.1,
            absoluteVelocityToleranceMetresPerSecond: numericalOptions.absoluteVelocityToleranceMetresPerSecond*0.1,
            absolutePropellantToleranceKilograms: numericalOptions.absolutePropellantToleranceKilograms*0.1,
            relativeTolerance: numericalOptions.relativeTolerance*0.1,
            initialStepSeconds: min(numericalOptions.initialStepSeconds,numericalOptions.maximumStepSeconds/2),
            minimumStepSeconds: min(numericalOptions.minimumStepSeconds,numericalOptions.maximumStepSeconds/2),
            maximumStepSeconds: numericalOptions.maximumStepSeconds/2,
            maximumStepAttempts: numericalOptions.maximumStepAttempts,
            maximumDerivativeEvaluations: numericalOptions.maximumDerivativeEvaluations)
        guard checkedOptions.absolutePositionToleranceMetres <= numericalOptions.absolutePositionToleranceMetres,
              checkedOptions.absoluteVelocityToleranceMetresPerSecond <= numericalOptions.absoluteVelocityToleranceMetresPerSecond,
              checkedOptions.absolutePropellantToleranceKilograms <= numericalOptions.absolutePropellantToleranceKilograms,
              checkedOptions.relativeTolerance <= numericalOptions.relativeTolerance,
              checkedOptions.maximumStepSeconds <= numericalOptions.maximumStepSeconds,
              checkedOptions != numericalOptions else { throw NumericalSequenceError.invalidInput }
        var evaluations = 0, rejected = 0, proposals = 0
        func evaluate(_ nodes: [DeepSpaceManeuverNode], _ requirements: [NumericalEncounterRequirement],
            options: NumericalPropagationOptions) throws -> NumericalSequenceEvaluation {
            evaluations += 1
            return try Self.evaluate(initialState: initialState,arrivalEpoch: arrivalEpoch,nodes: nodes,requirements: requirements,
                context: context,backgroundModels: backgroundModels,constraints: constraints,numericalOptions: options,
                eventOptions: eventOptions,maximumStateQueries: maximumStateQueriesPerEvaluation)
        }
        let initial = try evaluate(nodes,requirements,options: numericalOptions)
        var best = initial, step = 1.0
        func encountersSatisfied(_ evaluation: NumericalSequenceEvaluation) -> Bool {
            evaluation.encounters.allSatisfy(\.satisfied) && !evaluation.clearance.contains { $0.status == .violated }
        }
        func trial(_ commands: [DeepSpaceManeuverNode], _ targets: [NumericalEncounterRequirement]) throws -> Bool {
            guard evaluations < maximumEvaluations else { return false }
            do {
                let candidate = try evaluate(commands,targets,options: numericalOptions)
                if candidate.normalizedResidual < best.normalizedResidual { best = candidate; return true }
            } catch is NumericalPropagationError { rejected += 1 }
            catch is NumericalSequenceError { rejected += 1 }
            catch is EncounterPlanningError { rejected += 1 }
            catch is WaypointRoutePropagationError { rejected += 1 }
            catch is PropagationEventError { rejected += 1 }
            return false
        }
        while !encountersSatisfied(best), evaluations < maximumEvaluations, proposals < maximumProposalAttempts, step >= minimumRelativeStep {
            var improved = false
            for i in best.nodes.indices {
                for coordinate in 0..<5 {
                    for sign in [-1.0,1.0] {
                        guard evaluations < maximumEvaluations, proposals < maximumProposalAttempts else { break }
                        proposals += 1
                        var commands = best.nodes
                        let old = commands[i]
                        do {
                            let command: DeepSpaceManeuverCommand
                            switch old.command {
                            case .impulse(let impulse,let ve):
                                if coordinate == 4 { continue }
                                var epoch = impulse.epoch, dv = impulse.deltaVelocityMetresPerSecond
                                if coordinate == 0 { epoch = try epoch.advanced(by: Duration(seconds: sign*step*initialTimeStepSeconds)) }
                                else { dv[coordinate-1] += sign*step*initialDeltaVelocityStepMetresPerSecond }
                                command = .impulse(try ImpulsiveManeuver(epoch: epoch,frame: impulse.frame,deltaVelocityMetresPerSecond: dv),effectiveExhaustVelocityMetresPerSecond: ve)
                            case .finiteBurn(let burn):
                                var start = burn.start, duration = try burn.end.duration(since: burn.start).seconds, direction = burn.thrust.direction
                                if coordinate == 0 { start = try start.advanced(by: Duration(seconds: sign*step*initialTimeStepSeconds)) }
                                else if coordinate == 1 { duration *= 1+sign*step*0.25 }
                                else { direction[coordinate-2] += sign*step*0.25 }
                                command = .finiteBurn(try FiniteBurn(start: start,end: start.advanced(by: Duration(seconds: duration)),
                                    thrust: FiniteThrust(frame: burn.thrust.frame,thrustNewtons: burn.thrust.thrustNewtons,
                                        effectiveExhaustVelocityMetresPerSecond: burn.thrust.effectiveExhaustVelocityMetresPerSecond,direction: direction)))
                            }
                            commands[i] = try DeepSpaceManeuverNode(window: old.window,command: command)
                            improved = try trial(commands,best.requirements) || improved
                        } catch is EncounterPlanningError { rejected += 1 }
                        catch is FiniteThrustError { rejected += 1 }
                        if encountersSatisfied(best) { break }
                    }
                    if encountersSatisfied(best) { break }
                }
                if encountersSatisfied(best) { break }
            }
            if !encountersSatisfied(best) {
                for i in best.requirements.indices {
                    for sign in [-1.0,1.0] where evaluations < maximumEvaluations && proposals < maximumProposalAttempts {
                        proposals += 1
                        var targets = best.requirements
                        let old = targets[i], epoch = try old.epoch.advanced(by: Duration(seconds: sign*step*initialTimeStepSeconds))
                        guard old.window.contains(epoch), epoch > initialState.epoch, epoch <= arrivalEpoch,
                              (i == 0 || epoch > targets[i-1].epoch), (i == targets.count-1 || epoch < targets[i+1].epoch) else { rejected += 1; continue }
                        targets[i] = try old.shifted(to: epoch)
                        improved = try trial(best.nodes,targets) || improved
                        if encountersSatisfied(best) { break }
                    }
                    if encountersSatisfied(best) { break }
                }
            }
            if !improved { step *= 0.5 }
        }
        var verification: NumericalSequenceEvaluation?, tolerance: EncounterSequenceToleranceReport?
        let status: EncounterSequenceRefinementStatus
        if encountersSatisfied(best), evaluations < maximumEvaluations {
            let checked = try evaluate(best.nodes,best.requirements,options: checkedOptions)
            verification = checked
            var p = routeNorm(best.arrivalState.cartesian.position-checked.arrivalState.cartesian.position)
            var v = routeNorm(best.arrivalState.cartesian.velocity-checked.arrivalState.cartesian.velocity)
            var m = abs(best.arrivalState.mass.totalMassKilograms-checked.arrivalState.mass.totalMassKilograms)
            var eventTime = 0.0, eventDistance = 0.0
            for (a,b) in zip(best.encounters,checked.encounters) {
                p = max(p,routeNorm(a.state.cartesian.position-b.state.cartesian.position))
                v = max(v,routeNorm(a.state.cartesian.velocity-b.state.cartesian.velocity))
                m = max(m,abs(a.state.mass.totalMassKilograms-b.state.mass.totalMassKilograms))
                eventTime = max(eventTime,abs(try a.closestApproachEpoch.duration(since: b.closestApproachEpoch).seconds))
                eventDistance = max(eventDistance,abs(a.closestApproachDistanceMetres-b.closestApproachDistanceMetres))
            }
            let consistent = p <= convergencePositionMetres && v <= convergenceVelocityMetresPerSecond && m <= convergenceMassKilograms
            tolerance = EncounterSequenceToleranceReport(maximumPositionDifferenceMetres: p,
                maximumVelocityDifferenceMetresPerSecond: v,maximumMassDifferenceKilograms: m,
                maximumClosestApproachEpochDifferenceSeconds: eventTime,maximumClosestApproachDistanceDifferenceMetres: eventDistance,consistent: consistent)
            if !consistent || !encountersSatisfied(checked) { status = .toleranceUnconverged }
            else if !checked.constraintsResolved { status = .unresolvedConstraints }
            else { status = .converged }
        } else { status = evaluations >= maximumEvaluations ? .evaluationLimitReached :
            (proposals >= maximumProposalAttempts ? .proposalLimitReached : .stepLimitReached) }
        return EncounterSequenceRefinementResult(status: status,evaluations: evaluations,rejectedTrials: rejected,proposalAttempts: proposals,initial: initial,
            best: best,verification: verification,toleranceReport: tolerance,numericalOptions: numericalOptions,verificationOptions: checkedOptions)
    }
}

public enum NumericalSequenceCandidateStatus: Equatable {
    case converged, infeasible, unresolved, unconverged, workLimitReached
}
public struct NumericalSequenceSearchRow {
    public let seedIndex: Int
    public let status: NumericalSequenceCandidateStatus
    public let refinement: EncounterSequenceRefinementResult?
    public let failure: String?
}
public struct NumericalSequenceSearchResult {
    public let status: EncounterSequenceSearchStatus
    public let totalSeeds: Int
    public let rows: [NumericalSequenceSearchRow]
    public let rankedConvergedRows: [NumericalSequenceSearchRow]
}
extension NumericalEncounterSequence {
    /// Verifies/refines a bounded set of patched-conic seeds. Only numerically
    /// converged, constraint-resolved rows participate in the final ranking.
    public static func search(seeds: [EncounterSequenceSeed], context: ForceEvaluationContext,
        backgroundModels: [any ForceModel], constraints: NumericalSequenceConstraints = try! NumericalSequenceConstraints(),
        objective: EncounterSequenceObjective = .propellant, maximumCandidates: Int = 100,
        maximumEvaluationsPerCandidate: Int = 200, maximumProposalAttemptsPerCandidate: Int = 10_000,
        numericalOptions: NumericalPropagationOptions = .preset(for: .high),
        eventOptions: PropagationEventOptions = try! PropagationEventOptions(),
        maximumStateQueriesPerEvaluation: Int = 10_000) throws -> NumericalSequenceSearchResult {
        guard !seeds.isEmpty, seeds.count <= 100_000, (1...10_000).contains(maximumCandidates) else { throw NumericalSequenceError.invalidInput }
        var rows: [NumericalSequenceSearchRow] = []
        for (index,seed) in seeds.prefix(maximumCandidates).enumerated() {
            do {
                let result = try refine(initialState: seed.initialState,arrivalEpoch: seed.arrivalEpoch,nodes: seed.maneuvers,
                    requirements: seed.requirements,context: context,backgroundModels: backgroundModels,constraints: constraints,
                    maximumEvaluations: maximumEvaluationsPerCandidate,maximumProposalAttempts: maximumProposalAttemptsPerCandidate,
                    numericalOptions: numericalOptions,eventOptions: eventOptions,maximumStateQueriesPerEvaluation: maximumStateQueriesPerEvaluation)
                let status: NumericalSequenceCandidateStatus
                switch result.status {
                case .converged: status = .converged
                case .unresolvedConstraints: status = .unresolved
                case .evaluationLimitReached, .proposalLimitReached: status = .workLimitReached
                case .stepLimitReached: status = result.best.clearance.contains { $0.status == .violated } ? .infeasible : .unconverged
                case .toleranceUnconverged: status = .unconverged
                }
                rows.append(NumericalSequenceSearchRow(seedIndex: index,status: status,refinement: result,failure: nil))
            } catch let error as NumericalSequenceError {
                let status: NumericalSequenceCandidateStatus = error == .queryLimitReached ? .workLimitReached : .infeasible
                rows.append(NumericalSequenceSearchRow(seedIndex: index,status: status,refinement: nil,failure: String(describing: error)))
            } catch let error as NumericalPropagationError {
                let exhausted = error == .stepLimitExceeded || error == .evaluationLimitExceeded
                rows.append(NumericalSequenceSearchRow(seedIndex: index,status: exhausted ? .workLimitReached : .unconverged,
                    refinement: nil,failure: String(describing: error)))
            } catch let error as WaypointRoutePropagationError {
                rows.append(NumericalSequenceSearchRow(seedIndex: index,status: .infeasible,refinement: nil,failure: String(describing: error)))
            } catch let error as EncounterPlanningError {
                rows.append(NumericalSequenceSearchRow(seedIndex: index,status: .infeasible,refinement: nil,failure: String(describing: error)))
            } catch let error as PropagationEventError {
                rows.append(NumericalSequenceSearchRow(seedIndex: index,status: .workLimitReached,refinement: nil,failure: String(describing: error)))
            } catch GravityModelError.singularPosition {
                rows.append(NumericalSequenceSearchRow(seedIndex: index,status: .infeasible,refinement: nil,failure: "gravity singularity"))
            }
        }
        var scored: [(NumericalSequenceSearchRow,Double)] = []
        for row in rows where row.status == .converged {
            let result = row.refinement!, evaluation = result.verification ?? result.best
            let value: Double
            switch objective {
            case .propellant: value = evaluation.propellantUsedKilograms
            case .flightTime: value = try seeds[row.seedIndex].arrivalEpoch.duration(since: seeds[row.seedIndex].initialState.epoch).seconds
            case .deltaVelocity:
                // Includes ideal-equivalent rocket delta-v of finite burns, not
                // just impulse norms, using the physically carried mass inventory.
                let mass = seeds[row.seedIndex].initialState.mass.totalMassKilograms
                var total = 0.0, currentMass = mass
                for node in evaluation.nodes.sorted(by: { sequenceInterval($0).lowerBound < sequenceInterval($1).lowerBound }) {
                    switch node.command {
                    case .impulse(let impulse,let ve):
                        let dv = routeNorm(impulse.deltaVelocityMetresPerSecond)
                        total += dv; currentMass *= exp(-dv/ve)
                    case .finiteBurn(let burn):
                        let fuel = try burn.end.duration(since: burn.start).seconds*burn.thrust.propellantConsumptionKilogramsPerSecond
                        total += burn.thrust.effectiveExhaustVelocityMetresPerSecond*log(currentMass/(currentMass-fuel))
                        currentMass -= fuel
                    }
                }
                value = total
            }
            scored.append((row,value))
        }
        scored.sort { $0.1 == $1.1 ? $0.0.seedIndex < $1.0.seedIndex : $0.1 < $1.1 }
        return NumericalSequenceSearchResult(status: rows.count == seeds.count ? .candidatesExhausted : .evaluationLimitReached,
            totalSeeds: seeds.count,rows: rows,rankedConvergedRows: scored.map { $0.0 })
    }
}
