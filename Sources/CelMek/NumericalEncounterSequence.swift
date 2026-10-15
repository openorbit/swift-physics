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

public enum NumericalSequenceError: Error, Equatable {
    case invalidInput, queryLimitReached, constraintViolation(String), insufficientReserve
}
public struct NumericalEncounterRequirement {
    public let identifier: String
    public let target: MovingWaypointRouteTarget
    public let epoch: Epoch
    public let window: ClosedRange<Epoch>
    public let matchVelocity: Bool
    public let positionToleranceMetres: Double
    public let velocityToleranceMetresPerSecond: Double
    public let minimumPassageSpeedMetresPerSecond: Double
    public let closestApproachTarget: MovingWaypointRouteTarget
    public init(identifier: String, target: MovingWaypointRouteTarget, epoch: Epoch, window: ClosedRange<Epoch>,
        matchVelocity: Bool = false, positionToleranceMetres: Double = 1,
        velocityToleranceMetresPerSecond: Double = 0.001, minimumPassageSpeedMetresPerSecond: Double = 0,
        closestApproachTarget: MovingWaypointRouteTarget? = nil) throws {
        guard !identifier.isEmpty, window.contains(epoch), window.lowerBound >= target.coverage.lowerBound,
              window.upperBound <= target.coverage.upperBound,
              [positionToleranceMetres,velocityToleranceMetresPerSecond].allSatisfy({ $0.isFinite && $0 > 0 }),
              minimumPassageSpeedMetresPerSecond.isFinite, minimumPassageSpeedMetresPerSecond >= 0 else { throw NumericalSequenceError.invalidInput }
        if let closestApproachTarget {
            guard closestApproachTarget.frame === target.frame, closestApproachTarget.origin == target.origin,
                  closestApproachTarget.coverage.lowerBound <= window.lowerBound,
                  closestApproachTarget.coverage.upperBound >= window.upperBound else { throw NumericalSequenceError.invalidInput }
        }
        self.identifier = identifier; self.target = target; self.epoch = epoch; self.window = window
        self.matchVelocity = matchVelocity; self.positionToleranceMetres = positionToleranceMetres
        self.velocityToleranceMetresPerSecond = velocityToleranceMetresPerSecond
        self.minimumPassageSpeedMetresPerSecond = minimumPassageSpeedMetresPerSecond
        self.closestApproachTarget = closestApproachTarget ?? target
    }
    func shifted(to epoch: Epoch) throws -> Self {
        try Self(identifier: identifier,target: target,epoch: epoch,window: window,matchVelocity: matchVelocity,
            positionToleranceMetres: positionToleranceMetres,velocityToleranceMetresPerSecond: velocityToleranceMetresPerSecond,
            minimumPassageSpeedMetresPerSecond: minimumPassageSpeedMetresPerSecond,closestApproachTarget: closestApproachTarget)
    }
}
public struct NumericalSequenceClearance {
    public let identifier: String
    public let centre: MovingWaypointRouteTarget
    public let minimumDistanceMetres: Double
    /// Caller-guaranteed global relative-speed bound. Enables conservative interval
    /// certification; observed speeds alone cannot establish a between-sample bound.
    public let relativeSpeedBoundMetresPerSecond: Double?
    public init(identifier: String, centre: MovingWaypointRouteTarget, minimumDistanceMetres: Double,
        relativeSpeedBoundMetresPerSecond: Double? = nil) throws {
        guard !identifier.isEmpty, minimumDistanceMetres.isFinite, minimumDistanceMetres > 0,
              relativeSpeedBoundMetresPerSecond == nil || (relativeSpeedBoundMetresPerSecond!.isFinite && relativeSpeedBoundMetresPerSecond! >= 0) else {
            throw NumericalSequenceError.invalidInput
        }
        self.identifier = identifier; self.centre = centre; self.minimumDistanceMetres = minimumDistanceMetres
        self.relativeSpeedBoundMetresPerSecond = relativeSpeedBoundMetresPerSecond
    }
}
public enum NumericalSequenceCheckStatus: Equatable { case satisfied, violated, unresolved }
public struct NumericalSequenceClearanceReport {
    public let identifier: String
    public let status: NumericalSequenceCheckStatus
    public let minimumSampledDistanceMetres: Double
    public let usesCallerSpeedBound: Bool
}
public struct NumericalSequenceConstraints {
    public let clearances: [NumericalSequenceClearance]
    public let physical: WaypointRoutePhysicalConstraints
    public let minimumRemainingPropellantKilograms: Double
    public let maximumThrustNewtons: Double?
    public init(clearances: [NumericalSequenceClearance] = [], physical: WaypointRoutePhysicalConstraints = .unconstrained,
        minimumRemainingPropellantKilograms: Double = 0, maximumThrustNewtons: Double? = nil) throws {
        guard clearances.count <= 100, Set(clearances.map(\.identifier)).count == clearances.count,
              minimumRemainingPropellantKilograms.isFinite, minimumRemainingPropellantKilograms >= 0,
              maximumThrustNewtons == nil || (maximumThrustNewtons!.isFinite && maximumThrustNewtons! > 0) else { throw NumericalSequenceError.invalidInput }
        self.clearances = clearances; self.physical = physical
        self.minimumRemainingPropellantKilograms = minimumRemainingPropellantKilograms; self.maximumThrustNewtons = maximumThrustNewtons
    }
}
public struct NumericalEncounterResidual {
    public let identifier: String
    public let epoch: Epoch
    public let state: SpacecraftState
    public let positionErrorMetres: Double
    public let relativeSpeedMetresPerSecond: Double
    public let satisfied: Bool
    public let closestApproachEpoch: Epoch
    public let closestApproachDistanceMetres: Double
    public let localizationBracketSeconds: Double
}
public struct NumericalSequenceEvaluation {
    public let nodes: [DeepSpaceManeuverNode]
    public let requirements: [NumericalEncounterRequirement]
    public let arrivalState: SpacecraftState
    public let encounters: [NumericalEncounterResidual]
    public let clearance: [NumericalSequenceClearanceReport]
    public let finalHeatJoules: Double?
    public let stateQueries: Int
    public let normalizedResidual: Double
    public let propellantUsedKilograms: Double
    public var constraintsResolved: Bool { clearance.allSatisfy { $0.status == .satisfied } }
    public var converged: Bool { encounters.allSatisfy(\.satisfied) && constraintsResolved }
    public let provenance: String
    public let forceModelTypes: [String]
    public let contextOrigin: BodyIdentifier
    public let numericalOptions: NumericalPropagationOptions
    public let eventOptions: PropagationEventOptions
    public let maximumStateQueries: Int
}
/// No reset at encounters: all queries replay the actual command sequence from
/// one mass-aware seed in one physical origin, under the supplied force models.
public struct NumericalEncounterTrajectory: SpacecraftStateProvider {
    public let initialState: SpacecraftState
    public let arrivalEpoch: Epoch
    public let nodes: [DeepSpaceManeuverNode]
    public let context: ForceEvaluationContext
    public let backgroundModels: [any ForceModel]
    public let options: NumericalPropagationOptions
    public var frame: Frame { initialState.frame }
    public init(initialState: SpacecraftState, arrivalEpoch: Epoch, nodes: [DeepSpaceManeuverNode],
        context: ForceEvaluationContext, backgroundModels: [any ForceModel],
        options: NumericalPropagationOptions = .preset(for: .high)) throws {
        guard arrivalEpoch > initialState.epoch, nodes.count <= 100, context.frame === initialState.frame else { throw NumericalSequenceError.invalidInput }
        var last: Epoch?
        for node in nodes.sorted(by: { sequenceInterval($0).lowerBound < sequenceInterval($1).lowerBound }) {
            let interval = sequenceInterval(node)
            guard interval.lowerBound >= initialState.epoch, interval.upperBound <= arrivalEpoch,
                  last == nil || interval.lowerBound > last! else { throw NumericalSequenceError.invalidInput }
            switch node.command {
            case .impulse(let maneuver,_): guard maneuver.frame === initialState.frame else { throw SpacecraftDynamicsError.frameMismatch }
            case .finiteBurn(let burn): guard burn.thrust.frame === initialState.frame else { throw SpacecraftDynamicsError.frameMismatch }
            }
            last = interval.upperBound
        }
        self.initialState = initialState; self.arrivalEpoch = arrivalEpoch; self.nodes = nodes; self.context = context
        self.backgroundModels = backgroundModels; self.options = options
    }
    public func state(at epoch: Epoch) throws -> CartesianState { try spacecraftState(at: epoch).cartesian }
    public func spacecraftState(at epoch: Epoch) throws -> SpacecraftState {
        guard epoch >= initialState.epoch, epoch <= arrivalEpoch else { throw LambertTransferPlanningError.outsideFlight }
        var prefix: [DeepSpaceManeuverNode] = []
        for node in nodes {
            switch node.command {
            case .impulse(let maneuver,_): if maneuver.epoch <= epoch { prefix.append(node) }
            case .finiteBurn(let burn):
                if burn.start < epoch {
                    let truncated = try FiniteBurn(start: burn.start,end: min(epoch,burn.end),thrust: burn.thrust)
                    prefix.append(try DeepSpaceManeuverNode(window: node.window,command: .finiteBurn(truncated)))
                }
            }
        }
        return try DeepSpaceManeuverNode.execute(prefix,initialState: initialState,arrivalEpoch: epoch,
            context: context,backgroundModels: backgroundModels,options: options)
    }
}
func sequenceInterval(_ node: DeepSpaceManeuverNode) -> ClosedRange<Epoch> {
    switch node.command { case .impulse(let m,_): return m.epoch...m.epoch
    case .finiteBurn(let b): return b.start...b.end }
}
private struct SequenceReserveMonitor: ForceModel {
    let minimum: Double
    func contribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
        guard state.mass.propellantMassKilograms >= minimum else { throw NumericalSequenceError.insufficientReserve }
        return try ForceContribution(epoch: state.epoch,frame: state.frame,accelerationMetresPerSecondSquared: .zero)
    }
}
private struct SequenceAttitudeMonitor: ForceModel {
    let cone: WaypointThrustCone
    let nodes: [DeepSpaceManeuverNode]
    func contribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
        for (index,node) in nodes.enumerated() {
            if case .finiteBurn(let burn) = node.command, state.epoch >= burn.start, state.epoch < burn.end {
                try cone.check(burn.thrust.direction*burn.thrust.thrustNewtons,at: state.epoch,leg: index)
            }
        }
        return try ForceContribution(epoch: state.epoch,frame: state.frame,accelerationMetresPerSecondSquared: .zero)
    }
}
private struct MovingClosestApproachEvent: PropagationEventDetector {
    let identifier: String
    let target: MovingWaypointRouteTarget
    let direction: EventCrossingDirection = .increasing
    let action: PropagationEventAction = .record
    func value(at state: SpacecraftState) throws -> Double {
        let body = try target.state(at: state.epoch)
        return ((state.cartesian.position-body.position)*(state.cartesian.velocity-body.velocity)).sum()
    }
}
private final class SequenceQueryCounter { var count = 0 }
private struct CountedSequenceTrajectory: SpacecraftStateProvider {
    let base: NumericalEncounterTrajectory
    let counter: SequenceQueryCounter
    let limit: Int
    var frame: Frame { base.frame }
    func state(at epoch: Epoch) throws -> CartesianState { try spacecraftState(at: epoch).cartesian }
    func spacecraftState(at epoch: Epoch) throws -> SpacecraftState {
        guard counter.count < limit else { throw NumericalSequenceError.queryLimitReached }
        counter.count += 1
        return try base.spacecraftState(at: epoch)
    }
}
public enum NumericalEncounterSequence {
    public static func evaluate(initialState: SpacecraftState, arrivalEpoch: Epoch, nodes: [DeepSpaceManeuverNode],
        requirements: [NumericalEncounterRequirement], context: ForceEvaluationContext, backgroundModels: [any ForceModel],
        constraints: NumericalSequenceConstraints = try! NumericalSequenceConstraints(),
        numericalOptions: NumericalPropagationOptions = .preset(for: .high),
        eventOptions: PropagationEventOptions = try! PropagationEventOptions(), maximumStateQueries: Int = 10_000) throws -> NumericalSequenceEvaluation {
        guard !requirements.isEmpty, requirements.count <= 100, (1...100_000).contains(maximumStateQueries),
              Set(requirements.map(\.identifier)).count == requirements.count,
              zip(requirements,requirements.dropFirst()).allSatisfy({ $0.epoch < $1.epoch }),
              requirements.allSatisfy({ $0.epoch > initialState.epoch && $0.epoch <= arrivalEpoch &&
                  $0.target.origin == context.origin && $0.target.frame === context.frame }) else { throw NumericalSequenceError.invalidInput }
        try constraints.physical.validate(frame: context.frame)
        guard constraints.clearances.allSatisfy({ $0.centre.frame === context.frame && $0.centre.origin == context.origin &&
            $0.centre.coverage.lowerBound <= initialState.epoch && $0.centre.coverage.upperBound >= arrivalEpoch }) else { throw NumericalSequenceError.invalidInput }
        var models = backgroundModels
        if constraints.minimumRemainingPropellantKilograms > 0 { models.append(SequenceReserveMonitor(minimum: constraints.minimumRemainingPropellantKilograms)) }
        if let cone = constraints.physical.thrustCone {
            // Impulsive engine orientation/thermal feasibility cannot be established
            // from a delta-v alone: require finite commands for these constraints.
            guard !nodes.contains(where: { if case .impulse = $0.command { return true }; return false }) else { throw NumericalSequenceError.constraintViolation("attitude requires finite commands") }
            models.append(SequenceAttitudeMonitor(cone: cone,nodes: nodes))
        }
        if constraints.physical.thermalEnvelope != nil || constraints.maximumThrustNewtons != nil {
            guard !nodes.contains(where: { if case .impulse = $0.command { return true }; return false }) else { throw NumericalSequenceError.constraintViolation("thermal/engine constraints require finite commands") }
        }
        for node in nodes {
            if case .finiteBurn(let burn) = node.command, let ceiling = constraints.maximumThrustNewtons,
               burn.thrust.thrustNewtons > ceiling { throw NumericalSequenceError.constraintViolation("engine ceiling") }
        }
        let trajectory = try NumericalEncounterTrajectory(initialState: initialState,arrivalEpoch: arrivalEpoch,nodes: nodes,
            context: context,backgroundModels: models,options: numericalOptions)
        let counter = SequenceQueryCounter(), counted = CountedSequenceTrajectory(base: trajectory,counter: counter,limit: maximumStateQueries)
        let arrival = try counted.spacecraftState(at: arrivalEpoch)
        guard arrival.mass.propellantMassKilograms >= constraints.minimumRemainingPropellantKilograms else { throw NumericalSequenceError.insufficientReserve }
        var residuals: [NumericalEncounterResidual] = [], score = 0.0
        for requirement in requirements {
            let state = try counted.spacecraftState(at: requirement.epoch), target = try requirement.target.state(at: requirement.epoch)
            let position = routeNorm(state.cartesian.position-target.position), speed = routeNorm(state.cartesian.velocity-target.velocity)
            var residual = position/requirement.positionToleranceMetres
            if requirement.matchVelocity { residual = hypot(residual,speed/requirement.velocityToleranceMetresPerSecond) }
            else if speed <= requirement.minimumPassageSpeedMetresPerSecond { residual = hypot(residual,1+(requirement.minimumPassageSpeedMetresPerSecond-speed)/requirement.velocityToleranceMetresPerSecond) }
            score = hypot(score,residual)
            let start = max(initialState.epoch,requirement.window.lowerBound), end = min(arrivalEpoch,requirement.window.upperBound)
            let event = MovingClosestApproachEvent(identifier: requirement.identifier,target: requirement.closestApproachTarget)
            // Split at every impulse/engine boundary: crossing localization must
            // never span a discontinuous velocity command.
            var dates = [start,end]
            for node in nodes { for date in [sequenceInterval(node).lowerBound,sequenceInterval(node).upperBound] where date > start && date < end { dates.append(date) } }
            dates = Array(Set(dates)).sorted()
            var candidates: [(SpacecraftState,Double)] = try dates.map { (try counted.spacecraftState(at: $0),0) }
            for (a,b) in zip(dates,dates.dropFirst()) {
                // Endpoint impulses are right-continuous. Their epoch is already
                // included as a candidate; a bracket ending there is conservatively
                // excluded from event localization to avoid a false jump root.
                let upper: Epoch
                if nodes.contains(where: { if case .impulse(let m,_) = $0.command { return m.epoch == b }; return false }) {
                    upper = try b.advanced(by: Duration(seconds: -min(eventOptions.timeToleranceSeconds,try b.duration(since: a).seconds/2)))
                } else { upper = b }
                if upper > a {
                    let events = try EventDetectedPropagation(trajectory: counted,detectors: [event],options: eventOptions)
                        .propagation(from: a,to: upper)
                    candidates += events.events.map { ($0.state,$0.bracketDurationSeconds) }
                }
            }
            let closest = try candidates.map { candidate in
                (candidate.0,candidate.1,routeNorm(candidate.0.cartesian.position - (try requirement.closestApproachTarget.state(at: candidate.0.epoch)).position))
            }.min { $0.2 < $1.2 }!
            residuals.append(NumericalEncounterResidual(identifier: requirement.identifier,epoch: requirement.epoch,state: state,
                positionErrorMetres: position,relativeSpeedMetresPerSecond: speed,
                satisfied: position <= requirement.positionToleranceMetres && (requirement.matchVelocity ? speed <= requirement.velocityToleranceMetresPerSecond : speed > requirement.minimumPassageSpeedMetresPerSecond),
                closestApproachEpoch: closest.0.epoch,closestApproachDistanceMetres: closest.2,localizationBracketSeconds: closest.1))
        }
        var clearanceReports: [NumericalSequenceClearanceReport] = []
        var clearancePenalty = 0.0
        if !constraints.clearances.isEmpty {
            let seconds = try arrivalEpoch.duration(since: initialState.epoch).seconds
            let countDouble = ceil(seconds/eventOptions.maximumCheckIntervalSeconds)
            guard countDouble.isFinite, countDouble <= Double(maximumStateQueries) else { throw NumericalSequenceError.queryLimitReached }
            let count = max(1,Int(countDouble))
            let states = try (0...count).map { index in
                try counted.spacecraftState(at: index == count ? arrivalEpoch : initialState.epoch.advanced(by: Duration(seconds: seconds*Double(index)/Double(count))))
            }
            for clearance in constraints.clearances {
                var distances: [Double] = [], speedValid = true
                for state in states {
                    let centre = try clearance.centre.state(at: state.epoch)
                    distances.append(routeNorm(state.cartesian.position-centre.position))
                    if let bound = clearance.relativeSpeedBoundMetresPerSecond,
                       routeNorm(state.cartesian.velocity-centre.velocity) > bound*(1+1e-12) { speedValid = false }
                }
                var status: NumericalSequenceCheckStatus = .unresolved
                if distances.min()! < clearance.minimumDistanceMetres || !speedValid { status = .violated }
                else if let bound = clearance.relativeSpeedBoundMetresPerSecond {
                    let dt = seconds/Double(count)
                    if zip(distances,distances.dropFirst()).allSatisfy({ min($0,$1)-bound*dt/2 >= clearance.minimumDistanceMetres }) { status = .satisfied }
                }
                if status == .violated {
                    clearancePenalty = hypot(clearancePenalty,1+max(0,clearance.minimumDistanceMetres-distances.min()!)/clearance.minimumDistanceMetres)
                }
                clearanceReports.append(NumericalSequenceClearanceReport(identifier: clearance.identifier,status: status,
                    minimumSampledDistanceMetres: distances.min()!,usesCallerSpeedBound: clearance.relativeSpeedBoundMetresPerSecond != nil))
            }
        }
        var heat = constraints.physical.thermalEnvelope?.initialHeatJoules
        if let envelope = constraints.physical.thermalEnvelope {
            var dates = [initialState.epoch,arrivalEpoch]
            for node in nodes { dates += [sequenceInterval(node).lowerBound,sequenceInterval(node).upperBound] }
            dates = Array(Set(dates)).sorted()
            for (index,pair) in zip(dates,dates.dropFirst()).enumerated() {
                var thrust = 0.0
                for node in nodes { if case .finiteBurn(let burn) = node.command, pair.0 >= burn.start, pair.0 < burn.end { thrust += burn.thrust.thrustNewtons } }
                heat = try envelope.advance(heat: heat!,duration: pair.1.duration(since: pair.0).seconds,thrustCeiling: thrust,leg: index)
            }
        }
        score = hypot(score,clearancePenalty)
        guard !score.isNaN else { throw NumericalSequenceError.invalidInput }
        return NumericalSequenceEvaluation(nodes: nodes,requirements: requirements,arrivalState: arrival,encounters: residuals,
            clearance: clearanceReports,finalHeatJoules: heat,stateQueries: counter.count,normalizedResidual: score,
            propellantUsedKilograms: initialState.mass.propellantMassKilograms-arrival.mass.propellantMassKilograms,
            provenance: "Continuous numerical replay under caller-supplied forces; localized moving-target closest approaches; sampled clearance with optional caller-guaranteed speed envelopes. Local error tolerances are not global accuracy certificates.",
            forceModelTypes: backgroundModels.map { String(reflecting: type(of: $0)) },contextOrigin: context.origin,
            numericalOptions: numericalOptions,eventOptions: eventOptions,maximumStateQueries: maximumStateQueries)
    }
}
