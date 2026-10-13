// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum EncounterSequenceError: Error, Equatable {
    case invalidInput, outsideWindow(index: Int), overlappingEncounters(index: Int)
    case legFailed(index: Int, reason: String), invalidEndpoint, incompatibleOrigin
}
public enum EncounterSequenceCondition {
    /// Geometric or target-relative passage; no velocity matching impulse.
    case passage
    /// Explicit endpoint rest in the target's coordinates (capture/landing dynamics are not inferred).
    case velocityMatch
    /// Finite local arc around periapsis, connected by explicit patch impulses.
    /// The tree must map the request's body-centred frame into route coordinates.
    case unpoweredFlyby(request: FlybyEncounterRequest, incomingSeconds: Double, outgoingSeconds: Double, tree: FrameTree)
}
public struct EncounterSequenceVisit {
    public let identifier: String
    public let target: MovingWaypointRouteTarget
    public let arrivalWindow: ClosedRange<Epoch>
    public let condition: EncounterSequenceCondition
    public init(identifier: String, target: MovingWaypointRouteTarget, arrivalWindow: ClosedRange<Epoch>,
        condition: EncounterSequenceCondition = .passage) throws {
        guard !identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              arrivalWindow.lowerBound >= target.coverage.lowerBound,
              arrivalWindow.upperBound <= target.coverage.upperBound else { throw EncounterSequenceError.invalidInput }
        self.identifier = identifier; self.target = target; self.arrivalWindow = arrivalWindow; self.condition = condition
    }
}
public struct EncounterSequenceCandidate {
    public let epochs: [Epoch]
    /// One branch per Lambert connector, including connectors before/after local flybys.
    public let branches: [LambertWindowBranch]
    public let deepSpaceManeuvers: [DeepSpaceManeuverNode]
    public init(epochs: [Epoch], branches: [LambertWindowBranch], deepSpaceManeuvers: [DeepSpaceManeuverNode] = []) {
        self.epochs = epochs; self.branches = branches; self.deepSpaceManeuvers = deepSpaceManeuvers
    }
}
public struct EncounterSequencePatchResidual {
    public let identifier: String
    public let epoch: Epoch
    public let positionMetres: Double
    public let preImpulseVelocityMetresPerSecond: Double
}
public struct EncounterSequenceArc {
    public let coverage: ClosedRange<Epoch>
    public let provenance: String
    private let query: (Epoch) throws -> CartesianState
    init(coverage: ClosedRange<Epoch>, provenance: String, query: @escaping (Epoch) throws -> CartesianState) {
        self.coverage = coverage; self.provenance = provenance; self.query = query
    }
    public func state(at epoch: Epoch) throws -> CartesianState {
        guard coverage.contains(epoch) else { throw LambertTransferPlanningError.outsideFlight }
        return try query(epoch)
    }
}
public struct EncounterSequenceSeed {
    public let initialState: SpacecraftState
    public let candidate: EncounterSequenceCandidate
    public let arcs: [EncounterSequenceArc]
    public let maneuvers: [DeepSpaceManeuverNode]
    public let requirements: [NumericalEncounterRequirement]
    public let patchResiduals: [EncounterSequencePatchResidual]
    public let finalMass: SpacecraftMassProperties
    public let totalDeltaVelocityMetresPerSecond: Double
    public var arrivalEpoch: Epoch { candidate.epochs.last! }
    public var propellantUsedKilograms: Double { initialState.mass.propellantMassKilograms-finalMass.propellantMassKilograms }
    public var provenance: String { "Lambert central-body connectors; finite unpowered two-body flybys; explicit instantaneous patch corrections. Optional DSMs perturb these nominal arcs; their residuals and all numerical dynamics are verified separately." }
}
public enum EncounterSequenceObjective: Equatable { case flightTime, deltaVelocity, propellant }
public enum EncounterSequenceSearchStatus: Equatable { case candidatesExhausted, evaluationLimitReached }
public struct EncounterSequenceSearchRow {
    public let candidateIndex: Int
    public let seed: EncounterSequenceSeed?
    public let failure: String?
}
public struct EncounterSequenceSearchResult {
    public let status: EncounterSequenceSearchStatus
    public let totalCandidates: Int
    public let rows: [EncounterSequenceSearchRow]
    public let rankedSeeds: [EncounterSequenceSeed]
}

public enum EncounterSequencePlanner {
    /// Explicit Cartesian grid of date, Lambert-root and independent DSM choices.
    /// The bounded cardinality is checked before allocating combinations.
    public static func candidates(epochOptions: [[Epoch]], branchOptions: [[LambertWindowBranch]],
        deepSpaceOptions: [[DeepSpaceManeuverNode]] = [[]], maximumCandidates: Int = 100_000) throws -> [EncounterSequenceCandidate] {
        guard !epochOptions.isEmpty, epochOptions.count <= 100, branchOptions.count == epochOptions.count,
              !deepSpaceOptions.isEmpty, (1...100_000).contains(maximumCandidates),
              epochOptions.allSatisfy({ !$0.isEmpty && $0.count <= maximumCandidates }),
              branchOptions.allSatisfy({ !$0.isEmpty && $0.count <= 256 }),
              deepSpaceOptions.count <= maximumCandidates,
              deepSpaceOptions.allSatisfy({ $0.count <= 100 }) else { throw EncounterSequenceError.invalidInput }
        let epochs = epochOptions.map { Array(Set($0)).sorted() }
        var cardinality = deepSpaceOptions.count
        for i in epochs.indices {
            for count in [epochs[i].count,branchOptions[i].count] {
                guard cardinality <= maximumCandidates/count else { throw EncounterSequenceError.invalidInput }
                cardinality *= count
            }
        }
        var result: [EncounterSequenceCandidate] = []
        func append(_ index: Int, _ dates: [Epoch], _ branches: [LambertWindowBranch]) {
            if index == epochs.count {
                for commands in deepSpaceOptions { result.append(EncounterSequenceCandidate(epochs: dates,branches: branches,deepSpaceManeuvers: commands)) }
                return
            }
            for epoch in epochs[index] { for branch in branchOptions[index] { append(index+1,dates+[epoch],branches+[branch]) } }
        }
        append(0,[],[])
        return result
    }

    /// Date and branch search over an explicit bounded set, including independent
    /// DSM command candidates. Provider contract errors propagate; expected physical
    /// and solver failures remain visible as candidate rows.
    public static func search(initialState: SpacecraftState, visits: [EncounterSequenceVisit],
        candidates: [EncounterSequenceCandidate], context: ForceEvaluationContext,
        gravitationalParameter: GravitationalParameter, referencePole: SIMD3<Double>,
        effectiveExhaustVelocityMetresPerSecond: Double, additionalPassages: [NumericalEncounterRequirement] = [], objective: EncounterSequenceObjective = .deltaVelocity,
        maximumEvaluations: Int = 1000, positionToleranceMetres: Double = 0.01,
        velocityToleranceMetresPerSecond: Double = 1e-5, lambertOptions: LambertOptions = .default,
        propagationOptions: TwoBodyPropagationOptions = .default) throws -> EncounterSequenceSearchResult {
        guard !candidates.isEmpty, candidates.count <= 100_000, (1...10_000).contains(maximumEvaluations) else {
            throw EncounterSequenceError.invalidInput
        }
        try validate(initialState: initialState,visits: visits,context: context,exhaustVelocity: effectiveExhaustVelocityMetresPerSecond)
        var rows: [EncounterSequenceSearchRow] = []
        for (i,candidate) in candidates.prefix(maximumEvaluations).enumerated() {
            do {
                let result = try seed(initialState: initialState,visits: visits,candidate: candidate,context: context,
                    gravitationalParameter: gravitationalParameter,referencePole: referencePole,
                    effectiveExhaustVelocityMetresPerSecond: effectiveExhaustVelocityMetresPerSecond,additionalPassages: additionalPassages,
                    positionToleranceMetres: positionToleranceMetres,velocityToleranceMetresPerSecond: velocityToleranceMetresPerSecond,
                    lambertOptions: lambertOptions,propagationOptions: propagationOptions)
                rows.append(EncounterSequenceSearchRow(candidateIndex: i,seed: result,failure: nil))
            } catch let error as EncounterSequenceError { rows.append(EncounterSequenceSearchRow(candidateIndex: i,seed: nil,failure: String(describing: error))) }
            catch let error as EncounterPlanningError { rows.append(EncounterSequenceSearchRow(candidateIndex: i,seed: nil,failure: String(describing: error))) }
            catch let error as FlybyError { rows.append(EncounterSequenceSearchRow(candidateIndex: i,seed: nil,failure: String(describing: error))) }
            catch let error as NumericalPropagationError { rows.append(EncounterSequenceSearchRow(candidateIndex: i,seed: nil,failure: String(describing: error))) }
        }
        var scored: [(Int,EncounterSequenceSeed,Double)] = []
        for row in rows {
            guard let result = row.seed else { continue }
            let value: Double
            switch objective {
            case .flightTime:
                let duration = try result.arrivalEpoch.duration(since: initialState.epoch)
                value = duration.seconds
            case .deltaVelocity: value = result.totalDeltaVelocityMetresPerSecond
            case .propellant: value = result.propellantUsedKilograms
            }
            scored.append((row.candidateIndex,result,value))
        }
        scored.sort { $0.2 == $1.2 ? $0.0 < $1.0 : $0.2 < $1.2 }
        let ranked = scored.map { $0.1 }
        return EncounterSequenceSearchResult(status: rows.count == candidates.count ? .candidatesExhausted : .evaluationLimitReached,
            totalCandidates: candidates.count,rows: rows,rankedSeeds: ranked)
    }
    private static func validate(initialState: SpacecraftState, visits: [EncounterSequenceVisit],
        context: ForceEvaluationContext, exhaustVelocity: Double) throws {
        guard !visits.isEmpty, visits.count <= 100, Set(visits.map(\.identifier)).count == visits.count,
              exhaustVelocity.isFinite, exhaustVelocity > 0 else { throw EncounterSequenceError.invalidInput }
        guard context.frame === initialState.frame, visits.allSatisfy({ $0.target.frame === context.frame }) else { throw SpacecraftDynamicsError.frameMismatch }
        guard visits.allSatisfy({ $0.target.origin == context.origin }) else { throw EncounterSequenceError.incompatibleOrigin }
        guard initialState.frame.pseudoInertial else { throw PropagationError.nonInertialFrame }
        if case .unpoweredFlyby = visits.last!.condition { throw EncounterSequenceError.invalidEndpoint }
    }
    public static func seed(initialState: SpacecraftState, visits: [EncounterSequenceVisit], candidate: EncounterSequenceCandidate,
        context: ForceEvaluationContext, gravitationalParameter: GravitationalParameter, referencePole: SIMD3<Double>,
        effectiveExhaustVelocityMetresPerSecond: Double, additionalPassages: [NumericalEncounterRequirement] = [], positionToleranceMetres: Double = 0.01,
        velocityToleranceMetresPerSecond: Double = 1e-5, lambertOptions: LambertOptions = .default,
        propagationOptions: TwoBodyPropagationOptions = .default) throws -> EncounterSequenceSeed {
        try validate(initialState: initialState,visits: visits,context: context,exhaustVelocity: effectiveExhaustVelocityMetresPerSecond)
        guard candidate.epochs.count == visits.count, candidate.branches.count == visits.count,
              candidate.deepSpaceManeuvers.count <= 100,
              [positionToleranceMetres,velocityToleranceMetresPerSecond].allSatisfy({ $0.isFinite && $0 > 0 }) else { throw EncounterSequenceError.invalidInput }
        var current = initialState.cartesian, mass = initialState.mass
        var arcs: [EncounterSequenceArc] = [], nodes: [DeepSpaceManeuverNode] = []
        var residuals: [EncounterSequencePatchResidual] = [], requirements: [NumericalEncounterRequirement] = []
        var totalDV = 0.0
        func impulse(_ maneuver: ImpulsiveManeuver) throws {
            let post = try EncounterImpulse.apply(maneuver,to: SpacecraftState(cartesian: current,mass: mass),
                effectiveExhaustVelocityMetresPerSecond: effectiveExhaustVelocityMetresPerSecond)
            current = post.cartesian; mass = post.mass; totalDV += routeNorm(maneuver.deltaVelocityMetresPerSecond)
            if routeNorm(maneuver.deltaVelocityMetresPerSecond) > 0 {
                nodes.append(try DeepSpaceManeuverNode(window: maneuver.epoch...maneuver.epoch,
                    command: .impulse(maneuver,effectiveExhaustVelocityMetresPerSecond: effectiveExhaustVelocityMetresPerSecond)))
            }
        }
        func connect(_ target: CartesianState, matchVelocity: Bool, index: Int, identifier: String) throws {
            guard target.epoch > current.epoch else { throw EncounterSequenceError.overlappingEncounters(index: index) }
            let branch = candidate.branches[index]
            let plan: LambertTransferPlan
            do {
                plan = try LambertTransferPlanner.plan(departureState: current,arrivalTargetState: target,
                    centralBody: context.origin,gravitationalParameter: gravitationalParameter,referencePole: referencePole,
                    direction: branch.direction,arrivalCondition: matchVelocity ? .velocityMatchedRendezvous : .flyby,
                    revolutions: branch.revolutions,multiRevolutionBranch: branch.multiRevolutionBranch,
                    positionToleranceMetres: positionToleranceMetres,velocityToleranceMetresPerSecond: velocityToleranceMetresPerSecond,
                    lambertOptions: lambertOptions,propagationOptions: propagationOptions)
            } catch let error as LambertError { throw EncounterSequenceError.legFailed(index: index,reason: String(describing: error)) }
            catch let error as PropagationError { throw EncounterSequenceError.legFailed(index: index,reason: String(describing: error)) }
            guard plan.arrivalConditionSatisfied else { throw EncounterSequenceError.legFailed(index: index,reason: "endpoint residual") }
            try impulse(plan.departureBurn)
            arcs.append(EncounterSequenceArc(coverage: current.epoch...target.epoch,provenance: "Central-body Lambert") { try plan.trajectory.coast.state(at: $0) })
            current = plan.preArrivalState
            residuals.append(EncounterSequencePatchResidual(identifier: identifier,epoch: target.epoch,
                positionMetres: routeNorm(current.position-target.position),preImpulseVelocityMetresPerSecond: routeNorm(current.velocity-target.velocity)))
            if let burn = plan.arrivalBurn { try impulse(burn) }
        }
        for i in visits.indices {
            let visit = visits[i], epoch = candidate.epochs[i]
            guard visit.arrivalWindow.contains(epoch) else { throw EncounterSequenceError.outsideWindow(index: i) }
            switch visit.condition {
            case .passage, .velocityMatch:
                if case .velocityMatch = visit.condition, i != visits.count-1 { throw EncounterSequenceError.invalidEndpoint }
                let target = try visit.target.state(at: epoch)
                let match: Bool
                if case .velocityMatch = visit.condition { match = true } else { match = false }
                try connect(target,matchVelocity: match,index: i,identifier: visit.identifier)
                requirements.append(try NumericalEncounterRequirement(identifier: visit.identifier,target: visit.target,
                    epoch: epoch,window: visit.arrivalWindow,matchVelocity: match,positionToleranceMetres: positionToleranceMetres,
                    velocityToleranceMetresPerSecond: velocityToleranceMetresPerSecond))
            case .unpoweredFlyby(let request,let lead,let trail,let tree):
                guard lead.isFinite, lead > 0, trail.isFinite, trail > 0 else { throw EncounterSequenceError.invalidInput }
                let shifted = try FlybyEncounterRequest(body: request.body,periapsisEpoch: epoch,frame: request.frame,
                    gravitationalParameter: request.gravitationalParameter,bodyRadiusMetres: request.bodyRadiusMetres,
                    minimumPeriapsisAltitudeMetres: request.minimumPeriapsisAltitudeMetres,
                    incomingExcessVelocityMetresPerSecond: request.incomingExcessVelocityMetresPerSecond,
                    outgoingExcessVelocityMetresPerSecond: request.requestedOutgoingExcessVelocityMetresPerSecond,
                    relativeSpeedTolerance: request.relativeSpeedTolerance)
                let flyby = try UnpoweredFlyby.solve(shifted,propagationOptions: propagationOptions)
                let start = try epoch.advanced(by: Duration(seconds: -lead)), end = try epoch.advanced(by: Duration(seconds: trail))
                guard visit.target.coverage.contains(start), visit.target.coverage.contains(end) else { throw EncounterSequenceError.outsideWindow(index: i) }
                func localState(_ date: Epoch) throws -> CartesianState {
                    try tree.transform(from: request.frame,to: context.frame,at: date)
                        .applying(to: flyby.trajectory.state(at: date),destination: context.frame)
                }
                for date in [start,epoch,end] {
                    let centre = try tree.transform(from: request.frame,to: context.frame,at: date).applying(
                        to: CartesianState(epoch: date,frame: request.frame,position: .zero,velocity: .zero),destination: context.frame)
                    let body = try visit.target.state(at: date)
                    guard routeNorm(centre.position-body.position) <= positionToleranceMetres,
                          routeNorm(centre.velocity-body.velocity) <= velocityToleranceMetresPerSecond else { throw EncounterSequenceError.incompatibleOrigin }
                }
                try connect(localState(start),matchVelocity: true,index: i,identifier: visit.identifier+" entry")
                let entry = try localState(start)
                let offset = current.position-entry.position
                // Retain the propagated patch residual as a constant translation;
                // never reset a carried position to the encounter template.
                arcs.append(EncounterSequenceArc(coverage: start...end,provenance: "Finite local unpowered two-body flyby") { date in
                    let state = try localState(date)
                    return try CartesianState(epoch: date,frame: context.frame,position: state.position+offset,velocity: state.velocity)
                })
                current = try arcs.last!.state(at: end)
                let periTarget = try SequenceClosureProvider(frame: context.frame) { date in
                    let centre = try visit.target.state(at: date)
                    let peri = try localState(epoch), body = try visit.target.state(at: epoch)
                    return try CartesianState(epoch: date,frame: context.frame,position: centre.position+peri.position-body.position,velocity: centre.velocity)
                }
                let target = try MovingWaypointRouteTarget(provider: periTarget,origin: context.origin,coverage: visit.target.coverage)
                requirements.append(try NumericalEncounterRequirement(identifier: visit.identifier,target: target,
                    epoch: epoch,window: visit.arrivalWindow,positionToleranceMetres: positionToleranceMetres,
                    velocityToleranceMetresPerSecond: velocityToleranceMetresPerSecond,closestApproachTarget: visit.target))
                residuals.append(EncounterSequencePatchResidual(identifier: visit.identifier+" exit",epoch: end,
                    positionMetres: routeNorm(offset),preImpulseVelocityMetresPerSecond: 0))
            }
        }
        // DSM candidates are actual commands, not hidden state changes. Their cost
        // is accounted by numerical replay of the complete combined command set.
        let allNodes = nodes+candidate.deepSpaceManeuvers
        do {
            _ = try NumericalEncounterTrajectory(initialState: initialState,arrivalEpoch: current.epoch,nodes: allNodes,
                context: context,backgroundModels: [])
        } catch is NumericalSequenceError { throw EncounterSequenceError.invalidInput }
        if !candidate.deepSpaceManeuvers.isEmpty {
            // Rocket bookkeeping depends only on command order and mass, not on
            // the perturbed trajectory. Actual DSM residuals are checked numerically.
            mass = initialState.mass
            for node in allNodes.sorted(by: { sequenceInterval($0).lowerBound < sequenceInterval($1).lowerBound }) {
                guard case .impulse(let maneuver,let ve) = node.command else { throw EncounterSequenceError.invalidInput }
                let inventory = try SpacecraftState(cartesian: CartesianState(epoch: maneuver.epoch,frame: context.frame,
                    position: initialState.cartesian.position,velocity: .zero),mass: mass)
                mass = try EncounterImpulse.apply(maneuver,to: inventory,effectiveExhaustVelocityMetresPerSecond: ve).mass
            }
            for node in candidate.deepSpaceManeuvers {
                guard case .impulse(let maneuver,_) = node.command else { throw EncounterSequenceError.invalidInput }
                totalDV += routeNorm(maneuver.deltaVelocityMetresPerSecond)
            }
        }
        guard additionalPassages.count+requirements.count <= 100,
              additionalPassages.allSatisfy({ $0.epoch > initialState.epoch && $0.epoch <= current.epoch &&
                  $0.target.frame === context.frame && $0.target.origin == context.origin }) else { throw EncounterSequenceError.invalidInput }
        requirements += additionalPassages
        requirements.sort { $0.epoch < $1.epoch }
        guard Set(requirements.map(\.identifier)).count == requirements.count,
              zip(requirements,requirements.dropFirst()).allSatisfy({ $0.epoch < $1.epoch }) else { throw EncounterSequenceError.invalidInput }
        return EncounterSequenceSeed(initialState: initialState,candidate: candidate,arcs: arcs,maneuvers: allNodes,
            requirements: requirements,patchResiduals: residuals,finalMass: mass,totalDeltaVelocityMetresPerSecond: totalDV)
    }
}
struct SequenceClosureProvider: StateProvider {
    let frame: Frame
    let query: (Epoch) throws -> CartesianState
    init(frame: Frame, query: @escaping (Epoch) throws -> CartesianState) throws { self.frame = frame; self.query = query }
    func state(at epoch: Epoch) throws -> CartesianState { try query(epoch) }
}
