// SPDX-License-Identifier: Apache-2.0
import Math

public enum LambertWindowSearchError: Error, Equatable { case invalidGrid, invalidBudget, invalidConstraints, invalidPole }
public enum LambertWindowObjective: Equatable, Sendable { case minimumDeltaVelocity, minimumFlightTime }
public enum LambertWindowSearchStatus: Equatable, Sendable { case gridExhausted, evaluationLimitReached }
public struct LambertWindowBranch: Equatable {
    public let direction: LambertDirection
    public let revolutions: Int
    public let multiRevolutionBranch: LambertMultiRevolutionBranch?
    public init(direction: LambertDirection, revolutions: Int = 0,
                multiRevolutionBranch: LambertMultiRevolutionBranch? = nil) throws {
        guard (0...32).contains(revolutions) else { throw LambertError.unsupportedRevolutions }
        guard (revolutions == 0) == (multiRevolutionBranch == nil) else { throw LambertError.invalidBranch }
        self.direction = direction; self.revolutions = revolutions; self.multiRevolutionBranch = multiRevolutionBranch
    }
}
public struct LambertWindowCandidate {
    public let departureEpoch: Epoch
    public let flightTime: Duration
    public let branch: LambertWindowBranch
}
public enum LambertWindowFailure: Equatable {
    case departureOutsideCoverage, arrivalOutsideCoverage, arrivalEpochResolutionExceeded
    case timeArithmetic(TimeArithmeticError)
    case lambert(LambertError), propagation(PropagationError), planning(LambertTransferPlanningError)
    case endpointResidual, deltaVelocityLimit, arrivalRelativeSpeedLimit
}
/// One table row per visited date/time/branch combination, including failed cases.
public struct LambertWindowRow {
    public let candidate: LambertWindowCandidate
    public let arrivalEpoch: Epoch?
    public let plan: LambertTransferPlan?
    public let failure: LambertWindowFailure?
    public var totalDeltaVelocityMetresPerSecond: Double? { plan?.totalDeltaVelocityMetresPerSecond }
    public var departureRelativeSpeedMetresPerSecond: Double? { plan.map { routeNorm($0.departureRelativeVelocityMetresPerSecond) } }
    public var arrivalRelativeSpeedMetresPerSecond: Double? { plan.map { routeNorm($0.incomingRelativeVelocityMetresPerSecond) } }
    public var feasible: Bool { plan != nil && failure == nil }
}
public struct LambertWindowSearchResult {
    public let status: LambertWindowSearchStatus
    public let totalCandidates: Int
    public let rows: [LambertWindowRow]
    public let rankedFeasibleRows: [LambertWindowRow]
    public var evaluatedCandidates: Int { rows.count }
    public var best: LambertWindowRow? { rankedFeasibleRows.first }
}

/// Discrete grid exploration; ranking is over evaluated feasible samples only.
/// No global optimum or continuous transfer-window coverage is implied.
public enum LambertTransferWindowSearch {
    public static func search(departure: MovingWaypointRouteTarget, arrival: MovingWaypointRouteTarget,
        centralBody: BodyIdentifier, gravitationalParameter: GravitationalParameter, referencePole: SIMD3<Double>,
        departureEpochs: [Epoch], flightTimes: [Duration], branches: [LambertWindowBranch],
        arrivalCondition: LambertArrivalCondition, objective: LambertWindowObjective = .minimumDeltaVelocity,
        maximumCandidateEvaluations: Int = 1000, maximumDeltaVelocityMetresPerSecond: Double? = nil,
        maximumArrivalRelativeSpeedMetresPerSecond: Double? = nil,
        positionToleranceMetres: Double = 0.01, velocityToleranceMetresPerSecond: Double = 1e-5,
        lambertOptions: LambertOptions = .default, propagationOptions: TwoBodyPropagationOptions = .default) throws -> LambertWindowSearchResult {
        guard !departureEpochs.isEmpty, !flightTimes.isEmpty, !branches.isEmpty,
              departureEpochs.count <= 100_000, flightTimes.count <= 100_000, branches.count <= 256,
              departureEpochs.count <= 100_000 / flightTimes.count,
              departureEpochs.count * flightTimes.count <= 100_000 / branches.count,
              flightTimes.allSatisfy({ $0 > .zero && $0.seconds.isFinite }) else { throw LambertWindowSearchError.invalidGrid }
        guard (1...10_000).contains(maximumCandidateEvaluations) else { throw LambertWindowSearchError.invalidBudget }
        guard [positionToleranceMetres,velocityToleranceMetresPerSecond].allSatisfy({ $0.isFinite && $0 > 0 }),
              [maximumDeltaVelocityMetresPerSecond,maximumArrivalRelativeSpeedMetresPerSecond].allSatisfy({
                  $0 == nil || ($0!.isFinite && $0! >= 0)
              }) else { throw LambertWindowSearchError.invalidConstraints }
        guard [referencePole.x,referencePole.y,referencePole.z].allSatisfy(\.isFinite),
              routeNorm(referencePole).isFinite, routeNorm(referencePole) > 0 else { throw LambertWindowSearchError.invalidPole }
        guard departure.frame === arrival.frame else { throw SpacecraftDynamicsError.frameMismatch }
        guard departure.frame.pseudoInertial else { throw PropagationError.nonInertialFrame }
        guard departure.origin == centralBody, arrival.origin == centralBody else { throw LambertTransferPlanningError.originMismatch }
        try centralBody.validate()
        guard gravitationalParameter.convention == .newtonian || gravitationalParameter.convention == .ttCompatible else {
            throw PropagationError.incompatibleGMConvention
        }
        let epochs = Array(Set(departureEpochs)).sorted(), times = Array(Set(flightTimes)).sorted()
        var uniqueBranches: [LambertWindowBranch] = []
        for branch in branches where !uniqueBranches.contains(branch) { uniqueBranches.append(branch) }
        let total = epochs.count * times.count * uniqueBranches.count
        var rows: [LambertWindowRow] = []
        func finish(_ status: LambertWindowSearchStatus) -> LambertWindowSearchResult {
            let ranked = rows.enumerated().filter { $0.element.feasible }.sorted { a,b in
                let ac = a.element.totalDeltaVelocityMetresPerSecond!, bc = b.element.totalDeltaVelocityMetresPerSecond!
                let at = a.element.candidate.flightTime, bt = b.element.candidate.flightTime
                switch objective {
                case .minimumDeltaVelocity: if ac != bc { return ac < bc }; if at != bt { return at < bt }
                case .minimumFlightTime: if at != bt { return at < bt }; if ac != bc { return ac < bc }
                }
                if a.element.candidate.departureEpoch != b.element.candidate.departureEpoch {
                    return a.element.candidate.departureEpoch < b.element.candidate.departureEpoch
                }
                return a.offset < b.offset
            }.map(\.element)
            return LambertWindowSearchResult(status: status,totalCandidates: total,rows: rows,rankedFeasibleRows: ranked)
        }
        for epoch in epochs {
            for time in times {
                // Cache only this date/time pair's checked source states across roots.
                var departureState: CartesianState?, arrivalState: CartesianState?
                for branch in uniqueBranches {
                    if rows.count == maximumCandidateEvaluations { return finish(.evaluationLimitReached) }
                    let candidate = LambertWindowCandidate(departureEpoch: epoch,flightTime: time,branch: branch)
                    let end: Epoch
                    do { end = try epoch.advanced(by: time) }
                    catch let error as TimeArithmeticError {
                        rows.append(LambertWindowRow(candidate: candidate,arrivalEpoch: nil,plan: nil,failure: .timeArithmetic(error)))
                        continue
                    }
                    var failure: LambertWindowFailure?, plan: LambertTransferPlan?
                    if end <= epoch { failure = .arrivalEpochResolutionExceeded }
                    else if !departure.coverage.contains(epoch) { failure = .departureOutsideCoverage }
                    else if !arrival.coverage.contains(end) { failure = .arrivalOutsideCoverage }
                    else {
                        // Provider errors are contract/data failures, not Lambert candidate failures.
                        if departureState == nil { departureState = try departure.state(at: epoch) }
                        if arrivalState == nil { arrivalState = try arrival.state(at: end) }
                        do {
                            plan = try LambertTransferPlanner.plan(departureState: departureState!,arrivalTargetState: arrivalState!,
                                centralBody: centralBody,gravitationalParameter: gravitationalParameter,referencePole: referencePole,
                                direction: branch.direction,arrivalCondition: arrivalCondition,revolutions: branch.revolutions,
                                multiRevolutionBranch: branch.multiRevolutionBranch,positionToleranceMetres: positionToleranceMetres,
                                velocityToleranceMetresPerSecond: velocityToleranceMetresPerSecond,lambertOptions: lambertOptions,
                                propagationOptions: propagationOptions)
                            if !plan!.arrivalConditionSatisfied { failure = .endpointResidual }
                            else if let limit = maximumDeltaVelocityMetresPerSecond,
                                plan!.totalDeltaVelocityMetresPerSecond > limit { failure = .deltaVelocityLimit }
                            else if let limit = maximumArrivalRelativeSpeedMetresPerSecond,
                                routeNorm(plan!.incomingRelativeVelocityMetresPerSecond) > limit { failure = .arrivalRelativeSpeedLimit }
                        } catch let error as LambertError { failure = .lambert(error) }
                        catch let error as PropagationError { failure = .propagation(error) }
                        catch let error as LambertTransferPlanningError { failure = .planning(error) }
                    }
                    rows.append(LambertWindowRow(candidate: candidate,arrivalEpoch: end,plan: plan,failure: failure))
                }
            }
        }
        return finish(.gridExhausted)
    }

    public static func search(departure: any EphemerisProvider, arrival: any EphemerisProvider,
        centralBody: BodyIdentifier, gravitationalParameter: GravitationalParameter, referencePole: SIMD3<Double>,
        departureEpochs: [Epoch], flightTimes: [Duration], branches: [LambertWindowBranch],
        arrivalCondition: LambertArrivalCondition, objective: LambertWindowObjective = .minimumDeltaVelocity,
        maximumCandidateEvaluations: Int = 1000, maximumDeltaVelocityMetresPerSecond: Double? = nil,
        maximumArrivalRelativeSpeedMetresPerSecond: Double? = nil,
        positionToleranceMetres: Double = 0.01, velocityToleranceMetresPerSecond: Double = 1e-5,
        lambertOptions: LambertOptions = .default, propagationOptions: TwoBodyPropagationOptions = .default) throws -> LambertWindowSearchResult {
        try search(departure: MovingWaypointRouteTarget(ephemeris: departure),arrival: MovingWaypointRouteTarget(ephemeris: arrival),
            centralBody: centralBody,gravitationalParameter: gravitationalParameter,referencePole: referencePole,
            departureEpochs: departureEpochs,flightTimes: flightTimes,branches: branches,arrivalCondition: arrivalCondition,
            objective: objective,maximumCandidateEvaluations: maximumCandidateEvaluations,
            maximumDeltaVelocityMetresPerSecond: maximumDeltaVelocityMetresPerSecond,
            maximumArrivalRelativeSpeedMetresPerSecond: maximumArrivalRelativeSpeedMetresPerSecond,
            positionToleranceMetres: positionToleranceMetres,velocityToleranceMetresPerSecond: velocityToleranceMetresPerSecond,
            lambertOptions: lambertOptions,propagationOptions: propagationOptions)
    }
}
