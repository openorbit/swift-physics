// SPDX-License-Identifier: Apache-2.0
import Math

public enum WaypointTimingSearchError: Error, Equatable {
    case invalidCandidates, invalidBudget, invalidCandidateSpacing
    case candidateGenerationLimitExceeded, epochResolutionExceeded
}
public enum WaypointTimingSearchStatus: Equatable, Sendable {
    case feasible, candidateSetExhausted, evaluationLimitReached
}
public enum WaypointTimingCandidateFailure: Equatable {
    case movingWaypointOutsideCoverage(index: Int)
    case targetOutsideCoverage
    case invalidOrdering
    case geometry(BrachistochroneRouteError)
    case engine(WaypointRoutePropagationError)
    case integration(NumericalPropagationError)
}
public struct WaypointTimingAttempt {
    public let arrivalEpochs: [Epoch]
    public let failure: WaypointTimingCandidateFailure
}
public struct WaypointTimingSearchResult {
    public let status: WaypointTimingSearchStatus
    public let evaluatedCandidates: Int
    public let route: WaypointBrachistochroneRoute?
    public let propagation: WaypointRoutePropagationResult?
    public let failedAttempts: [WaypointTimingAttempt]
}

/// Bounded lexicographic search of explicit arrival candidates. All spatial
/// passages retain their supplied order, while their order relative to engine
/// phase boundaries is determined by each candidate. No continuous optimum or
/// global infeasibility claim follows from exhausting this discrete set.
public enum WaypointRouteTimingPlanner {
    /// Generate bounded candidates for each ordered spatial passage. Endpoints
    /// and engine-phase boundaries supplement the regular grid; exact windows
    /// yield one candidate. This does not enumerate their Cartesian product.
    public static func arrivalCandidates(in windows: [ClosedRange<Epoch>], spacing: Duration,
        including timedEnginePhases: [WaypointRouteEnginePhase] = [],
        maximumTotalCandidates: Int = 100_000) throws -> [[Epoch]] {
        guard spacing > .zero, spacing.seconds.isFinite, spacing.seconds > 0 else {
            throw WaypointTimingSearchError.invalidCandidateSpacing
        }
        guard (1...100_000).contains(maximumTotalCandidates) else { throw WaypointTimingSearchError.invalidBudget }
        guard !windows.isEmpty, windows.count <= 10_000, timedEnginePhases.count <= 10_000 else {
            throw WaypointTimingSearchError.invalidCandidates
        }
        let sortedBoundaries = timedEnginePhases.flatMap { [$0.startEpoch, $0.endEpoch] }.sorted()
        let boundaries = sortedBoundaries.enumerated().compactMap { index, epoch in
            index == 0 || sortedBoundaries[index-1] != epoch ? epoch : nil
        }
        var result: [[Epoch]] = [], count = 0
        for window in windows {
            var epochs = [window.lowerBound], current = window.lowerBound
            guard count+epochs.count <= maximumTotalCandidates else {
                throw WaypointTimingSearchError.candidateGenerationLimitExceeded
            }
            while current < window.upperBound {
                let remaining = try window.upperBound.duration(since: current)
                let next = spacing >= remaining ? window.upperBound : try current.advanced(by: spacing)
                guard next > current else { throw WaypointTimingSearchError.epochResolutionExceeded }
                epochs.append(next); current = next
                guard count+epochs.count <= maximumTotalCandidates else {
                    throw WaypointTimingSearchError.candidateGenerationLimitExceeded
                }
            }
            // Find only the temporal boundaries in this window, avoiding a full
            // scan of every engine phase for every spatial waypoint.
            var low = 0, high = boundaries.count
            while low < high {
                let middle = low+(high-low)/2
                if boundaries[middle] < window.lowerBound { low = middle+1 } else { high = middle }
            }
            var boundaryIndex = low
            while boundaryIndex < boundaries.count, boundaries[boundaryIndex] <= window.upperBound {
                epochs.append(boundaries[boundaryIndex]); boundaryIndex += 1
            }
            let sorted = epochs.sorted()
            let unique = sorted.enumerated().compactMap { index, epoch in
                index == 0 || sorted[index-1] != epoch ? epoch : nil
            }
            count += unique.count
            guard count <= maximumTotalCandidates else {
                throw WaypointTimingSearchError.candidateGenerationLimitExceeded
            }
            result.append(unique)
        }
        return result
    }

    public static func search(departureEpoch: Epoch, frame: Frame, origin: BodyIdentifier,
        waypoints: [BrachistochroneWaypoint], legAccelerationsMetresPerSecondSquared: [Double],
        candidateArrivalEpochs: [[Epoch]], initialMass: SpacecraftMassProperties,
        backgroundDynamics: SpacecraftDynamics, legEngineLimits: [WaypointRouteEngineLimits],
        timedEnginePhases: [WaypointRouteEnginePhase] = [], obstacles: [SphericalExclusionRegion] = [],
        maximumCandidateEvaluations: Int = 1000,
        numericalOptions: NumericalPropagationOptions = try! NumericalPropagationOptions(),
        endpointPositionToleranceMetres: Double = 0.01,
        endpointVelocityToleranceMetresPerSecond: Double = 1e-5,
        movingTarget: MovingWaypointRouteTarget? = nil,
        movingWaypoints: [Int: MovingWaypointPassage] = [:],
        movingObstacles: [MovingSphericalExclusionRegion] = []) throws -> WaypointTimingSearchResult {
        try searchImpl(departureEpoch: departureEpoch, frame: frame, origin: origin, waypoints: waypoints,
            legAccelerationsMetresPerSecondSquared: legAccelerationsMetresPerSecondSquared,
            candidateArrivalEpochs: candidateArrivalEpochs, initialMass: initialMass,
            backgroundDynamics: backgroundDynamics, legEngineLimits: legEngineLimits,
            timedEnginePhases: timedEnginePhases, obstacles: obstacles,
            maximumCandidateEvaluations: maximumCandidateEvaluations, numericalOptions: numericalOptions,
            endpointPositionToleranceMetres: endpointPositionToleranceMetres,
            endpointVelocityToleranceMetresPerSecond: endpointVelocityToleranceMetresPerSecond,
            excluding: [], movingTarget: movingTarget, movingWaypoints: movingWaypoints, movingObstacles: movingObstacles)
    }

    static func searchImpl(departureEpoch: Epoch, frame: Frame, origin: BodyIdentifier,
        waypoints: [BrachistochroneWaypoint], legAccelerationsMetresPerSecondSquared: [Double],
        candidateArrivalEpochs: [[Epoch]], initialMass: SpacecraftMassProperties,
        backgroundDynamics: SpacecraftDynamics, legEngineLimits: [WaypointRouteEngineLimits],
        timedEnginePhases: [WaypointRouteEnginePhase], obstacles: [SphericalExclusionRegion],
        maximumCandidateEvaluations: Int, numericalOptions: NumericalPropagationOptions,
        endpointPositionToleranceMetres: Double, endpointVelocityToleranceMetresPerSecond: Double,
        excluding: Set<[Epoch]>, movingTarget: MovingWaypointRouteTarget? = nil,
        movingWaypoints: [Int: MovingWaypointPassage] = [:],
        movingObstacles: [MovingSphericalExclusionRegion] = []) throws -> WaypointTimingSearchResult {
        guard (1...10_000).contains(maximumCandidateEvaluations) else { throw WaypointTimingSearchError.invalidBudget }
        guard (2...10_001).contains(waypoints.count), candidateArrivalEpochs.count == waypoints.count-1,
              candidateArrivalEpochs.allSatisfy({ !$0.isEmpty }),
              candidateArrivalEpochs.reduce(0, { $0+$1.count }) <= 100_000 else {
            throw WaypointTimingSearchError.invalidCandidates
        }
        let sortedPhases = try validateTimingConfiguration(departureEpoch: departureEpoch, frame: frame, origin: origin,
            waypoints: waypoints, legAccelerationsMetresPerSecondSquared: legAccelerationsMetresPerSecondSquared,
            backgroundDynamics: backgroundDynamics, legEngineLimits: legEngineLimits,
            timedEnginePhases: timedEnginePhases, obstacles: obstacles,
            endpointPositionToleranceMetres: endpointPositionToleranceMetres,
            endpointVelocityToleranceMetresPerSecond: endpointVelocityToleranceMetresPerSecond, movingTarget: movingTarget, movingWaypoints: movingWaypoints, movingObstacles: movingObstacles)
        let candidates = candidateArrivalEpochs.map { values -> [Epoch] in
            let sorted = values.sorted()
            return sorted.enumerated().compactMap { index, epoch in
                index == 0 || sorted[index-1] != epoch ? epoch : nil
            }
        }
        var indices = Array(repeating: 0, count: candidates.count)
        var failures: [WaypointTimingAttempt] = []
        var evaluated = 0
        while true {
            let arrivals = candidates.indices.map { candidates[$0][indices[$0]] }
            if !excluding.contains(arrivals) {
                evaluated += 1
                do {
                    var previous = departureEpoch, durations: [Duration] = []
                    var ordered = true
                    for arrival in arrivals {
                        if arrival <= previous { ordered = false; break }
                        durations.append(try arrival.duration(since: previous)); previous = arrival
                    }
                    if !ordered {
                        failures.append(WaypointTimingAttempt(arrivalEpochs: arrivals, failure: .invalidOrdering))
                    } else {
                        var resolved = waypoints
                        for index in movingWaypoints.keys.sorted() {
                            let passage = movingWaypoints[index]!
                            let epoch = arrivals[index-1]
                            guard passage.provider.coverage.contains(epoch) else {
                                throw MovingWaypointTargetError.waypointOutsideCoverage(index: index)
                            }
                            resolved[index] = try passage.waypoint(at: epoch, template: waypoints[index])
                        }
                        if let movingTarget {
                            let target = try movingTarget.state(at: arrivals.last!)
                            resolved[resolved.count-1] = try BrachistochroneWaypoint(positionMetres: target.position,
                                arrivalWindow: waypoints.last!.arrivalWindow, velocityMetresPerSecond: target.velocity)
                        }
                        let route = try WaypointBrachistochroneRoute(departureEpoch: departureEpoch, frame: frame, origin: origin,
                            waypoints: resolved, legAccelerationsMetresPerSecondSquared: legAccelerationsMetresPerSecondSquared,
                            obstacles: obstacles, legDurations: durations, movingObstacles: movingObstacles)
                        let execution = try WaypointRoutePropagator(route: route, initialMass: initialMass,
                            backgroundDynamics: backgroundDynamics, legEngineLimits: legEngineLimits,
                            timedEnginePhases: sortedPhases, numericalOptions: numericalOptions,
                            endpointPositionToleranceMetres: endpointPositionToleranceMetres,
                            endpointVelocityToleranceMetresPerSecond: endpointVelocityToleranceMetresPerSecond)
                        let result = try execution.propagation(at: route.arrivalEpoch)
                        return WaypointTimingSearchResult(status: .feasible, evaluatedCandidates: evaluated,
                            route: route, propagation: result, failedAttempts: failures)
                    }
                } catch MovingWaypointTargetError.waypointOutsideCoverage(let index) {
                    failures.append(WaypointTimingAttempt(arrivalEpochs: arrivals, failure: .movingWaypointOutsideCoverage(index: index)))
                } catch MovingWaypointTargetError.outsideCoverage {
                    failures.append(WaypointTimingAttempt(arrivalEpochs: arrivals, failure: .targetOutsideCoverage))
                } catch let error as BrachistochroneRouteError {
                    failures.append(WaypointTimingAttempt(arrivalEpochs: arrivals, failure: .geometry(error)))
                } catch let error as WaypointRoutePropagationError {
                    // A shorter candidate can fail to cover a scheduled phase. Other
                    // static errors and background mass flow must remain caller errors.
                    switch error {
                    case .accelerationLimitExceeded, .thrustLimitExceeded, .endpointToleranceExceeded, .timedPhaseOutsideCoverage:
                        failures.append(WaypointTimingAttempt(arrivalEpochs: arrivals, failure: .engine(error)))
                    default: throw error
                    }
                } catch let error as NumericalPropagationError {
                    failures.append(WaypointTimingAttempt(arrivalEpochs: arrivals, failure: .integration(error)))
                }
            }
            var position = indices.count-1
            while indices[position]+1 == candidates[position].count {
                indices[position] = 0
                if position == 0 {
                    return WaypointTimingSearchResult(status: .candidateSetExhausted, evaluatedCandidates: evaluated,
                        route: nil, propagation: nil, failedAttempts: failures)
                }
                position -= 1
            }
            indices[position] += 1
            if evaluated == maximumCandidateEvaluations {
                return WaypointTimingSearchResult(status: .evaluationLimitReached, evaluatedCandidates: evaluated,
                    route: nil, propagation: nil, failedAttempts: failures)
            }
        }
    }
    static func validateTimingConfiguration(departureEpoch: Epoch, frame: Frame, origin: BodyIdentifier,
        waypoints: [BrachistochroneWaypoint], legAccelerationsMetresPerSecondSquared: [Double],
        backgroundDynamics: SpacecraftDynamics, legEngineLimits: [WaypointRouteEngineLimits],
        timedEnginePhases: [WaypointRouteEnginePhase], obstacles: [SphericalExclusionRegion],
        endpointPositionToleranceMetres: Double, endpointVelocityToleranceMetresPerSecond: Double,
        movingTarget: MovingWaypointRouteTarget? = nil,
        movingWaypoints: [Int: MovingWaypointPassage] = [:],
        movingObstacles: [MovingSphericalExclusionRegion] = []
    ) throws -> [WaypointRouteEnginePhase] {
        try validateMovingObstacles(movingObstacles, staticObstacles: obstacles, frame: frame, origin: origin)
        // Validate frame, budget and engine configuration before searching.
        guard frame === backgroundDynamics.context.frame else { throw SpacecraftDynamicsError.frameMismatch }
        guard origin == backgroundDynamics.context.origin else { throw BrachistochroneRouteError.originMismatch }
        guard frame.pseudoInertial else { throw SpacecraftDynamicsError.nonInertialFrame }
        for index in movingWaypoints.keys.sorted() {
            guard index > 0, index < waypoints.count-1 else { throw MovingWaypointTargetError.invalidWaypointIndex }
            let passage = movingWaypoints[index]!
            guard passage.provider.frame === frame else { throw MovingWaypointTargetError.frameMismatch }
            guard passage.provider.origin == origin else { throw MovingWaypointTargetError.originMismatch }
            if passage.velocityMode == .matchProviderVelocity, waypoints[index].velocityMetresPerSecond != nil {
                throw MovingWaypointTargetError.conflictingWaypointVelocity(index: index)
            }
        }
        if let movingTarget {
            guard movingTarget.frame === frame else { throw MovingWaypointTargetError.frameMismatch }
            guard movingTarget.origin == origin else { throw MovingWaypointTargetError.originMismatch }
            guard waypoints.last!.velocityMetresPerSecond == nil else { throw MovingWaypointTargetError.conflictingEndpointVelocity }
        }
        guard legEngineLimits.count == waypoints.count-1,
              legAccelerationsMetresPerSecondSquared.count == waypoints.count-1 else {
            throw BrachistochroneRouteError.legCountMismatch
        }
        guard legAccelerationsMetresPerSecondSquared.allSatisfy({ $0.isFinite && $0 > 0 }),
              [endpointPositionToleranceMetres, endpointVelocityToleranceMetresPerSecond]
                .allSatisfy({ $0.isFinite && $0 > 0 }) else { throw WaypointRoutePropagationError.invalidLimits }
        guard waypoints.allSatisfy({ $0.dwellDuration == .zero }) else { throw WaypointRoutePropagationError.unsupportedDwell }
        guard obstacles.count <= 10_000, Set(obstacles.map(\.identifier)).count == obstacles.count else {
            throw BrachistochroneRouteError.invalidObstacle
        }
        guard obstacles.allSatisfy({ $0.frame === frame }) else { throw SpacecraftDynamicsError.frameMismatch }
        guard obstacles.allSatisfy({ $0.origin == origin }) else { throw BrachistochroneRouteError.originMismatch }
        if let window = waypoints[0].arrivalWindow, !window.contains(departureEpoch) {
            throw BrachistochroneRouteError.departureOutsideArrivalWindow
        }
        let sortedPhases = timedEnginePhases.sorted { $0.startEpoch < $1.startEpoch }
        guard sortedPhases.count <= 10_000 else { throw WaypointRoutePropagationError.invalidTimedPhase }
        for (index, phase) in sortedPhases.enumerated() {
            guard phase.startEpoch >= departureEpoch else { throw WaypointRoutePropagationError.timedPhaseOutsideCoverage }
            guard index == 0 || sortedPhases[index-1].endEpoch <= phase.startEpoch else {
                throw WaypointRoutePropagationError.overlappingTimedPhases
            }
        }
        return sortedPhases
    }

}
