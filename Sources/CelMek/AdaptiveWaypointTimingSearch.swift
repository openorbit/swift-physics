// SPDX-License-Identifier: Apache-2.0
import Math

public struct WaypointTimingRefinementOptions: Equatable, Sendable {
    public let initialSpacingSeconds, minimumSpacingSeconds: Double
    public let maximumRounds, maximumCandidateEvaluations, maximumGeneratedCandidates: Int
    public init(initialSpacingSeconds: Double, minimumSpacingSeconds: Double,
                maximumRounds: Int = 8, maximumCandidateEvaluations: Int = 1000,
                maximumGeneratedCandidates: Int = 100_000) throws {
        guard initialSpacingSeconds.isFinite, minimumSpacingSeconds.isFinite,
              minimumSpacingSeconds > 0, initialSpacingSeconds >= minimumSpacingSeconds else {
            throw WaypointTimingSearchError.invalidCandidateSpacing
        }
        guard (1...32).contains(maximumRounds), (1...10_000).contains(maximumCandidateEvaluations),
              (1...100_000).contains(maximumGeneratedCandidates) else { throw WaypointTimingSearchError.invalidBudget }
        self.initialSpacingSeconds = initialSpacingSeconds; self.minimumSpacingSeconds = minimumSpacingSeconds
        self.maximumRounds = maximumRounds; self.maximumCandidateEvaluations = maximumCandidateEvaluations
        self.maximumGeneratedCandidates = maximumGeneratedCandidates
    }
}
public enum WaypointTimingRefinementStatus: Equatable, Sendable {
    case feasible, resolutionLimitReached, evaluationLimitReached
    case refinementLimitReached, candidateGenerationLimitReached
}
public struct WaypointTimingRefinementRound {
    public let spacingSeconds: Double
    public let generatedCandidates: Int
    public let evaluatedCandidates: Int
    public let status: WaypointTimingSearchStatus
}
public struct WaypointTimingRefinementResult {
    public let status: WaypointTimingRefinementStatus
    public let evaluatedCandidates: Int
    public let rounds: [WaypointTimingRefinementRound]
    public let route: WaypointBrachistochroneRoute?
    public let propagation: WaypointRoutePropagationResult?
    public let failedAttempts: [WaypointTimingAttempt]
}

extension WaypointRouteTimingPlanner {
    /// Refine the grid only after its candidates fail. Halving spacing retains
    /// phase boundaries and avoids replaying previously rejected combinations.
    /// This searches for feasibility, not a minimum-time or continuous optimum.
    public static func adaptiveSearch(departureEpoch: Epoch, frame: Frame, origin: BodyIdentifier,
        waypoints: [BrachistochroneWaypoint], legAccelerationsMetresPerSecondSquared: [Double],
        arrivalWindows: [ClosedRange<Epoch>], refinementOptions: WaypointTimingRefinementOptions,
        initialMass: SpacecraftMassProperties, backgroundDynamics: SpacecraftDynamics,
        legEngineLimits: [WaypointRouteEngineLimits], timedEnginePhases: [WaypointRouteEnginePhase] = [],
        obstacles: [SphericalExclusionRegion] = [],
        numericalOptions: NumericalPropagationOptions = try! NumericalPropagationOptions(),
        endpointPositionToleranceMetres: Double = 0.01,
        endpointVelocityToleranceMetresPerSecond: Double = 1e-5,
        movingTarget: MovingWaypointRouteTarget? = nil,
        movingWaypoints: [Int: MovingWaypointPassage] = [:],
        movingObstacles: [MovingSphericalExclusionRegion] = [],
        physicalConstraints: WaypointRoutePhysicalConstraints = .unconstrained) throws -> WaypointTimingRefinementResult {
        guard (2...10_001).contains(waypoints.count), arrivalWindows.count == waypoints.count-1 else {
            throw WaypointTimingSearchError.invalidCandidates
        }
        _ = try validateTimingConfiguration(departureEpoch: departureEpoch, frame: frame, origin: origin,
            waypoints: waypoints, legAccelerationsMetresPerSecondSquared: legAccelerationsMetresPerSecondSquared,
            backgroundDynamics: backgroundDynamics, legEngineLimits: legEngineLimits,
            timedEnginePhases: timedEnginePhases, obstacles: obstacles,
            endpointPositionToleranceMetres: endpointPositionToleranceMetres,
            endpointVelocityToleranceMetresPerSecond: endpointVelocityToleranceMetresPerSecond, movingTarget: movingTarget, movingWaypoints: movingWaypoints, movingObstacles: movingObstacles, physicalConstraints: physicalConstraints)
        var windows = try intersectPassageCoverage(arrivalWindows, movingWaypoints: movingWaypoints)
        for index in windows.indices {
            if let constraint = waypoints[index+1].arrivalWindow {
                let lower = max(windows[index].lowerBound, constraint.lowerBound)
                let upper = min(windows[index].upperBound, constraint.upperBound)
                guard lower <= upper else { throw WaypointTimingSearchError.invalidCandidates }
                windows[index] = lower...upper
            }
        }
        if let movingTarget {
            let last = windows.count-1
            let lower = max(windows[last].lowerBound, movingTarget.coverage.lowerBound)
            let upper = min(windows[last].upperBound, movingTarget.coverage.upperBound)
            guard lower <= upper else { throw WaypointTimingSearchError.invalidCandidates }
            windows[last] = lower...upper
        }
        var spacing = refinementOptions.initialSpacingSeconds
        var failures: [WaypointTimingAttempt] = [], rounds: [WaypointTimingRefinementRound] = []
        var evaluated = 0, excluded = Set<[Epoch]>()
        func finish(_ status: WaypointTimingRefinementStatus,
                    _ success: WaypointTimingSearchResult? = nil) -> WaypointTimingRefinementResult {
            WaypointTimingRefinementResult(status: status, evaluatedCandidates: evaluated, rounds: rounds,
                route: success?.route, propagation: success?.propagation, failedAttempts: failures)
        }
        for _ in 0..<refinementOptions.maximumRounds {
            let candidates: [[Epoch]]
            do {
                candidates = try arrivalCandidates(in: windows, spacing: Duration(seconds: spacing),
                    including: timedEnginePhases, maximumTotalCandidates: refinementOptions.maximumGeneratedCandidates)
            } catch WaypointTimingSearchError.candidateGenerationLimitExceeded {
                return finish(.candidateGenerationLimitReached)
            }
            let result = try searchImpl(departureEpoch: departureEpoch, frame: frame, origin: origin,
                waypoints: waypoints, legAccelerationsMetresPerSecondSquared: legAccelerationsMetresPerSecondSquared,
                candidateArrivalEpochs: candidates, initialMass: initialMass,
                backgroundDynamics: backgroundDynamics, legEngineLimits: legEngineLimits,
                timedEnginePhases: timedEnginePhases, obstacles: obstacles,
                maximumCandidateEvaluations: refinementOptions.maximumCandidateEvaluations-evaluated,
                numericalOptions: numericalOptions, endpointPositionToleranceMetres: endpointPositionToleranceMetres,
                endpointVelocityToleranceMetresPerSecond: endpointVelocityToleranceMetresPerSecond, excluding: excluded, movingTarget: movingTarget, movingWaypoints: movingWaypoints, movingObstacles: movingObstacles, physicalConstraints: physicalConstraints)
            evaluated += result.evaluatedCandidates
            failures.append(contentsOf: result.failedAttempts)
            for attempt in result.failedAttempts { excluded.insert(attempt.arrivalEpochs) }
            rounds.append(WaypointTimingRefinementRound(spacingSeconds: spacing,
                generatedCandidates: candidates.reduce(0) { $0+$1.count },
                evaluatedCandidates: result.evaluatedCandidates, status: result.status))
            if result.status == .feasible { return finish(.feasible, result) }
            if result.status == .evaluationLimitReached { return finish(.evaluationLimitReached) }
            if spacing <= refinementOptions.minimumSpacingSeconds || windows.allSatisfy({ $0.lowerBound == $0.upperBound }) {
                return finish(.resolutionLimitReached)
            }
            if evaluated == refinementOptions.maximumCandidateEvaluations { return finish(.evaluationLimitReached) }
            spacing = max(refinementOptions.minimumSpacingSeconds, spacing/2)
        }
        return finish(.refinementLimitReached)
    }
}
