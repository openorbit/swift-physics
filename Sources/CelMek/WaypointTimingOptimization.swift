// SPDX-License-Identifier: Apache-2.0
import Math

public struct WaypointTimingOptimizationOptions: Equatable, Sendable {
    public let initialStepSeconds, minimumStepSeconds: Double
    public let maximumIterations, maximumCandidateEvaluations: Int
    public init(initialStepSeconds: Double, minimumStepSeconds: Double,
                maximumIterations: Int = 256, maximumCandidateEvaluations: Int = 1000) throws {
        guard initialStepSeconds.isFinite, minimumStepSeconds.isFinite,
              minimumStepSeconds > 0, initialStepSeconds >= minimumStepSeconds,
              (try? Duration(seconds: initialStepSeconds)) != nil,
              (try? Duration(seconds: minimumStepSeconds)) != nil else {
            throw WaypointTimingSearchError.invalidCandidateSpacing
        }
        guard (1...1000).contains(maximumIterations), (1...10_000).contains(maximumCandidateEvaluations) else {
            throw WaypointTimingSearchError.invalidBudget
        }
        self.initialStepSeconds = initialStepSeconds; self.minimumStepSeconds = minimumStepSeconds
        self.maximumIterations = maximumIterations; self.maximumCandidateEvaluations = maximumCandidateEvaluations
    }
}
public enum WaypointTimingOptimizationStatus: Equatable, Sendable {
    case seedRejected, stepLimitReached, evaluationLimitReached, iterationLimitReached
}
public struct WaypointTimingOptimizationIteration {
    public let stepSeconds: Double
    public let evaluatedCandidates: Int
    public let improved: Bool
    public let bestArrivalEpoch: Epoch
}
public struct WaypointTimingOptimizationResult {
    public let status: WaypointTimingOptimizationStatus
    public let evaluatedCandidates: Int
    public let iterations: [WaypointTimingOptimizationIteration]
    public let route: WaypointBrachistochroneRoute?
    public let propagation: WaypointRoutePropagationResult?
    public let failedAttempts: [WaypointTimingAttempt]
}

extension WaypointRouteTimingPlanner {
    /// Feasible-seed direct search for an earlier final arrival. Only strict
    /// improvements are accepted. The tested pattern and step size do not prove
    /// stationarity, global optimality or continuous constraint feasibility.
    public static func optimizeTiming(departureEpoch: Epoch, frame: Frame, origin: BodyIdentifier,
        waypoints: [BrachistochroneWaypoint], legAccelerationsMetresPerSecondSquared: [Double],
        arrivalWindows: [ClosedRange<Epoch>], seedArrivalEpochs: [Epoch],
        optimizationOptions: WaypointTimingOptimizationOptions, initialMass: SpacecraftMassProperties,
        backgroundDynamics: SpacecraftDynamics, legEngineLimits: [WaypointRouteEngineLimits],
        timedEnginePhases: [WaypointRouteEnginePhase] = [], obstacles: [SphericalExclusionRegion] = [],
        numericalOptions: NumericalPropagationOptions = try! NumericalPropagationOptions(),
        endpointPositionToleranceMetres: Double = 0.01,
        endpointVelocityToleranceMetresPerSecond: Double = 1e-5,
        movingTarget: MovingWaypointRouteTarget? = nil,
        movingWaypoints: [Int: MovingWaypointPassage] = [:],
        movingObstacles: [MovingSphericalExclusionRegion] = [],
        physicalConstraints: WaypointRoutePhysicalConstraints = .unconstrained) throws -> WaypointTimingOptimizationResult {
        guard (2...10_001).contains(waypoints.count), arrivalWindows.count == waypoints.count-1,
              seedArrivalEpochs.count == arrivalWindows.count else { throw WaypointTimingSearchError.invalidCandidates }
        _ = try validateTimingConfiguration(departureEpoch: departureEpoch, frame: frame, origin: origin,
            waypoints: waypoints, legAccelerationsMetresPerSecondSquared: legAccelerationsMetresPerSecondSquared,
            backgroundDynamics: backgroundDynamics, legEngineLimits: legEngineLimits,
            timedEnginePhases: timedEnginePhases, obstacles: obstacles,
            endpointPositionToleranceMetres: endpointPositionToleranceMetres,
            endpointVelocityToleranceMetresPerSecond: endpointVelocityToleranceMetresPerSecond, movingTarget: movingTarget, movingWaypoints: movingWaypoints, movingObstacles: movingObstacles, physicalConstraints: physicalConstraints)
        var windows = try intersectPassageCoverage(arrivalWindows, movingWaypoints: movingWaypoints)
        for index in windows.indices {
            if let constraint = waypoints[index+1].arrivalWindow {
                let low = max(windows[index].lowerBound, constraint.lowerBound)
                let high = min(windows[index].upperBound, constraint.upperBound)
                guard low <= high else { throw WaypointTimingSearchError.invalidCandidates }
                windows[index] = low...high
            }
            if index == windows.count-1, let movingTarget {
                let lower = max(windows[index].lowerBound, movingTarget.coverage.lowerBound)
                let upper = min(windows[index].upperBound, movingTarget.coverage.upperBound)
                guard lower <= upper else { throw WaypointTimingSearchError.invalidCandidates }
                windows[index] = lower...upper
            }
            guard windows[index].contains(seedArrivalEpochs[index]) else { throw WaypointTimingSearchError.invalidCandidates }
        }
        var evaluated = 0, tried = Set<[Epoch]>()
        var failures: [WaypointTimingAttempt] = [], iterations: [WaypointTimingOptimizationIteration] = []
        func evaluate(_ epochs: [Epoch]) throws -> WaypointTimingSearchResult {
            tried.insert(epochs); evaluated += 1
            let result = try searchImpl(departureEpoch: departureEpoch, frame: frame, origin: origin, waypoints: waypoints,
                legAccelerationsMetresPerSecondSquared: legAccelerationsMetresPerSecondSquared,
                candidateArrivalEpochs: epochs.map { [$0] }, initialMass: initialMass,
                backgroundDynamics: backgroundDynamics, legEngineLimits: legEngineLimits,
                timedEnginePhases: timedEnginePhases, obstacles: obstacles, maximumCandidateEvaluations: 1,
                numericalOptions: numericalOptions, endpointPositionToleranceMetres: endpointPositionToleranceMetres,
                endpointVelocityToleranceMetresPerSecond: endpointVelocityToleranceMetresPerSecond, excluding: [], movingTarget: movingTarget, movingWaypoints: movingWaypoints, movingObstacles: movingObstacles, physicalConstraints: physicalConstraints)
            failures.append(contentsOf: result.failedAttempts)
            return result
        }
        var best = try evaluate(seedArrivalEpochs)
        var bestEpochs = seedArrivalEpochs
        func finish(_ status: WaypointTimingOptimizationStatus) -> WaypointTimingOptimizationResult {
            WaypointTimingOptimizationResult(status: status, evaluatedCandidates: evaluated, iterations: iterations,
                route: best.route, propagation: best.propagation, failedAttempts: failures)
        }
        guard best.status == .feasible else { return finish(.seedRejected) }
        var step = optimizationOptions.initialStepSeconds
        let last = bestEpochs.count-1
        for _ in 0..<optimizationOptions.maximumIterations {
            let before = evaluated
            var improved = false, exhaustedBudget = false
            let delta = try Duration(seconds: step)
            let available = try bestEpochs[last].duration(since: windows[last].lowerBound)
            let final = delta >= available ? windows[last].lowerBound : try bestEpochs[last].advanced(by: Duration(seconds: -step))
            if final < bestEpochs[last], final > departureEpoch {
                // Final-only, proportional compression, then +/- one intermediate
                // passage paired with an earlier final arrival. Generate lazily.
                for pattern in 0..<(4+6*last) {
                    var proposal = bestEpochs
                    proposal[last] = final
                    if pattern == 1 {
                        let ratio = try final.duration(since: departureEpoch).seconds /
                            bestEpochs[last].duration(since: departureEpoch).seconds
                        for index in 0..<last {
                            proposal[index] = try departureEpoch.advanced(by: Duration(seconds:
                                bestEpochs[index].duration(since: departureEpoch).seconds*ratio))
                        }
                    } else if pattern >= 2 && pattern < 2+2*last {
                        let index = (pattern-2)/2
                        if pattern.isMultiple(of: 2) {
                            let available = try proposal[index].duration(since: windows[index].lowerBound)
                            proposal[index] = delta >= available ? windows[index].lowerBound :
                                try proposal[index].advanced(by: Duration(seconds: -step))
                        } else {
                            let available = try windows[index].upperBound.duration(since: proposal[index])
                            proposal[index] = delta >= available ? windows[index].upperBound :
                                try proposal[index].advanced(by: delta)
                        }
                    }
                    if pattern >= 2+2*last {
                        // Joint shifts of all intermediates, then each prefix/suffix.
                        let joint = pattern-(2+2*last)
                        let direction = joint.isMultiple(of: 2) ? -step : step
                        let selection: Range<Int>
                        if joint < 2 { selection = 0..<last }
                        else {
                            let group = (joint-2)/4
                            selection = ((joint-2)%4 < 2) ? 0..<(group+1) : group..<last
                        }
                        for index in selection {
                            let available = try direction < 0 ?
                                proposal[index].duration(since: windows[index].lowerBound) :
                                windows[index].upperBound.duration(since: proposal[index])
                            proposal[index] = delta >= available ?
                                (direction < 0 ? windows[index].lowerBound : windows[index].upperBound) :
                                try proposal[index].advanced(by: Duration(seconds: direction))
                        }
                    }
                    for index in proposal.indices {
                        proposal[index] = max(windows[index].lowerBound, min(proposal[index], windows[index].upperBound))
                    }
                    guard !tried.contains(proposal) else { continue }
                    if evaluated == optimizationOptions.maximumCandidateEvaluations { exhaustedBudget = true; break }
                    let candidate = try evaluate(proposal)
                    if candidate.status == .feasible {
                        best = candidate; bestEpochs = proposal; improved = true; break
                    }
                }
            }
            iterations.append(WaypointTimingOptimizationIteration(stepSeconds: step,
                evaluatedCandidates: evaluated-before, improved: improved, bestArrivalEpoch: bestEpochs[last]))
            if exhaustedBudget { return finish(.evaluationLimitReached) }
            if !improved {
                if step <= optimizationOptions.minimumStepSeconds { return finish(.stepLimitReached) }
                step = max(optimizationOptions.minimumStepSeconds, step/2)
            }
        }
        return finish(.iterationLimitReached)
    }
}


public enum WaypointTimingMultiStartStatus: Equatable, Sendable {
    case seedsExhausted, evaluationLimitReached
}
public struct WaypointTimingMultiStartResult {
    public let status: WaypointTimingMultiStartStatus
    public let evaluatedCandidates: Int
    public let starts: [WaypointTimingOptimizationResult]
    public let route: WaypointBrachistochroneRoute?
    public let propagation: WaypointRoutePropagationResult?
}
extension WaypointRouteTimingPlanner {
    /// Explore bounded, caller-selected seeds, retaining the earliest feasible
    /// local result. Seeds may span the full windows; no global optimum is implied.
    public static func optimizeTimingFromSeeds(departureEpoch: Epoch, frame: Frame, origin: BodyIdentifier,
        waypoints: [BrachistochroneWaypoint], legAccelerationsMetresPerSecondSquared: [Double],
        arrivalWindows: [ClosedRange<Epoch>], seedArrivalEpochs: [[Epoch]], maximumTotalCandidateEvaluations: Int = 1000,
        optimizationOptions: WaypointTimingOptimizationOptions, initialMass: SpacecraftMassProperties,
        backgroundDynamics: SpacecraftDynamics, legEngineLimits: [WaypointRouteEngineLimits],
        timedEnginePhases: [WaypointRouteEnginePhase] = [], obstacles: [SphericalExclusionRegion] = [],
        numericalOptions: NumericalPropagationOptions = try! NumericalPropagationOptions(),
        endpointPositionToleranceMetres: Double = 0.01,
        endpointVelocityToleranceMetresPerSecond: Double = 1e-5,
        movingTarget: MovingWaypointRouteTarget? = nil,
        movingWaypoints: [Int: MovingWaypointPassage] = [:],
        movingObstacles: [MovingSphericalExclusionRegion] = [],
        physicalConstraints: WaypointRoutePhysicalConstraints = .unconstrained) throws -> WaypointTimingMultiStartResult {
        guard !seedArrivalEpochs.isEmpty, seedArrivalEpochs.count <= 256,
              (1...10_000).contains(maximumTotalCandidateEvaluations) else {
            throw WaypointTimingSearchError.invalidBudget
        }
        var starts: [WaypointTimingOptimizationResult] = []
        var seen = Set<[Epoch]>(), evaluated = 0
        var best: WaypointTimingOptimizationResult?
        for seed in seedArrivalEpochs where seen.insert(seed).inserted {
            if evaluated == maximumTotalCandidateEvaluations { break }
            let options = try WaypointTimingOptimizationOptions(
                initialStepSeconds: optimizationOptions.initialStepSeconds,
                minimumStepSeconds: optimizationOptions.minimumStepSeconds,
                maximumIterations: optimizationOptions.maximumIterations,
                maximumCandidateEvaluations: min(optimizationOptions.maximumCandidateEvaluations,
                    maximumTotalCandidateEvaluations-evaluated))
            let result = try optimizeTiming(departureEpoch: departureEpoch, frame: frame, origin: origin,
                waypoints: waypoints, legAccelerationsMetresPerSecondSquared: legAccelerationsMetresPerSecondSquared,
                arrivalWindows: arrivalWindows, seedArrivalEpochs: seed, optimizationOptions: options,
                initialMass: initialMass, backgroundDynamics: backgroundDynamics, legEngineLimits: legEngineLimits,
                timedEnginePhases: timedEnginePhases, obstacles: obstacles, numericalOptions: numericalOptions,
                endpointPositionToleranceMetres: endpointPositionToleranceMetres,
                endpointVelocityToleranceMetresPerSecond: endpointVelocityToleranceMetresPerSecond,
                movingTarget: movingTarget, movingWaypoints: movingWaypoints, movingObstacles: movingObstacles, physicalConstraints: physicalConstraints)
            evaluated += result.evaluatedCandidates; starts.append(result)
            if let route = result.route, best?.route == nil || route.arrivalEpoch < best!.route!.arrivalEpoch {
                best = result
            }
        }
        return WaypointTimingMultiStartResult(
            status: evaluated == maximumTotalCandidateEvaluations ? .evaluationLimitReached : .seedsExhausted,
            evaluatedCandidates: evaluated, starts: starts, route: best?.route, propagation: best?.propagation)
    }
}
