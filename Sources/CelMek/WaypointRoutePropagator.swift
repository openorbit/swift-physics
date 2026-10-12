// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum WaypointRoutePropagationError: Error, Equatable {
    case invalidLimits, unsupportedDwell, backgroundMassFlow
    case accelerationLimitExceeded(leg: Int), thrustLimitExceeded(leg: Int)
    case endpointToleranceExceeded(leg: Int)
    case invalidTimedPhase, overlappingTimedPhases, timedPhaseOutsideCoverage
}

/// Engine limits apply to thrust acceleration, including compensation for gravity.
public struct WaypointRouteEngineLimits: Equatable, Sendable {
    /// Standard gravity, for converting an explicitly requested limit in g to SI.
    public static let standardGravityMetresPerSecondSquared = 9.80665
    public let maximumAccelerationMetresPerSecondSquared: Double
    public let maximumThrustNewtons: Double
    public let effectiveExhaustVelocityMetresPerSecond: Double
    public init(maximumAccelerationMetresPerSecondSquared: Double, maximumThrustNewtons: Double,
                effectiveExhaustVelocityMetresPerSecond: Double) throws {
        guard [maximumAccelerationMetresPerSecondSquared, maximumThrustNewtons,
               effectiveExhaustVelocityMetresPerSecond].allSatisfy({ $0.isFinite && $0 > 0 }) else {
            throw WaypointRoutePropagationError.invalidLimits
        }
        self.maximumAccelerationMetresPerSecondSquared = maximumAccelerationMetresPerSecondSquared
        self.maximumThrustNewtons = maximumThrustNewtons
        self.effectiveExhaustVelocityMetresPerSecond = effectiveExhaustVelocityMetresPerSecond
    }
}

/// Engine-limit override on [startEpoch, endEpoch), independent of spatial legs.
/// At final arrival the incoming phase supplies the left-limit engine evaluation.
public struct WaypointRouteEnginePhase: Equatable {
    public let startEpoch: Epoch
    public let endEpoch: Epoch
    public let limits: WaypointRouteEngineLimits
    public init(startEpoch: Epoch, endEpoch: Epoch, limits: WaypointRouteEngineLimits) throws {
        guard endEpoch > startEpoch else { throw WaypointRoutePropagationError.invalidTimedPhase }
        self.startEpoch = startEpoch; self.endEpoch = endEpoch; self.limits = limits
    }
    public init(departureEpoch: Epoch, startOffset: Duration, duration: Duration,
                limits: WaypointRouteEngineLimits) throws {
        guard startOffset >= .zero, duration > .zero else { throw WaypointRoutePropagationError.invalidTimedPhase }
        let start = try departureEpoch.advanced(by: startOffset)
        try self.init(startEpoch: start, endEpoch: start.advanced(by: duration), limits: limits)
    }
}

public struct WaypointRouteExecutionSegment {
    public let legIndex: Int
    /// Index in the propagator's chronologically sorted phases; nil uses leg limits.
    public let timedEnginePhaseIndex: Int?
    public let initialState: SpacecraftState
    public let finalState: SpacecraftState
    public let diagnostics: NumericalPropagationDiagnostics
}

public struct WaypointRouteCheckpoint {
    public let legIndex: Int
    public let state: SpacecraftState
    public let positionErrorMetres: Double
    public let velocityErrorMetresPerSecond: Double
    public let diagnostics: NumericalPropagationDiagnostics
}

public struct WaypointRoutePropagationResult {
    public let state: SpacecraftState
    /// Completed legs only; checkpoint mass is the next leg's initial mass.
    public let checkpoints: [WaypointRouteCheckpoint]
    /// Positive-duration intervals split at both spatial and temporal boundaries.
    public let segments: [WaypointRouteExecutionSegment]
}

/// Numerically executes a prescribed cubic route using ideal, instantaneous
/// thrust steering. This is inverse dynamics, not a trajectory/time optimizer.
/// Queries replay forward from departure and split exactly at every waypoint.
public struct WaypointRoutePropagator: SpacecraftStateProvider {
    public let route: WaypointBrachistochroneRoute
    public let initialMass: SpacecraftMassProperties
    public let backgroundDynamics: SpacecraftDynamics
    public let legEngineLimits: [WaypointRouteEngineLimits]
    public let timedEnginePhases: [WaypointRouteEnginePhase]
    private let temporalBoundaries: [Epoch]
    public let numericalOptions: NumericalPropagationOptions
    public let endpointPositionToleranceMetres: Double
    public let endpointVelocityToleranceMetresPerSecond: Double
    public var frame: Frame { route.frame }
    public var coverage: ClosedRange<Epoch> { route.coverage }

    public init(route: WaypointBrachistochroneRoute, initialMass: SpacecraftMassProperties,
                backgroundDynamics: SpacecraftDynamics, legEngineLimits: [WaypointRouteEngineLimits],
                timedEnginePhases: [WaypointRouteEnginePhase] = [],
                numericalOptions: NumericalPropagationOptions = try! NumericalPropagationOptions(),
                endpointPositionToleranceMetres: Double = 0.01,
                endpointVelocityToleranceMetresPerSecond: Double = 1e-5) throws {
        guard route.frame === backgroundDynamics.context.frame else { throw SpacecraftDynamicsError.frameMismatch }
        guard route.origin == backgroundDynamics.context.origin else { throw BrachistochroneRouteError.originMismatch }
        guard legEngineLimits.count == route.legs.count,
              [endpointPositionToleranceMetres, endpointVelocityToleranceMetresPerSecond]
                .allSatisfy({ $0.isFinite && $0 > 0 }) else { throw WaypointRoutePropagationError.invalidLimits }
        guard route.waypoints.allSatisfy({ $0.dwellDuration == .zero }) else {
            throw WaypointRoutePropagationError.unsupportedDwell
        }
        guard timedEnginePhases.count <= 10_000 else { throw WaypointRoutePropagationError.invalidTimedPhase }
        let phases = timedEnginePhases.sorted { $0.startEpoch < $1.startEpoch }
        for (index, phase) in phases.enumerated() {
            guard phase.startEpoch >= route.departureEpoch, phase.endEpoch <= route.arrivalEpoch else {
                throw WaypointRoutePropagationError.timedPhaseOutsideCoverage
            }
            guard index == 0 || phases[index-1].endEpoch <= phase.startEpoch else {
                throw WaypointRoutePropagationError.overlappingTimedPhases
            }
        }
        self.timedEnginePhases = phases
        let boundaries = phases.flatMap { [$0.startEpoch, $0.endEpoch] }.sorted()
        self.temporalBoundaries = boundaries.enumerated().compactMap { index, epoch in
            index == 0 || boundaries[index-1] != epoch ? epoch : nil
        }
        self.route = route; self.initialMass = initialMass; self.backgroundDynamics = backgroundDynamics
        self.legEngineLimits = legEngineLimits; self.numericalOptions = numericalOptions
        self.endpointPositionToleranceMetres = endpointPositionToleranceMetres
        self.endpointVelocityToleranceMetresPerSecond = endpointVelocityToleranceMetresPerSecond
    }

    public func state(at epoch: Epoch) throws -> CartesianState { try spacecraftState(at: epoch).cartesian }
    public func spacecraftState(at epoch: Epoch) throws -> SpacecraftState { try propagation(at: epoch).state }

    public func propagation(at epoch: Epoch) throws -> WaypointRoutePropagationResult {
        guard coverage.contains(epoch) else { throw BrachistochroneError.outsideCoverage }
        var current = SpacecraftState(cartesian: try route.legs[0].state(at: route.departureEpoch), mass: initialMass)
        var checkpoints: [WaypointRouteCheckpoint] = [], segments: [WaypointRouteExecutionSegment] = []
        var boundaryIndex = 0
        for (index, leg) in route.legs.enumerated() {
            guard epoch >= leg.departureEpoch else { break }
            var diagnostics: NumericalPropagationDiagnostics?
            repeat {
                while boundaryIndex < temporalBoundaries.count && temporalBoundaries[boundaryIndex] <= current.epoch {
                    boundaryIndex += 1
                }
                let phaseIndex = enginePhaseIndex(at: current.epoch)
                let model = RouteExecutionModel(leg: leg, index: index, background: backgroundDynamics,
                    limits: phaseIndex.map { timedEnginePhases[$0].limits } ?? legEngineLimits[index])
                let dynamics = try SpacecraftDynamics(models: [model], context: backgroundDynamics.context)
                var end = min(epoch, leg.arrivalEpoch)
                if boundaryIndex < temporalBoundaries.count { end = min(end, temporalBoundaries[boundaryIndex]) }
                let initial = current
                let result = try NumericalSpacecraftPropagator(initialState: initial, dynamics: dynamics,
                    options: numericalOptions).propagation(at: end)
                current = result.state
                diagnostics = mergeRouteDiagnostics(diagnostics, result.diagnostics)
                if end > initial.epoch {
                    segments.append(WaypointRouteExecutionSegment(legIndex: index, timedEnginePhaseIndex: phaseIndex,
                        initialState: initial, finalState: current, diagnostics: result.diagnostics))
                }
                if end == epoch || end == leg.arrivalEpoch { break }
            } while true
            if current.epoch == leg.arrivalEpoch {
                let end = current.epoch
                let reference = try leg.state(at: end)
                let positionError = executionNorm(current.cartesian.position-reference.position)
                let velocityError = executionNorm(current.cartesian.velocity-reference.velocity)
                guard positionError <= endpointPositionToleranceMetres,
                      velocityError <= endpointVelocityToleranceMetresPerSecond else {
                    throw WaypointRoutePropagationError.endpointToleranceExceeded(leg: index)
                }
                checkpoints.append(WaypointRouteCheckpoint(legIndex: index, state: current,
                    positionErrorMetres: positionError, velocityErrorMetresPerSecond: velocityError,
                    diagnostics: diagnostics!))
            }
            if epoch < leg.arrivalEpoch {
                // A phase change at the query epoch takes effect immediately.
                // The incoming segment used its own left limit for all RK stages.
                let phaseIndex = enginePhaseIndex(at: epoch)
                let model = RouteExecutionModel(leg: leg, index: index, background: backgroundDynamics,
                    limits: phaseIndex.map { timedEnginePhases[$0].limits } ?? legEngineLimits[index])
                let force = try model.contribution(for: current, in: backgroundDynamics.context)
                current = SpacecraftState(cartesian: try CartesianState(epoch: epoch, frame: frame,
                    position: current.cartesian.position, velocity: current.cartesian.velocity,
                    acceleration: force.accelerationMetresPerSecondSquared), mass: current.mass)
            }
            // At a shared boundary evaluate the outgoing leg too, for right-continuous
            // acceleration and engine constraints. Never reset position/velocity/mass.
            if epoch < leg.arrivalEpoch || index == route.legs.count-1 { break }
        }
        return WaypointRoutePropagationResult(state: current, checkpoints: checkpoints, segments: segments)
    }

    private func enginePhaseIndex(at epoch: Epoch) -> Int? {
        var low = 0, high = timedEnginePhases.count
        while low < high {
            let middle = low+(high-low)/2
            if timedEnginePhases[middle].startEpoch <= epoch { low = middle+1 } else { high = middle }
        }
        guard low > 0, epoch < timedEnginePhases[low-1].endEpoch else { return nil }
        return low-1
    }
}

private func mergeRouteDiagnostics(_ previous: NumericalPropagationDiagnostics?,
                                   _ next: NumericalPropagationDiagnostics) -> NumericalPropagationDiagnostics {
    guard let previous else { return next }
    return NumericalPropagationDiagnostics(acceptedSteps: previous.acceptedSteps+next.acceptedSteps,
        rejectedSteps: previous.rejectedSteps+next.rejectedSteps,
        derivativeEvaluations: previous.derivativeEvaluations+next.derivativeEvaluations,
        smallestAcceptedStepSeconds: [previous.smallestAcceptedStepSeconds, next.smallestAcceptedStepSeconds].compactMap { $0 }.min(),
        largestAcceptedStepSeconds: [previous.largestAcceptedStepSeconds, next.largestAcceptedStepSeconds].compactMap { $0 }.max(),
        maximumAcceptedNormalizedError: max(previous.maximumAcceptedNormalizedError, next.maximumAcceptedNormalizedError))
}

private struct RouteExecutionModel: ForceModel {
    let leg: WaypointBrachistochroneLeg
    let index: Int
    let background: SpacecraftDynamics
    let limits: WaypointRouteEngineLimits
    func contribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
        let environment = try background.derivative(at: state)
        guard environment.propellantMassRateKilogramsPerSecond == 0 else {
            throw WaypointRoutePropagationError.backgroundMassFlow
        }
        let reference = try leg.state(at: state.epoch)
        let acceleration = reference.acceleration!
        let thrustAcceleration = executionNorm(acceleration-environment.velocityRateMetresPerSecondSquared)
        guard thrustAcceleration.isFinite else { throw SpacecraftDynamicsError.numericalRangeExceeded }
        guard thrustAcceleration <= limits.maximumAccelerationMetresPerSecondSquared else {
            throw WaypointRoutePropagationError.accelerationLimitExceeded(leg: index)
        }
        let thrust = thrustAcceleration*state.mass.totalMassKilograms
        guard thrust.isFinite else { throw SpacecraftDynamicsError.numericalRangeExceeded }
        guard thrust <= limits.maximumThrustNewtons else {
            throw WaypointRoutePropagationError.thrustLimitExceeded(leg: index)
        }
        guard thrust == 0 || state.mass.propellantMassKilograms > 0 else {
            throw NumericalPropagationError.propellantDepleted
        }
        // Sum the environment and its compensating engine analytically to avoid
        // losing a small prescribed acceleration in cancellation of large forces.
        return try ForceContribution(epoch: state.epoch, frame: context.frame,
            accelerationMetresPerSecondSquared: acceleration,
            propellantMassRateKilogramsPerSecond: -thrust/limits.effectiveExhaustVelocityMetresPerSecond)
    }
}
private func executionNorm(_ value: SIMD3<Double>) -> Double { hypot(hypot(value.x,value.y),value.z) }
