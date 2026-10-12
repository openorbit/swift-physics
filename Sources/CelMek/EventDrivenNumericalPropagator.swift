// SPDX-License-Identifier: Apache-2.0
import Math

/// Explicit event reset. State epoch/frame must match the localized occurrence;
/// dynamics must retain the same frame and physical origin. The caller implements
/// forward and inverse resets when backward replay is required.
public struct EventStateTransition {
    public let state: SpacecraftState
    public let dynamics: SpacecraftDynamics
    public init(state: SpacecraftState, dynamics: SpacecraftDynamics) {
        self.state = state; self.dynamics = dynamics
    }
}

/// Accepted-step crossing detection with explicit restarts and optional rearming.
/// Detectors are one-shot unless given a positive rearming threshold. Queries restart from the seed;
/// no closure may mutate hidden simulation state. Forces must remain evaluable
/// throughout a crossing bracket, including its trial stages.
public struct EventDrivenNumericalPropagator {
    public let initialState: SpacecraftState
    public let dynamics: SpacecraftDynamics
    public let detectors: [any PropagationEventDetector]
    public let numericalOptions: NumericalPropagationOptions
    public let eventOptions: PropagationEventOptions
    public let rearmingThresholds: [String: Double]
    private let transition: ((PropagationEventOccurrence, Bool, SpacecraftDynamics) throws -> EventStateTransition)?

    public init(initialState: SpacecraftState, dynamics: SpacecraftDynamics,
                detectors: [any PropagationEventDetector],
                numericalOptions: NumericalPropagationOptions = try! NumericalPropagationOptions(),
                eventOptions: PropagationEventOptions = try! PropagationEventOptions(),
                rearmingThresholds: [String: Double] = [:],
                transition: ((PropagationEventOccurrence, Bool, SpacecraftDynamics) throws -> EventStateTransition)? = nil) throws {
        _ = try NumericalSpacecraftPropagator(initialState: initialState, dynamics: dynamics, options: numericalOptions)
        _ = try EventDetectedPropagation(trajectory: NumericalSpacecraftPropagator(initialState: initialState,
            dynamics: dynamics, options: numericalOptions), detectors: detectors, options: eventOptions)
        guard rearmingThresholds.allSatisfy({ entry in
            detectors.contains(where: { $0.identifier == entry.key }) && entry.value.isFinite && entry.value > 0
        }) else { throw PropagationEventError.invalidOptions }
        self.rearmingThresholds = rearmingThresholds
        self.initialState = initialState; self.dynamics = dynamics; self.detectors = detectors
        self.numericalOptions = numericalOptions; self.eventOptions = eventOptions; self.transition = transition
    }

    public func propagation(at target: Epoch) throws -> EventPropagationResult {
        var current = initialState, activeDynamics = dynamics
        var fired = Set<String>(), events: [PropagationEventOccurrence] = []
        var queries = 0, acceptedChecks = 0
        let forward = target >= initialState.epoch
        // Bound accepted-step spacing independently of the numerical accuracy preset.
        let options = try NumericalPropagationOptions(
            absolutePositionToleranceMetres: numericalOptions.absolutePositionToleranceMetres,
            absoluteVelocityToleranceMetresPerSecond: numericalOptions.absoluteVelocityToleranceMetresPerSecond,
            absolutePropellantToleranceKilograms: numericalOptions.absolutePropellantToleranceKilograms,
            relativeTolerance: numericalOptions.relativeTolerance,
            initialStepSeconds: min(numericalOptions.initialStepSeconds,eventOptions.maximumCheckIntervalSeconds),
            minimumStepSeconds: min(numericalOptions.minimumStepSeconds,eventOptions.maximumCheckIntervalSeconds),
            maximumStepSeconds: min(numericalOptions.maximumStepSeconds,eventOptions.maximumCheckIntervalSeconds),
            maximumStepAttempts: numericalOptions.maximumStepAttempts,
            maximumDerivativeEvaluations: numericalOptions.maximumDerivativeEvaluations)
        func search(_ seed: SpacecraftState, to end: Epoch, endpoint: SpacecraftState? = nil) throws -> PropagationEventOccurrence? {
            if endpoint != nil {
                for detector in detectors where fired.contains(detector.identifier) {
                    if let threshold = rearmingThresholds[detector.identifier] {
                        let value = try detector.value(at: seed)
                        guard value.isFinite else { throw PropagationEventError.nonFiniteValue }
                        if abs(value) > threshold { fired.remove(detector.identifier) }
                    }
                }
            }
            var remaining = detectors.filter { !fired.contains($0.identifier) }
            if let endpoint {
                // Reuse accepted states. Only an actual crossing triggers local
                // reintegration; ordinary steps require no additional state queries.
                remaining = try remaining.filter { detector in
                    let a = try detector.value(at: seed), b = try detector.value(at: endpoint)
                    guard a.isFinite, b.isFinite else { throw PropagationEventError.nonFiniteValue }
                    if a == 0 { return true }
                    guard b == 0 || (a < 0) != (b < 0) else { return false }
                    let increasing = forward ? b > a : a > b
                    return detector.direction == .any || (detector.direction == .increasing) == increasing
                }
            } else {
                remaining = try remaining.filter { detector in
                    let value = try detector.value(at: seed)
                    guard value.isFinite else { throw PropagationEventError.nonFiniteValue }
                    return value == 0
                }
            }
            if remaining.isEmpty { return nil }
            guard queries < eventOptions.maximumStateEvaluations else { throw PropagationEventError.evaluationLimitExceeded }
            let bounded = try PropagationEventOptions(maximumCheckIntervalSeconds: eventOptions.maximumCheckIntervalSeconds,
                timeToleranceSeconds: eventOptions.timeToleranceSeconds,
                maximumLocalizationIterations: eventOptions.maximumLocalizationIterations,
                maximumStateEvaluations: eventOptions.maximumStateEvaluations-queries,
                maximumEvents: eventOptions.maximumEvents)
            let local = try NumericalSpacecraftPropagator(initialState: seed, dynamics: activeDynamics, options: options)
            let firstCrossings: [any PropagationEventDetector] = try remaining.map {
                if let timed = $0 as? EpochEvent {
                    return try EpochEvent(identifier: timed.identifier, epoch: timed.epoch, action: .stop)
                }
                return FirstCrossing(detector: $0)
            }
            let result = try EventDetectedPropagation(trajectory: local,
                detectors: firstCrossings, options: bounded)
                .propagation(from: seed.epoch, to: end)
            queries += result.stateEvaluations
            guard let found = result.events.first,
                  let original = remaining.first(where: { $0.identifier == found.identifier }) else { return nil }
            return PropagationEventOccurrence(identifier: found.identifier, action: original.action,
                state: found.state, bracketDurationSeconds: found.bracketDurationSeconds)
        }
        var initial = true
        while true {
            var occurrence: PropagationEventOccurrence?
            if initial { occurrence = try search(current,to: current.epoch); initial = false }
            if occurrence == nil && current.epoch != target {
                let integrator = try NumericalSpacecraftPropagator(initialState: current, dynamics: activeDynamics, options: options)
                let result = try integrator.propagation(at: target) { left, right in
                    guard acceptedChecks < eventOptions.maximumStateEvaluations else { throw PropagationEventError.evaluationLimitExceeded }
                    acceptedChecks += 1
                    occurrence = try search(left,to: right.epoch,endpoint: right)
                    return occurrence == nil
                }
                if occurrence == nil { current = result.state }
            }
            guard let event = occurrence else { break }
            guard events.count < eventOptions.maximumEvents else { throw PropagationEventError.eventLimitExceeded }
            events.append(event); fired.insert(event.identifier); current = event.state
            if event.action == .stop {
                return EventPropagationResult(state: current, events: events,
                    reachedTarget: current.epoch == target, stateEvaluations: queries)
            }
            if let transition {
                let reset = try transition(event,forward,activeDynamics)
                guard reset.state.epoch == current.epoch else { throw SpacecraftDynamicsError.epochMismatch }
                guard reset.state.frame === initialState.frame, reset.dynamics.context.frame === initialState.frame else {
                    throw SpacecraftDynamicsError.frameMismatch
                }
                guard reset.dynamics.context.origin == dynamics.context.origin else { throw EnvironmentalForceError.originMismatch }
                current = reset.state; activeDynamics = reset.dynamics
            }
            // Recheck remaining detectors at the reset epoch, including simultaneous events.
            initial = true
        }
        // Recompute endpoint acceleration after an event reset at the requested epoch.
        let derivative = try activeDynamics.derivative(at: current)
        let state = SpacecraftState(cartesian: try CartesianState(epoch: current.epoch, frame: current.frame,
            position: current.cartesian.position, velocity: current.cartesian.velocity,
            acceleration: derivative.velocityRateMetresPerSecondSquared), mass: current.mass)
        return EventPropagationResult(state: state, events: events, reachedTarget: state.epoch == target, stateEvaluations: queries)
    }
}

private struct FirstCrossing: PropagationEventDetector {
    let detector: any PropagationEventDetector
    var identifier: String { detector.identifier }
    var direction: EventCrossingDirection { detector.direction }
    var action: PropagationEventAction { .stop }
    func value(at state: SpacecraftState) throws -> Double { try detector.value(at: state) }
}
