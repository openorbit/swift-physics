// SPDX-License-Identifier: Apache-2.0
import Math

public enum FuelDepletionPolicy: Sendable { case fail, cutoff }
public enum ScheduledThrustPropagationError: Error, Equatable {
    case backgroundMassFlowUnsupported, depletionEpochResolutionExceeded
}
public struct ScheduledThrustPropagationResult {
    public let state: SpacecraftState
    public let depletionEpoch: Epoch?
}

/// Immutable-seed numerical replay with exact scheduled ignition/cutoff boundaries.
/// Background forces must be continuous between boundaries; their own events are
/// not detected. Options/budgets apply to each segment. Optional automatic fuel cutoff uses
/// the exact constant scheduled mass flow; background mass flow is unsupported in that mode.
public struct ScheduledThrustPropagator: StateProvider {
    public let initialState: SpacecraftState
    public let schedule: ScheduledThrust
    public let context: ForceEvaluationContext
    public let options: NumericalPropagationOptions
    public let fuelDepletionPolicy: FuelDepletionPolicy
    private let backgroundModels: [any ForceModel]
    public var frame: Frame { initialState.frame }
    public var referenceEpoch: Epoch { initialState.epoch }

    public init(initialState: SpacecraftState, schedule: ScheduledThrust,
                backgroundModels: [any ForceModel] = [], context: ForceEvaluationContext,
                options: NumericalPropagationOptions = try! NumericalPropagationOptions(),
                fuelDepletionPolicy: FuelDepletionPolicy = .fail) throws {
        guard initialState.frame === schedule.frame, initialState.frame === context.frame else {
            throw SpacecraftDynamicsError.frameMismatch
        }
        self.initialState = initialState; self.schedule = schedule
        self.backgroundModels = backgroundModels; self.context = context; self.options = options
        self.fuelDepletionPolicy = fuelDepletionPolicy
    }

    public func state(at epoch: Epoch) throws -> CartesianState { try spacecraftState(at: epoch).cartesian }
    public func spacecraftState(at target: Epoch) throws -> SpacecraftState { try propagation(at: target).state }
    public func propagation(at target: Epoch) throws -> ScheduledThrustPropagationResult {
        let forward = target >= referenceEpoch
        let lower = min(target, referenceEpoch), upper = max(target, referenceEpoch)
        var boundaries = schedule.burns.flatMap { [$0.start,$0.end] }
            .filter { $0 > lower && $0 < upper }
            .sorted { forward ? $0 < $1 : $0 > $1 }
        // Shared boundaries are integrated only once, preserving exact epoch identity.
        boundaries = boundaries.enumerated().compactMap { index, epoch in
            index == 0 || epoch != boundaries[index-1] ? epoch : nil
        }
        boundaries.append(target)
        var current = initialState
        var depletionEpoch: Epoch?
        for end in boundaries where end != current.epoch {
            let start = current.epoch
            // Segment interior has a constant set of active engines. Evaluate its
            // smooth extension at both RK endpoints, irrespective of replay direction.
            let active = schedule.burns.filter { $0.start <= min(start,end) && $0.end >= max(start,end) }
            if fuelDepletionPolicy == .cutoff {
                var flow = active.reduce(0.0) { $0 + $1.thrust.propellantConsumptionKilogramsPerSecond }
                guard flow.isFinite else { throw NumericalPropagationError.numericalRangeExceeded }
                let fuel = current.mass.propellantMassKilograms
                if fuel == 0 && forward {
                    if flow > 0 && depletionEpoch == nil { depletionEpoch = start }
                    flow = 0
                }
                // Backward replay from an empty independent seed lacks the cutoff
                // history needed to distinguish coast from a historical firing.
                if fuel == 0 && !forward && flow > 0 { throw NumericalPropagationError.propellantDepleted }
                let interval = try end.duration(since: start)
                let duration = interval.seconds
                let depletionSeconds = flow > 0 ? fuel / flow : Double.infinity
                let depletes = try forward && depletionSeconds.isFinite && depletionSeconds <= duration
                    && interval >= Duration(seconds: depletionSeconds)
                let burnEnd = depletes ? try start.advanced(by: Duration(seconds: fuel / flow)) : end
                if depletes && fuel > 0 && burnEnd <= start {
                    throw ScheduledThrustPropagationError.depletionEpochResolutionExceeded
                }
                if burnEnd != start {
                    let segment = ScheduledFuelSegment(seed: current, flow: flow, engines: flow > 0 ? active.map(\.thrust) : [],
                        background: backgroundModels, end: burnEnd, exhaustsFuel: depletes)
                    let dynamics = try SpacecraftDynamics(models: [segment], context: context)
                    let propagated = try NumericalSpacecraftPropagator(initialState: current, dynamics: dynamics, options: options)
                        .spacecraftState(at: burnEnd)
                    let remainingFuel = depletes ? 0 : fuel - flow * duration
                    current = SpacecraftState(cartesian: propagated.cartesian,
                        mass: try SpacecraftMassProperties(dryMassKilograms: current.mass.dryMassKilograms,
                            propellantMassKilograms: remainingFuel))
                }
                if depletes {
                    depletionEpoch = burnEnd
                    if burnEnd != end {
                        let coast = ScheduledFuelSegment(seed: current, flow: 0, engines: [], background: backgroundModels,
                            end: end, exhaustsFuel: false)
                        current = try NumericalSpacecraftPropagator(initialState: current,
                            dynamics: SpacecraftDynamics(models: [coast], context: context), options: options).spacecraftState(at: end)
                    }
                }
                continue
            }
            let models: [any ForceModel] = backgroundModels + active.map { $0.thrust as any ForceModel }
            let dynamics = try SpacecraftDynamics(models: models, context: context)
            current = try NumericalSpacecraftPropagator(initialState: current, dynamics: dynamics, options: options)
                .spacecraftState(at: end)
        }
        // Returned acceleration follows the right-continuous [start,end) schedule,
        // not the smooth extension used for integration up to a boundary.
        let endpointModels: [any ForceModel]
        if fuelDepletionPolicy == .cutoff {
            let engines = current.mass.propellantMassKilograms > 0 ? schedule.burns.filter {
                target >= $0.start && target < $0.end
            }.map(\.thrust) : []
            endpointModels = [ScheduledFuelSegment(seed: current, flow: 0, engines: engines,
                background: backgroundModels, end: target, exhaustsFuel: false)]
        } else { endpointModels = backgroundModels + [schedule] }
        let endpoint = try SpacecraftDynamics(models: endpointModels, context: context).derivative(at: current)
        let result = SpacecraftState(cartesian: try CartesianState(epoch: target, frame: frame,
            position: current.cartesian.position, velocity: current.cartesian.velocity,
            acceleration: endpoint.velocityRateMetresPerSecondSquared), mass: current.mass)
        return ScheduledThrustPropagationResult(state: result, depletionEpoch: depletionEpoch)
    }
}


/// Integrate Cartesian motion while supplying analytically evolving fuel to every
/// force. RK's reservoir coordinate is held constant, avoiding negative trial fuel.
private struct ScheduledFuelSegment: ForceModel {
    let seed: SpacecraftState
    let flow: Double
    let engines: [FiniteThrust]
    let background: [any ForceModel]
    let end: Epoch
    let exhaustsFuel: Bool
    func contribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
        let elapsed = try state.epoch.duration(since: seed.epoch).seconds
        let fuel = exhaustsFuel && state.epoch == end ? 0 : seed.mass.propellantMassKilograms - flow * elapsed
        guard fuel.isFinite, fuel >= 0 else { throw NumericalPropagationError.propellantDepleted }
        let corrected = SpacecraftState(cartesian: state.cartesian,
            mass: try SpacecraftMassProperties(dryMassKilograms: seed.mass.dryMassKilograms, propellantMassKilograms: fuel))
        var acceleration = SIMD3<Double>.zero
        for model in background {
            let result = try model.contribution(for: corrected, in: context)
            guard result.epoch == state.epoch else { throw SpacecraftDynamicsError.epochMismatch }
            guard result.frame === state.frame else { throw SpacecraftDynamicsError.frameMismatch }
            guard result.propellantMassRateKilogramsPerSecond == 0 else {
                throw ScheduledThrustPropagationError.backgroundMassFlowUnsupported
            }
            acceleration += result.accelerationMetresPerSecondSquared
        }
        // Smooth left-limit thrust at depletion, evaluated with dry mass. The
        // public FiniteThrust model still rejects actual firing with zero fuel.
        for engine in engines {
            acceleration += engine.direction * (engine.thrustNewtons / corrected.mass.totalMassKilograms)
        }
        return try ForceContribution(epoch: state.epoch, frame: state.frame, accelerationMetresPerSecondSquared: acceleration)
    }
}
