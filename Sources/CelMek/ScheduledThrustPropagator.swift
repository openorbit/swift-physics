// SPDX-License-Identifier: Apache-2.0
import Math

/// Immutable-seed numerical replay with exact scheduled ignition/cutoff boundaries.
/// Background forces must be continuous between boundaries; their own events are
/// not detected. Options/budgets apply to each segment. No automatic fuel cutoff.
public struct ScheduledThrustPropagator: StateProvider {
    public let initialState: SpacecraftState
    public let schedule: ScheduledThrust
    public let context: ForceEvaluationContext
    public let options: NumericalPropagationOptions
    private let backgroundModels: [any ForceModel]
    public var frame: Frame { initialState.frame }
    public var referenceEpoch: Epoch { initialState.epoch }

    public init(initialState: SpacecraftState, schedule: ScheduledThrust,
                backgroundModels: [any ForceModel] = [], context: ForceEvaluationContext,
                options: NumericalPropagationOptions = try! NumericalPropagationOptions()) throws {
        guard initialState.frame === schedule.frame, initialState.frame === context.frame else {
            throw SpacecraftDynamicsError.frameMismatch
        }
        self.initialState = initialState; self.schedule = schedule
        self.backgroundModels = backgroundModels; self.context = context; self.options = options
    }

    public func state(at epoch: Epoch) throws -> CartesianState { try spacecraftState(at: epoch).cartesian }
    public func spacecraftState(at target: Epoch) throws -> SpacecraftState {
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
        for end in boundaries where end != current.epoch {
            let start = current.epoch
            // Segment interior has a constant set of active engines. Evaluate its
            // smooth extension at both RK endpoints, irrespective of replay direction.
            let active = schedule.burns.filter { $0.start <= min(start,end) && $0.end >= max(start,end) }
            let models: [any ForceModel] = backgroundModels + active.map { $0.thrust as any ForceModel }
            let dynamics = try SpacecraftDynamics(models: models, context: context)
            current = try NumericalSpacecraftPropagator(initialState: current, dynamics: dynamics, options: options)
                .spacecraftState(at: end)
        }
        // Returned acceleration follows the right-continuous [start,end) schedule,
        // not the smooth extension used for integration up to a boundary.
        let endpoint = try SpacecraftDynamics(models: backgroundModels + [schedule], context: context).derivative(at: current)
        return SpacecraftState(cartesian: try CartesianState(epoch: target, frame: frame,
            position: current.cartesian.position, velocity: current.cartesian.velocity,
            acceleration: endpoint.velocityRateMetresPerSecondSquared), mass: current.mass)
    }
}
