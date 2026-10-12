// SPDX-License-Identifier: Apache-2.0
import Math

/// Pure state/epoch-dependent engine command, including throttle and direction.
/// The command must be defined at dry mass for depletion-boundary localization.
/// It must remain continuous within integration segments; schedule discontinuities
/// through explicit event resets. Normal evaluation still rejects firing dry.
public struct VariableThrust: FuelConsumingForceModel {
    public let frame: Frame
    private let command: (SpacecraftState, ForceEvaluationContext) throws -> FiniteThrust
    public init(frame: Frame, command: @escaping (SpacecraftState, ForceEvaluationContext) throws -> FiniteThrust) throws {
        guard frame.pseudoInertial else { throw SpacecraftDynamicsError.nonInertialFrame }
        self.frame = frame; self.command = command
    }
    public func contribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
        try engine(state, context).contribution(for: state, in: context)
    }
    public func depletionBoundaryContribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
        try engine(state, context).depletionBoundaryContribution(for: state, in: context)
    }
    private func engine(_ state: SpacecraftState, _ context: ForceEvaluationContext) throws -> FiniteThrust {
        guard state.frame === frame, context.frame === frame else { throw SpacecraftDynamicsError.frameMismatch }
        let result = try command(state, context)
        guard result.frame === frame else { throw SpacecraftDynamicsError.frameMismatch }
        return result
    }
}
