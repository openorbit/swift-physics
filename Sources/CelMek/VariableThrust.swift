//
// SPDX-License-Identifier: Apache-2.0
//
// Copyright 2026 Mattias Holm
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
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
