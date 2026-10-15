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

public struct FuelLimitedPropagationResult {
    public let state: SpacecraftState
    public let depletionEpoch: Epoch?
    public let depletionBracketSeconds: Double?
}

/// Variable-flow depletion search. Consumers explicitly supply a dry-mass limit;
/// other models must have zero mass rate. All consumers switch off together.
public struct FuelLimitedNumericalPropagator: SpacecraftStateProvider {
    public let initialState: SpacecraftState
    public let dynamics: SpacecraftDynamics
    public let numericalOptions: NumericalPropagationOptions
    public let eventOptions: PropagationEventOptions
    public var frame: Frame { initialState.frame }
    public init(initialState: SpacecraftState, dynamics: SpacecraftDynamics,
                numericalOptions: NumericalPropagationOptions = try! NumericalPropagationOptions(),
                eventOptions: PropagationEventOptions = try! PropagationEventOptions()) throws {
        _ = try NumericalSpacecraftPropagator(initialState: initialState, dynamics: dynamics, options: numericalOptions)
        self.initialState = initialState; self.dynamics = dynamics
        self.numericalOptions = numericalOptions; self.eventOptions = eventOptions
    }
    public func state(at epoch: Epoch) throws -> CartesianState { try spacecraftState(at: epoch).cartesian }
    public func spacecraftState(at epoch: Epoch) throws -> SpacecraftState { try propagation(at: epoch).state }
    public func propagation(at target: Epoch) throws -> FuelLimitedPropagationResult {
        var queries = 0
        func consumeQueryBudget() throws {
            guard queries < eventOptions.maximumStateEvaluations else { throw PropagationEventError.evaluationLimitExceeded }
            queries += 1
        }
        if target < initialState.epoch {
            try consumeQueryBudget()
            return FuelLimitedPropagationResult(state: try NumericalSpacecraftPropagator(initialState: initialState,
                dynamics: dynamics, options: numericalOptions).spacecraftState(at: target), depletionEpoch: nil, depletionBracketSeconds: nil)
        }
        let coast = try dynamics.withoutFuelConsumers()
        if initialState.mass.propellantMassKilograms == 0 {
            try consumeQueryBudget()
            return FuelLimitedPropagationResult(state: try NumericalSpacecraftPropagator(initialState: initialState,
                dynamics: coast, options: numericalOptions).spacecraftState(at: target), depletionEpoch: initialState.epoch, depletionBracketSeconds: 0)
        }
        var bracket: (SpacecraftState,SpacecraftState)?
        try consumeQueryBudget()
        let integrated = try NumericalSpacecraftPropagator(initialState: initialState, dynamics: dynamics, options: numericalOptions)
            .propagation(at: target, acceptedStep: nil, extendFuelBoundary: true) { left, leftFuel, right, rightFuel in
                if leftFuel > 0 && rightFuel <= 0 { bracket = (left,right); return false }
                return true
            }
        guard var (left,right) = bracket else {
            return FuelLimitedPropagationResult(state: integrated.state, depletionEpoch: nil, depletionBracketSeconds: nil)
        }
        var iterations = 0
        var width = try right.epoch.duration(since: left.epoch).seconds
        while width > eventOptions.timeToleranceSeconds {
            guard iterations < eventOptions.maximumLocalizationIterations else { throw PropagationEventError.localizationLimitExceeded }
            try consumeQueryBudget()
            iterations += 1
            let date = try left.epoch.advanced(by: Duration(seconds: width/2))
            guard date > left.epoch && date < right.epoch else { throw PropagationEventError.epochResolutionExceeded }
            let trial = try NumericalSpacecraftPropagator(initialState: left, dynamics: dynamics, options: numericalOptions)
                .propagation(at: date, acceptedStep: nil, extendFuelBoundary: true)
            if trial.rawPropellantMassKilograms > 0 { left = trial.state } else { right = trial.state }
            width = try right.epoch.duration(since: left.epoch).seconds
        }
        // Select the last valid side, at most timeTolerance before exhaustion;
        // replace only its bounded residual fuel with the known boundary value.
        let boundary = SpacecraftState(cartesian: left.cartesian,
            mass: try SpacecraftMassProperties(dryMassKilograms: left.mass.dryMassKilograms))
        try consumeQueryBudget()
        let final = try NumericalSpacecraftPropagator(initialState: boundary, dynamics: coast, options: numericalOptions)
            .spacecraftState(at: target)
        return FuelLimitedPropagationResult(state: final, depletionEpoch: boundary.epoch, depletionBracketSeconds: width)
    }
}
