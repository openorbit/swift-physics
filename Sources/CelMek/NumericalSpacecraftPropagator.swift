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
import Foundation
import Math

public enum NumericalPropagationError: Error, Equatable {
    case invalidOptions, stepLimitExceeded, evaluationLimitExceeded
    case minimumStepReached, epochResolutionExceeded, propellantDepleted, numericalRangeExceeded
}

/// Local embedded-error tolerances, not certified global/physical accuracy bounds.
public struct NumericalPropagationOptions: Equatable, Sendable {
    public let absolutePositionToleranceMetres: Double
    public let absoluteVelocityToleranceMetresPerSecond: Double
    public let absolutePropellantToleranceKilograms: Double
    public let relativeTolerance: Double
    public let initialStepSeconds, minimumStepSeconds, maximumStepSeconds: Double
    public let maximumStepAttempts, maximumDerivativeEvaluations: Int
    public init(absolutePositionToleranceMetres: Double = 1e-3,
                absoluteVelocityToleranceMetresPerSecond: Double = 1e-6,
                absolutePropellantToleranceKilograms: Double = 1e-9,
                relativeTolerance: Double = 1e-9, initialStepSeconds: Double = 10,
                minimumStepSeconds: Double = 1e-6, maximumStepSeconds: Double = 300,
                maximumStepAttempts: Int = 100_000, maximumDerivativeEvaluations: Int = 700_001) throws {
        guard [absolutePositionToleranceMetres,absoluteVelocityToleranceMetresPerSecond,
               absolutePropellantToleranceKilograms,initialStepSeconds,minimumStepSeconds,maximumStepSeconds]
            .allSatisfy({ $0.isFinite && $0 > 0 }), relativeTolerance.isFinite, relativeTolerance >= 0,
              minimumStepSeconds <= initialStepSeconds, initialStepSeconds <= maximumStepSeconds,
              maximumStepAttempts > 0, maximumDerivativeEvaluations > 0 else {
            throw NumericalPropagationError.invalidOptions
        }
        self.absolutePositionToleranceMetres = absolutePositionToleranceMetres
        self.absoluteVelocityToleranceMetresPerSecond = absoluteVelocityToleranceMetresPerSecond
        self.absolutePropellantToleranceKilograms = absolutePropellantToleranceKilograms
        self.relativeTolerance = relativeTolerance; self.initialStepSeconds = initialStepSeconds
        self.minimumStepSeconds = minimumStepSeconds; self.maximumStepSeconds = maximumStepSeconds
        self.maximumStepAttempts = maximumStepAttempts; self.maximumDerivativeEvaluations = maximumDerivativeEvaluations
    }
}

extension NumericalPropagationOptions {
    /// Version 1 local-error presets. These change numerical tolerances only;
    /// every level evaluates all supplied force models with the same RK method.
    public static func preset(for accuracy: AccuracyLevel) -> Self {
        let position, velocity, mass, relative: Double
        switch accuracy {
        case .low: (position,velocity,mass,relative) = (10,1e-2,1e-4,1e-5)
        case .balanced: (position,velocity,mass,relative) = (1e-1,1e-4,1e-6,1e-7)
        case .high: (position,velocity,mass,relative) = (1e-3,1e-6,1e-9,1e-9)
        case .full: (position,velocity,mass,relative) = (1e-5,1e-8,1e-11,1e-11)
        }
        return try! Self(absolutePositionToleranceMetres: position,
            absoluteVelocityToleranceMetresPerSecond: velocity,
            absolutePropellantToleranceKilograms: mass,relativeTolerance: relative)
    }
}

public struct NumericalPropagationDiagnostics {
    public let acceptedSteps, rejectedSteps, derivativeEvaluations: Int
    public let smallestAcceptedStepSeconds, largestAcceptedStepSeconds: Double?
    public let maximumAcceptedNormalizedError: Double
}

public struct NumericalSpacecraftPropagationResult {
    public let state: SpacecraftState
    public let diagnostics: NumericalPropagationDiagnostics
    let rawPropellantMassKilograms: Double
}

/// Native Swift Dormand–Prince 5(4), integrating position/velocity/propellant.
/// Every query restarts from the immutable seed, including backward queries.
/// Internal adaptive steps are independent of renderer cadence; no dense output,
/// force discontinuity detection or propellant-depletion event localization yet.
public struct NumericalSpacecraftPropagator: StateProvider {
    public let initialState: SpacecraftState
    public let dynamics: SpacecraftDynamics
    public let options: NumericalPropagationOptions
    /// Requested hint, nil for the direct-options initializer.
    public let accuracy: AccuracyLevel?
    /// Selected preset, nil when explicit settings override it. This reports a
    /// configuration, not achieved trajectory accuracy or a physical error bound.
    public let effectiveAccuracy: AccuracyLevel?
    public var frame: Frame { initialState.frame }
    public var referenceEpoch: Epoch { initialState.epoch }
    public init(initialState: SpacecraftState, dynamics: SpacecraftDynamics,
                options: NumericalPropagationOptions = try! NumericalPropagationOptions()) throws {
        guard initialState.frame === dynamics.context.frame else { throw SpacecraftDynamicsError.frameMismatch }
        self.initialState = initialState; self.dynamics = dynamics; self.options = options
        accuracy = nil; effectiveAccuracy = nil
    }
    /// Explicit options are authoritative; the requested hint remains visible,
    /// while effectiveAccuracy is unspecified for a custom configuration.
    public init(initialState: SpacecraftState, dynamics: SpacecraftDynamics,
                accuracy: AccuracyLevel, options: NumericalPropagationOptions? = nil) throws {
        guard initialState.frame === dynamics.context.frame else { throw SpacecraftDynamicsError.frameMismatch }
        self.initialState = initialState; self.dynamics = dynamics
        self.options = options ?? .preset(for: accuracy)
        self.accuracy = accuracy; effectiveAccuracy = options == nil ? accuracy : nil
    }
    public func state(at epoch: Epoch) throws -> CartesianState { try propagation(at: epoch).state.cartesian }
    public func spacecraftState(at epoch: Epoch) throws -> SpacecraftState { try propagation(at: epoch).state }

    public func propagation(at target: Epoch) throws -> NumericalSpacecraftPropagationResult {
        try propagation(at: target, acceptedStep: nil)
    }
    /// Internal observer receives only error-controlled accepted steps. Returning
    /// false ends this integration at the accepted endpoint for event restart.
    func propagation(at target: Epoch,
                     acceptedStep: ((SpacecraftState, SpacecraftState) throws -> Bool)?,
                     extendFuelBoundary: Bool = false,
                     fuelStep: ((SpacecraftState, Double, SpacecraftState, Double) throws -> Bool)? = nil) throws -> NumericalSpacecraftPropagationResult {
        var current = referenceEpoch
        let forward = target >= current
        var y = [initialState.cartesian.position.x,initialState.cartesian.position.y,initialState.cartesian.position.z,
                 initialState.cartesian.velocity.x,initialState.cartesian.velocity.y,initialState.cartesian.velocity.z,
                 initialState.mass.propellantMassKilograms]
        var accepted = 0, rejected = 0, evaluations = 0
        var smallest: Double?, largest: Double?, maxError = 0.0
        var nextStep = options.initialStepSeconds
        func state(_ values: [Double], at epoch: Epoch, acceleration: SIMD3<Double>? = nil) throws -> SpacecraftState {
            guard values.allSatisfy(\.isFinite) else { throw NumericalPropagationError.numericalRangeExceeded }
            return SpacecraftState(cartesian: try CartesianState(epoch: epoch,frame: frame,
                position: SIMD3(values[0],values[1],values[2]),velocity: SIMD3(values[3],values[4],values[5]),acceleration: acceleration),
                mass: try SpacecraftMassProperties(dryMassKilograms: initialState.mass.dryMassKilograms,propellantMassKilograms: extendFuelBoundary ? max(0,values[6]) : values[6]))
        }
        func derivative(_ values: [Double], at epoch: Epoch) throws -> [Double] {
            guard evaluations < options.maximumDerivativeEvaluations else { throw NumericalPropagationError.evaluationLimitExceeded }
            evaluations += 1
            let local = try state(values,at: epoch)
            let d = try extendFuelBoundary ? dynamics.fuelBoundaryDerivative(at: local) : dynamics.derivative(at: local)
            return [d.positionRateMetresPerSecond.x,d.positionRateMetresPerSecond.y,d.positionRateMetresPerSecond.z,
                    d.velocityRateMetresPerSecondSquared.x,d.velocityRateMetresPerSecondSquared.y,
                    d.velocityRateMetresPerSecondSquared.z,d.propellantMassRateKilogramsPerSecond]
        }
        let c: [Double] = [0.0,1.0/5,3.0/10,4.0/5,8.0/9,1,1]
        let a: [[Double]] = [[],[1.0/5],[3.0/40,9.0/40],[44.0/45,-56.0/15,32.0/9],
            [19372.0/6561,-25360.0/2187,64448.0/6561,-212.0/729],
            [9017.0/3168,-355.0/33,46732.0/5247,49.0/176,-5103.0/18656],
            [35.0/384,0,500.0/1113,125.0/192,-2187.0/6784,11.0/84]]
        let errorWeights: [Double] = [71.0/57600,0,-71.0/16695,71.0/1920,-17253.0/339200,22.0/525,-1.0/40]
        while current != target {
            guard accepted+rejected < options.maximumStepAttempts else { throw NumericalPropagationError.stepLimitExceeded }
            let remaining = try forward ? target.duration(since: current) : current.duration(since: target)
            let magnitude = min(nextStep,remaining.seconds)
            let final = magnitude >= remaining.seconds
            var end = final ? target : try current.advanced(by: Duration(seconds: forward ? magnitude : -magnitude))
            // Split epochs canonicalize their fractional seconds in 0..<1. A
            // normal step can land a few rounding units before the endpoint
            // (e.g. 0.1 + 0.5 on a 0.6 s query), leaving an unresolvable RK
            // remainder. Integrate that rounding sliver with this step instead
            // of dropping elapsed time or fabricating a state at a later epoch.
            // Do not merge genuinely tiny requested/attempted steps: their stage
            // resolution checks below must still reject collapsed timestamps.
            let endpointRoundoffSeconds = 8*Double.ulpOfOne
            if !final, magnitude > 16*endpointRoundoffSeconds,
               abs(try target.duration(since:end).seconds) <= endpointRoundoffSeconds {
                end = target
            }
            guard forward ? end > current : end < current else { throw NumericalPropagationError.epochResolutionExceeded }
            let h = try end.duration(since: current).seconds
            var stages: [[Double]] = [], candidate = y, depleted = false
            for stage in 0..<7 {
                var values = y
                for component in 0..<7 {
                    var sum = 0.0
                    for j in 0..<stage { sum += a[stage][j]*stages[j][component] }
                    values[component] += h*sum
                }
                guard values.allSatisfy(\.isFinite) else { throw NumericalPropagationError.numericalRangeExceeded }
                if values[6] < 0 && !extendFuelBoundary { depleted = true; break }
                let date = stage == 0 ? current : c[stage] == 1 ? end : try current.advanced(by: Duration(seconds: c[stage]*h))
                if c[stage] > 0 && c[stage] < 1 && !(forward ? date > current && date < end : date < current && date > end) {
                    throw NumericalPropagationError.epochResolutionExceeded
                }
                stages.append(try derivative(values,at: date))
                if stage == 6 { candidate = values }
            }
            var error = Double.infinity
            if !depleted {
                error = 0
                for i in 0..<7 {
                    let absolute = i < 3 ? options.absolutePositionToleranceMetres : i < 6
                        ? options.absoluteVelocityToleranceMetresPerSecond : options.absolutePropellantToleranceKilograms
                    let scale = absolute + options.relativeTolerance*max(abs(y[i]),abs(candidate[i]))
                    var weighted = 0.0
                    for stage in 0..<7 { weighted += errorWeights[stage]*stages[stage][i] }
                    let value = abs(h*weighted)/scale
                    guard scale.isFinite, value.isFinite else { throw NumericalPropagationError.numericalRangeExceeded }
                    error = max(error,value)
                }
            }
            if error <= 1 {
                let previous = acceptedStep == nil && fuelStep == nil ? nil : try state(y, at: current)
                let previousFuel = y[6]
                y = candidate; current = end; accepted += 1
                smallest = min(smallest ?? abs(h),abs(h)); largest = max(largest ?? abs(h),abs(h))
                maxError = max(maxError,error)
                if let fuelStep, let previous,
                   try !fuelStep(previous, previousFuel, state(y, at: current), y[6]) { break }
                if let acceptedStep, let previous,
                   try !acceptedStep(previous, state(y, at: current)) { break }
                let factor = error == 0 ? 5 : min(5,max(0.2,0.9*pow(error,-0.2)))
                nextStep = max(options.minimumStepSeconds,min(options.maximumStepSeconds,abs(h)*factor))
            } else {
                rejected += 1
                guard magnitude > options.minimumStepSeconds else {
                    throw depleted ? NumericalPropagationError.propellantDepleted : NumericalPropagationError.minimumStepReached
                }
                let factor = depleted ? 0.5 : max(0.1,min(0.5,0.9*pow(error,-0.2)))
                nextStep = max(options.minimumStepSeconds,magnitude*factor)
            }
        }
        let finalDerivative = try derivative(y,at: current)
        let result = try state(y,at: current,acceleration: SIMD3(finalDerivative[3],finalDerivative[4],finalDerivative[5]))
        return NumericalSpacecraftPropagationResult(state: result,diagnostics: NumericalPropagationDiagnostics(
            acceptedSteps: accepted,rejectedSteps: rejected,derivativeEvaluations: evaluations,
            smallestAcceptedStepSeconds: smallest,largestAcceptedStepSeconds: largest,maximumAcceptedNormalizedError: maxError), rawPropellantMassKilograms: y[6])
    }
}
