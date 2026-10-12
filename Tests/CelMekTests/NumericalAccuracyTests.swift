// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct NumericalAccuracyTests {
    private struct Exponential: ForceModel {
        func contribution(for state: SpacecraftState,in context: ForceEvaluationContext) throws -> ForceContribution {
            try ForceContribution(epoch: state.epoch,frame: state.frame,accelerationMetresPerSecondSquared: state.cartesian.velocity)
        }
    }
    @Test func presetsTightenAndKeepWorkControlsExplicit() {
        let levels = AccuracyLevel.allCases, presets = levels.map { NumericalPropagationOptions.preset(for: $0) }
        #expect(presets[2] == (try? NumericalPropagationOptions()))
        for (loose,tight) in zip(presets,presets.dropFirst()) {
            #expect(tight.absolutePositionToleranceMetres < loose.absolutePositionToleranceMetres)
            #expect(tight.absoluteVelocityToleranceMetresPerSecond < loose.absoluteVelocityToleranceMetresPerSecond)
            #expect(tight.absolutePropellantToleranceKilograms < loose.absolutePropellantToleranceKilograms)
            #expect(tight.relativeTolerance < loose.relativeTolerance)
            #expect(tight.minimumStepSeconds == loose.minimumStepSeconds && tight.maximumStepSeconds == loose.maximumStepSeconds)
            #expect(tight.maximumStepAttempts == loose.maximumStepAttempts)
        }
    }
    @Test func hintsAndAuthoritativeOverrides() throws {
        let frame = Frame(name: "inertial",pseudoInertial: true)
        let seed = SpacecraftState(cartesian: try CartesianState(epoch: .j2000,frame: frame,position: .zero,velocity: SIMD3(1,0,0)),
            mass: try SpacecraftMassProperties(dryMassKilograms: 1))
        let dynamics = try SpacecraftDynamics(models: [Exponential()],context: ForceEvaluationContext(frame: frame,origin: .earth))
        let target = try Epoch.j2000.advanced(by: Duration(seconds: 3))
        var errors: [Double] = [], counts: [Int] = []
        for level in AccuracyLevel.allCases {
            let provider = try NumericalSpacecraftPropagator(initialState: seed,dynamics: dynamics,accuracy: level)
            #expect(provider.accuracy == level && provider.effectiveAccuracy == level)
            #expect(provider.options == .preset(for: level))
            let result = try provider.propagation(at: target)
            errors.append(abs(result.state.cartesian.velocity.x-exp(3)))
            counts.append(result.diagnostics.derivativeEvaluations)
            #expect(result.diagnostics.maximumAcceptedNormalizedError <= 1)
        }
        #expect(errors[3] < errors[0]*0.01)
        #expect(counts[3] > counts[0])
        let custom = try NumericalPropagationOptions(absolutePositionToleranceMetres: 1e-8,
            absoluteVelocityToleranceMetresPerSecond: 1e-10,relativeTolerance: 1e-12,maximumStepSeconds: 20)
        let override = try NumericalSpacecraftPropagator(initialState: seed,dynamics: dynamics,accuracy: .low,options: custom)
        let direct = try NumericalSpacecraftPropagator(initialState: seed,dynamics: dynamics,options: custom)
        #expect(override.accuracy == .low && override.effectiveAccuracy == nil)
        #expect(direct.accuracy == nil && direct.effectiveAccuracy == nil)
        #expect(override.options == custom)
        #expect(try override.state(at: target) == direct.state(at: target))
    }
}
