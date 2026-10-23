// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct NumericalPropagationTests {
    private struct Constant: ForceModel {
        let acceleration: SIMD3<Double>
        var massRate = 0.0
        func contribution(for state: SpacecraftState,in context: ForceEvaluationContext) throws -> ForceContribution {
            try ForceContribution(epoch: state.epoch,frame: state.frame,accelerationMetresPerSecondSquared: acceleration,
                propellantMassRateKilogramsPerSecond: massRate)
        }
    }
    private struct Exponential: ForceModel {
        func contribution(for state: SpacecraftState,in context: ForceEvaluationContext) throws -> ForceContribution {
            try ForceContribution(epoch: state.epoch,frame: state.frame,
                accelerationMetresPerSecondSquared: state.cartesian.velocity)
        }
    }
    private func norm(_ v: SIMD3<Double>) -> Double { hypot(hypot(v.x,v.y),v.z) }
    private func epoch(_ t: Double) throws -> Epoch { try .j2000.advanced(by: Duration(seconds: t)) }
    private func setup(_ models: [any ForceModel], options: NumericalPropagationOptions = try! NumericalPropagationOptions(),
                       propellant: Double = 10) throws -> NumericalSpacecraftPropagator {
        let frame = Frame(name: "test inertial",pseudoInertial: true)
        let seed = SpacecraftState(cartesian: try CartesianState(epoch: .j2000,frame: frame,
            position: SIMD3(1,2,3),velocity: SIMD3(1,0,0)),
            mass: try SpacecraftMassProperties(dryMassKilograms: 100,propellantMassKilograms: propellant))
        return try NumericalSpacecraftPropagator(initialState: seed,
            dynamics: SpacecraftDynamics(models: models,context: ForceEvaluationContext(frame: frame,origin: .earth)),options: options)
    }
    @Test func polynomialMotionMassAndExactEpochs() throws {
        let acceleration = SIMD3<Double>(2,-3,4)
        let model = try setup([Constant(acceleration: acceleration,massRate: -0.01)])
        for seconds in [100.0,-100,0,15.125,100] {
            let date = try epoch(seconds), result = try model.propagation(at: date)
            let p = model.initialState.cartesian.position + model.initialState.cartesian.velocity*seconds + acceleration*(0.5*seconds*seconds)
            #expect(norm(result.state.cartesian.position-p) < 1e-8)
            #expect(norm(result.state.cartesian.velocity-model.initialState.cartesian.velocity-acceleration*seconds) < 1e-10)
            #expect(abs(result.state.mass.propellantMassKilograms-(10-0.01*seconds)) < 1e-12)
            #expect(result.state.cartesian.acceleration == acceleration)
            #expect(result.state.epoch == date && result.state.frame === model.frame)
            #expect(result.diagnostics.maximumAcceptedNormalizedError <= 1)
            #expect(result.diagnostics.derivativeEvaluations == 7*(result.diagnostics.acceptedSteps+result.diagnostics.rejectedSteps)+1)
        }
        let zero = try model.propagation(at: .j2000)
        #expect(zero.diagnostics.acceptedSteps == 0 && zero.diagnostics.derivativeEvaluations == 1)
        #expect(zero.diagnostics.smallestAcceptedStepSeconds == nil)
        let free = try setup([]), date = try epoch(50)
        #expect(norm(try free.state(at: date).position-SIMD3(51,2,3)) < 1e-12)
        #expect(try free.spacecraftState(at: date).mass == free.initialState.mass)
        let samples = try TrajectorySampler.sample(model,at: [date,.j2000,date])
        #expect(samples[0] == samples[2])
    }
    @Test func orbitalReferenceAndToleranceConvergence() throws {
        let frame = Frame(name: "Earth inertial",pseudoInertial: true), mu = 3.986004418e14
        let seed = SpacecraftState(cartesian: try CartesianState(epoch: .j2000,frame: frame,
            position: SIMD3(7e6,0,0),velocity: SIMD3(0,7500,0)),mass: try SpacecraftMassProperties(dryMassKilograms: 100))
        let gm = try GravitationalParameter(cubicMetresPerSecondSquared: mu,convention: .newtonian,
            provenance: PhysicalDataProvenance(source: "Test Earth GM",version: "1"))
        let dynamics = try SpacecraftDynamics(models: [CentralGravity(body: .earth,gravitationalParameter: gm)],
            context: ForceEvaluationContext(frame: frame,origin: .earth))
        let analytic = try TwoBodyPropagator(initialState: seed.cartesian,centralBody: .earth,gravitationalParameter: gm)
        let loose = try NumericalSpacecraftPropagator(initialState: seed,dynamics: dynamics,
            options: NumericalPropagationOptions(absolutePositionToleranceMetres: 100,absoluteVelocityToleranceMetresPerSecond: 0.1,
                relativeTolerance: 1e-4,initialStepSeconds: 1000,maximumStepSeconds: 1000))
        let fine = try NumericalSpacecraftPropagator(initialState: seed,dynamics: dynamics,
            options: NumericalPropagationOptions(absolutePositionToleranceMetres: 1e-4,absoluteVelocityToleranceMetresPerSecond: 1e-7,
                relativeTolerance: 1e-11,initialStepSeconds: 1000,maximumStepSeconds: 1000))
        for seconds in [3600.0,-3600] {
            let date = try epoch(seconds), reference = try analytic.state(at: date)
            let coarse = try loose.propagation(at: date), accurate = try fine.propagation(at: date)
            let coarseError = norm(coarse.state.cartesian.position-reference.position)
            let fineError = norm(accurate.state.cartesian.position-reference.position)
            #expect(fineError < 0.1 && fineError < coarseError*0.01)
            #expect(norm(accurate.state.cartesian.velocity-reference.velocity) < 1e-4)
            #expect(accurate.diagnostics.acceptedSteps > coarse.diagnostics.acceptedSteps)
            #expect(accurate.diagnostics.rejectedSteps > 0)
        }
        let end = try fine.spacecraftState(at: epoch(3600))
        let reverse = try NumericalSpacecraftPropagator(initialState: end,dynamics: dynamics,options: fine.options)
        #expect(norm(try reverse.state(at: .j2000).position-seed.cartesian.position) < 0.2)
    }
    @Test func adaptiveErrorAndBoundedWork() throws {
        let options = try NumericalPropagationOptions(absolutePositionToleranceMetres: 1e-12,
            absoluteVelocityToleranceMetresPerSecond: 1e-12,relativeTolerance: 1e-12,
            initialStepSeconds: 1,maximumStepSeconds: 1)
        let model = try setup([Exponential()],options: options)
        let result = try model.propagation(at: epoch(1))
        #expect(abs(result.state.cartesian.velocity.x-exp(1)) < 1e-10)
        #expect(abs(result.state.cartesian.position.x-exp(1)) < 1e-10)
        #expect(result.diagnostics.rejectedSteps > 0)
        #expect(result.diagnostics.largestAcceptedStepSeconds! <= options.maximumStepSeconds)
        let limited = try setup([],options: NumericalPropagationOptions(initialStepSeconds: 1,maximumStepSeconds: 1,maximumStepAttempts: 1))
        #expect(throws: NumericalPropagationError.stepLimitExceeded) { try limited.state(at: epoch(2)) }
        let evaluations = try setup([],options: NumericalPropagationOptions(maximumDerivativeEvaluations: 1))
        #expect(throws: NumericalPropagationError.evaluationLimitExceeded) { try evaluations.state(at: epoch(1)) }
        let minimum = try setup([Exponential()],options: NumericalPropagationOptions(absolutePositionToleranceMetres: 1e-14,
            absoluteVelocityToleranceMetresPerSecond: 1e-14,relativeTolerance: 0,initialStepSeconds: 1,minimumStepSeconds: 1,maximumStepSeconds: 1))
        #expect(throws: NumericalPropagationError.minimumStepReached) { try minimum.state(at: epoch(1)) }
        let short = try setup([],options: NumericalPropagationOptions(initialStepSeconds: 1,minimumStepSeconds: 1,maximumStepSeconds: 1))
        #expect(try short.propagation(at: epoch(0.125)).diagnostics.smallestAcceptedStepSeconds == 0.125)
    }
    private enum Failure: Error { case environment }
    private struct Broken: ForceModel {
        func contribution(for state: SpacecraftState,in context: ForceEvaluationContext) throws -> ForceContribution { throw Failure.environment }
    }
    @Test func invalidOptionsDepletionAndProviderFailure() throws {
        #expect(throws: NumericalPropagationError.invalidOptions) { try NumericalPropagationOptions(relativeTolerance: -.infinity) }
        #expect(throws: NumericalPropagationError.invalidOptions) { try NumericalPropagationOptions(absolutePositionToleranceMetres: 0) }
        #expect(throws: NumericalPropagationError.invalidOptions) { try NumericalPropagationOptions(initialStepSeconds: 1000,maximumStepSeconds: 10) }
        #expect(throws: NumericalPropagationError.invalidOptions) { try NumericalPropagationOptions(maximumStepAttempts: 0) }
        let depleted = try setup([Constant(acceleration: .zero,massRate: -1)],propellant: 0)
        #expect(throws: NumericalPropagationError.propellantDepleted) { try depleted.state(at: epoch(1)) }
        let consuming = try setup([Constant(acceleration: .zero,massRate: -1)],propellant: 1)
        #expect(throws: NumericalPropagationError.propellantDepleted) { try consuming.state(at: epoch(2)) }
        #expect(consuming.initialState.mass.propellantMassKilograms == 1)
        let broken = try setup([Broken()])
        #expect(throws: Failure.environment) { try broken.state(at: .j2000) }
        #expect(throws: Failure.environment) { try broken.state(at: epoch(1)) }
    }
    @Test func collapsedStageEpochsAndSeedFramesAreRejected() throws {
        let frame = Frame(name: "inertial",pseudoInertial: true)
        let start = try epoch(0.5)
        let seed = SpacecraftState(cartesian: try CartesianState(epoch: start,frame: frame,position: .zero,velocity: .zero),
            mass: try SpacecraftMassProperties(dryMassKilograms: 1))
        let dynamics = try SpacecraftDynamics(models: [],context: ForceEvaluationContext(frame: frame,origin: .earth))
        let model = try NumericalSpacecraftPropagator(initialState: seed,dynamics: dynamics)
        let target = try start.advanced(by: Duration(seconds: 1e-16))
        #expect(target > start)
        #expect(throws: NumericalPropagationError.epochResolutionExceeded) { try model.state(at: target) }
        let other = try SpacecraftDynamics(models: [],context: ForceEvaluationContext(
            frame: Frame(name: "other",pseudoInertial: true),origin: .earth))
        #expect(throws: SpacecraftDynamicsError.frameMismatch) {
            try NumericalSpacecraftPropagator(initialState: seed,dynamics: other)
        }
    }

    @Test func fractionalStepRemaindersReachExactEndpointWithoutCollapsedStages() throws {
        let frame=Frame(name:"fractional endpoint",pseudoInertial:true)
        let dynamics=try SpacecraftDynamics(models:[Constant(acceleration:SIMD3(2,0,0))],context:ForceEvaluationContext(frame:frame,origin:.earth))
        for fraction in [0.5,0.2693127683618036,0.99] {
            for sign in [1.0,-1.0] {
                let start=Epoch(offsetFromJ2000:try Duration(wholeSeconds:800000000,fractionalSecond:fraction))
                let seed=SpacecraftState(cartesian:try CartesianState(epoch:start,frame:frame,position:.zero,velocity:SIMD3(3,0,0)),mass:try SpacecraftMassProperties(dryMassKilograms:1))
                let options=try NumericalPropagationOptions(initialStepSeconds:0.1,maximumStepSeconds:2)
                let model=try NumericalSpacecraftPropagator(initialState:seed,dynamics:dynamics,options:options)
                let target=try start.advanced(by:Duration(seconds:sign*0.6))
                let result=try model.propagation(at:target)
                #expect(result.state.epoch==target)
                #expect(abs(result.state.cartesian.position.x-(3*sign*0.6+0.36))<1e-12)
                #expect(abs(result.state.cartesian.velocity.x-(3+2*sign*0.6))<1e-12)
            }
        }
    }

}
