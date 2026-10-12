// SPDX-License-Identifier: Apache-2.0
import Math
import Testing
@testable import CelMek

@Suite struct SpacecraftDynamicsTests {
    private struct Acceleration: ForceModel {
        let value: SIMD3<Double>
        var massRate = 0.0
        func contribution(for state: SpacecraftState,in context: ForceEvaluationContext) throws -> ForceContribution {
            try ForceContribution(epoch: state.epoch,frame: state.frame,
                accelerationMetresPerSecondSquared: value,propellantMassRateKilogramsPerSecond: massRate)
        }
    }
    private struct Thrust: ForceModel {
        func contribution(for state: SpacecraftState,in context: ForceEvaluationContext) throws -> ForceContribution {
            try ForceContribution(epoch: state.epoch,frame: state.frame,
                accelerationMetresPerSecondSquared: SIMD3(100,0,0)/state.mass.totalMassKilograms,
                propellantMassRateKilogramsPerSecond: -0.1)
        }
    }
    private func state(_ frame: Frame, propellant: Double = 50) throws -> SpacecraftState {
        SpacecraftState(cartesian: try CartesianState(epoch: .j2000,frame: frame,
            position: SIMD3(7e6,0,0),velocity: SIMD3(0,7500,0),acceleration: SIMD3(999,999,999)),
            mass: try SpacecraftMassProperties(dryMassKilograms: 50,propellantMassKilograms: propellant))
    }
    @Test func compositionIsStateDependentAndDoesNotDoubleCountAcceleration() throws {
        let frame = Frame(name: "Earth inertial",pseudoInertial: true)
        let context = try ForceEvaluationContext(frame: frame,origin: .earth)
        let dynamics = try SpacecraftDynamics(models: [Acceleration(value: SIMD3(0,-8,0)),Thrust()],context: context)
        let heavy = try state(frame), lighter = try state(frame,propellant: 0)
        let derivative = try dynamics.derivative(at: heavy)
        #expect(derivative.positionRateMetresPerSecond == SIMD3(0,7500,0))
        #expect(derivative.velocityRateMetresPerSecondSquared == SIMD3(1,-8,0))
        #expect(derivative.propellantMassRateKilogramsPerSecond == -0.1)
        #expect(derivative.frame === frame && derivative.epoch == heavy.epoch)
        #expect(try dynamics.derivative(at: lighter).velocityRateMetresPerSecondSquared == SIMD3(2,-8,0))
        #expect(try dynamics.derivative(at: heavy).velocityRateMetresPerSecondSquared == derivative.velocityRateMetresPerSecondSquared)
        #expect(heavy.mass.propellantMassKilograms == 50)
        let empty = try SpacecraftDynamics(models: [],context: context).derivative(at: heavy)
        #expect(empty.velocityRateMetresPerSecondSquared == .zero && empty.propellantMassRateKilogramsPerSecond == 0)
    }
    @Test func physicalInputValidation() throws {
        for dry in [0.0,-1,.nan,.infinity] {
            #expect(throws: SpacecraftDynamicsError.invalidMass) { try SpacecraftMassProperties(dryMassKilograms: dry) }
        }
        for propellant in [-1.0,.nan,.infinity] {
            #expect(throws: SpacecraftDynamicsError.invalidMass) {
                try SpacecraftMassProperties(dryMassKilograms: 1,propellantMassKilograms: propellant)
            }
        }
        #expect(throws: SpacecraftDynamicsError.invalidMass) {
            try SpacecraftMassProperties(dryMassKilograms: Double.greatestFiniteMagnitude,
                propellantMassKilograms: Double.greatestFiniteMagnitude)
        }
        #expect(try SpacecraftMassProperties(dryMassKilograms: 10).totalMassKilograms == 10)
        let frame = Frame(name: "test",pseudoInertial: true)
        #expect(throws: SpacecraftDynamicsError.nonFiniteAcceleration) {
            try ForceContribution(epoch: .j2000,frame: frame,accelerationMetresPerSecondSquared: SIMD3(.nan,0,0))
        }
        #expect(throws: SpacecraftDynamicsError.invalidMassRate) {
            try ForceContribution(epoch: .j2000,frame: frame,accelerationMetresPerSecondSquared: .zero,
                propellantMassRateKilogramsPerSecond: .infinity)
        }
    }
    private enum Failure: Error { case unavailableEnvironment }
    private struct Broken: ForceModel {
        enum Mode { case frame, epoch, environment }
        let mode: Mode
        func contribution(for state: SpacecraftState,in context: ForceEvaluationContext) throws -> ForceContribution {
            if mode == .environment { throw Failure.unavailableEnvironment }
            return try ForceContribution(epoch: mode == .epoch ? state.epoch.advanced(by: Duration(seconds: 1)) : state.epoch,
                frame: mode == .frame ? Frame(name: "wrong",pseudoInertial: true) : state.frame,
                accelerationMetresPerSecondSquared: .zero)
        }
    }
    @Test func contractFailuresAndOverflowRemainVisible() throws {
        let frame = Frame(name: "test",pseudoInertial: true), context = try ForceEvaluationContext(frame: frame,origin: .earth)
        let input = try state(frame)
        #expect(throws: SpacecraftDynamicsError.frameMismatch) {
            try SpacecraftDynamics(models: [Broken(mode: .frame)],context: context).derivative(at: input)
        }
        #expect(throws: SpacecraftDynamicsError.epochMismatch) {
            try SpacecraftDynamics(models: [Broken(mode: .epoch)],context: context).derivative(at: input)
        }
        #expect(throws: Failure.unavailableEnvironment) {
            try SpacecraftDynamics(models: [Broken(mode: .environment)],context: context).derivative(at: input)
        }
        #expect(throws: SpacecraftDynamicsError.frameMismatch) {
            try SpacecraftDynamics(models: [],context: context).derivative(at: state(Frame(name: "other",pseudoInertial: true)))
        }
        #expect(throws: SpacecraftDynamicsError.nonInertialFrame) {
            try SpacecraftDynamics(models: [],context: ForceEvaluationContext(frame: Frame(name: "rotating",pseudoInertial: false),origin: .earth))
        }
        let huge = Acceleration(value: SIMD3(Double.greatestFiniteMagnitude,0,0))
        #expect(throws: SpacecraftDynamicsError.numericalRangeExceeded) {
            try SpacecraftDynamics(models: [huge,huge],context: context).derivative(at: input)
        }
        let mass = Acceleration(value: .zero,massRate: Double.greatestFiniteMagnitude)
        #expect(throws: SpacecraftDynamicsError.numericalRangeExceeded) {
            try SpacecraftDynamics(models: [mass,mass],context: context).derivative(at: input)
        }
    }
}
