// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct WaypointRoutePropagationTests {
    let frame = Frame(name: "Execution", pseudoInertial: true)
    func route(_ points: [SIMD3<Double>]) throws -> WaypointBrachistochroneRoute {
        try WaypointBrachistochroneRoute(departureEpoch: .j2000, frame: frame, origin: .sun,
            waypoints: points.map { try BrachistochroneWaypoint(positionMetres: $0) },
            legAccelerationsMetresPerSecondSquared: Array(repeating: 1, count: points.count-1))
    }
    func limits(_ acceleration: Double = 100, thrust: Double = 10_000) throws -> WaypointRouteEngineLimits {
        try WaypointRouteEngineLimits(maximumAccelerationMetresPerSecondSquared: acceleration,
            maximumThrustNewtons: thrust, effectiveExhaustVelocityMetresPerSecond: 100)
    }
    @Test func analyticFuelAndFlyThroughContinuity() throws {
        let plan = try route([.zero, SIMD3(4,0,0)])
        let background = try SpacecraftDynamics(models: [], context: ForceEvaluationContext(frame: frame, origin: .sun))
        let mass = try SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50)
        let execution = try WaypointRoutePropagator(route: plan, initialMass: mass, backgroundDynamics: background,
            legEngineLimits: [limits()], numericalOptions: .preset(for: .full))
        let result = try execution.propagation(at: plan.arrivalEpoch)
        // Cubic rest-to-rest speed peaks at 1.5*d/T; total |a| integral is 3*d/T.
        let expectedMass = 100*exp(-12/plan.duration.seconds/100)
        #expect(abs(result.state.mass.totalMassKilograms-expectedMass) < 1e-7)
        #expect(result.checkpoints.count == 1)
        #expect(result.checkpoints[0].positionErrorMetres < 1e-9)
        #expect(result.checkpoints[0].velocityErrorMetresPerSecond < 1e-9)

        let flyThrough = try route([.zero, SIMD3(4,0,0), SIMD3(8,2,0)])
        let multiple = try WaypointRoutePropagator(route: flyThrough, initialMass: mass,
            backgroundDynamics: background, legEngineLimits: [limits(), limits()], numericalOptions: .preset(for: .full))
        let boundary = try multiple.spacecraftState(at: flyThrough.legs[0].arrivalEpoch)
        let completed = try multiple.propagation(at: flyThrough.arrivalEpoch)
        #expect(boundary.cartesian.velocity.x > 0)
        #expect(boundary.mass == completed.checkpoints[0].state.mass)
        #expect(completed.state.mass.propellantMassKilograms < boundary.mass.propellantMassKilograms)
        #expect(completed.checkpoints.count == 2)
        #expect(try multiple.spacecraftState(at: flyThrough.arrivalEpoch) == completed.state)
    }
    @Test func compensationAndExplicitFailures() throws {
        let plan = try route([SIMD3(10,0,0), SIMD3(14,0,0)])
        let context = try ForceEvaluationContext(frame: frame, origin: .sun)
        let mass = try SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50)
        let empty = try SpacecraftDynamics(models: [], context: context)
        let environment = try SpacecraftDynamics(models: [ConstantAcceleration()], context: context)
        let execution = try WaypointRoutePropagator(route: plan, initialMass: mass, backgroundDynamics: environment,
            legEngineLimits: [limits()], numericalOptions: .preset(for: .full))
        let result = try execution.propagation(at: plan.arrivalEpoch)
        #expect(result.checkpoints[0].positionErrorMetres < 1e-9)
        // Environment is +2 m/s²; cubic |a| < 1, so compensation is always negative.
        // Integral |a-2| dt = 2*T since endpoint velocities are equal.
        #expect(abs(result.state.mass.totalMassKilograms-100*exp(-2*plan.duration.seconds/100)) < 1e-7)
        #expect(throws: WaypointRoutePropagationError.accelerationLimitExceeded(leg: 0)) {
            try WaypointRoutePropagator(route: plan, initialMass: mass, backgroundDynamics: environment,
                legEngineLimits: [limits(0.01)]).spacecraftState(at: plan.arrivalEpoch)
        }
        #expect(throws: WaypointRoutePropagationError.thrustLimitExceeded(leg: 0)) {
            try WaypointRoutePropagator(route: plan, initialMass: mass, backgroundDynamics: empty,
                legEngineLimits: [limits(thrust: 0.01)]).spacecraftState(at: plan.arrivalEpoch)
        }
        #expect(throws: NumericalPropagationError.propellantDepleted) {
            try WaypointRoutePropagator(route: plan, initialMass: SpacecraftMassProperties(dryMassKilograms: 50),
                backgroundDynamics: empty, legEngineLimits: [limits()]).spacecraftState(at: plan.departureEpoch)
        }
        #expect(throws: WaypointRoutePropagationError.backgroundMassFlow) {
            let inflow = try SpacecraftDynamics(models: [ConstantAcceleration(massRate: 1)], context: context)
            _ = try WaypointRoutePropagator(route: plan, initialMass: mass, backgroundDynamics: inflow,
                legEngineLimits: [limits()]).spacecraftState(at: plan.arrivalEpoch)
        }
    }
    @Test func centralGravityFuelAgainstIndependentQuadrature() throws {
        let plan = try route([SIMD3(10,0,0), SIMD3(14,0,0)])
        let gravity = try CentralGravity(body: .sun, gravitationalParameter: GravitationalParameter(
            cubicMetresPerSecondSquared: 1000, convention: .newtonian,
            provenance: PhysicalDataProvenance(source: "Synthetic test GM", version: "1")))
        let background = try SpacecraftDynamics(models: [gravity],
            context: ForceEvaluationContext(frame: frame, origin: .sun))
        let execution = try WaypointRoutePropagator(route: plan,
            initialMass: SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50),
            backgroundDynamics: background, legEngineLimits: [limits()], numericalOptions: NumericalPropagationOptions(
                absolutePositionToleranceMetres: 1e-8, absoluteVelocityToleranceMetresPerSecond: 1e-10,
                absolutePropellantToleranceKilograms: 1e-11, relativeTolerance: 1e-11,
                initialStepSeconds: 0.01, maximumStepSeconds: 0.1))
        let result = try execution.propagation(at: plan.arrivalEpoch)
        // Gravity dominates cubic acceleration throughout, so |a_thrust| = a_ref + GM/r².
        // Integral a_ref dt = 0. Independently integrate GM/r² with Simpson's rule.
        let count = 10_000
        var integral = 0.0
        for i in 0...count {
            let u = Double(i)/Double(count)
            let radius = 10+4*(3*u*u-2*u*u*u)
            let weight = i == 0 || i == count ? 1.0 : (i.isMultiple(of: 2) ? 2 : 4)
            integral += weight*1000/(radius*radius)
        }
        integral *= plan.duration.seconds/(3*Double(count))
        #expect(abs(result.state.mass.totalMassKilograms-100*exp(-integral/100)) < 1e-7)
        #expect(result.checkpoints[0].positionErrorMetres < 1e-9)
        #expect(result.checkpoints[0].velocityErrorMetresPerSecond < 1e-9)
    }
    private struct ConstantAcceleration: ForceModel {
        var massRate = 0.0
        func contribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
            try ForceContribution(epoch: state.epoch, frame: state.frame,
                accelerationMetresPerSecondSquared: SIMD3(2,0,0), propellantMassRateKilogramsPerSecond: massRate)
        }
    }
}
