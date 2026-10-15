// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Math
@testable import CelMek

@Suite struct FiniteTransferPlannerTests {
    let frame = Frame(name: "Finite transfer", pseudoInertial: true)
    func epoch(_ seconds: Double) throws -> Epoch { try Epoch.j2000.advanced(by: CelMek.Duration(seconds: seconds)) }
    func seed() throws -> SpacecraftState {
        try SpacecraftState(cartesian: CartesianState(epoch: .j2000, frame: frame,
            position: SIMD3(1,0,0), velocity: SIMD3(0,1,0)),
            mass: SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50))
    }
    @Test func rocketEquationAndEndpointPlacement() throws {
        let maneuvers = try [ImpulsiveManeuver(epoch: .j2000, frame: frame, deltaVelocityMetresPerSecond: SIMD3(1,0,0)),
            ImpulsiveManeuver(epoch: epoch(10), frame: frame, deltaVelocityMetresPerSecond: SIMD3(0,2,0))]
        let schedule = try FiniteTransferPlanner.initialGuess(initialState: seed(), arrivalEpoch: epoch(10),
            maneuvers: maneuvers, thrustNewtons: 100, effectiveExhaustVelocityMetresPerSecond: 100)
        #expect(schedule.burns.first?.start == .j2000)
        let end = try epoch(10)
        #expect(schedule.burns.last?.end == end)
        let used = try schedule.burns.reduce(0.0) { try $0 + $1.end.duration(since: $1.start).seconds }
        #expect(abs(used - 100 * -expm1(-0.03)) < 1e-12)
        #expect(throws: FiniteTransferPlanningError.insufficientPropellant) {
            try FiniteTransferPlanner.initialGuess(initialState: seed(), arrivalEpoch: epoch(10),
                maneuvers: maneuvers, thrustNewtons: 100, effectiveExhaustVelocityMetresPerSecond: 1)
        }
        #expect(throws: FiniteTransferPlanningError.burnsDoNotFit) {
            try FiniteTransferPlanner.initialGuess(initialState: seed(), arrivalEpoch: epoch(10),
                maneuvers: maneuvers, thrustNewtons: 0.01, effectiveExhaustVelocityMetresPerSecond: 100)
        }
    }
    @Test func shootingImprovesRealEndpointAndReportsBudget() throws {
        let context = try ForceEvaluationContext(frame: frame, origin: .sun)
        let engine = try FiniteThrust(frame: frame, thrustNewtons: 10,
            effectiveExhaustVelocityMetresPerSecond: 100, direction: SIMD3(1,0,0))
        let exact = try ScheduledThrust(frame: frame, burns: [FiniteBurn(start: epoch(2), end: epoch(4), thrust: engine)])
        let target = try ScheduledThrustPropagator(initialState: seed(), schedule: exact, context: context).state(at: epoch(10))
        let perturbed = try ScheduledThrust(frame: frame, burns: [FiniteBurn(start: epoch(2), end: epoch(3.6), thrust: engine)])
        let result = try FiniteTransferPlanner.refine(initialState: seed(), initialSchedule: perturbed,
            target: target, context: context, positionToleranceMetres: 1e-6, velocityToleranceMetresPerSecond: 1e-6)
        #expect(result.status == .converged)
        #expect(result.best.positionErrorMetres < result.initial.positionErrorMetres)
        #expect(abs(result.best.propellantUsedKilograms - 0.2) < 1e-6)
        let limited = try FiniteTransferPlanner.refine(initialState: seed(), initialSchedule: perturbed,
            target: target, context: context, maximumEvaluations: 1)
        #expect(limited.status == .evaluationLimitReached)
        #expect(limited.evaluations == 1)
        #expect(throws: FiniteTransferPlanningError.insufficientPropellant) {
            try FiniteTransferPlanner.refine(initialState: seed(), initialSchedule: perturbed,
                target: target, context: context, minimumRemainingPropellantKilograms: 50)
        }
    }
    @Test func shorterBurnsApproachImpulsiveHohmannUnderGravity() throws {
        let gm = try GravitationalParameter(cubicMetresPerSecondSquared: 1, convention: .newtonian, provenance: PhysicalDataProvenance(source: "Test", version: "1"))
        let transfer = try HohmannTransfer(departureRadiusMetres: 1, arrivalRadiusMetres: 2, gravitationalParameter: gm)
        let plan = try transfer.plan(departureEpoch: .j2000, frame: frame, centralBody: .sun,
            radialDirection: SIMD3(1,0,0), progradeDirection: SIMD3(0,1,0))
        let arrival = plan.arrivalEpoch
        let target = try plan.trajectory.state(at: arrival)
        let context = try ForceEvaluationContext(frame: frame, origin: .sun)
        let gravity = try CentralGravity(body: .sun, gravitationalParameter: gm)
        var errors: [Double] = []
        for thrust in [1_000.0, 10_000.0] {
            let schedule = try FiniteTransferPlanner.initialGuess(initialState: seed(), arrivalEpoch: arrival,
                maneuvers: [plan.departureBurn, plan.arrivalBurn], thrustNewtons: thrust,
                effectiveExhaustVelocityMetresPerSecond: 100)
            let actual = try ScheduledThrustPropagator(initialState: seed(), schedule: schedule,
                backgroundModels: [gravity], context: context, options: .preset(for: .full)).state(at: arrival)
            errors.append(routeNorm(actual.position - target.position) + routeNorm(actual.velocity - target.velocity))
            if thrust == 1_000 {
                let refined = try FiniteTransferPlanner.refine(initialState: seed(), initialSchedule: schedule,
                    target: target, context: context, backgroundModels: [gravity],
                    positionToleranceMetres: 1e-5, velocityToleranceMetresPerSecond: 1e-5,
                    maximumEvaluations: 100, numericalOptions: .preset(for: .full))
                #expect(hypot(refined.best.positionErrorMetres, refined.best.velocityErrorMetresPerSecond) <
                    hypot(refined.initial.positionErrorMetres, refined.initial.velocityErrorMetresPerSecond))
                #expect(refined.evaluations <= 100)
            }
        }
        #expect(errors[1] < errors[0] * 0.2)
        #expect(errors[1] > 0)
    }
}
