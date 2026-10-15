// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct CommandedBurnPlannerTests {
    let frame = Frame(name: "Commanded burns", pseudoInertial: true)
    func epoch(_ t: Double) throws -> Epoch { try Epoch.j2000.advanced(by: CelMek.Duration(seconds: t)) }
    struct LinearGoal: StateProvider {
        let frame: Frame
        let p, v: SIMD3<Double>
        func state(at epoch: Epoch) throws -> CartesianState {
            try CartesianState(epoch: epoch, frame: frame,
                position: p + v*epoch.duration(since: .j2000).seconds, velocity: v, acceleration: .zero)
        }
    }
    struct Gravity: ForceModel {
        func contribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
            try ForceContribution(epoch: state.epoch, frame: context.frame,
                accelerationMetresPerSecondSquared: SIMD3(0.1,0,0))
        }
    }
    func goal(position: Double, velocity: Double, at t: Double) throws -> MovingWaypointRouteTarget {
        try MovingWaypointRouteTarget(provider: LinearGoal(frame: frame,
            p: SIMD3(position-velocity*t,0,0), v: SIMD3(velocity,0,0)),
            origin: .sun, coverage: epoch(0)...epoch(20))
    }
    func candidates() throws -> [CommandedBurnCandidate] {
        let forward = try FiniteThrust(frame: frame, thrustNewtons: 10,
            effectiveExhaustVelocityMetresPerSecond: 100, direction: SIMD3(1,0,0))
        let reverse = try FiniteThrust(frame: frame, thrustNewtons: 10,
            effectiveExhaustVelocityMetresPerSecond: 100, direction: SIMD3(-1,0,0))
        let schedule = try ScheduledThrust(frame: frame, burns: [
            FiniteBurn(start: .j2000, end: epoch(5), thrust: forward),
            FiniteBurn(start: epoch(5), end: epoch(10), thrust: reverse)])
        return try [CommandedBurnCandidate(arrivalEpoch: epoch(10), schedule: schedule),
            CommandedBurnCandidate(arrivalEpoch: epoch(8), schedule: schedule),
            CommandedBurnCandidate(arrivalEpoch: epoch(9), schedule: ScheduledThrust(frame: frame, burns: []))]
    }
    func reference(_ t: Double) -> (p: Double, v: Double) {
        let firstTime = min(t,5), m1 = 100-0.1*firstTime
        let v1 = 100*log(100/m1)
        let p1 = 100*(firstTime-m1/0.1*log(100/m1))
        if t <= 5 { return (p1+0.05*t*t, v1+0.1*t) }
        let dt = t-5, m2 = m1-0.1*dt
        return (p1+v1*dt-100*(dt-m2/0.1*log(m1/m2))+0.05*t*t,
            v1-100*log(m1/m2)+0.1*t)
    }
    func search(budget: Int = 1000, reserve: Double = 0) throws -> CommandedBurnSearchResult {
        let arrival = reference(10), passage = reference(2.5)
        return try CommandedBurnPlanner.search(initialState:
            SpacecraftState(cartesian: CartesianState(epoch: .j2000, frame: frame, position: .zero, velocity: .zero),
                mass: SpacecraftMassProperties(dryMassKilograms: 50, propellantMassKilograms: 50)),
            context: ForceEvaluationContext(frame: frame, origin: .sun), backgroundModels: [Gravity()],
            candidates: candidates(), arrivalTarget: goal(position: arrival.p, velocity: arrival.v, at: 10),
            passages: [CommandedBurnPassage(epoch: epoch(2.5),
                target: goal(position: passage.p, velocity: passage.v, at: 2.5))],
            maximumCandidateEvaluations: budget, minimumRemainingPropellantKilograms: reserve,
            positionToleranceMetres: 1e-7, velocityToleranceMetresPerSecond: 1e-8,
            numericalOptions: .preset(for: .full))
    }
    @Test func commandedDirectionsCarryFuelThroughIndependentPassage() throws {
        let found = try search()
        #expect(found.status == .feasible); #expect(found.candidateIndex == 0)
        #expect(found.evaluatedCandidates == 3)
        #expect(found.failedAttempts[0].failure == .scheduleOutsideFlight)
        #expect(found.failedAttempts[1].failure == .passageResidual(index: 0))
        #expect(found.checkpoints.count == 2)
        #expect(found.checkpoints[0].state.cartesian.velocity.x > 0)
        #expect(abs(found.checkpoints[0].state.mass.propellantMassKilograms-49.75) < 1e-10)
        #expect(abs(found.checkpoints[1].state.mass.propellantMassKilograms-49) < 1e-10)
        #expect(found.checkpoints.allSatisfy { $0.positionErrorMetres < 1e-7 && $0.velocityErrorMetresPerSecond < 1e-8 })
    }
    @Test func budgetAndFuelReserveAreExplicitFailures() throws {
        let limited = try search(budget: 1)
        #expect(limited.status == .evaluationLimitReached); #expect(limited.evaluatedCandidates == 1)
        #expect(limited.candidate == nil)
        let reserve = try search(reserve: 50)
        #expect(reserve.status == .candidateSetExhausted)
        #expect(reserve.failedAttempts.last?.failure == .propellantReserve)
    }
}
