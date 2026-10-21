// SPDX-License-Identifier: Apache-2.0
import Testing
import Math
@testable import CelMek

@Suite struct FiniteBurnPositionTargetingTests {
    let frame = Frame(name:"Position corrector test",pseudoInertial:true)
    func epoch(_ seconds: Double) throws -> Epoch { try Epoch.j2000.advanced(by:CelMek.Duration(seconds:seconds)) }
    func state() throws -> SpacecraftState {
        try SpacecraftState(cartesian:CartesianState(epoch:.j2000,frame:frame,position:SIMD3(10,0,0),velocity:SIMD3(0,1,0)),
            mass:SpacecraftMassProperties(dryMassKilograms:50,propellantMassKilograms:50))
    }
    func burn(_ seconds: Double, _ direction: SIMD3<Double>) throws -> FiniteBurn {
        try FiniteBurn(start:epoch(1),end:epoch(1+seconds),thrust:FiniteThrust(frame:frame,thrustNewtons:10,
            effectiveExhaustVelocityMetresPerSecond:100,direction:direction))
    }
    @Test func correctsThreeDimensionalPositionWithFixedIgnitionAndFiniteFuel() throws {
        let initial = try state(), context = try ForceEvaluationContext(frame:frame,origin:.earth)
        let exact = try ScheduledThrust(frame:frame,burns:[burn(2,SIMD3(1,0.3,0.2))])
        let target = try ScheduledThrustPropagator(initialState:initial,schedule:exact,context:context).state(at:epoch(10))
        let result = try FiniteBurnPositionTargeting.refine(initialState:initial,initialBurn:burn(1.6,SIMD3(1,0.1,0.3)),
            target:target,context:context,positionToleranceMetres:1e-7,maximumEvaluations:100,burnWindow:epoch(1)...epoch(6))
        #expect(result.status == .positionConverged)
        #expect(result.best.positionErrorMetres < 1e-7)
        #expect(result.best.positionErrorMetres < result.initial.positionErrorMetres)
        #expect(result.best.schedule.burns[0].start == (try epoch(1)))
        #expect(result.evaluations <= 100)
        #expect(abs(result.best.propellantUsedKilograms-0.2) < 1e-5)
    }
    @Test func respectsEvaluationReserveAndTimingLimits() throws {
        let initial = try state(), context = try ForceEvaluationContext(frame:frame,origin:.earth)
        let target = try CartesianState(epoch:epoch(10),frame:frame,position:SIMD3(12,10,0),velocity:.zero)
        let limited = try FiniteBurnPositionTargeting.refine(initialState:initial,initialBurn:burn(1,SIMD3(1,0,0)),target:target,
            context:context,positionToleranceMetres:1e-12,maximumEvaluations:1)
        #expect(limited.status == .evaluationLimitReached && limited.evaluations == 1)
        #expect(throws:FiniteTransferPlanningError.invalidInput) {
            try FiniteBurnPositionTargeting.refine(initialState:initial,initialBurn:burn(1,SIMD3(1,0,0)),target:target,
                context:context,minimumRemainingPropellantKilograms:50)
        }
        #expect(throws:FiniteTransferPlanningError.invalidInput) {
            try FiniteBurnPositionTargeting.refine(initialState:initial,initialBurn:burn(1,SIMD3(1,0,0)),target:target,
                context:context,burnWindow:epoch(3)...epoch(5))
        }
    }
}

@Suite struct FiniteBurnStateTargetingTests {
    @Test func correctsTwoBurnsAgainstPositionAndVelocityAndRejectsInvalidInputs() throws {
        let frame = Frame(name:"State target",pseudoInertial:true)
        func epoch(_ s: Double) throws -> Epoch { try Epoch.j2000.advanced(by:CelMek.Duration(seconds:s)) }
        let initial = try SpacecraftState(cartesian:CartesianState(epoch:.j2000,frame:frame,position:.zero,velocity:SIMD3(0,1,0)),mass:SpacecraftMassProperties(dryMassKilograms:50,propellantMassKilograms:50))
        let context = try ForceEvaluationContext(frame:frame,origin:.earth)
        func schedule(_ duration: Double, _ bias: Double) throws -> ScheduledThrust {
            try ScheduledThrust(frame:frame,burns:[1.0,8.0].map { start in
                try FiniteBurn(start:epoch(start),end:epoch(start+duration),thrust:FiniteThrust(frame:frame,thrustNewtons:10,effectiveExhaustVelocityMetresPerSecond:100,direction:SIMD3(1,bias,start/10)))
            })
        }
        let exact = try schedule(2,0.3), seed = try schedule(1.8,0.2)
        let target = try ScheduledThrustPropagator(initialState:initial,schedule:exact,context:context).state(at:epoch(15))
        let result = try FiniteBurnStateTargeting.refine(initialState:initial,initialSchedule:seed,target:target,context:context,backgroundModels:[],positionToleranceMetres:1e-6,velocityToleranceMetresPerSecond:1e-7)
        #expect(result.positionErrorMetres < 1e-6 && result.velocityErrorMetresPerSecond < 1e-7)
        #expect(abs(result.propellantUsedKilograms-0.4) < 1e-5)
        #expect(throws:FiniteTransferPlanningError.invalidInput) {
            try FiniteBurnStateTargeting.refine(initialState:initial,initialSchedule:seed,target:target,context:context,backgroundModels:[],minimumRemainingPropellantKilograms:50)
        }
        #expect(throws:FiniteTransferPlanningError.invalidInput) {
            try FiniteBurnStateTargeting.refine(initialState:initial,initialSchedule:seed,target:target,context:context,backgroundModels:[],positionToleranceMetres:.nan)
        }
    }
}
