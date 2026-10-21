import Testing
import CelMek
import Math

@Test func finiteBurnStateTargetingAcceptsFineControlStepsAndRejectsInvalidOnes() throws {
    let frame=CelestialFrames.icrf,epoch=Epoch.j2000
    let initial=try SpacecraftState(cartesian:CartesianState(epoch:epoch,frame:frame,position:SIMD3(1e8,0,0),velocity:.zero),
        mass:SpacecraftMassProperties(dryMassKilograms:1000,propellantMassKilograms:1000))
    let context=try ForceEvaluationContext(frame:frame,origin:.solarSystemBarycentre)
    func burn(_ start: Double,_ duration: Double,_ direction: SIMD3<Double>) throws -> FiniteBurn {
        try FiniteBurn(start:epoch.advanced(by:CelMek.Duration(seconds:start)),end:epoch.advanced(by:CelMek.Duration(seconds:start+duration)),
            thrust:FiniteThrust(frame:frame,thrustNewtons:100,effectiveExhaustVelocityMetresPerSecond:10000,direction:direction))
    }
    let reference=try ScheduledThrust(frame:frame,burns:[burn(10,100,SIMD3(1,0,0)),burn(700,100,SIMD3(0,1,0))])
    let arrival=try epoch.advanced(by:CelMek.Duration(seconds:1000))
    let target=try ScheduledThrustPropagator(initialState:initial,schedule:reference,context:context).state(at:arrival)
    let seed=try ScheduledThrust(frame:frame,burns:[burn(10,105,SIMD3(1,0,0)),burn(700,95,SIMD3(0,1,0))])
    let result=try FiniteBurnStateTargeting.refine(initialState:initial,initialSchedule:seed,target:target,context:context,backgroundModels:[],
        positionToleranceMetres:0.01,velocityToleranceMetresPerSecond:1e-5,relativeControlStep:1e-6)
    #expect(result.positionErrorMetres<0.01 && result.velocityErrorMetresPerSecond<1e-5)
    for value in [0.0,-0.1,Double.nan,Double.infinity,0.1] {
        #expect(throws:FiniteTransferPlanningError.invalidInput) {
            try FiniteBurnStateTargeting.refine(initialState:initial,initialSchedule:seed,target:target,context:context,backgroundModels:[],relativeControlStep:value)
        }
    }
}
