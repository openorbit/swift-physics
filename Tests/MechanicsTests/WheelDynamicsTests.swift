import Foundation
import Testing
import CelMek
import Mechanics
import Math

private func wheelBody(_ z:Double,_ v:SIMD3<Double> = .zero) throws -> MechanicalBody {
    try MechanicalBody(id:"rover",massKilograms:100,inertia:RigidBodyInertia(xx:100,yy:100,zz:100),positionMetres:SIMD3(0,0,z),velocityMetresPerSecond:v)
}
@Test func wheelGroundTractionCircleAirborneAndFrictionlessSpin() throws {
    let ground=try MechanicalGround(pointMetres:.zero,normal:SIMD3(0,0,1)),body=try wheelBody(0.48)
    let wheel=try SuspendedWheel(id:"w",mountBodyMetres:.zero)
    let command=try WheelCommand(motorTorqueNewtonMetres:3),spin=try WheelState(spinRadiansPerSecond:10)
    let hit=try wheel.evaluate(body:body,state:spin,ground:ground,command:command,stepSeconds:0.001)
    #expect(hit.normalLoad>0);#expect(hit.wrench.forceReferenceNewtons.x>0)
    #expect(hypot(hit.wrench.forceReferenceNewtons.x,hit.wrench.forceReferenceNewtons.y)<=0.8*hit.normalLoad+1e-10)
    let air=try wheel.evaluate(body:wheelBody(3),state:spin,ground:ground,command:command,stepSeconds:0.001)
    #expect(air.normalLoad==0);#expect(air.wrench.forceReferenceNewtons == .zero)
    #expect(air.spinAcceleration>0);#expect(air.wrench.torqueBodyNewtonMetres.y<0)
    let ice=try SuspendedWheel(id:"ice",mountBodyMetres:.zero,friction:0)
    let iceSample=try ice.evaluate(body:body,state:spin,ground:ground,command:command,stepSeconds:0.001)
    #expect(iceSample.wrench.forceReferenceNewtons.x==0)
    #expect(throws:LandingMechanicsError.self) {try WheelCommand(brakeFraction:1.1)}
}
@Test func wheelAirborneMotorConservesTotalAngularMomentumAndStepRollback() throws {
    let wheel=try SuspendedWheel(id:"w",mountBodyMetres:.zero)
    var set=try SuspendedWheelSet(bodyID:"rover",wheels:[wheel]),system=try MechanicalSystem(frame:CelestialFrames.icrf,epoch:.j2000,bodies:[wheelBody(3)])
    let command=try WheelCommand(motorTorqueNewtonMetres:3),ground=try MechanicalGround(pointMetres:.zero,normal:SIMD3(0,0,1))
    for _ in 0..<1000 {
        try set.step(system:&system,seconds:0.001,commands:[command],ground:{_,_,_ in ground},externalLoads:{_,_ in [MechanicalWrench()]})
    }
    let omega=system.bodies[0].angularVelocityBody.y,spin=set.states[0].spinRadiansPerSecond
    #expect(abs(100*omega+0.25*spin)<1e-9)
    #expect(spin>10);#expect(omega<0)
    let epoch=system.epoch,previousSpin=spin
    #expect(throws:LandingMechanicsError.self) {
        try set.step(system:&system,seconds:0.001,commands:[command],ground:{_,_,_ in throw LandingMechanicsError.invalidParameter},externalLoads:{_,_ in []})
    }
    #expect(system.epoch==epoch);#expect(set.states[0].spinRadiansPerSecond==previousSpin)
}
