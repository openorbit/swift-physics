import Foundation
import Testing
import CelMek
import Mechanics
import Math

private func mechanicalBody(_ id:String,_ p:SIMD3<Double>,_ v:SIMD3<Double> = .zero) throws -> MechanicalBody {
    try MechanicalBody(id:id,massKilograms:10,inertia:RigidBodyInertia(xx:4,yy:4,zz:4),positionMetres:p,velocityMetresPerSecond:v)
}
@Test func landingTetherSlackPayoutAndBreakAreExplicit() throws {
    let a=try mechanicalBody("a",.zero), b=try mechanicalBody("b",SIMD3(3,0,0),SIMD3(1,0,0))
    var cable=try MechanicalTether(id:"c",firstBody:"a",secondBody:"b",lengthMetres:2,stiffnessNewtonsPerMetre:100,dampingNewtonSecondsPerMetre:20,breakingTensionNewtons:110)
    #expect(try cable.load(first:a,second:b).tensionNewtons == 120)
    #expect(try cable.load(first:a,second:b,payoutMetresPerSecond:1).tensionNewtons == 100)
    var system=try MechanicalSystem(frame:CelestialFrames.icrf,epoch:.j2000,bodies:[a,b],tethers:[cable])
    #expect(try system.step(seconds:0.001){_,states in states.map{_ in MechanicalWrench()} } == ["c"])
    #expect(system.tethers[0].connected == false)
    #expect(system.bodies[0].velocityMetresPerSecond == .zero)
    try cable.setLength(metres:4)
    #expect(try cable.load(first:a,second:b).tensionNewtons == 0)
    cable.cut();try cable.setLength(metres:1)
    #expect(try cable.load(first:a,second:b).tensionNewtons == 0)
}
@Test func landingCableGraphConservesLinearMomentumAndCutHasNoImpulse() throws {
    let a=try mechanicalBody("a",.zero), b=try mechanicalBody("b",SIMD3(0,0,3))
    let cables=try [-1.0,0,1].map{ x in try MechanicalTether(id:"c\(x)",firstBody:"a",secondBody:"b",firstMountMetres:SIMD3(x,0,0),secondMountMetres:SIMD3(x,0,0),lengthMetres:2,stiffnessNewtonsPerMetre:100,dampingNewtonSecondsPerMetre:10) }
    var system=try MechanicalSystem(frame:CelestialFrames.icrf,epoch:.j2000,bodies:[a,b],tethers:cables)
    for _ in 0..<1000 { try system.step(seconds:0.001){_,states in states.map{_ in MechanicalWrench()} } }
    let momentum=(system.bodies[0].velocityMetresPerSecond+system.bodies[1].velocityMetresPerSecond)*10
    #expect((momentum*momentum).sum()<1e-20)
    let before=system.bodies
    try system.cut(tether:cables[0].id)
    #expect(system.bodies[0].positionMetres == before[0].positionMetres)
    #expect(system.bodies[1].velocityMetresPerSecond == before[1].velocityMetresPerSecond)
    let epoch=system.epoch
    #expect(throws:LandingMechanicsError.self) { try system.step(seconds:0.001,payout:["missing":1]){_,_ in []} }
    #expect(system.epoch == epoch)
}
@Test func landingContactSupportsBouncesAndDissipatesEnergy() throws {
    let ground=try MechanicalGround(pointMetres:.zero,normal:SIMD3(0,0,1))
    let contact=try MechanicalContact(body:"lander",mountMetres:.zero,radiusMetres:1,stiffnessNewtonsPerMetre:10000,dampingNewtonSecondsPerMetre:150,frictionCoefficient:0.7)
    var system=try MechanicalSystem(frame:CelestialFrames.icrf,epoch:.j2000,bodies:[mechanicalBody("lander",SIMD3(0,0,4),SIMD3(2,0,0))])
    var rebounded=false,previousVelocity=0.0,lowest=4.0
    for _ in 0..<20000 {
        try system.step(seconds:0.001){_,states in
            var load=try contact.load(on:states[0],ground:ground)
            load.forceReferenceNewtons += SIMD3(0,0,-37.1)
            return [load]
        }
        let b=system.bodies[0];lowest=min(lowest,b.positionMetres.z)
        if previousVelocity < 0 && b.velocityMetresPerSecond.z > 0 { rebounded=true }
        previousVelocity=b.velocityMetresPerSecond.z
    }
    let b=system.bodies[0]
    #expect(rebounded);#expect(lowest > 0.8)
    #expect(abs(b.positionMetres.z-(1-37.1/10000))<0.001)
    #expect(abs(b.velocityMetresPerSecond.z)<0.01)
    #expect(abs(b.velocityMetresPerSecond.x)<2)
    // A single sphere rolls without slip; friction alone must not halt rolling.
    let local=b.bodyToReference.invMatrix*SIMD3<Double>(0,0,-1)
    #expect(abs(b.pointVelocity(local).x)<0.01)
    let airborne=try mechanicalBody("lander",SIMD3(0,0,2),SIMD3(0,0,-10))
    #expect(try contact.load(on:airborne,ground:ground).forceReferenceNewtons == .zero)
}
@Test func landingGasCompressionVentingAndMountTorque() throws {
    var bag=try GasCushion(volumeCubicMetres:1,areaSquareMetres:1,absolutePressurePascals:1000)
    let initial=try bag.forceNewtons(compressionMetres:0,ambientPressurePascals:600)
    let compressed=try bag.forceNewtons(compressionMetres:0.5,ambientPressurePascals:600)
    #expect(compressed > initial)
    try bag.vent(retaining:0)
    #expect(try bag.forceNewtons(compressionMetres:0.5,ambientPressurePascals:600) == 0)
    #expect(throws:LandingMechanicsError.self) { try bag.vent(retaining:1.1) }
    let foot=try MechanicalContact(body:"b",mountMetres:SIMD3(1,0,0),radiusMetres:1,stiffnessNewtonsPerMetre:100,dampingNewtonSecondsPerMetre:0)
    let load=try foot.load(on:mechanicalBody("b",SIMD3(0,0,0.5)),ground:MechanicalGround(pointMetres:.zero,normal:SIMD3(0,0,1)))
    #expect(load.forceReferenceNewtons.z == 50)
    #expect(load.torqueBodyNewtonMetres.y == -50)
}
