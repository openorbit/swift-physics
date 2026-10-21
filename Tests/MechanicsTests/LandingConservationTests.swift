import Foundation
import Testing
import CelMek
import Mechanics
import Math

private func lcCross(_ a:SIMD3<Double>,_ b:SIMD3<Double>) -> SIMD3<Double> { SIMD3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x) }
@Test func landingOffCentreBridleConservesAngularMomentumAndConverges() throws {
    func run(_ dt:Double,_ offset:SIMD3<Double>) throws -> (MechanicalSystem,SIMD3<Double>) {
        let tensor=try RigidBodyInertia(xx:10,yy:10,zz:10)
        let a=try MechanicalBody(id:"a",massKilograms:10,inertia:tensor,positionMetres:offset)
        let b=try MechanicalBody(id:"b",massKilograms:20,inertia:tensor,positionMetres:offset+SIMD3(0,0,3))
        let cable=try MechanicalTether(id:"c",firstBody:"a",secondBody:"b",firstMountMetres:SIMD3(1,0,0),secondMountMetres:SIMD3(-1,0,0),lengthMetres:3,stiffnessNewtonsPerMetre:100,dampingNewtonSecondsPerMetre:5)
        var s=try MechanicalSystem(frame:CelestialFrames.icrf,epoch:.j2000,bodies:[a,b],tethers:[cable])
        for _ in 0..<Int((1/dt).rounded()) { try s.step(seconds:dt){_,states in states.map{_ in MechanicalWrench()} } }
        var angular=SIMD3<Double>.zero
        for body in s.bodies {
            angular += lcCross(body.positionMetres-offset,body.velocityMetresPerSecond*body.massKilograms)
                + body.bodyToReference.matrix*(try body.inertia.angularMomentumBody(for:body.angularVelocityBody))
        }
        return (s,angular)
    }
    let coarse=try run(0.002,.zero),fine=try run(0.001,.zero),shift=try run(0.001,SIMD3(1e8,-1e8,1e8))
    #expect(sqrt((fine.1*fine.1).sum())<0.001)
    #expect(abs(fine.0.bodies[0].angularVelocityBody.y)>0.01)
    let delta=coarse.0.bodies[0].positionMetres-fine.0.bodies[0].positionMetres
    #expect(sqrt((delta*delta).sum())<0.001)
    let shifted=shift.0.bodies[0].positionMetres-SIMD3<Double>(1e8,-1e8,1e8)-fine.0.bodies[0].positionMetres
    #expect(sqrt((shifted*shifted).sum())<0.001)
}
