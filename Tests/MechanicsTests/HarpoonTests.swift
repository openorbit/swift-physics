import Foundation
import Testing
import CelMek
import Mechanics
import Math

private func cometNorm(_ v:SIMD3<Double>)->Double {sqrt((v*v).sum())}
private func cometCross(_ a:SIMD3<Double>,_ b:SIMD3<Double>)->SIMD3<Double> {SIMD3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x)}

@Test func harpoonInternalLaunchConservesMassAndMomentumIncludingRecoilTorque() throws {
    let parent=try MechanicalBody(id:"probe",massKilograms:100,inertia:RigidBodyInertia(xx:90,yy:90,zz:80),positionMetres:SIMD3(1000,20,30),
        velocityMetresPerSecond:SIMD3(1,2,3),angularVelocityBody:SIMD3(0.1,0.2,0.3))
    let launch=try parent.launching(id:"dart",massKilograms:0.05,inertia:RigidBodyInertia(xx:0.0001,yy:0.0001,zz:0.00005),centreBodyMetres:SIMD3(0.3,0,-0.2),impulseBodyNewtonSeconds:SIMD3(0,0,-0.75))
    #expect(abs(launch.retained.massKilograms+launch.projectile.massKilograms-parent.massKilograms)<1e-12)
    let delta=launch.retained.velocityMetresPerSecond*launch.retained.massKilograms+launch.projectile.velocityMetresPerSecond*launch.projectile.massKilograms-parent.velocityMetresPerSecond*parent.massKilograms
    #expect(cometNorm(delta)<1e-10)
    var momentum=SIMD3<Double>.zero
    for b in [launch.retained,launch.projectile] {
        momentum += cometCross(b.positionMetres-parent.positionMetres,b.velocityMetresPerSecond*b.massKilograms)+b.bodyToReference.matrix*(try b.inertia.angularMomentumBody(for:b.angularVelocityBody))
    }
    #expect(cometNorm(momentum-(try parent.inertia.angularMomentumBody(for:parent.angularVelocityBody)))<1e-9)
    #expect(launch.retained.velocityMetresPerSecond.z>parent.pointVelocity(launch.retainedCentreBodyMetres).z)
    #expect(launch.retained.angularVelocityBody.y != parent.angularVelocityBody.y)
}
@Test func harpoonAttachmentDepthFramePulloutCutAndCableBreak() throws {
    let fixed=FrameTree(name:"prescribed surface").root
    let substrate=try HarpoonSubstrate(penetrationResistanceNewtons:25,minimumDepthMetres:0.04,maximumDepthMetres:0.8,holdingForceNewtons:20)
    var anchor=try HarpoonAnchor(id:"a",surfaceFrame:fixed)
    try anchor.fire()
    #expect(throws:LandingMechanicsError.self) {try anchor.impact(pointMetres:.zero,in:CelestialFrames.icrf,massKilograms:0.05,normalSpeedMetresPerSecond:15,incidenceCosine:1,substrate:substrate)}
    #expect(anchor.phase == .flying)
    #expect(try anchor.impact(pointMetres:SIMD3(1,2,3),in:fixed,massKilograms:0.05,normalSpeedMetresPerSecond:15,incidenceCosine:1,substrate:substrate))
    #expect(abs(anchor.penetrationMetres-0.225)<1e-12)
    var pullout=anchor,broken=anchor
    try pullout.commitLoad(tensionNewtons:21,breakingTensionNewtons:30);#expect(pullout.phase == .pulledOut)
    try broken.commitLoad(tensionNewtons:31,breakingTensionNewtons:30);#expect(broken.cableState == .broken);#expect(broken.phase == .anchored)
    anchor.cut();try anchor.commitLoad(tensionNewtons:100,breakingTensionNewtons:30)
    #expect(anchor.phase == .anchored);#expect(anchor.cableState == .cut)
    var rejected=try HarpoonAnchor(id:"r",surfaceFrame:fixed);try rejected.fire()
    #expect(try !rejected.impact(pointMetres:.zero,in:fixed,massKilograms:0.05,normalSpeedMetresPerSecond:0.1,incidenceCosine:1,substrate:substrate))
    #expect(rejected.phase == .rejected)
    #expect(throws:LandingMechanicsError.self) {try HarpoonSubstrate(penetrationResistanceNewtons:0,minimumDepthMetres:1,maximumDepthMetres:2,holdingForceNewtons:1)}
}
