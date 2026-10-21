import Foundation
import Testing
import CelMek
import Math

@Test func nestedMassHierarchyDetachmentPreservesMassAndMomentum() throws {
    let inertia=try RigidBodyInertia(xx:2,yy:2,zz:2)
    let rotation=Quaternion(x:0,y:0,z:sin(.pi/4),w:cos(.pi/4))
    let a=try MassComponent(id:"a",massKilograms:10,inertiaAtCentre:inertia)
    let b=try MassComponent(id:"b",massKilograms:20,inertiaAtCentre:inertia,translationParentMetres:SIMD3(2,0,0))
    let group=try MassComponent(id:"group",translationParentMetres:SIMD3(0,3,0),localToParent:rotation,children:[b])
    let root=try MassComponent(id:"root",children:[a,group])
    let whole=try root.massGeometry(), split=try root.detaching(id:"b")
    let kept=try split.retained.massGeometry(), detached=try split.detached.massGeometry()
    #expect(abs(whole.massKilograms-kept.massKilograms-detached.massKilograms)<1e-12)
    #expect(abs(detached.centreMetres.y-5)<1e-12)
    #expect(abs(whole.centreMetres.y-100/30)<1e-12)
    let w=SIMD3<Double>(0,0,2)
    let k=kept.offsetKinematics(from:whole,angularVelocity:w)
    let d=detached.offsetKinematics(from:whole,angularVelocity:w)
    let momentum=k.velocity*kept.massKilograms+d.velocity*detached.massKilograms
    #expect((momentum*momentum).sum()<1e-18)
    #expect(abs((whole.inertia.matrix*SIMD3<Double>(1,0,0)).x-(4+10*pow(whole.centreMetres.y,2)+20*pow(5-whole.centreMetres.y,2)))<1e-9)
    #expect(throws:MassGeometryError.self) { try root.detaching(id:"missing") }
    #expect(throws:MassGeometryError.self) { try MassComponent(id:"duplicate",children:[a,a]) }
}

@Test func rotatedComponentInertiaAndEmptyStructure() throws {
    let q=Quaternion(x:0,y:0,z:sin(.pi/4),w:cos(.pi/4))
    let node=try MassComponent(id:"rotated",massKilograms:1,
        inertiaAtCentre:RigidBodyInertia(xx:2,yy:3,zz:4),localToParent:q)
    let result=try node.massGeometry()
    #expect(abs((result.inertia.matrix*SIMD3<Double>(1,0,0)).x-3)<1e-12)
    #expect(abs((result.inertia.matrix*SIMD3<Double>(0,1,0)).y-2)<1e-12)
    #expect(throws:MassGeometryError.self) { try MassComponent(id:"empty").massGeometry() }
}
