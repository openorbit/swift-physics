// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Mattias Holm
import Foundation
import Math

public enum MassGeometryError: Error { case invalidComponent, duplicateIdentifier, missingComponent, cannotDetachRoot, emptyAssembly }

/// Rigid, instantaneous mass hierarchy. Translation and rotation map each local
/// structural frame into its parent's frame. Mass and intrinsic inertia belong
/// to the node itself; children are additional mass, never included twice.
/// Empty structural nodes have zero mass and nil inertia. No slosh is implied.
public struct MassComponent {
    public let id: String
    public let massKilograms: Double
    public let inertiaAtCentre: RigidBodyInertia?
    public let centreLocalMetres: SIMD3<Double>
    public let translationParentMetres: SIMD3<Double>
    public let localToParent: Quaternion
    public private(set) var children: [MassComponent]

    public init(id: String, massKilograms: Double = 0, inertiaAtCentre: RigidBodyInertia? = nil,
                centreLocalMetres: SIMD3<Double> = .zero, translationParentMetres: SIMD3<Double> = .zero,
                localToParent: Quaternion = Quaternion(x:0,y:0,z:0,w:1), children: [MassComponent] = []) throws {
        let q=[localToParent.x,localToParent.y,localToParent.z,localToParent.w]
        guard !id.isEmpty, massKilograms.isFinite, massKilograms >= 0,
              (massKilograms == 0) == (inertiaAtCentre == nil),
              (0..<3).allSatisfy({ centreLocalMetres[$0].isFinite && translationParentMetres[$0].isFinite }),
              q.allSatisfy(\.isFinite), abs(q.reduce(0) { $0+$1*$1 }-1)<1e-10 else { throw MassGeometryError.invalidComponent }
        self.id=id; self.massKilograms=massKilograms; self.inertiaAtCentre=inertiaAtCentre
        self.centreLocalMetres=centreLocalMetres; self.translationParentMetres=translationParentMetres
        self.localToParent=localToParent; self.children=children
        var ids=Set<String>()
        func validate(_ node: MassComponent) throws {
            guard ids.insert(node.id).inserted else { throw MassGeometryError.duplicateIdentifier }
            for child in node.children { try validate(child) }
        }
        try validate(self)
    }

    public static func cylinderInertia(mass: Double, length: Double, radius: Double) throws -> RigidBodyInertia {
        guard [mass,length,radius].allSatisfy({ $0.isFinite && $0>0 }) else { throw MassGeometryError.invalidComponent }
        return try RigidBodyInertia(xx:mass*radius*radius/2,
            yy:mass*(3*radius*radius+length*length)/12,zz:mass*(3*radius*radius+length*length)/12)
    }

    /// Nonmutating partition. Both results retain the original parent's axes and
    /// origin, so their centres can be mapped directly into the same rigid state.
    /// Failed lookup leaves the original tree untouched. Separation impulses are
    /// a dynamics concern and must be applied as equal/opposite impulses there.
    public func detaching(id target: String) throws -> (retained: MassComponent, detached: MassComponent) {
        guard target != id else { throw MassGeometryError.cannotDetachRoot }
        func remove(_ node: MassComponent) throws -> (MassComponent, MassComponent?) {
            var kept=node
            for i in node.children.indices {
                let child=node.children[i]
                if child.id==target {
                    kept.children.remove(at:i)
                    return (kept,child)
                }
                let (replacement,found)=try remove(child)
                if let found {
                    kept.children[i]=replacement
                    // Carry the detached subtree through this intervening frame.
                    let transformed=try MassComponent(id:found.id,massKilograms:found.massKilograms,
                        inertiaAtCentre:found.inertiaAtCentre,centreLocalMetres:found.centreLocalMetres,
                        translationParentMetres:child.translationParentMetres+child.localToParent.matrix*found.translationParentMetres,
                        localToParent:child.localToParent*found.localToParent,children:found.children)
                    return (kept,transformed)
                }
            }
            return (kept,nil)
        }
        let (retained,found)=try remove(self)
        guard let found else { throw MassGeometryError.missingComponent }
        // Root's placement also applies to the detached result.
        let detached=try MassComponent(id:found.id,massKilograms:found.massKilograms,
            inertiaAtCentre:found.inertiaAtCentre,centreLocalMetres:found.centreLocalMetres,
            translationParentMetres:translationParentMetres+localToParent.matrix*found.translationParentMetres,
            localToParent:localToParent*found.localToParent,children:found.children)
        return (retained,detached)
    }

    public func massGeometry() throws -> CompositeMassGeometry {
        var elements:[(Double,SIMD3<Double>,Quaternion,RigidBodyInertia)]=[]
        func visit(_ node: MassComponent,_ origin: SIMD3<Double>,_ rotation: Quaternion) {
            let p=origin+rotation.matrix*node.translationParentMetres, q=rotation*node.localToParent
            if let inertia=node.inertiaAtCentre { elements.append((node.massKilograms,p+q.matrix*node.centreLocalMetres,q,inertia)) }
            for child in node.children { visit(child,p,q) }
        }
        visit(self,.zero,Quaternion(x:0,y:0,z:0,w:1))
        let mass=elements.reduce(0) { $0+$1.0 }
        guard mass>0, mass.isFinite else { throw MassGeometryError.emptyAssembly }
        let centre=elements.reduce(SIMD3<Double>.zero) { $0+$1.0*$1.1 }/mass
        let axes=[SIMD3<Double>(1,0,0),SIMD3<Double>(0,1,0),SIMD3<Double>(0,0,1)]
        var tensor=Array(repeating:0.0,count:9)
        for (m,p,q,inertia) in elements {
            let r=p-centre, r2=(r*r).sum()
            for i in 0..<3 { for j in 0..<3 {
                let a=q.conj.matrix*axes[i], b=q.conj.matrix*axes[j]
                let rotated=(a*(inertia.matrix*b)).sum()
                let diagonal: Double = i==j ? r2 : 0
                let shifted=m*(diagonal-r[i]*r[j])
                tensor[i*3+j] += rotated+shifted
            } }
        }
        return CompositeMassGeometry(massKilograms:mass,centreMetres:centre,
            inertia:try RigidBodyInertia(xx:tensor[0],yy:tensor[4],zz:tensor[8],xy:tensor[1],xz:tensor[2],yz:tensor[5]))
    }
}

public struct CompositeMassGeometry {
    public let massKilograms: Double
    public let centreMetres: SIMD3<Double>
    public let inertia: RigidBodyInertia

    /// Instantaneous rigid point kinematics in the assembly's axes. Includes the
    /// rotational velocity at this COM, without inventing a separation impulse.
    public func offsetKinematics(from assembly: CompositeMassGeometry, angularVelocity: SIMD3<Double>) -> (position: SIMD3<Double>, velocity: SIMD3<Double>) {
        let r=centreMetres-assembly.centreMetres
        return (r,SIMD3(angularVelocity.y*r.z-angularVelocity.z*r.y,
            angularVelocity.z*r.x-angularVelocity.x*r.z,angularVelocity.x*r.y-angularVelocity.y*r.x))
    }
}
