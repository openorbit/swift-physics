// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Mattias Holm
import Foundation
import Math

public enum LandingMechanicsError: Error { case invalidParameter, invalidTopology, numericalRange }

/// Centre-of-mass state, in the ensemble's shared inertial reference frame.
/// Connections are a graph; they do not change the frame tree or combine inertia.
public struct MechanicalBody {
    public let id: String
    public let massKilograms: Double
    public let inertia: RigidBodyInertia
    public private(set) var positionMetres: SIMD3<Double>
    public private(set) var velocityMetresPerSecond: SIMD3<Double>
    public private(set) var bodyToReference: Quaternion
    public private(set) var angularVelocityBody: SIMD3<Double>
    public init(id: String, massKilograms: Double, inertia: RigidBodyInertia,
        positionMetres: SIMD3<Double>, velocityMetresPerSecond: SIMD3<Double> = .zero,
        bodyToReference: Quaternion = .identity, angularVelocityBody: SIMD3<Double> = .zero) throws {
        guard !id.isEmpty, massKilograms.isFinite, massKilograms > 0,
              attitudeFinite(positionMetres), attitudeFinite(velocityMetresPerSecond), attitudeFinite(angularVelocityBody) else {
            throw LandingMechanicsError.invalidParameter
        }
        self.id=id; self.massKilograms=massKilograms; self.inertia=inertia
        self.positionMetres=positionMetres; self.velocityMetresPerSecond=velocityMetresPerSecond
        self.bodyToReference=try mechanicalUnit(bodyToReference); self.angularVelocityBody=angularVelocityBody
    }
    public func point(_ local: SIMD3<Double>) -> SIMD3<Double> { positionMetres + bodyToReference.matrix*local }
    public func pointVelocity(_ local: SIMD3<Double>) -> SIMD3<Double> {
        velocityMetresPerSecond + bodyToReference.matrix*attitudeCross(angularVelocityBody,local)
    }
}

public struct MechanicalWrench {
    public var forceReferenceNewtons: SIMD3<Double> = .zero
    public var torqueBodyNewtonMetres: SIMD3<Double> = .zero
    public init(forceReferenceNewtons: SIMD3<Double> = .zero, torqueBodyNewtonMetres: SIMD3<Double> = .zero) {
        self.forceReferenceNewtons=forceReferenceNewtons; self.torqueBodyNewtonMetres=torqueBodyNewtonMetres
    }
    public mutating func apply(force: SIMD3<Double>, at local: SIMD3<Double>, to body: MechanicalBody) {
        forceReferenceNewtons += force
        torqueBodyNewtonMetres += attitudeCross(local,body.bodyToReference.invMatrix*force)
    }
}

/// Tension-only Kelvin–Voigt cable. Payout is positive, retrieval negative.
/// The caller advances commanded length and commits break/cut events between steps.
public struct MechanicalTether {
    public let id: String
    public let firstBody: String
    public let secondBody: String
    public let firstMountMetres: SIMD3<Double>
    public let secondMountMetres: SIMD3<Double>
    public let stiffnessNewtonsPerMetre: Double
    public let dampingNewtonSecondsPerMetre: Double
    public let breakingTensionNewtons: Double
    public private(set) var lengthMetres: Double
    public private(set) var connected: Bool = true
    public init(id: String, firstBody: String, secondBody: String,
        firstMountMetres: SIMD3<Double> = .zero, secondMountMetres: SIMD3<Double> = .zero,
        lengthMetres: Double, stiffnessNewtonsPerMetre: Double, dampingNewtonSecondsPerMetre: Double,
        breakingTensionNewtons: Double = .greatestFiniteMagnitude) throws {
        guard !id.isEmpty, !firstBody.isEmpty, !secondBody.isEmpty, firstBody != secondBody,
              attitudeFinite(firstMountMetres), attitudeFinite(secondMountMetres),
              [lengthMetres,stiffnessNewtonsPerMetre,breakingTensionNewtons].allSatisfy({$0.isFinite && $0 > 0}),
              dampingNewtonSecondsPerMetre.isFinite, dampingNewtonSecondsPerMetre >= 0 else { throw LandingMechanicsError.invalidParameter }
        self.id=id; self.firstBody=firstBody; self.secondBody=secondBody
        self.firstMountMetres=firstMountMetres; self.secondMountMetres=secondMountMetres
        self.lengthMetres=lengthMetres; self.stiffnessNewtonsPerMetre=stiffnessNewtonsPerMetre
        self.dampingNewtonSecondsPerMetre=dampingNewtonSecondsPerMetre; self.breakingTensionNewtons=breakingTensionNewtons
    }
    public mutating func cut() { connected=false }
    public mutating func setLength(metres: Double) throws {
        guard metres.isFinite, metres > 0 else { throw LandingMechanicsError.invalidParameter }
        lengthMetres=metres
    }
    public func load(first: MechanicalBody, second: MechanicalBody, payoutMetresPerSecond: Double = 0) throws -> TetherLoad {
        guard first.id == firstBody, second.id == secondBody, payoutMetresPerSecond.isFinite else { throw LandingMechanicsError.invalidTopology }
        let displacement=second.point(secondMountMetres)-first.point(firstMountMetres)
        let distance=sqrt((displacement*displacement).sum())
        guard distance.isFinite else { throw LandingMechanicsError.numericalRange }
        guard connected, distance > lengthMetres else { return TetherLoad(forceOnFirst: .zero, tensionNewtons:0, exceedsBreakingLoad:false) }
        let direction=displacement/distance
        let speed=((second.pointVelocity(secondMountMetres)-first.pointVelocity(firstMountMetres))*direction).sum()-payoutMetresPerSecond
        let tension=max(0,stiffnessNewtonsPerMetre*(distance-lengthMetres)+dampingNewtonSecondsPerMetre*speed)
        guard tension.isFinite else { throw LandingMechanicsError.numericalRange }
        return TetherLoad(forceOnFirst:direction*tension,tensionNewtons:tension,exceedsBreakingLoad:tension > breakingTensionNewtons)
    }
}
public struct TetherLoad {
    public let forceOnFirst: SIMD3<Double>
    public let tensionNewtons: Double
    public let exceedsBreakingLoad: Bool
}

/// Adiabatic ideal-gas cushion, with atmospheric pressure subtracted. Inflation
/// and venting are explicit gas-inventory changes (pressure-volume^gamma units).
/// Effective area and uncompressed volume are engineering inputs, not a soft mesh.
public struct GasCushion {
    public let uncompressedVolumeCubicMetres: Double
    public let minimumVolumeFraction: Double
    public let areaSquareMetres: Double
    public let heatCapacityRatio: Double
    public private(set) var gasInventory: Double
    public init(volumeCubicMetres: Double, areaSquareMetres: Double, absolutePressurePascals: Double,
        heatCapacityRatio: Double = 1.4, minimumVolumeFraction: Double = 0.15) throws {
        guard [volumeCubicMetres,areaSquareMetres,absolutePressurePascals,heatCapacityRatio].allSatisfy({$0.isFinite && $0 > 0}),
              heatCapacityRatio > 1, minimumVolumeFraction.isFinite, minimumVolumeFraction > 0, minimumVolumeFraction < 1 else {
            throw LandingMechanicsError.invalidParameter
        }
        uncompressedVolumeCubicMetres=volumeCubicMetres; self.areaSquareMetres=areaSquareMetres
        self.heatCapacityRatio=heatCapacityRatio; self.minimumVolumeFraction=minimumVolumeFraction
        gasInventory=absolutePressurePascals*pow(volumeCubicMetres,heatCapacityRatio)
        guard gasInventory.isFinite else { throw LandingMechanicsError.numericalRange }
    }
    /// Fraction of current gas retained. Venting never adds gas or stored energy.
    public mutating func vent(retaining fraction: Double) throws {
        guard fraction.isFinite, (0...1).contains(fraction) else { throw LandingMechanicsError.invalidParameter }
        gasInventory *= fraction
    }
    public func forceNewtons(compressionMetres: Double, ambientPressurePascals: Double) throws -> Double {
        guard compressionMetres.isFinite, compressionMetres >= 0, ambientPressurePascals.isFinite, ambientPressurePascals >= 0 else {
            throw LandingMechanicsError.invalidParameter
        }
        let volume=max(uncompressedVolumeCubicMetres*minimumVolumeFraction,uncompressedVolumeCubicMetres-areaSquareMetres*compressionMetres)
        let force=max(0,gasInventory/pow(volume,heatCapacityRatio)-ambientPressurePascals)*areaSquareMetres
        guard force.isFinite else { throw LandingMechanicsError.numericalRange }
        return force
    }
}

/// Local terrain plane and its instantaneous surface velocity. The terrain
/// provider must refresh these in the same reference frame and epoch as bodies.
public struct MechanicalGround {
    public let pointMetres: SIMD3<Double>
    public let normal: SIMD3<Double>
    public let velocityMetresPerSecond: SIMD3<Double>
    public init(pointMetres: SIMD3<Double>, normal: SIMD3<Double>, velocityMetresPerSecond: SIMD3<Double> = .zero) throws {
        guard attitudeFinite(pointMetres), attitudeFinite(velocityMetresPerSecond), attitudeFinite(normal),
              (normal*normal).sum().isFinite, (normal*normal).sum() > 1e-24 else { throw LandingMechanicsError.invalidParameter }
        self.pointMetres=pointMetres; self.normal=normal/sqrt((normal*normal).sum()); self.velocityMetresPerSecond=velocityMetresPerSecond
    }
}

/// Spherical contact lobe/foot at a body-local mount. Spring plus compression-
/// only damping; regularized Coulomb friction cannot add tangential energy.
/// A gas cushion adds pressure support; the spring remains the bottom-out guard.
public struct MechanicalContact {
    public let body: String
    public let mountMetres: SIMD3<Double>
    public let radiusMetres: Double
    public let stiffnessNewtonsPerMetre: Double
    public let dampingNewtonSecondsPerMetre: Double
    public let frictionCoefficient: Double
    public let frictionTransitionMetresPerSecond: Double
    public let rollingResistanceMetres: Double
    public var cushion: GasCushion?
    public init(body: String, mountMetres: SIMD3<Double>, radiusMetres: Double,
        stiffnessNewtonsPerMetre: Double, dampingNewtonSecondsPerMetre: Double,
        frictionCoefficient: Double = 0.6, frictionTransitionMetresPerSecond: Double = 0.05, rollingResistanceMetres: Double = 0, cushion: GasCushion? = nil) throws {
        guard !body.isEmpty, attitudeFinite(mountMetres),
              [radiusMetres,stiffnessNewtonsPerMetre,dampingNewtonSecondsPerMetre,frictionCoefficient,rollingResistanceMetres].allSatisfy({$0.isFinite && $0 >= 0}),
              frictionTransitionMetresPerSecond.isFinite, frictionTransitionMetresPerSecond > 0 else { throw LandingMechanicsError.invalidParameter }
        self.body=body; self.mountMetres=mountMetres; self.radiusMetres=radiusMetres
        self.stiffnessNewtonsPerMetre=stiffnessNewtonsPerMetre; self.dampingNewtonSecondsPerMetre=dampingNewtonSecondsPerMetre
        self.frictionCoefficient=frictionCoefficient; self.frictionTransitionMetresPerSecond=frictionTransitionMetresPerSecond; self.cushion=cushion; self.rollingResistanceMetres=rollingResistanceMetres
    }
    public func load(on state: MechanicalBody, ground: MechanicalGround, ambientPressurePascals: Double = 0) throws -> MechanicalWrench {
        guard state.id == body else { throw LandingMechanicsError.invalidTopology }
        let penetration=radiusMetres-((state.point(mountMetres)-ground.pointMetres)*ground.normal).sum()
        guard penetration > 0 else { return MechanicalWrench() }
        let local=mountMetres-state.bodyToReference.invMatrix*(ground.normal*radiusMetres)
        let velocity=state.pointVelocity(local)-ground.velocityMetresPerSecond
        let normalSpeed=(velocity*ground.normal).sum()
        let normalForce=stiffnessNewtonsPerMetre*penetration+dampingNewtonSecondsPerMetre*max(0,-normalSpeed)
            + (try cushion?.forceNewtons(compressionMetres:penetration,ambientPressurePascals:ambientPressurePascals) ?? 0)*min(1,penetration/max(0.001,radiusMetres*0.5))
        let tangent=velocity-ground.normal*normalSpeed, speed=sqrt((tangent*tangent).sum())
        let force=ground.normal*normalForce-tangent*(frictionCoefficient*normalForce/max(frictionTransitionMetresPerSecond,speed))
        guard attitudeFinite(force) else { throw LandingMechanicsError.numericalRange }
        var result=MechanicalWrench(); result.apply(force:force,at:local,to:state)
        let nBody=state.bodyToReference.invMatrix*ground.normal
        let roll=state.angularVelocityBody-nBody*(state.angularVelocityBody*nBody).sum()
        result.torqueBodyNewtonMetres -= roll*(normalForce*rollingResistanceMetres/max(0.1,sqrt((roll*roll).sum())))
        return result
    }
}

/// Shared inertial frame and epoch, independent 6-DOF bodies, arbitrary cable graph.
/// Explicit midpoint integration with bounded user-selected substeps. This is a
/// compliant model: select timestep by stiffness/mass and verify convergence.
/// External loads must include gravity/controls and be pure at each stage.
public struct MechanicalSystem {
    public let frame: Frame
    public private(set) var epoch: Epoch
    public private(set) var bodies: [MechanicalBody]
    public private(set) var tethers: [MechanicalTether]
    public init(frame: Frame, epoch: Epoch, bodies: [MechanicalBody], tethers: [MechanicalTether] = []) throws {
        let ids=Set(bodies.map(\.id))
        guard frame.pseudoInertial, !bodies.isEmpty, ids.count == bodies.count,
              Set(tethers.map(\.id)).count == tethers.count,
              tethers.allSatisfy({ids.contains($0.firstBody) && ids.contains($0.secondBody)}) else { throw LandingMechanicsError.invalidTopology }
        self.frame=frame; self.epoch=epoch; self.bodies=bodies; self.tethers=tethers
    }
    public mutating func cut(tether id: String) throws {
        guard let i=tethers.firstIndex(where:{$0.id == id}) else { throw LandingMechanicsError.invalidTopology }
        tethers[i].cut()
    }
    public func cableLoads(payout: [String:Double] = [:]) throws -> [TetherLoad] {
        try tethers.map { t in try t.load(first:bodies.first{$0.id == t.firstBody}!,second:bodies.first{$0.id == t.secondBody}!,payoutMetresPerSecond:payout[t.id] ?? 0) }
    }
    /// Break events are committed at accepted step boundaries, never during load
    /// evaluation. Returns IDs of newly broken cables. Entire step is transactional.
    @discardableResult public mutating func step(seconds: Double, payout: [String:Double] = [:],
        externalLoads: (Epoch,[MechanicalBody]) throws -> [MechanicalWrench]) throws -> [String] {
        guard seconds.isFinite, seconds > 0, seconds <= 0.01,
              payout.allSatisfy({ key,value in value.isFinite && tethers.contains{$0.id == key} }) else { throw LandingMechanicsError.invalidParameter }
        var next=self
        var broken:[String]=[]
        func breakOverloaded() throws {
            let loads=try next.cableLoads(payout:payout)
            for i in loads.indices where next.tethers[i].connected && loads[i].exceedsBreakingLoad {
                broken.append(next.tethers[i].id); next.tethers[i].cut()
            }
        }
        try breakOverloaded()
        let start=next.bodies
        func rates(_ states:[MechanicalBody], at time:Epoch, cableFraction:Double) throws -> [MechanicalWrench] {
            var loads=try externalLoads(time,states)
            guard loads.count == states.count else { throw LandingMechanicsError.invalidTopology }
            for t in next.tethers {
                let i=states.firstIndex{$0.id == t.firstBody}!, j=states.firstIndex{$0.id == t.secondBody}!
                var cable=t; try cable.setLength(metres:t.lengthMetres+(payout[t.id] ?? 0)*seconds*cableFraction)
                let f=try cable.load(first:states[i],second:states[j],payoutMetresPerSecond:payout[t.id] ?? 0).forceOnFirst
                loads[i].apply(force:f,at:t.firstMountMetres,to:states[i]); loads[j].apply(force:-f,at:t.secondMountMetres,to:states[j])
            }
            return loads
        }
        func advanced(_ base:MechanicalBody, using sample:MechanicalBody, load:MechanicalWrench, dt:Double) throws -> MechanicalBody {
            let momentum=try sample.inertia.angularMomentumBody(for:sample.angularVelocityBody)
            let angularAcceleration=try sample.inertia.solve(load.torqueBodyNewtonMetres-attitudeCross(sample.angularVelocityBody,momentum))
            let qdot=0.5*(sample.bodyToReference*Quaternion(x:sample.angularVelocityBody.x,y:sample.angularVelocityBody.y,z:sample.angularVelocityBody.z,w:0))
            return try MechanicalBody(id:base.id,massKilograms:base.massKilograms,inertia:base.inertia,
                positionMetres:base.positionMetres+sample.velocityMetresPerSecond*dt,
                velocityMetresPerSecond:base.velocityMetresPerSecond+load.forceReferenceNewtons*(dt/base.massKilograms),
                bodyToReference:base.bodyToReference+dt*qdot,angularVelocityBody:base.angularVelocityBody+angularAcceleration*dt)
        }
        let first=try rates(start,at:epoch,cableFraction:0)
        let half=try start.indices.map{try advanced(start[$0],using:start[$0],load:first[$0],dt:seconds/2)}
        let middle=try rates(half,at:epoch.advanced(by:Duration(seconds:seconds/2)),cableFraction:0.5)
        next.bodies=try start.indices.map{try advanced(start[$0],using:half[$0],load:middle[$0],dt:seconds)}
        for i in next.tethers.indices {
            try next.tethers[i].setLength(metres:next.tethers[i].lengthMetres+(payout[next.tethers[i].id] ?? 0)*seconds)
        }
        next.epoch=try epoch.advanced(by:Duration(seconds:seconds))
        try breakOverloaded()
        self=next; return broken
    }
}
private func mechanicalUnit(_ q:Quaternion) throws -> Quaternion {
    let norm=sqrt(q.x*q.x+q.y*q.y+q.z*q.z+q.w*q.w)
    guard norm.isFinite, norm > 1e-12 else { throw LandingMechanicsError.invalidParameter }
    return Quaternion(x:q.x/norm,y:q.y/norm,z:q.z/norm,w:q.w/norm)
}
