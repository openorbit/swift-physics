// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Mattias Holm
import Foundation
import Math
import CelMek

public enum HarpoonPhase: String, Codable { case stowed, flying, anchored, rejected, pulledOut, misfire }
public enum AnchorCableState: String, Codable { case connected, cut, broken }
/// Explicit engineering assumptions, not an inference of subsurface strength from
/// the visible terrain. Normal kinetic energy is spent against constant resistance.
public struct HarpoonSubstrate {
    public let penetrationResistanceNewtons, minimumDepthMetres, maximumDepthMetres, holdingForceNewtons, minimumIncidenceCosine: Double
    public init(penetrationResistanceNewtons: Double, minimumDepthMetres: Double, maximumDepthMetres: Double,
        holdingForceNewtons: Double, minimumIncidenceCosine: Double = 0.5) throws {
        guard [penetrationResistanceNewtons,minimumDepthMetres,maximumDepthMetres,holdingForceNewtons].allSatisfy({$0.isFinite && $0 > 0}),
              maximumDepthMetres >= minimumDepthMetres, minimumIncidenceCosine.isFinite, (0...1).contains(minimumIncidenceCosine) else { throw LandingMechanicsError.invalidParameter }
        self.penetrationResistanceNewtons=penetrationResistanceNewtons; self.minimumDepthMetres=minimumDepthMetres
        self.maximumDepthMetres=maximumDepthMetres; self.holdingForceNewtons=holdingForceNewtons; self.minimumIncidenceCosine=minimumIncidenceCosine
    }
}
/// Terrain anchor state is in the declared rotating surface frame. Its tether
/// couples to an independent inertial body; it never reparents or grounds it.
public struct HarpoonAnchor {
    public let id:String
    public let surfaceFrame:Frame
    public private(set) var phase:HarpoonPhase = .stowed
    public private(set) var cableState:AnchorCableState = .connected
    public private(set) var surfacePointMetres:SIMD3<Double>?
    public private(set) var penetrationMetres=0.0
    public private(set) var holdingForceNewtons=0.0
    public init(id:String,surfaceFrame:Frame) throws {
        guard !id.isEmpty else { throw LandingMechanicsError.invalidParameter }
        self.id=id; self.surfaceFrame=surfaceFrame
    }
    public mutating func fire(succeeded:Bool = true) throws {
        guard phase == .stowed else { throw LandingMechanicsError.invalidParameter }
        phase=succeeded ? .flying : .misfire
    }
    /// Returns true only for a captured, incoming impact within the supported
    /// burial layer. Rejected impacts remain dynamic and can bounce away.
    @discardableResult public mutating func impact(pointMetres:SIMD3<Double>,in frame:Frame,
        massKilograms:Double,normalSpeedMetresPerSecond:Double,incidenceCosine:Double,substrate:HarpoonSubstrate) throws -> Bool {
        guard phase == .flying,frame === surfaceFrame,mechanicsFinite(pointMetres),massKilograms.isFinite,massKilograms > 0,
              normalSpeedMetresPerSecond.isFinite,normalSpeedMetresPerSecond > 0,incidenceCosine.isFinite,(0...1).contains(incidenceCosine) else { throw LandingMechanicsError.invalidParameter }
        let depth=0.5*massKilograms*normalSpeedMetresPerSecond*normalSpeedMetresPerSecond/substrate.penetrationResistanceNewtons
        guard depth.isFinite else { throw LandingMechanicsError.numericalRange }
        penetrationMetres=depth
        if incidenceCosine < substrate.minimumIncidenceCosine || depth < substrate.minimumDepthMetres || depth > substrate.maximumDepthMetres {
            phase = .rejected; return false
        }
        surfacePointMetres=pointMetres; holdingForceNewtons=substrate.holdingForceNewtons;phase = .anchored;return true
    }
    /// Commit only at accepted integration boundaries. Break precedence is
    /// explicit when both thresholds are crossed. No change to vehicle velocity.
    public mutating func commitLoad(tensionNewtons:Double,breakingTensionNewtons:Double) throws {
        guard tensionNewtons.isFinite,tensionNewtons >= 0,breakingTensionNewtons.isFinite,breakingTensionNewtons > 0 else { throw LandingMechanicsError.invalidParameter }
        guard cableState == .connected else { return }
        if tensionNewtons > breakingTensionNewtons { cableState = .broken }
        else if phase == .anchored && tensionNewtons > holdingForceNewtons { phase = .pulledOut }
    }
    public mutating func cut() { cableState = .cut }
}

/// A prescribed terrain endpoint uses the same tension law as two-body cables.
/// The equal and opposite force is a reaction on the prescribed celestial body,
/// not a force to apply a second time to the lander or an integrated dummy mass.
public struct TerrainTether {
    public private(set) var cable:MechanicalTether
    public init(id:String,body:String,mountMetres:SIMD3<Double>,lengthMetres:Double,
        stiffnessNewtonsPerMetre:Double,dampingNewtonSecondsPerMetre:Double,breakingTensionNewtons:Double) throws {
        cable=try MechanicalTether(id:id,firstBody:body,secondBody:"terrain:\(id)",firstMountMetres:mountMetres,
            lengthMetres:lengthMetres,stiffnessNewtonsPerMetre:stiffnessNewtonsPerMetre,
            dampingNewtonSecondsPerMetre:dampingNewtonSecondsPerMetre,breakingTensionNewtons:breakingTensionNewtons)
    }
    public mutating func setLength(metres:Double) throws { try cable.setLength(metres:metres) }
    public func load(body:MechanicalBody,pointReferenceMetres:SIMD3<Double>,velocityReferenceMetresPerSecond:SIMD3<Double>,payoutMetresPerSecond:Double=0) throws -> TetherLoad {
        let endpoint=try MechanicalBody(id:cable.secondBody,massKilograms:1,inertia:RigidBodyInertia(xx:1,yy:1,zz:1),positionMetres:pointReferenceMetres,velocityMetresPerSecond:velocityReferenceMetresPerSecond)
        return try cable.load(first:body,second:endpoint,payoutMetresPerSecond:payoutMetresPerSecond)
    }
}

public struct MechanicalLaunch {
    public let retained,projectile:MechanicalBody
    public let retainedCentreBodyMetres:SIMD3<Double>
}
extension MechanicalBody {
    /// Separate a small rigid component, then apply an internal firing impulse
    /// at the component centre. Mass, linear and angular momentum are conserved;
    /// launch energy is supplied by the firing mechanism. Inertia is explicitly
    /// removed and translated to the remaining body's new centre of mass.
    public func launching(id projectileID:String,massKilograms m:Double,inertia projectileInertia:RigidBodyInertia,
        centreBodyMetres r:SIMD3<Double>,impulseBodyNewtonSeconds j:SIMD3<Double>) throws -> MechanicalLaunch {
        guard m.isFinite,m > 0,m < massKilograms,projectileID != id,mechanicsFinite(r),mechanicsFinite(j) else { throw LandingMechanicsError.invalidParameter }
        let remaining=massKilograms-m,s = -r*(m/remaining)
        func parallel(_ v:SIMD3<Double>,_ mass:Double,_ i:Int,_ k:Int)->Double { mass*((i == k ? (v*v).sum() : 0)-v[i]*v[k]) }
        func entry(_ i:Int,_ k:Int)->Double { inertia.matrix[i,k]-projectileInertia.matrix[i,k]-parallel(r,m,i,k)-parallel(s,remaining,i,k) }
        let retainedInertia=try RigidBodyInertia(xx:entry(0,0),yy:entry(1,1),zz:entry(2,2),xy:entry(0,1),xz:entry(0,2),yz:entry(1,2))
        let impulse=bodyToReference.matrix*j
        let retained=try MechanicalBody(id:id,massKilograms:remaining,inertia:retainedInertia,positionMetres:point(s),
            velocityMetresPerSecond:pointVelocity(s)-impulse/remaining,bodyToReference:bodyToReference,
            angularVelocityBody:angularVelocityBody+retainedInertia.solve(mechanicsCross(r-s,-j)))
        let projectile=try MechanicalBody(id:projectileID,massKilograms:m,inertia:projectileInertia,positionMetres:point(r),
            velocityMetresPerSecond:pointVelocity(r)+impulse/m,bodyToReference:bodyToReference,angularVelocityBody:angularVelocityBody)
        return MechanicalLaunch(retained:retained,projectile:projectile,retainedCentreBodyMetres:s)
    }
}
