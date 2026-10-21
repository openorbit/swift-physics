// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Mattias Holm
import Foundation
import Math
import CelMek

public struct WheelCommand {
    public let motorTorqueNewtonMetres,brakeFraction,steeringRadians:Double
    public init(motorTorqueNewtonMetres:Double=0,brakeFraction:Double=0,steeringRadians:Double=0) throws {
        guard motorTorqueNewtonMetres.isFinite,brakeFraction.isFinite,(0...1).contains(brakeFraction),
              steeringRadians.isFinite,abs(steeringRadians)<=0.7 else {throw LandingMechanicsError.invalidParameter}
        self.motorTorqueNewtonMetres=motorTorqueNewtonMetres;self.brakeFraction=brakeFraction;self.steeringRadians=steeringRadians
    }
}
public struct WheelState {
    public let spinRadiansPerSecond,angleRadians,steeringRadians:Double
    public init(spinRadiansPerSecond:Double=0,angleRadians:Double=0,steeringRadians:Double=0) throws {
        guard [spinRadiansPerSecond,angleRadians,steeringRadians].allSatisfy(\.isFinite),abs(steeringRadians)<=0.7 else {throw LandingMechanicsError.invalidParameter}
        self.spinRadiansPerSecond=spinRadiansPerSecond;self.angleRadians=angleRadians;self.steeringRadians=steeringRadians
    }
    func advanced(_ sample:WheelSample,seconds:Double) throws -> WheelState {
        try WheelState(spinRadiansPerSecond:spinRadiansPerSecond+sample.spinAcceleration*seconds,
            angleRadians:(angleRadians+sample.spinRate*seconds).truncatingRemainder(dividingBy:2 * .pi),
            steeringRadians:steeringRadians+sample.steeringRate*seconds)
    }
}
/// Body axes: +X forward, +Y left, +Z suspension up. Positive spin is forward
/// rolling about +Y. Geometry/mass/inertia of wheels are already in the chassis;
/// only relative rotor momentum is added, never another translational mass.
public struct SuspendedWheel {
    public let id:String,mountBodyMetres:SIMD3<Double>
    public let radiusMetres,rotorInertia,restLengthMetres,travelMetres,spring,damping,friction,slipStiffness,maximumBrakeTorque,rollingResistance:Double
    public init(id:String,mountBodyMetres:SIMD3<Double>,radiusMetres:Double=0.28,rotorInertia:Double=0.25,
        restLengthMetres:Double=0.24,travelMetres:Double=0.16,springNewtonsPerMetre:Double=3500,dampingNewtonSecondsPerMetre:Double=400,
        friction:Double=0.8,slipStiffnessNewtonSecondsPerMetre:Double=250,maximumBrakeTorque:Double=15,rollingResistance:Double=0.02) throws {
        guard !id.isEmpty,mechanicsFinite(mountBodyMetres),[radiusMetres,rotorInertia,restLengthMetres,travelMetres,springNewtonsPerMetre,slipStiffnessNewtonSecondsPerMetre].allSatisfy({$0.isFinite && $0>0}),
              [dampingNewtonSecondsPerMetre,friction,maximumBrakeTorque,rollingResistance].allSatisfy({$0.isFinite && $0>=0}),travelMetres<restLengthMetres else {throw LandingMechanicsError.invalidParameter}
        self.id=id;self.mountBodyMetres=mountBodyMetres;self.radiusMetres=radiusMetres;self.rotorInertia=rotorInertia
        self.restLengthMetres=restLengthMetres;self.travelMetres=travelMetres;spring=springNewtonsPerMetre;damping=dampingNewtonSecondsPerMetre
        self.friction=friction;slipStiffness=slipStiffnessNewtonSecondsPerMetre;self.maximumBrakeTorque=maximumBrakeTorque;self.rollingResistance=rollingResistance
    }
    /// Massless suspension strut and planar contact sample. Rotor acceleration,
    /// steering angular momentum transport and gyroscopic reactions all enter
    /// the chassis wrench. No contact means no driving or braking ground force.
    public func evaluate(body:MechanicalBody,state:WheelState,ground:MechanicalGround,command:WheelCommand,stepSeconds:Double) throws -> WheelSample {
        guard stepSeconds.isFinite,stepSeconds>0,stepSeconds<=0.01 else {throw LandingMechanicsError.invalidParameter}
        let q=body.bodyToReference,up=q.matrix*SIMD3<Double>(0,0,1),alignment=(up*ground.normal).sum()
        let steeringRate=max(-0.7,min(0.7,(command.steeringRadians-state.steeringRadians)/0.15))
        let forwardBody=SIMD3(cos(state.steeringRadians),sin(state.steeringRadians),0.0),axleBody=SIMD3(-forwardBody.y,forwardBody.x,0)
        let axleWorld=q.matrix*axleBody
        let height=((body.point(mountBodyMetres)-ground.pointMetres)*ground.normal).sum()
        let requiredLength=alignment>0.2 ? (height-radiusMetres)/alignment : restLengthMetres+travelMetres+1
        let length=max(restLengthMetres-travelMetres,min(restLengthMetres+travelMetres,requiredLength))
        let centre=mountBodyMetres-SIMD3(0,0,length),point=centre-q.invMatrix*(ground.normal*radiusMetres)
        let contact=alignment>0.2 && requiredLength<=restLengthMetres+travelMetres
        var force=SIMD3<Double>.zero,normalLoad=0.0,slip=0.0,lateralSpeed=0.0,groundTorque=0.0
        if contact {
            let velocity=body.pointVelocity(point)-ground.velocityMetresPerSecond,vn=(velocity*ground.normal).sum()
            normalLoad=max(0,spring*max(0,restLengthMetres-requiredLength)-damping*vn)
                + spring*8*max(0,restLengthMetres-travelMetres-requiredLength)
            var forward=q.matrix*forwardBody;forward -= ground.normal*(forward*ground.normal).sum()
            forward /= sqrt((forward*forward).sum())
            let side=mechanicsCross(ground.normal,forward)
            let effectiveRadius=radiusMetres*(axleWorld*side).sum()
            slip=effectiveRadius*state.spinRadiansPerSecond-(velocity*forward).sum();lateralSpeed=(velocity*side).sum()
            var longitudinal=slipStiffness*slip,lateral = -slipStiffness*lateralSpeed
            let demand=hypot(longitudinal,lateral),limit=friction*normalLoad
            if demand>limit && demand>0 {longitudinal *= limit/demand;lateral *= limit/demand}
            force=ground.normal*normalLoad+forward*longitudinal+side*lateral
            groundTorque=(mechanicsCross(-ground.normal*radiusMetres,force)*axleWorld).sum()
        }
        let rollingTorque=rollingResistance*normalLoad*radiusMetres*tanh(state.spinRadiansPerSecond/0.1)
        let freeTorque=command.motorTorqueNewtonMetres+groundTorque-rollingTorque
        let brakeLimit=maximumBrakeTorque*command.brakeFraction
        let brake=max(-brakeLimit,min(brakeLimit,freeTorque+rotorInertia*state.spinRadiansPerSecond/stepSeconds))
        let acceleration=(freeTorque-brake)/rotorInertia
        var wrench=MechanicalWrench();wrench.apply(force:force,at:point,to:body)
        let momentum=axleBody*(rotorInertia*state.spinRadiansPerSecond)
        let momentumRate=axleBody*(rotorInertia*acceleration)-forwardBody*(steeringRate*rotorInertia*state.spinRadiansPerSecond)
        wrench.torqueBodyNewtonMetres -= momentumRate+mechanicsCross(body.angularVelocityBody,momentum)
        guard mechanicsFinite(wrench.forceReferenceNewtons),mechanicsFinite(wrench.torqueBodyNewtonMetres),[acceleration,normalLoad,slip,lateralSpeed].allSatisfy(\.isFinite) else {throw LandingMechanicsError.numericalRange}
        return WheelSample(wrench:wrench,spinAcceleration:acceleration,spinRate:state.spinRadiansPerSecond,steeringRate:steeringRate,
            centreBodyMetres:centre,normalLoad:normalLoad,compression:restLengthMetres-length,longitudinalSlip:slip,lateralSpeed:lateralSpeed)
    }
}
public struct WheelSample {
    public let wrench:MechanicalWrench
    public let spinAcceleration,spinRate,steeringRate:Double
    public let centreBodyMetres:SIMD3<Double>
    public let normalLoad,compression,longitudinalSlip,lateralSpeed:Double
}
/// Couples relative wheel spin and the existing multi-body midpoint integrator at
/// the same two stages. External loads/terrain must be pure stage evaluations.
/// Failed steps leave both chassis and wheel state unchanged.
public struct SuspendedWheelSet {
    public let bodyID:String,wheels:[SuspendedWheel]
    public private(set) var states:[WheelState]
    public init(bodyID:String,wheels:[SuspendedWheel],states:[WheelState]?=nil) throws {
        guard !bodyID.isEmpty,!wheels.isEmpty,Set(wheels.map(\.id)).count==wheels.count,states == nil || states!.count==wheels.count else {throw LandingMechanicsError.invalidTopology}
        self.bodyID=bodyID;self.wheels=wheels;self.states=try states ?? wheels.map{_ in try WheelState()}
    }
    public func samples(body:MechanicalBody,commands:[WheelCommand],stepSeconds:Double,
        ground:(MechanicalBody,SuspendedWheel) throws -> MechanicalGround) throws -> [WheelSample] {
        guard body.id==bodyID,commands.count==wheels.count else {throw LandingMechanicsError.invalidTopology}
        return try wheels.indices.map{try wheels[$0].evaluate(body:body,state:states[$0],ground:ground(body,wheels[$0]),command:commands[$0],stepSeconds:stepSeconds)}
    }
    @discardableResult public mutating func step(system:inout MechanicalSystem,seconds:Double,commands:[WheelCommand],payout:[String:Double]=[:],
        ground:(Epoch,MechanicalBody,SuspendedWheel) throws -> MechanicalGround,
        externalLoads:(Epoch,[MechanicalBody]) throws -> [MechanicalWrench]) throws -> [String] {
        guard commands.count==wheels.count,let bodyIndex=system.bodies.firstIndex(where:{$0.id==bodyID}) else {throw LandingMechanicsError.invalidTopology}
        let initial=states,config=wheels,startEpoch=system.epoch
        var half=states,midpoint:[WheelSample]=[],next=system
        let broken=try next.step(seconds:seconds,payout:payout){epoch,bodies in
            let second=epoch>startEpoch,samples=try config.indices.map { i in
                try config[i].evaluate(body:bodies[bodyIndex],state:second ? half[i] : initial[i],ground:ground(epoch,bodies[bodyIndex],config[i]),command:commands[i],stepSeconds:seconds)
            }
            if second {midpoint=samples}
            else {half=try initial.indices.map{try initial[$0].advanced(samples[$0],seconds:seconds/2)}}
            var loads=try externalLoads(epoch,bodies)
            guard loads.count==bodies.count else {throw LandingMechanicsError.invalidTopology}
            for sample in samples {
                loads[bodyIndex].forceReferenceNewtons += sample.wrench.forceReferenceNewtons
                loads[bodyIndex].torqueBodyNewtonMetres += sample.wrench.torqueBodyNewtonMetres
            }
            return loads
        }
        guard midpoint.count==initial.count else {throw LandingMechanicsError.numericalRange}
        let final=try initial.indices.map{try initial[$0].advanced(midpoint[$0],seconds:seconds)}
        system=next;states=final;return broken
    }
}
