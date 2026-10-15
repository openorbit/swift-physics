//
// SPDX-License-Identifier: Apache-2.0
//
// Copyright 2026 Mattias Holm
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
import Foundation
import Math

private func encounterCross(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> SIMD3<Double> {
    SIMD3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x)
}
private func encounterDot(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> Double {
    a.x*b.x+a.y*b.y+a.z*b.z
}
public enum EncounterPlanningError: Error, Equatable {
    case degenerateAxes, invalidGeometry, invalidWindow, insufficientPropellant, invalidExhaustVelocity
}

/// S follows incoming velocity; T = unit(S × pole), R = S × T.
/// B is the incoming asymptote offset: B × S follows angular momentum.
public struct FlybyBPlane {
    public let s: SIMD3<Double>
    public let t: SIMD3<Double>
    public let r: SIMD3<Double>
    public let bTMetres: Double
    public let bRMetres: Double
    public var vectorMetres: SIMD3<Double> { t*bTMetres+r*bRMetres }
    public init(incomingVelocity: SIMD3<Double>, referencePole: SIMD3<Double>, bTMetres: Double, bRMetres: Double) throws {
        let speed = routeNorm(incomingVelocity), poleLength = routeNorm(referencePole)
        guard speed.isFinite, speed > 0, poleLength.isFinite, poleLength > 0,
              bTMetres.isFinite, bRMetres.isFinite, hypot(bTMetres,bRMetres).isFinite,
              hypot(bTMetres,bRMetres) > 0 else { throw EncounterPlanningError.invalidGeometry }
        s = incomingVelocity/speed
        let axis = encounterCross(s,referencePole/poleLength), length = routeNorm(axis)
        guard length > 1e-12 else { throw EncounterPlanningError.degenerateAxes }
        t = axis/length; r = encounterCross(s,t)
        self.bTMetres = bTMetres; self.bRMetres = bRMetres
    }
    public init(flyby: UnpoweredFlybyResult, referencePole: SIMD3<Double>) throws {
        let b = flyby.periapsisRadiusMetres * (flyby.periapsisSpeedMetresPerSecond/flyby.excessSpeedMetresPerSecond)
        let vector = encounterCross(flyby.request.incomingExcessVelocityMetresPerSecond/flyby.excessSpeedMetresPerSecond,flyby.orbitalNormal)*b
        let axes = try Self(incomingVelocity: flyby.request.incomingExcessVelocityMetresPerSecond,
            referencePole: referencePole,bTMetres: 1,bRMetres: 0)
        try self.init(incomingVelocity: flyby.request.incomingExcessVelocityMetresPerSecond,
            referencePole: referencePole,bTMetres: encounterDot(vector,axes.t),bRMetres: encounterDot(vector,axes.r))
    }
    /// The request supplies body, epoch, axes, incoming speed and safety limits;
    /// its outgoing vector is replaced by the B-plane prediction.
    public func solve(_ request: FlybyEncounterRequest) throws -> UnpoweredFlybyResult {
        let speed = routeNorm(request.incomingExcessVelocityMetresPerSecond)
        guard routeNorm(request.incomingExcessVelocityMetresPerSecond/speed-s) < 1e-12 else {
            throw EncounterPlanningError.invalidGeometry
        }
        let b = hypot(bTMetres,bRMetres), length = request.gravitationalParameter.cubicMetresPerSecondSquared/speed/speed
        let radius = b*(b/(hypot(length,b)+length))
        let c = 1/(1+radius/length), angle = 2*asin(c)
        let normal = encounterCross(vectorMetres/b,s)
        let outgoing = (s*cos(angle)+encounterCross(normal,s)*sin(angle))*speed
        return try UnpoweredFlyby.solve(FlybyEncounterRequest(body: request.body,periapsisEpoch: request.periapsisEpoch,
            frame: request.frame,gravitationalParameter: request.gravitationalParameter,bodyRadiusMetres: request.bodyRadiusMetres,
            minimumPeriapsisAltitudeMetres: request.minimumPeriapsisAltitudeMetres,
            incomingExcessVelocityMetresPerSecond: request.incomingExcessVelocityMetresPerSecond,
            outgoingExcessVelocityMetresPerSecond: outgoing))
    }
}

/// Ideal instantaneous rocket equation, with an explicit exhaust velocity.
public enum EncounterImpulse {
    public static func apply(_ maneuver: ImpulsiveManeuver, to state: SpacecraftState,
        effectiveExhaustVelocityMetresPerSecond: Double) throws -> SpacecraftState {
        guard maneuver.frame === state.frame else { throw SpacecraftDynamicsError.frameMismatch }
        guard maneuver.epoch == state.epoch else { throw SpacecraftDynamicsError.epochMismatch }
        let ve = effectiveExhaustVelocityMetresPerSecond
        guard ve.isFinite, ve > 0 else { throw EncounterPlanningError.invalidExhaustVelocity }
        let cost = state.mass.totalMassKilograms * -expm1(-routeNorm(maneuver.deltaVelocityMetresPerSecond)/ve)
        guard cost.isFinite, cost <= state.mass.propellantMassKilograms else { throw EncounterPlanningError.insufficientPropellant }
        return SpacecraftState(cartesian: try CartesianState(epoch: state.epoch,frame: state.frame,
            position: state.cartesian.position,velocity: state.cartesian.velocity+maneuver.deltaVelocityMetresPerSecond),
            mass: try SpacecraftMassProperties(dryMassKilograms: state.mass.dryMassKilograms,
                propellantMassKilograms: state.mass.propellantMassKilograms-cost))
    }
}
public struct PoweredFlybyResult {
    public let incomingTrajectory: TwoBodyPropagator
    public let outgoingTrajectory: TwoBodyPropagator
    public let maneuver: ImpulsiveManeuver
    public let preBurnState: SpacecraftState
    public let postBurnState: SpacecraftState
    public var propellantUsedKilograms: Double { preBurnState.mass.propellantMassKilograms-postBurnState.mass.propellantMassKilograms }
}
public enum PoweredFlyby {
    /// Explicit common periapsis geometry. Supports unequal excess speeds and
    /// plane-changing impulses; inconsistent radial/asymptote geometry is rejected.
    public static func solve(_ request: FlybyEncounterRequest, periapsisRadiusMetres radius: Double,
        radialDirection: SIMD3<Double>, mass: SpacecraftMassProperties,
        effectiveExhaustVelocityMetresPerSecond: Double) throws -> PoweredFlybyResult {
        guard radius.isFinite, radius > request.bodyRadiusMetres,
              radius >= request.bodyRadiusMetres+request.minimumPeriapsisAltitudeMetres else { throw FlybyError.unsafePeriapsis }
        let norm = routeNorm(radialDirection)
        guard norm.isFinite, norm > 0 else { throw EncounterPlanningError.invalidGeometry }
        let radial = radialDirection/norm, mu = request.gravitationalParameter.cubicMetresPerSecondSquared
        func branch(_ velocity: SIMD3<Double>, incoming: Bool) throws -> TwoBodyPropagator {
            let speed = routeNorm(velocity), u = velocity/speed
            let c = 1/(1+radius*speed*speed/mu), expected = incoming ? c : -c
            guard c > 0, c < 1, abs(encounterDot(u,radial)-expected) < 1e-12 else { throw EncounterPlanningError.invalidGeometry }
            let tangent = u-radial*expected, tangentLength = routeNorm(tangent)
            guard tangentLength > 1e-12 else { throw EncounterPlanningError.invalidGeometry }
            let seed = try CartesianState(epoch: request.periapsisEpoch,frame: request.frame,position: radial*radius,
                velocity: tangent/tangentLength*hypot(speed,sqrt(2)*sqrt(mu)/sqrt(radius)))
            return try TwoBodyPropagator(initialState: seed,centralBody: request.body,gravitationalParameter: request.gravitationalParameter)
        }
        let before = try branch(request.incomingExcessVelocityMetresPerSecond,incoming: true)
        let after = try branch(request.requestedOutgoingExcessVelocityMetresPerSecond,incoming: false)
        let pre = SpacecraftState(cartesian: try before.state(at: request.periapsisEpoch),mass: mass)
        let maneuver = try ImpulsiveManeuver(epoch: request.periapsisEpoch,frame: request.frame,
            deltaVelocityMetresPerSecond: after.initialState.velocity-before.initialState.velocity)
        let post = try EncounterImpulse.apply(maneuver,to: pre,effectiveExhaustVelocityMetresPerSecond: effectiveExhaustVelocityMetresPerSecond)
        return PoweredFlybyResult(incomingTrajectory: before,outgoingTrajectory: after,maneuver: maneuver,preBurnState: pre,postBurnState: post)
    }
}

public enum DeepSpaceManeuverCommand {
    case impulse(ImpulsiveManeuver, effectiveExhaustVelocityMetresPerSecond: Double)
    case finiteBurn(FiniteBurn)
}
public struct DeepSpaceManeuverNode {
    public let window: ClosedRange<Epoch>
    public let command: DeepSpaceManeuverCommand
    public init(window: ClosedRange<Epoch>, command: DeepSpaceManeuverCommand) throws {
        switch command {
        case .impulse(let impulse, let ve):
            guard window.contains(impulse.epoch), ve.isFinite, ve > 0 else { throw EncounterPlanningError.invalidWindow }
        case .finiteBurn(let burn):
            guard window.contains(burn.start), window.contains(burn.end) else { throw EncounterPlanningError.invalidWindow }
        }
        self.window = window; self.command = command
    }
    /// Converts one impulse to an engine-duration seed, placed wholly inside its
    /// window. Call sequentially with the carried mass for multiple nodes.
    public func finiteBurnSeed(initialState: SpacecraftState, thrustNewtons: Double) throws -> Self {
        guard case .impulse(let maneuver,let ve) = command else { throw EncounterPlanningError.invalidGeometry }
        guard initialState.epoch <= window.lowerBound else { throw EncounterPlanningError.invalidWindow }
        let schedule = try FiniteTransferPlanner.initialGuess(initialState:
            SpacecraftState(cartesian: CartesianState(epoch: window.lowerBound,frame: initialState.frame,
                position: initialState.cartesian.position,velocity: initialState.cartesian.velocity),mass: initialState.mass),
            arrivalEpoch: window.upperBound,maneuvers: [maneuver],thrustNewtons: thrustNewtons,
            effectiveExhaustVelocityMetresPerSecond: ve)
        guard let burn = schedule.burns.first else { throw EncounterPlanningError.invalidGeometry }
        return try Self(window: window,command: .finiteBurn(burn))
    }
    /// Refines finite nodes in one propagation context. Returned burn order maps
    /// to chronological node order; every trial retains the respective window.
    public static func refine(_ nodes: [Self], initialState: SpacecraftState, target: CartesianState,
        context: ForceEvaluationContext, backgroundModels: [any ForceModel] = [],
        maximumEvaluations: Int = 500) throws -> FiniteTransferRefinementResult {
        var pairs: [(FiniteBurn,ClosedRange<Epoch>)] = []
        for node in nodes {
            guard case .finiteBurn(let burn) = node.command else { throw EncounterPlanningError.invalidGeometry }
            pairs.append((burn,node.window))
        }
        pairs.sort { $0.0.start < $1.0.start }
        return try FiniteTransferPlanner.refine(initialState: initialState,
            initialSchedule: ScheduledThrust(frame: initialState.frame,burns: pairs.map { $0.0 }),
            target: target,context: context,backgroundModels: backgroundModels,
            maximumEvaluations: maximumEvaluations,burnWindows: pairs.map { $0.1 })
    }
    private var interval: ClosedRange<Epoch> {
        switch command { case .impulse(let m,_): return m.epoch...m.epoch
        case .finiteBurn(let b): return b.start...b.end }
    }
    /// Executes chronologically with caller-supplied gravity; no implicit forces.
    public static func execute(_ nodes: [Self], initialState: SpacecraftState, arrivalEpoch: Epoch,
        context: ForceEvaluationContext, backgroundModels: [any ForceModel] = [],
        options: NumericalPropagationOptions = try! NumericalPropagationOptions()) throws -> SpacecraftState {
        guard nodes.count <= 100, arrivalEpoch >= initialState.epoch else { throw EncounterPlanningError.invalidWindow }
        let ordered = nodes.sorted { $0.interval.lowerBound < $1.interval.lowerBound }
        guard ordered.allSatisfy({ $0.interval.lowerBound >= initialState.epoch && $0.interval.upperBound <= arrivalEpoch }),
              zip(ordered,ordered.dropFirst()).allSatisfy({ $0.interval.upperBound < $1.interval.lowerBound }) else {
            throw EncounterPlanningError.invalidWindow
        }
        var state = initialState
        func propagate(_ to: Epoch, burns: [FiniteBurn] = []) throws -> SpacecraftState {
            try ScheduledThrustPropagator(initialState: state,schedule: ScheduledThrust(frame: state.frame,burns: burns),
                backgroundModels: backgroundModels,context: context,options: options).spacecraftState(at: to)
        }
        for node in ordered {
            state = try propagate(node.interval.lowerBound)
            switch node.command {
            case .impulse(let maneuver,let ve): state = try EncounterImpulse.apply(maneuver,to: state,effectiveExhaustVelocityMetresPerSecond: ve)
            case .finiteBurn(let burn): state = try propagate(burn.end,burns: [burn])
            }
        }
        return try propagate(arrivalEpoch)
    }
}

/// An explicitly chosen epoch is the patch surface. Frame providers must encode
/// the moving origins and axes. This converts the carried state, never a template.
/// Switching force models after this handover is the caller's patched-conic approximation.
public struct EncounterPatchBoundary {
    public let epoch: Epoch
    public let source: ForceEvaluationContext
    public let destination: ForceEvaluationContext
    public init(epoch: Epoch, source: ForceEvaluationContext, destination: ForceEvaluationContext) {
        self.epoch = epoch; self.source = source; self.destination = destination
    }
    public func handover(_ state: SpacecraftState, using tree: FrameTree) throws -> SpacecraftState {
        guard state.epoch == epoch else { throw SpacecraftDynamicsError.epochMismatch }
        guard state.frame === source.frame else { throw SpacecraftDynamicsError.frameMismatch }
        let transform = try tree.transform(from: source.frame,to: destination.frame,at: epoch)
        return SpacecraftState(cartesian: try transform.applying(to: state.cartesian,destination: destination.frame),mass: state.mass)
    }
}
