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

public enum HohmannRendezvousError: Error, Equatable {
    case incompatibleCentralBody, incompatibleGM, nonCoplanarOrbits, opposingOrbitalDirections
    case invalidTolerance, invalidBudget, relativeMotionNotRepresentable
}

/// An explicit circular two-body model. State adapters validate and report
/// the limits of a local circular approximation; they do not fit a full ephemeris.
/// The caller guarantees the frame's physical centre is centralBody.
public struct CircularHohmannOrbit: StateProvider {
    public let referenceEpoch: Epoch
    public let frame: Frame
    public let centralBody: BodyIdentifier
    public let radiusMetres: Double
    public let gravitationalParameter: GravitationalParameter
    public let angularRateRadiansPerSecond: Double
    public let circularSpeedMetresPerSecond: Double
    private let radial: SIMD3<Double>
    private let tangent: SIMD3<Double>
    var normal: SIMD3<Double> { hohmannCross(radial, tangent) }

    public init(referenceEpoch: Epoch, frame: Frame, centralBody: BodyIdentifier,
                radiusMetres: Double, radialDirection: SIMD3<Double>, progradeDirection: SIMD3<Double>,
                gravitationalParameter: GravitationalParameter) throws {
        try centralBody.validate()
        guard frame.pseudoInertial else { throw PropagationError.nonInertialFrame }
        let circular = try HohmannTransfer(departureRadiusMetres: radiusMetres,
            arrivalRadiusMetres: radiusMetres, gravitationalParameter: gravitationalParameter)
        let r = try hohmannUnit(radialDirection), t = try hohmannUnit(progradeDirection)
        let dot = (r*t).sum()
        guard abs(dot) <= 1e-12 else { throw HohmannTransferError.invalidDirections }
        radial = r; tangent = try hohmannUnit(t-r*dot)
        self.referenceEpoch = referenceEpoch; self.frame = frame; self.centralBody = centralBody
        self.radiusMetres = radiusMetres; self.gravitationalParameter = gravitationalParameter
        circularSpeedMetresPerSecond = circular.departureCircularSpeedMetresPerSecond
        angularRateRadiansPerSecond = circularSpeedMetresPerSecond/radiusMetres
        guard angularRateRadiansPerSecond.isFinite, angularRateRadiansPerSecond > 0 else {
            throw HohmannTransferError.numericalRangeExceeded
        }
    }

    func directions(at epoch: Epoch) throws -> (radial: SIMD3<Double>, tangent: SIMD3<Double>) {
        let angle = angularRateRadiansPerSecond * (try epoch.duration(since: referenceEpoch)).seconds
        guard angle.isFinite else { throw HohmannTransferError.numericalRangeExceeded }
        let reduced = angle.truncatingRemainder(dividingBy: 2 * .pi)
        let c = cos(reduced), s = sin(reduced)
        return (radial*c+tangent*s, tangent*c-radial*s)
    }
    public func state(at epoch: Epoch) throws -> CartesianState {
        let axes = try directions(at: epoch)
        return try CartesianState(epoch: epoch, frame: frame,
            position: axes.radial*radiusMetres, velocity: axes.tangent*circularSpeedMetresPerSecond,
            acceleration: -axes.radial*(circularSpeedMetresPerSecond*angularRateRadiansPerSecond))
    }
}

public enum HohmannRendezvousStatus: Equatable, Sendable {
    case feasible, departureWindowExhausted, evaluationLimitReached
    case alreadyRendezvoused, sameOrbitPhasingRequired
}
public struct HohmannRendezvousAttempt {
    public let departureEpoch: Epoch
    public let phaseErrorRadians: Double
    public let positionErrorMetres: Double
    public let velocityErrorMetresPerSecond: Double
}
public struct HohmannRendezvousResult {
    public let status: HohmannRendezvousStatus
    /// Target angle ahead of the departure radius, positive in orbital motion, [0, 2π).
    public let requiredPhaseAngleRadians: Double
    public let relativeAngularRateRadiansPerSecond: Double
    public let attempts: [HohmannRendezvousAttempt]
    /// Nil for already-coincident targets and searches without a verified transfer.
    public let plan: HohmannTransferPlan?
    public let departureEpoch: Epoch?
    public let targetArrivalState: CartesianState?
}

/// Analytic launch opportunities for explicit circular, coplanar, co-directional
/// models, followed by propagated rendezvous checks. No arbitrary epoch grid.
public enum HohmannRendezvousPlanner {
    public static func search(departureOrbit: CircularHohmannOrbit, targetOrbit: CircularHohmannOrbit,
        departureWindow: ClosedRange<Epoch>, maximumOpportunityEvaluations: Int = 1000,
        positionToleranceMetres: Double = 0.01, velocityToleranceMetresPerSecond: Double = 1e-5,
        propagationOptions: TwoBodyPropagationOptions = try! TwoBodyPropagationOptions()) throws -> HohmannRendezvousResult {
        guard (1...10_000).contains(maximumOpportunityEvaluations) else { throw HohmannRendezvousError.invalidBudget }
        guard [positionToleranceMetres, velocityToleranceMetresPerSecond].allSatisfy({ $0.isFinite && $0 > 0 }) else {
            throw HohmannRendezvousError.invalidTolerance
        }
        guard departureOrbit.frame === targetOrbit.frame else { throw SpacecraftDynamicsError.frameMismatch }
        guard departureOrbit.centralBody == targetOrbit.centralBody else { throw HohmannRendezvousError.incompatibleCentralBody }
        let gm = departureOrbit.gravitationalParameter, targetGM = targetOrbit.gravitationalParameter
        guard gm.cubicMetresPerSecondSquared == targetGM.cubicMetresPerSecondSquared,
              gm.convention == targetGM.convention else { throw HohmannRendezvousError.incompatibleGM }
        let n1 = departureOrbit.normal, n2 = targetOrbit.normal
        if hohmannNorm(n1+n2) <= 1e-12 { throw HohmannRendezvousError.opposingOrbitalDirections }
        guard hohmannNorm(n1-n2) <= 1e-12 else { throw HohmannRendezvousError.nonCoplanarOrbits }
        let transfer = try HohmannTransfer(departureRadiusMetres: departureOrbit.radiusMetres,
            arrivalRadiusMetres: targetOrbit.radiusMetres, gravitationalParameter: gm)
        let angle = targetOrbit.angularRateRadiansPerSecond*transfer.duration.seconds
        guard angle.isFinite else { throw HohmannTransferError.numericalRangeExceeded }
        let required = departureOrbit.radiusMetres == targetOrbit.radiusMetres ? 0 : hohmannWrap(.pi-angle)
        // Stable mean-motion difference for nearly equal radii.
        let ratioDifference = (departureOrbit.radiusMetres-targetOrbit.radiusMetres)/targetOrbit.radiusMetres
        let rate = abs(ratioDifference) < 0.5 ?
            departureOrbit.angularRateRadiansPerSecond*expm1(1.5*log1p(ratioDifference)) :
            targetOrbit.angularRateRadiansPerSecond-departureOrbit.angularRateRadiansPerSecond
        guard rate.isFinite else { throw HohmannRendezvousError.relativeMotionNotRepresentable }
        var attempts: [HohmannRendezvousAttempt] = []
        func finish(_ status: HohmannRendezvousStatus, plan: HohmannTransferPlan? = nil,
                    departure: Epoch? = nil, target: CartesianState? = nil) -> HohmannRendezvousResult {
            HohmannRendezvousResult(status: status, requiredPhaseAngleRadians: required,
                relativeAngularRateRadiansPerSecond: rate, attempts: attempts,
                plan: plan, departureEpoch: departure, targetArrivalState: target)
        }
        let lower = departureWindow.lowerBound
        let axes = try departureOrbit.directions(at: lower)
        let targetAxes = try targetOrbit.directions(at: lower)
        let phase = atan2((targetAxes.radial*axes.tangent).sum(), (targetAxes.radial*axes.radial).sum())
        if departureOrbit.radiusMetres == targetOrbit.radiusMetres {
            let departure = try departureOrbit.state(at: lower), target = try targetOrbit.state(at: lower)
            let p = hohmannNorm(departure.position-target.position), v = hohmannNorm(departure.velocity-target.velocity)
            attempts.append(HohmannRendezvousAttempt(departureEpoch: lower, phaseErrorRadians: phase,
                positionErrorMetres: p, velocityErrorMetresPerSecond: v))
            if p <= positionToleranceMetres && v <= velocityToleranceMetresPerSecond {
                return finish(.alreadyRendezvoused, departure: lower, target: target)
            }
            return finish(.sameOrbitPhasingRequired)
        }
        guard rate != 0 else { throw HohmannRendezvousError.relativeMotionNotRepresentable }
        let gap = rate > 0 ? hohmannWrap(required-phase) : hohmannWrap(phase-required)
        let angularTolerance = min(64*Double.ulpOfOne,
            min(positionToleranceMetres/targetOrbit.radiusMetres,
                velocityToleranceMetresPerSecond/targetOrbit.circularSpeedMetresPerSecond))
        let firstOffset = min(gap, 2 * .pi-gap) <= angularTolerance ? 0 : gap/abs(rate)
        let period = 2 * .pi/abs(rate)
        let span = try departureWindow.upperBound.duration(since: lower).seconds
        var previous: Epoch?
        var opportunity = 0
        while true {
            // Multiplying infinity by zero must not hide a valid first opportunity.
            let offset = opportunity == 0 ? firstOffset : firstOffset+Double(opportunity)*period
            if !offset.isFinite || offset > span { return finish(.departureWindowExhausted) }
            if attempts.count == maximumOpportunityEvaluations { return finish(.evaluationLimitReached) }
            let epoch = try lower.advanced(by: Duration(seconds: offset))
            if epoch > departureWindow.upperBound { return finish(.departureWindowExhausted) }
            guard previous == nil || epoch > previous! else { throw HohmannRendezvousError.relativeMotionNotRepresentable }
            previous = epoch
            let departureAxes = try departureOrbit.directions(at: epoch)
            let plan = try transfer.plan(departureEpoch: epoch, frame: departureOrbit.frame,
                centralBody: departureOrbit.centralBody, radialDirection: departureAxes.radial,
                progradeDirection: departureAxes.tangent, options: propagationOptions)
            let arrival = try plan.trajectory.state(at: plan.arrivalEpoch)
            let target = try targetOrbit.state(at: plan.arrivalEpoch)
            let targetDeparture = try targetOrbit.directions(at: epoch)
            let actualPhase = atan2((targetDeparture.radial*departureAxes.tangent).sum(),
                (targetDeparture.radial*departureAxes.radial).sum())
            let phaseError = atan2(sin(actualPhase-required), cos(actualPhase-required))
            let p = hohmannNorm(arrival.position-target.position), v = hohmannNorm(arrival.velocity-target.velocity)
            attempts.append(HohmannRendezvousAttempt(departureEpoch: epoch, phaseErrorRadians: phaseError,
                positionErrorMetres: p, velocityErrorMetresPerSecond: v))
            if p <= positionToleranceMetres && v <= velocityToleranceMetresPerSecond {
                return finish(.feasible, plan: plan, departure: epoch, target: target)
            }
            opportunity += 1
        }
    }
}
private func hohmannWrap(_ angle: Double) -> Double {
    let remainder = angle.truncatingRemainder(dividingBy: 2 * .pi)
    return remainder < 0 ? remainder+2 * .pi : remainder
}
private func hohmannUnit(_ v: SIMD3<Double>) throws -> SIMD3<Double> {
    let n = hohmannNorm(v)
    guard n.isFinite, n > 0 else { throw HohmannTransferError.invalidDirections }
    return v/n
}
private func hohmannCross(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> SIMD3<Double> {
    SIMD3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x)
}

private func hohmannNorm(_ v: SIMD3<Double>) -> Double { hypot(hypot(v.x,v.y),v.z) }
