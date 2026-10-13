// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum BiEllipticTransferError: Error, Equatable {
    case invalidIntermediateRadius, numericalRangeExceeded, invalidComparison
}

/// Three-impulse transfer between circular orbits via a caller-selected apoapsis.
/// An optional plane change is combined with the impulse at that apoapsis.
public struct BiEllipticTransfer {
    public let departureRadiusMetres, arrivalRadiusMetres, intermediateRadiusMetres: Double
    public let gravitationalParameter: GravitationalParameter
    public let firstSemiMajorAxisMetres, secondSemiMajorAxisMetres: Double
    public let firstLegDuration, secondLegDuration, duration: Duration
    public let departureDeltaVelocityMetresPerSecond, arrivalDeltaVelocityMetresPerSecond: Double
    /// Signed coplanar speed difference; the combined middle impulse may be larger.
    public let intermediateSpeedChangeMetresPerSecond: Double
    public let intermediateImpulse: PlaneChangeImpulse
    public let intermediateDeltaVelocityMetresPerSecond: Double
    public let totalDeltaVelocityMetresPerSecond: Double
    private let departureCircularSpeed, arrivalCircularSpeed: Double

    public init(departureRadiusMetres: Double, arrivalRadiusMetres: Double, intermediateRadiusMetres: Double,
                gravitationalParameter: GravitationalParameter, planeChangeRadians: Double = 0) throws {
        guard intermediateRadiusMetres.isFinite, intermediateRadiusMetres > max(departureRadiusMetres, arrivalRadiusMetres) else {
            throw BiEllipticTransferError.invalidIntermediateRadius
        }
        let first = try HohmannTransfer(departureRadiusMetres: departureRadiusMetres,
            arrivalRadiusMetres: intermediateRadiusMetres, gravitationalParameter: gravitationalParameter)
        let second = try HohmannTransfer(departureRadiusMetres: intermediateRadiusMetres,
            arrivalRadiusMetres: arrivalRadiusMetres, gravitationalParameter: gravitationalParameter)
        let x = departureRadiusMetres/intermediateRadiusMetres, y = arrivalRadiusMetres/intermediateRadiusMetres
        let apoCircular = first.arrivalCircularSpeedMetresPerSecond
        let before = apoCircular*sqrt(2*x/(1+x)), after = apoCircular*sqrt(2*y/(1+y))
        let speedChange = apoCircular*(2*((arrivalRadiusMetres-departureRadiusMetres)/intermediateRadiusMetres) /
            ((1+x)*(1+y)))/(sqrt(2*x/(1+x))+sqrt(2*y/(1+y)))
        let impulse = try PlaneChangeImpulse(beforeSpeedMetresPerSecond: before,
            afterSpeedMetresPerSecond: after, angleRadians: planeChangeRadians)
        // Keep the rationalized coplanar difference for nearly equal endpoint radii.
        let middleCost = planeChangeRadians == 0 ? abs(speedChange) : impulse.deltaVelocityMetresPerSecond
        let total = abs(first.departureDeltaVelocityMetresPerSecond)+middleCost+abs(second.arrivalDeltaVelocityMetresPerSecond)
        guard [before, after, speedChange, total].allSatisfy(\.isFinite), before > 0, after > 0 else {
            throw BiEllipticTransferError.numericalRangeExceeded
        }
        self.departureRadiusMetres = departureRadiusMetres; self.arrivalRadiusMetres = arrivalRadiusMetres
        self.intermediateRadiusMetres = intermediateRadiusMetres; self.gravitationalParameter = gravitationalParameter
        firstSemiMajorAxisMetres = first.semiMajorAxisMetres; secondSemiMajorAxisMetres = second.semiMajorAxisMetres
        firstLegDuration = first.duration; secondLegDuration = second.duration
        duration = try first.duration.adding(second.duration)
        departureDeltaVelocityMetresPerSecond = first.departureDeltaVelocityMetresPerSecond
        arrivalDeltaVelocityMetresPerSecond = second.arrivalDeltaVelocityMetresPerSecond
        intermediateSpeedChangeMetresPerSecond = speedChange; intermediateImpulse = impulse
        totalDeltaVelocityMetresPerSecond = total; intermediateDeltaVelocityMetresPerSecond = middleCost
        departureCircularSpeed = first.departureCircularSpeedMetresPerSecond
        arrivalCircularSpeed = second.arrivalCircularSpeedMetresPerSecond
    }

    public func plan(departureEpoch: Epoch, frame: Frame, centralBody: BodyIdentifier,
        radialDirection: SIMD3<Double>, progradeDirection: SIMD3<Double>, planeChangeSense: PlaneChangeSense = .positive,
        propagationOptions: TwoBodyPropagationOptions = try! TwoBodyPropagationOptions()) throws -> BiEllipticTransferPlan {
        guard frame.pseudoInertial else { throw PropagationError.nonInertialFrame }
        let radial = try biUnit(radialDirection), supplied = try biUnit(progradeDirection)
        let dot = (radial*supplied).sum()
        guard abs(dot) <= 1e-12 else { throw HohmannTransferError.invalidDirections }
        let tangent = try biUnit(supplied-radial*dot)
        let intermediateEpoch = try departureEpoch.advanced(by: firstLegDuration)
        let arrivalEpoch = try intermediateEpoch.advanced(by: secondLegDuration)
        guard intermediateEpoch > departureEpoch, arrivalEpoch > intermediateEpoch else {
            throw BiEllipticTransferError.numericalRangeExceeded
        }
        let initial = try CartesianState(epoch: departureEpoch, frame: frame,
            position: radial*departureRadiusMetres, velocity: tangent*departureCircularSpeed)
        let base = try TwoBodyPropagator(initialState: initial, centralBody: centralBody,
            gravitationalParameter: gravitationalParameter, options: propagationOptions)
        let departure = try ImpulsiveManeuver(epoch: departureEpoch, frame: frame,
            deltaVelocityMetresPerSecond: tangent*departureDeltaVelocityMetresPerSecond)
        let middle: ImpulsiveManeuver
        if intermediateImpulse.angleRadians == 0 {
            middle = try ImpulsiveManeuver(epoch: intermediateEpoch, frame: frame,
                deltaVelocityMetresPerSecond: -tangent*intermediateSpeedChangeMetresPerSecond)
        } else {
            middle = try intermediateImpulse.maneuver(epoch: intermediateEpoch, frame: frame,
                radialDirection: -radial, progradeDirection: -tangent, sense: planeChangeSense)
        }
        let finalTangent = Quaternion(angle: planeChangeSense == .positive ?
            intermediateImpulse.angleRadians : -intermediateImpulse.angleRadians, axis: -radial).matrix*tangent
        let arrival = try ImpulsiveManeuver(epoch: arrivalEpoch, frame: frame,
            deltaVelocityMetresPerSecond: finalTangent*arrivalDeltaVelocityMetresPerSecond)
        let trajectory = try ManeuveredTwoBodyPropagator(base: base, maneuvers: [departure,middle,arrival])
        let expected = try CartesianState(epoch: arrivalEpoch, frame: frame,
            position: radial*arrivalRadiusMetres, velocity: finalTangent*arrivalCircularSpeed)
        return BiEllipticTransferPlan(transfer: self, intermediateEpoch: intermediateEpoch, arrivalEpoch: arrivalEpoch,
            departureBurn: departure, intermediateBurn: middle, arrivalBurn: arrival,
            expectedArrivalState: expected, trajectory: trajectory)
    }
}
public struct BiEllipticTransferPlan {
    public let transfer: BiEllipticTransfer
    public let intermediateEpoch, arrivalEpoch: Epoch
    public let departureBurn, intermediateBurn, arrivalBurn: ImpulsiveManeuver
    public let expectedArrivalState: CartesianState
    public let trajectory: ManeuveredTwoBodyPropagator
}
private func biUnit(_ v: SIMD3<Double>) throws -> SIMD3<Double> {
    let n = hypot(hypot(v.x,v.y),v.z)
    guard n.isFinite, n > 0 else { throw HohmannTransferError.invalidDirections }
    return v/n
}
