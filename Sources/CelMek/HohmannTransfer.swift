// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum HohmannTransferError: Error, Equatable {
    case invalidRadius, invalidDirections, numericalRangeExceeded
}

/// Newtonian two-impulse transfer between coplanar circular orbits about one
/// central body. Radii are measured from its centre, not surface altitudes.
public struct HohmannTransfer {
    public let departureRadiusMetres: Double
    public let arrivalRadiusMetres: Double
    public let gravitationalParameter: GravitationalParameter
    public let semiMajorAxisMetres: Double
    public let eccentricity: Double
    public let departureCircularSpeedMetresPerSecond: Double
    public let arrivalCircularSpeedMetresPerSecond: Double
    /// Signed along the local direction of orbital motion at each burn.
    public let departureDeltaVelocityMetresPerSecond: Double
    public let arrivalDeltaVelocityMetresPerSecond: Double
    public let totalDeltaVelocityMetresPerSecond: Double
    public let duration: Duration

    public init(departureRadiusMetres r1: Double, arrivalRadiusMetres r2: Double,
                gravitationalParameter gm: GravitationalParameter) throws {
        guard r1.isFinite, r2.isFinite, r1 > 0, r2 > 0 else { throw HohmannTransferError.invalidRadius }
        guard gm.convention == .newtonian || gm.convention == .ttCompatible else {
            throw PropagationError.incompatibleGMConvention
        }
        let scale = max(r1,r2), x = r1/scale, y = r2/scale, sum = x+y
        let a = r1/2+r2/2
        let v1 = sqrt(gm.cubicMetresPerSecondSquared)/sqrt(r1)
        let v2 = sqrt(gm.cubicMetresPerSecondSquared)/sqrt(r2)
        // Rationalized speed differences preserve small radius changes.
        let difference = (r2-r1)/scale
        let dv1 = v1 * (difference/sum)/(sqrt(2*y/sum)+1)
        let dv2 = v2 * (difference/sum)/(sqrt(2*x/sum)+1)
        let seconds = Double.pi * (a/(sqrt(gm.cubicMetresPerSecondSquared)/sqrt(a)))
        let total = abs(dv1)+abs(dv2), e = abs(difference)/sum
        guard [a,v1,v2,dv1,dv2,seconds,total,e].allSatisfy(\.isFinite),
              a > 0, v1 > 0, v2 > 0, seconds > 0, e < 1 else {
            throw HohmannTransferError.numericalRangeExceeded
        }
        departureRadiusMetres = r1; arrivalRadiusMetres = r2; gravitationalParameter = gm
        semiMajorAxisMetres = a; eccentricity = e
        departureCircularSpeedMetresPerSecond = v1; arrivalCircularSpeedMetresPerSecond = v2
        departureDeltaVelocityMetresPerSecond = dv1; arrivalDeltaVelocityMetresPerSecond = dv2
        totalDeltaVelocityMetresPerSecond = total; duration = try Duration(seconds: seconds)
    }

    /// Construct the initial circular state and both burns in an explicit inertial
    /// frame centred on centralBody. Directions are normalized and must be perpendicular.
    /// This reaches the destination orbit; it does not solve target rendezvous phasing.
    public func plan(departureEpoch: Epoch, frame: Frame, centralBody: BodyIdentifier,
                     radialDirection: SIMD3<Double>, progradeDirection: SIMD3<Double>,
                     options: TwoBodyPropagationOptions = try! TwoBodyPropagationOptions()) throws -> HohmannTransferPlan {
        guard frame.pseudoInertial else { throw PropagationError.nonInertialFrame }
        func unit(_ v: SIMD3<Double>) throws -> SIMD3<Double> {
            let n = hypot(hypot(v.x,v.y),v.z)
            guard n.isFinite, n > 0 else { throw HohmannTransferError.invalidDirections }
            return v/n
        }
        let radial = try unit(radialDirection), tangent = try unit(progradeDirection)
        guard abs((radial*tangent).sum()) <= 1e-12 else { throw HohmannTransferError.invalidDirections }
        let initial = try CartesianState(epoch: departureEpoch, frame: frame,
            position: radial*departureRadiusMetres, velocity: tangent*departureCircularSpeedMetresPerSecond)
        let base = try TwoBodyPropagator(initialState: initial, centralBody: centralBody,
            gravitationalParameter: gravitationalParameter, options: options)
        let arrival = try departureEpoch.advanced(by: duration)
        guard arrival > departureEpoch else { throw HohmannTransferError.numericalRangeExceeded }
        let departureBurn = try ImpulsiveManeuver(epoch: departureEpoch, frame: frame,
            deltaVelocityMetresPerSecond: tangent*departureDeltaVelocityMetresPerSecond)
        // Arrival is diametrically opposite departure; local prograde is reversed.
        let arrivalBurn = try ImpulsiveManeuver(epoch: arrival, frame: frame,
            deltaVelocityMetresPerSecond: -tangent*arrivalDeltaVelocityMetresPerSecond)
        let trajectory = try ManeuveredTwoBodyPropagator(base: base, maneuvers: [departureBurn,arrivalBurn])
        return HohmannTransferPlan(transfer: self, arrivalEpoch: arrival, departureBurn: departureBurn,
            arrivalBurn: arrivalBurn, trajectory: trajectory)
    }
}

public struct HohmannTransferPlan {
    public let transfer: HohmannTransfer
    public let arrivalEpoch: Epoch
    public let departureBurn: ImpulsiveManeuver
    public let arrivalBurn: ImpulsiveManeuver
    /// Query before/after burns or sample with TrajectorySampler.
    public let trajectory: ManeuveredTwoBodyPropagator
}
