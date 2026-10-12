// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum GravityModelError: Error, Equatable {
    case incompatibleGMConvention, originMismatch, identicalBodies, missingEphemerides
    case singularPosition, numericalRangeExceeded
}

private func validateGravityGM(_ gm: GravitationalParameter) throws {
    guard gm.convention == .newtonian || gm.convention == .ttCompatible else {
        throw GravityModelError.incompatibleGMConvention
    }
}
private func validateGravityState(_ state: SpacecraftState, _ context: ForceEvaluationContext) throws {
    guard state.frame === context.frame else { throw SpacecraftDynamicsError.frameMismatch }
    guard state.frame.pseudoInertial else { throw SpacecraftDynamicsError.nonInertialFrame }
}
private func gravityNorm(_ v: SIMD3<Double>) -> Double { hypot(hypot(v.x,v.y),v.z) }
private func attraction(_ displacement: SIMD3<Double>, mu: Double) throws -> SIMD3<Double> {
    let distance = gravityNorm(displacement)
    guard distance.isFinite else { throw GravityModelError.numericalRangeExceeded }
    guard distance > 0 else { throw GravityModelError.singularPosition }
    let result = (displacement/distance)*((mu/distance)/distance)
    guard [result.x,result.y,result.z].allSatisfy(\.isFinite) else { throw GravityModelError.numericalRangeExceeded }
    return result
}

/// Newtonian point-mass gravity for a body located at the propagation origin.
/// No radius/softening, collision handling or implicit body GM lookup.
public struct CentralGravity: ForceModel {
    public let body: BodyIdentifier
    public let gravitationalParameter: GravitationalParameter
    public init(body: BodyIdentifier, gravitationalParameter: GravitationalParameter) throws {
        try body.validate(); try validateGravityGM(gravitationalParameter)
        self.body = body; self.gravitationalParameter = gravitationalParameter
    }
    public func contribution(for state: SpacecraftState,in context: ForceEvaluationContext) throws -> ForceContribution {
        try validateGravityState(state,context)
        guard context.origin == body else { throw GravityModelError.originMismatch }
        return try ForceContribution(epoch: state.epoch,frame: state.frame,
            accelerationMetresPerSecondSquared: attraction(-state.cartesian.position,mu: gravitationalParameter.cubicMetresPerSecondSquared))
    }
}

/// Differential perturbation relative to a central body's nonrotating origin:
/// μ[(R-r)/|R-r|³ - R/|R|³]. The indirect term removes that body's acceleration
/// of the origin; this is not absolute attraction in a barycentric inertial frame.
public struct ThirdBodyGravity: ForceModel {
    public let body: BodyIdentifier
    public let centralBody: BodyIdentifier
    public let gravitationalParameter: GravitationalParameter
    public init(body: BodyIdentifier, centralBody: BodyIdentifier, gravitationalParameter: GravitationalParameter) throws {
        try body.validate(); try centralBody.validate(); try validateGravityGM(gravitationalParameter)
        guard body != centralBody else { throw GravityModelError.identicalBodies }
        self.body = body; self.centralBody = centralBody; self.gravitationalParameter = gravitationalParameter
    }
    public func contribution(for state: SpacecraftState,in context: ForceEvaluationContext) throws -> ForceContribution {
        try validateGravityState(state,context)
        guard context.origin == centralBody else { throw GravityModelError.originMismatch }
        guard let collection = context.ephemerides else { throw GravityModelError.missingEphemerides }
        let provider = try collection.provider(for: body)
        guard provider.metadata.origin == centralBody else { throw GravityModelError.originMismatch }
        guard provider.frame === state.frame else { throw SpacecraftDynamicsError.frameMismatch }
        let location = try provider.checkedState(at: state.epoch).position
        let radius = gravityNorm(location)
        guard radius.isFinite else { throw GravityModelError.numericalRangeExceeded }
        guard radius > 0 else { throw GravityModelError.singularPosition }
        let position = state.cartesian.position, mu = gravitationalParameter.cubicMetresPerSecondSquared
        let acceleration: SIMD3<Double>
        if position == .zero {
            acceleration = .zero
        } else if gravityNorm(position)/radius <= 0.1 {
            // Rationalize s³-1 using q = s²-1, avoiding cancellation of
            // nearly equal accelerations when the spacecraft is near the origin.
            let u = position/radius, n = location/radius, s = gravityNorm(n-u)
            let q = (u*(u-2*n)).sum()
            let f = q*(s*s+s+1)/(s+1)
            acceleration = -(u+f*n)*(((mu/radius)/radius)/(s*s*s))
        } else {
            acceleration = try attraction(location-position,mu: mu)-attraction(location,mu: mu)
        }
        guard [acceleration.x,acceleration.y,acceleration.z].allSatisfy(\.isFinite) else {
            throw GravityModelError.numericalRangeExceeded
        }
        return try ForceContribution(epoch: state.epoch,frame: state.frame,accelerationMetresPerSecondSquared: acceleration)
    }
}
