// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum J2GravityError: Error, Equatable {
    case invalidCoefficient
}

/// Axisymmetric degree-two gravity correction, excluding central attraction.
/// J2 is dimensionless and equals minus the unnormalized C20 coefficient.
/// The orientation maps propagation axes into body axes, whose z axis is the pole.
public struct J2Gravity: ForceModel {
    public let body: BodyIdentifier
    public let gravitationalParameter: GravitationalParameter
    public let equatorialRadius: SourcedPhysicalValue
    public let j2: Double
    public let coefficientProvenance: PhysicalDataProvenance
    public let orientation: any OrientationProvider

    public init(body: BodyIdentifier, gravitationalParameter: GravitationalParameter,
                equatorialRadius: SourcedPhysicalValue, j2: Double,
                coefficientProvenance: PhysicalDataProvenance,
                orientation: any OrientationProvider) throws {
        try body.validate()
        guard gravitationalParameter.convention == .newtonian || gravitationalParameter.convention == .ttCompatible else {
            throw GravityModelError.incompatibleGMConvention
        }
        guard j2.isFinite else { throw J2GravityError.invalidCoefficient }
        self.body = body
        self.gravitationalParameter = gravitationalParameter
        self.equatorialRadius = equatorialRadius
        self.j2 = j2
        self.coefficientProvenance = coefficientProvenance
        self.orientation = orientation
    }

    public func contribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
        guard state.frame === context.frame, orientation.referenceFrame === state.frame else {
            throw SpacecraftDynamicsError.frameMismatch
        }
        guard state.frame.pseudoInertial else { throw SpacecraftDynamicsError.nonInertialFrame }
        guard context.origin == body else { throw GravityModelError.originMismatch }
        let axes = try orientation.orientation(at: state.epoch)
        guard axes.referenceFrame === state.frame else { throw SpacecraftDynamicsError.frameMismatch }
        guard axes.epoch == state.epoch else { throw SpacecraftDynamicsError.epochMismatch }
        let p = axes.rotation.matrix * state.cartesian.position
        let r = hypot(hypot(p.x, p.y), p.z)
        guard r.isFinite else { throw GravityModelError.numericalRangeExceeded }
        guard r > 0 else { throw GravityModelError.singularPosition }
        let acceleration: SIMD3<Double>
        if j2 == 0 {
            acceleration = .zero
        } else {
            let n = p / r
            let radiusRatio = equatorialRadius.value / r
            let scale = 1.5 * j2 * ((gravitationalParameter.cubicMetresPerSecondSquared / r) / r)
                * radiusRatio * radiusRatio
            let zTerm = 5 * n.z * n.z
            let bodyAcceleration = scale * SIMD3(n.x * (zTerm - 1), n.y * (zTerm - 1), n.z * (zTerm - 3))
            // Transform a physical force vector only: no rotating-frame kinematic terms.
            acceleration = axes.rotation.invMatrix * bodyAcceleration
        }
        guard [acceleration.x, acceleration.y, acceleration.z].allSatisfy(\.isFinite) else {
            throw GravityModelError.numericalRangeExceeded
        }
        return try ForceContribution(epoch: state.epoch, frame: state.frame,
            accelerationMetresPerSecondSquared: acceleration)
    }
}
