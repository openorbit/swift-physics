// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

/// Constant centre-of-mass inertia tensor (kg m²), in body axes. Off-diagonal
/// parameters are tensor entries, not unsigned products of inertia. Validation
/// requires positive definiteness and physical principal-moment triangle bounds.
public struct RigidBodyInertia {
    public let matrix: Matrix3x3
    private let scale: Double
    private let factor: [Double]
    public init(xx: Double, yy: Double, zz: Double, xy: Double = 0, xz: Double = 0, yz: Double = 0) throws {
        let values = [xx,yy,zz,xy,xz,yz]
        guard values.allSatisfy(\.isFinite) else { throw AttitudeError.invalidInertia }
        scale = values.map(abs).max()!
        guard scale > 0 else { throw AttitudeError.invalidInertia }
        let a = xx/scale, b = yy/scale, c = zz/scale, d = xy/scale, e = xz/scale, f = yz/scale
        guard a > 0 else { throw AttitudeError.invalidInertia }
        let l00 = sqrt(a), l10 = d/l00, l20 = e/l00, pivot1 = b-l10*l10
        guard pivot1 > 0 else { throw AttitudeError.invalidInertia }
        let l11 = sqrt(pivot1), l21 = (f-l20*l10)/l11, pivot2 = c-l20*l20-l21*l21
        guard pivot2 > 0 else { throw AttitudeError.invalidInertia }
        let traceHalf = (a+b+c)/2, s0 = traceHalf-a, s1 = traceHalf-b, s2 = traceHalf-c
        let tolerance = 128*Double.ulpOfOne
        let determinant = s0*s1*s2-2*d*e*f-s0*f*f-s1*e*e-s2*d*d
        guard [s0,s1,s2,s0*s1-d*d,s0*s2-e*e,s1*s2-f*f,determinant].allSatisfy({ $0 >= -tolerance }) else { throw AttitudeError.invalidInertia }
        factor = [l00,l10,l20,l11,l21,sqrt(pivot2)]
        matrix = Matrix3x3(elements: [xx,xy,xz,xy,yy,yz,xz,yz,zz])
    }
    public init(principalMomentsKilogramSquareMetres: SIMD3<Double>) throws {
        try self.init(xx: principalMomentsKilogramSquareMetres.x,yy: principalMomentsKilogramSquareMetres.y,zz: principalMomentsKilogramSquareMetres.z)
    }
    public func angularMomentumBody(for angularVelocity: SIMD3<Double>) throws -> SIMD3<Double> {
        guard attitudeFinite(angularVelocity) else { throw AttitudeError.invalidVector }
        let result = matrix*angularVelocity
        guard attitudeFinite(result) else { throw AttitudeError.numericalRangeExceeded }
        return result
    }
    public func kineticEnergyJoules(for angularVelocity: SIMD3<Double>) throws -> Double {
        let result = 0.5*(angularVelocity * (try angularMomentumBody(for: angularVelocity))).sum()
        guard result.isFinite else { throw AttitudeError.numericalRangeExceeded }
        return result
    }
    public func solve(_ torque: SIMD3<Double>) throws -> SIMD3<Double> {
        guard attitudeFinite(torque) else { throw AttitudeError.invalidVector }
        let rhs = torque/scale
        let y0 = rhs.x/factor[0], y1 = (rhs.y-factor[1]*y0)/factor[3]
        let y2 = (rhs.z-factor[2]*y0-factor[4]*y1)/factor[5]
        let z = y2/factor[5], y = (y1-factor[4]*z)/factor[3], x = (y0-factor[1]*y-factor[2]*z)/factor[0]
        let result = SIMD3(x,y,z)
        guard attitudeFinite(result) else { throw AttitudeError.numericalRangeExceeded }
        return result
    }
}
public struct RotationalEvaluationContext {
    public let inertia: RigidBodyInertia
    public let bodyFrame: Frame
    public let referenceContext: ForceEvaluationContext
    public var referenceFrame: Frame { referenceContext.frame }
    public init(inertia: RigidBodyInertia, bodyFrame: Frame, referenceContext: ForceEvaluationContext) throws {
        guard referenceContext.frame.pseudoInertial else { throw AttitudeError.nonInertialReference }
        guard bodyFrame !== referenceContext.frame else { throw AttitudeError.frameMismatch }
        self.inertia = inertia; self.bodyFrame = bodyFrame; self.referenceContext = referenceContext
    }
}
public struct TorqueContribution {
    public let epoch: Epoch
    public let bodyFrame: Frame
    public let torqueBodyNewtonMetres: SIMD3<Double>
    public init(epoch: Epoch, bodyFrame: Frame, torqueBodyNewtonMetres: SIMD3<Double>) throws {
        guard attitudeFinite(torqueBodyNewtonMetres) else { throw AttitudeError.invalidVector }
        self.epoch = epoch; self.bodyFrame = bodyFrame; self.torqueBodyNewtonMetres = torqueBodyNewtonMetres
    }
}
/// Pure, repeatable evaluations. Torque acts about the spacecraft centre of mass.
/// The optional translation is prescribed for standalone rotation, carried at
/// every RK stage for coupled dynamics. Models must declare missing inputs as errors.
public protocol TorqueModel {
    func contribution(for attitude: AttitudeState, spacecraft: SpacecraftState?,
        in context: RotationalEvaluationContext) throws -> TorqueContribution
}
func validateRotationalState(_ attitude: AttitudeState, spacecraft: SpacecraftState?, context: RotationalEvaluationContext) throws {
    guard attitude.referenceFrame === context.referenceFrame, attitude.bodyFrame === context.bodyFrame else { throw AttitudeError.frameMismatch }
    if let spacecraft {
        guard spacecraft.frame === context.referenceFrame else { throw AttitudeError.frameMismatch }
        guard spacecraft.epoch == attitude.epoch else { throw AttitudeError.epochMismatch }
    }
}
public struct RigidBodyDerivative {
    public let quaternionRate: Quaternion
    public let angularAccelerationBodyRadiansPerSecondSquared: SIMD3<Double>
}
/// Euler equation I ωdot = τ − ω × Iω; fixed inertia and no internal wheel
/// momentum. Actuator torques are external/ideal torques unless modelled as wrenches.
public struct RigidBodyDynamics {
    public let context: RotationalEvaluationContext
    public let models: [any TorqueModel]
    public init(context: RotationalEvaluationContext, models: [any TorqueModel] = []) { self.context = context; self.models = models }
    public func derivative(at attitude: AttitudeState, spacecraft: SpacecraftState? = nil,
        additionalTorqueBodyNewtonMetres: SIMD3<Double> = .zero) throws -> RigidBodyDerivative {
        try validateRotationalState(attitude,spacecraft: spacecraft,context: context)
        var torque = additionalTorqueBodyNewtonMetres
        guard attitudeFinite(torque) else { throw AttitudeError.invalidVector }
        for model in models {
            let contribution = try model.contribution(for: attitude,spacecraft: spacecraft,in: context)
            guard contribution.epoch == attitude.epoch else { throw AttitudeError.epochMismatch }
            guard contribution.bodyFrame === context.bodyFrame else { throw AttitudeError.frameMismatch }
            torque += contribution.torqueBodyNewtonMetres
            guard attitudeFinite(torque) else { throw AttitudeError.numericalRangeExceeded }
        }
        let omega = attitude.angularVelocityBodyRadiansPerSecond
        let momentum = try context.inertia.angularMomentumBody(for: omega)
        let alpha = try context.inertia.solve(torque-attitudeCross(omega,momentum))
        return RigidBodyDerivative(quaternionRate: 0.5*(attitude.bodyToReference*Quaternion(x: omega.x,y: omega.y,z: omega.z,w: 0)),
            angularAccelerationBodyRadiansPerSecondSquared: alpha)
    }
}
public struct ConstantBodyTorque: TorqueModel {
    public let bodyFrame: Frame
    public let torqueBodyNewtonMetres: SIMD3<Double>
    public init(bodyFrame: Frame, torqueBodyNewtonMetres: SIMD3<Double>) throws {
        guard attitudeFinite(torqueBodyNewtonMetres) else { throw AttitudeError.invalidVector }
        self.bodyFrame = bodyFrame; self.torqueBodyNewtonMetres = torqueBodyNewtonMetres
    }
    public func contribution(for attitude: AttitudeState, spacecraft: SpacecraftState?, in context: RotationalEvaluationContext) throws -> TorqueContribution {
        try validateRotationalState(attitude,spacecraft: spacecraft,context: context)
        guard bodyFrame === context.bodyFrame else { throw AttitudeError.frameMismatch }
        return try TorqueContribution(epoch: attitude.epoch,bodyFrame: bodyFrame,torqueBodyNewtonMetres: torqueBodyNewtonMetres)
    }
}
/// Ideal bounded external torque controller; no reaction-wheel momentum,
/// desaturation, thruster allocation, sensor delay or power model is inferred.
public struct AttitudePDController: TorqueModel {
    public let desired: any AttitudeProvider
    public let proportionalGainNewtonMetres: Double
    public let derivativeGainNewtonMetreSeconds: Double
    public let maximumTorqueNewtonMetres: Double
    public init(desired: any AttitudeProvider, proportionalGainNewtonMetres: Double,
        derivativeGainNewtonMetreSeconds: Double, maximumTorqueNewtonMetres: Double) throws {
        guard [proportionalGainNewtonMetres,derivativeGainNewtonMetreSeconds].allSatisfy({ $0.isFinite && $0 >= 0 }),
              maximumTorqueNewtonMetres.isFinite, maximumTorqueNewtonMetres > 0 else { throw AttitudeError.invalidConfiguration }
        self.desired = desired; self.proportionalGainNewtonMetres = proportionalGainNewtonMetres
        self.derivativeGainNewtonMetreSeconds = derivativeGainNewtonMetreSeconds; self.maximumTorqueNewtonMetres = maximumTorqueNewtonMetres
    }
    public func contribution(for attitude: AttitudeState, spacecraft: SpacecraftState?, in context: RotationalEvaluationContext) throws -> TorqueContribution {
        try validateRotationalState(attitude,spacecraft: spacecraft,context: context)
        guard desired.referenceFrame === context.referenceFrame, desired.bodyFrame === context.bodyFrame else { throw AttitudeError.frameMismatch }
        let goal = try desired.checkedAttitude(at: attitude.epoch)
        var error = attitude.bodyToReference.conj*goal.bodyToReference
        let tieAxis = [error.x,error.y,error.z].first { $0 != 0 } ?? 0
        if error.w < 0 || (error.w == 0 && tieAxis < 0) { error = -error }
        let goalOmega = try attitude.bodyVector(fromReference: goal.referenceVector(fromBody: goal.angularVelocityBodyRadiansPerSecond))
        var torque = 2*proportionalGainNewtonMetres*error.vector + derivativeGainNewtonMetreSeconds*(goalOmega-attitude.angularVelocityBodyRadiansPerSecond)
        let length = routeNorm(torque)
        guard length.isFinite else { throw AttitudeError.numericalRangeExceeded }
        if length > maximumTorqueNewtonMetres { torque *= maximumTorqueNewtonMetres/length }
        return try TorqueContribution(epoch: attitude.epoch,bodyFrame: context.bodyFrame,torqueBodyNewtonMetres: torque)
    }
}
/// Leading-order point-mass gravity gradient about the centre of mass:
/// τ = 3μ/r³ rHatBody × (I rHatBody). The spacecraft must be small relative to r.
public struct GravityGradientTorque: TorqueModel {
    public let centre: MovingWaypointRouteTarget
    public let gravitationalParameter: GravitationalParameter
    public init(centre: MovingWaypointRouteTarget, gravitationalParameter: GravitationalParameter) throws {
        guard gravitationalParameter.convention == .newtonian || gravitationalParameter.convention == .ttCompatible else { throw GravityModelError.incompatibleGMConvention }
        self.centre = centre; self.gravitationalParameter = gravitationalParameter
    }
    public func contribution(for attitude: AttitudeState, spacecraft: SpacecraftState?, in context: RotationalEvaluationContext) throws -> TorqueContribution {
        try validateRotationalState(attitude,spacecraft: spacecraft,context: context)
        guard centre.frame === context.referenceFrame, centre.origin == context.referenceContext.origin else { throw GravityModelError.originMismatch }
        guard let spacecraft else { throw AttitudeError.invalidConfiguration }
        let radial = spacecraft.cartesian.position - (try centre.state(at: attitude.epoch)).position, radius = routeNorm(radial)
        guard radius.isFinite, radius > 0 else { throw GravityModelError.singularPosition }
        let unit = try attitude.bodyVector(fromReference: radial/radius), inertiaUnit = context.inertia.matrix*unit
        let strength = 3*((gravitationalParameter.cubicMetresPerSecondSquared/radius)/radius)/radius
        let torque = attitudeCross(unit,inertiaUnit)*strength
        return try TorqueContribution(epoch: attitude.epoch,bodyFrame: context.bodyFrame,torqueBodyNewtonMetres: torque)
    }
}
