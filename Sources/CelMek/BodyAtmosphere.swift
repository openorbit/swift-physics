// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Atmospheres

public extension BodyIdentifier {
    /// Nominal built-in atmosphere for a body centre; nil means no supplied model.
    var atmosphericModel: BodyAtmosphere? {
        guard case .naif(let code) = self else { return nil }
        return BodyAtmosphere.model(naifBodyCode: code)
    }
}

/// Spherical altitude adapter for drag. Caller supplies the radius of the profile's
/// datum and angular velocity in the state's inertial axes. Atmospheric co-rotation
/// is included as Ω×r; no additional meteorological winds are inferred.
/// Origin must be the body's centre; no implicit ephemeris or orientation transform.
public struct SphericalBodyAtmosphere: AtmosphereProvider {
    public let bodyIdentifier: BodyIdentifier
    public let profile: any AtmosphericProfile
    public let referenceRadiusMetres: Double
    public let angularVelocityRadiansPerSecond: SIMD3<Double>
    public let frame: Frame
    public init(bodyIdentifier: BodyIdentifier, profile: any AtmosphericProfile,
                referenceRadiusMetres: Double, frame: Frame,
                angularVelocityRadiansPerSecond: SIMD3<Double> = .zero) throws {
        try bodyIdentifier.validate()
        guard referenceRadiusMetres.isFinite, referenceRadiusMetres > 0,
              [angularVelocityRadiansPerSecond.x, angularVelocityRadiansPerSecond.y,
               angularVelocityRadiansPerSecond.z].allSatisfy(\.isFinite) else {
            throw EnvironmentalForceError.invalidParameter
        }
        if let identified = profile as? BodyAtmosphere,
           bodyIdentifier != .naif(identified.body.rawValue) { throw BodyModelError.bodyMismatch }
        guard frame.pseudoInertial else { throw SpacecraftDynamicsError.nonInertialFrame }
        self.bodyIdentifier = bodyIdentifier; self.profile = profile; self.referenceRadiusMetres = referenceRadiusMetres
        self.frame = frame; self.angularVelocityRadiansPerSecond = angularVelocityRadiansPerSecond
    }
    public func sample(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> AtmosphereSample {
        guard context.origin == bodyIdentifier else { throw EnvironmentalForceError.originMismatch }
        guard state.frame === frame, context.frame === frame else { throw SpacecraftDynamicsError.frameMismatch }
        let p = state.cartesian.position, w = angularVelocityRadiansPerSecond
        let altitude = hypot(hypot(p.x,p.y),p.z)-referenceRadiusMetres
        let properties = try profile.properties(at: altitude)
        return try AtmosphereSample(epoch: state.epoch, frame: frame, origin: bodyIdentifier,
            densityKilogramsPerCubicMetre: properties.densityKilogramsPerCubicMetre,
            velocityMetresPerSecond: SIMD3(w.y*p.z-w.z*p.y, w.z*p.x-w.x*p.z, w.x*p.y-w.y*p.x))
    }
}
