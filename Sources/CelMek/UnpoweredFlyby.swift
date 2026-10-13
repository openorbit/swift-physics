// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum FlybyError: Error, Equatable {
    case invalidBodyRadius, invalidAltitude, invalidVelocity, invalidTolerance
    case excessSpeedMismatch, degenerateTurn, unsafePeriapsis, numericalRangeExceeded
}
/// Incoming/outgoing vectors are actual velocity directions at infinity, not
/// position-asymptote directions. The caller guarantees body-centred inertial axes.
public struct FlybyEncounterRequest {
    public let body: BodyIdentifier
    /// Interpreted as the epoch of closest approach in the local two-body model.
    public let periapsisEpoch: Epoch
    public let frame: Frame
    public let gravitationalParameter: GravitationalParameter
    public let bodyRadiusMetres: Double
    public let minimumPeriapsisAltitudeMetres: Double
    public let incomingExcessVelocityMetresPerSecond: SIMD3<Double>
    public let requestedOutgoingExcessVelocityMetresPerSecond: SIMD3<Double>
    public let relativeSpeedTolerance: Double
    public init(body: BodyIdentifier, periapsisEpoch: Epoch, frame: Frame,
        gravitationalParameter: GravitationalParameter, bodyRadiusMetres: Double,
        minimumPeriapsisAltitudeMetres: Double = 0,
        incomingExcessVelocityMetresPerSecond: SIMD3<Double>,
        outgoingExcessVelocityMetresPerSecond: SIMD3<Double>, relativeSpeedTolerance: Double = 1e-12) throws {
        try body.validate()
        guard frame.pseudoInertial else { throw PropagationError.nonInertialFrame }
        guard gravitationalParameter.convention == .newtonian || gravitationalParameter.convention == .ttCompatible else {
            throw PropagationError.incompatibleGMConvention
        }
        guard bodyRadiusMetres.isFinite, bodyRadiusMetres > 0 else { throw FlybyError.invalidBodyRadius }
        guard minimumPeriapsisAltitudeMetres.isFinite, minimumPeriapsisAltitudeMetres >= 0,
              (bodyRadiusMetres+minimumPeriapsisAltitudeMetres).isFinite else { throw FlybyError.invalidAltitude }
        guard relativeSpeedTolerance.isFinite, relativeSpeedTolerance >= 0, relativeSpeedTolerance < 1 else { throw FlybyError.invalidTolerance }
        for velocity in [incomingExcessVelocityMetresPerSecond,outgoingExcessVelocityMetresPerSecond] {
            guard [velocity.x,velocity.y,velocity.z].allSatisfy(\.isFinite), routeNorm(velocity).isFinite,
                  routeNorm(velocity) > 0 else { throw FlybyError.invalidVelocity }
        }
        self.body = body; self.periapsisEpoch = periapsisEpoch; self.frame = frame
        self.gravitationalParameter = gravitationalParameter; self.bodyRadiusMetres = bodyRadiusMetres
        self.minimumPeriapsisAltitudeMetres = minimumPeriapsisAltitudeMetres
        self.incomingExcessVelocityMetresPerSecond = incomingExcessVelocityMetresPerSecond
        requestedOutgoingExcessVelocityMetresPerSecond = outgoingExcessVelocityMetresPerSecond
        self.relativeSpeedTolerance = relativeSpeedTolerance
    }
}
public struct UnpoweredFlybyResult {
    public let request: FlybyEncounterRequest
    public let excessSpeedMetresPerSecond: Double
    public let outgoingExcessVelocityMetresPerSecond: SIMD3<Double>
    /// Residual between the speed-conserving model and the original requested vector.
    public let outgoingVelocityResidualMetresPerSecond: Double
    public let requestedOutgoingSpeedDifferenceMetresPerSecond: Double
    public let turningAngleRadians: Double
    public let periapsisRadiusMetres: Double
    public let periapsisAltitudeMetres: Double
    public let eccentricity: Double
    public let periapsisSpeedMetresPerSecond: Double
    public let orbitalNormal: SIMD3<Double>
    public let periapsisRadialDirection: SIMD3<Double>
    public let periapsisProgradeDirection: SIMD3<Double>
    public let periapsisState: CartesianState
    public let trajectory: TwoBodyPropagator
}

/// Ideal spherical point-mass patched-conic geometry, with no impulse/fuel cost.
/// Does not patch a heliocentric arc or certify atmospheric/terrain clearance.
public enum UnpoweredFlyby {
    public static func solve(_ request: FlybyEncounterRequest,
        propagationOptions: TwoBodyPropagationOptions = .default) throws -> UnpoweredFlybyResult {
        let incoming = request.incomingExcessVelocityMetresPerSecond
        let requested = request.requestedOutgoingExcessVelocityMetresPerSecond
        let speed = routeNorm(incoming), requestedSpeed = routeNorm(requested)
        guard abs(requestedSpeed/speed-1) <= request.relativeSpeedTolerance else { throw FlybyError.excessSpeedMismatch }
        let u = incoming/speed, w = requested/requestedSpeed
        let difference = u-w, sum = u+w, dn = routeNorm(difference), sn = routeNorm(sum)
        // No finite flyby for zero turn; opposite vectors have no unique plane
        // and would require a zero-radius collision. Never choose an arbitrary plane.
        guard dn > 1e-12, sn > 1e-12 else { throw FlybyError.degenerateTurn }
        func cross(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> SIMD3<Double> {
            SIMD3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x)
        }
        let plane = cross(u,w), normal = plane/routeNorm(plane)
        let radial: SIMD3<Double>, tangent: SIMD3<Double>
        // Use the better-conditioned bisector, deriving its orthogonal partner
        // rather than retaining subtraction noise as a radial periapsis velocity.
        if dn < sn {
            tangent = sum/sn
            let direction = cross(tangent,normal); radial = direction/routeNorm(direction)
        } else {
            radial = difference/dn
            let direction = cross(normal,radial); tangent = direction/routeNorm(direction)
        }
        let halfNorm = hypot(dn,sn), c = dn/halfNorm, s = sn/halfNorm
        let turn = 2*atan2(c,s), length = sqrt(request.gravitationalParameter.cubicMetresPerSecondSquared)/speed
        // (1/sin(turn/2)-1), rationalized near a pi turn.
        let radius = (length*s)*(length*s)/(c*(1+c))
        guard radius.isFinite, radius > 0 else { throw FlybyError.numericalRangeExceeded }
        guard radius > request.bodyRadiusMetres,
              radius >= request.bodyRadiusMetres+request.minimumPeriapsisAltitudeMetres else { throw FlybyError.unsafePeriapsis }
        let periSpeed = hypot(speed,sqrt(2)*sqrt(request.gravitationalParameter.cubicMetresPerSecondSquared)/sqrt(radius))
        let outgoing = w*speed, eccentricity = 1/c
        let residual = routeNorm(outgoing-requested)
        guard [periSpeed,eccentricity,residual].allSatisfy(\.isFinite), periSpeed > speed,
              eccentricity > 1 else { throw FlybyError.numericalRangeExceeded }
        let seed = try CartesianState(epoch: request.periapsisEpoch,frame: request.frame,
            position: radial*radius,velocity: tangent*periSpeed)
        let trajectory = try TwoBodyPropagator(initialState: seed,centralBody: request.body,
            gravitationalParameter: request.gravitationalParameter,options: propagationOptions)
        return UnpoweredFlybyResult(request: request,excessSpeedMetresPerSecond: speed,
            outgoingExcessVelocityMetresPerSecond: outgoing,outgoingVelocityResidualMetresPerSecond: residual,
            requestedOutgoingSpeedDifferenceMetresPerSecond: requestedSpeed-speed,turningAngleRadians: turn,
            periapsisRadiusMetres: radius,periapsisAltitudeMetres: radius-request.bodyRadiusMetres,
            eccentricity: eccentricity,periapsisSpeedMetresPerSecond: periSpeed,orbitalNormal: normal,
            periapsisRadialDirection: radial,periapsisProgradeDirection: tangent,
            periapsisState: try trajectory.state(at: request.periapsisEpoch),trajectory: trajectory)
    }
}
