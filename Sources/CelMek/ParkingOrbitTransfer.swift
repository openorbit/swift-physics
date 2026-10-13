// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum ParkingOrbitTransferError: Error, Equatable {
    case invalidRadius, unsafePeriapsis, invalidExcessSpeed, invalidDirections
    case excessSpeedMismatch, zeroExcessDirection, numericalRangeExceeded
}
public enum ParkingOrbitManeuverKind: Equatable, Sendable { case escape, capture }

/// Tangential periapsis impulse between a circular parking orbit and a local
/// hyperbolic (or zero-excess parabolic) orbit. All radii are centre-relative.
public struct ParkingOrbitTransfer {
    public let gravitationalParameter: GravitationalParameter
    public let parkingRadiusMetres, bodyRadiusMetres, minimumPeriapsisAltitudeMetres: Double
    public let hyperbolicExcessSpeedMetresPerSecond: Double
    public let circularSpeedMetresPerSecond, escapeSpeedMetresPerSecond, periapsisSpeedMetresPerSecond: Double
    public let deltaVelocityMetresPerSecond, characteristicEnergySquareMetresPerSquareSecond: Double
    public let hyperbolicEccentricity: Double
    /// Acute angle between the asymptotic line and the periapsis apse line.
    public let asymptoteAngleRadians: Double
    private let radialAsymptoteComponent, tangentialAsymptoteComponent: Double

    public init(parkingRadiusMetres: Double, bodyRadiusMetres: Double,
                minimumPeriapsisAltitudeMetres: Double = 0,
                hyperbolicExcessSpeedMetresPerSecond: Double,
                gravitationalParameter: GravitationalParameter) throws {
        guard parkingRadiusMetres.isFinite, parkingRadiusMetres > 0,
              bodyRadiusMetres.isFinite, bodyRadiusMetres > 0,
              minimumPeriapsisAltitudeMetres.isFinite, minimumPeriapsisAltitudeMetres >= 0,
              (bodyRadiusMetres+minimumPeriapsisAltitudeMetres).isFinite else {
            throw ParkingOrbitTransferError.invalidRadius
        }
        guard parkingRadiusMetres > bodyRadiusMetres,
              parkingRadiusMetres >= bodyRadiusMetres+minimumPeriapsisAltitudeMetres else {
            throw ParkingOrbitTransferError.unsafePeriapsis
        }
        let excess = hyperbolicExcessSpeedMetresPerSecond
        guard excess.isFinite, excess >= 0 else { throw ParkingOrbitTransferError.invalidExcessSpeed }
        guard gravitationalParameter.convention == .newtonian || gravitationalParameter.convention == .ttCompatible else {
            throw PropagationError.incompatibleGMConvention
        }
        let circular = sqrt(gravitationalParameter.cubicMetresPerSecondSquared)/sqrt(parkingRadiusMetres)
        let escape = sqrt(2)*circular, periapsis = hypot(excess, escape)
        let ratio = excess/circular, x = ratio*ratio, eccentricity = 1+x
        let radial = 1/eccentricity
        // Avoid cancellation in 1-1/e² near the parabolic limit and e² overflow.
        let tangent = sqrt(x)/eccentricity*sqrt(2+x)
        let c3 = excess*excess, delta = periapsis-circular
        guard [circular, escape, periapsis, eccentricity, radial, tangent, c3, delta].allSatisfy(\.isFinite),
              circular > 0, escape > 0, periapsis > 0, delta > 0 else {
            throw ParkingOrbitTransferError.numericalRangeExceeded
        }
        self.gravitationalParameter = gravitationalParameter
        self.parkingRadiusMetres = parkingRadiusMetres; self.bodyRadiusMetres = bodyRadiusMetres
        self.minimumPeriapsisAltitudeMetres = minimumPeriapsisAltitudeMetres
        self.hyperbolicExcessSpeedMetresPerSecond = excess
        circularSpeedMetresPerSecond = circular; escapeSpeedMetresPerSecond = escape
        periapsisSpeedMetresPerSecond = periapsis; deltaVelocityMetresPerSecond = delta
        characteristicEnergySquareMetresPerSquareSecond = c3; hyperbolicEccentricity = eccentricity
        radialAsymptoteComponent = radial; tangentialAsymptoteComponent = tangent
        asymptoteAngleRadians = atan2(tangent, radial)
    }

    /// The excess vector is the actual incoming/outgoing velocity at infinity in
    /// body-centred native axes. The normal fixes the right-handed motion plane.
    /// A zero-excess scalar calculation has no vector-defined asymptote plan.
    public func plan(kind: ParkingOrbitManeuverKind, periapsisEpoch: Epoch, frame: Frame,
        centralBody: BodyIdentifier, hyperbolicExcessVelocityMetresPerSecond: SIMD3<Double>,
        orbitalNormal: SIMD3<Double>, relativeSpeedTolerance: Double = 1e-12,
        propagationOptions: TwoBodyPropagationOptions = try! TwoBodyPropagationOptions()) throws -> ParkingOrbitManeuverPlan {
        guard frame.pseudoInertial else { throw PropagationError.nonInertialFrame }
        guard hyperbolicExcessSpeedMetresPerSecond > 0, tangentialAsymptoteComponent > 0 else {
            throw ParkingOrbitTransferError.zeroExcessDirection
        }
        let magnitude = parkingNorm(hyperbolicExcessVelocityMetresPerSecond)
        guard magnitude.isFinite, magnitude > 0, relativeSpeedTolerance.isFinite,
              relativeSpeedTolerance >= 0, relativeSpeedTolerance < 1 else { throw ParkingOrbitTransferError.invalidDirections }
        guard abs(magnitude/hyperbolicExcessSpeedMetresPerSecond-1) <= relativeSpeedTolerance else {
            throw ParkingOrbitTransferError.excessSpeedMismatch
        }
        let direction = hyperbolicExcessVelocityMetresPerSecond/magnitude
        let normal = try parkingUnit(orbitalNormal)
        guard abs((direction*normal).sum()) <= 1e-12 else { throw ParkingOrbitTransferError.invalidDirections }
        let perpendicular = try parkingUnit(parkingCross(normal, direction))
        let c = radialAsymptoteComponent, s = tangentialAsymptoteComponent
        let radial = kind == .escape ? -c*direction-s*perpendicular : c*direction-s*perpendicular
        let tangent = kind == .escape ? s*direction-c*perpendicular : s*direction+c*perpendicular
        let position = radial*parkingRadiusMetres
        let beforeSpeed = kind == .escape ? circularSpeedMetresPerSecond : periapsisSpeedMetresPerSecond
        let afterSpeed = kind == .escape ? periapsisSpeedMetresPerSecond : circularSpeedMetresPerSecond
        let before = try CartesianState(epoch: periapsisEpoch, frame: frame, position: position,
            velocity: tangent*beforeSpeed)
        let after = try CartesianState(epoch: periapsisEpoch, frame: frame, position: position,
            velocity: tangent*afterSpeed)
        let impulse = try ImpulsiveManeuver(epoch: periapsisEpoch, frame: frame,
            deltaVelocityMetresPerSecond: tangent*(kind == .escape ? deltaVelocityMetresPerSecond : -deltaVelocityMetresPerSecond))
        let base = try TwoBodyPropagator(initialState: before, centralBody: centralBody,
            gravitationalParameter: gravitationalParameter, options: propagationOptions)
        let trajectory = try ManeuveredTwoBodyPropagator(base: base, maneuvers: [impulse])
        return ParkingOrbitManeuverPlan(kind: kind, transfer: self,
            hyperbolicExcessVelocityMetresPerSecond: hyperbolicExcessVelocityMetresPerSecond,
            periapsisRadialDirection: radial, periapsisProgradeDirection: tangent,
            preBurnState: before, postBurnState: after, maneuver: impulse, trajectory: trajectory)
    }
}
public struct ParkingOrbitManeuverPlan {
    public let kind: ParkingOrbitManeuverKind
    public let transfer: ParkingOrbitTransfer
    public let hyperbolicExcessVelocityMetresPerSecond: SIMD3<Double>
    public let periapsisRadialDirection, periapsisProgradeDirection: SIMD3<Double>
    public let preBurnState, postBurnState: CartesianState
    public let maneuver: ImpulsiveManeuver
    public let trajectory: ManeuveredTwoBodyPropagator
}
private func parkingNorm(_ v: SIMD3<Double>) -> Double { hypot(hypot(v.x,v.y),v.z) }
private func parkingUnit(_ v: SIMD3<Double>) throws -> SIMD3<Double> {
    let n = parkingNorm(v)
    guard n.isFinite, n > 0 else { throw ParkingOrbitTransferError.invalidDirections }
    return v/n
}
private func parkingCross(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> SIMD3<Double> {
    SIMD3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x)
}
