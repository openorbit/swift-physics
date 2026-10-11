// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum OrbitalElementsError: Error, Equatable {
    case invalidElements, invalidOptions, nonInertialFrame
    case zeroPosition, radialTrajectory, numericalRangeExceeded, anomalyOutsideConic
}

/// Thresholds for canonicalizing effectively circular/equatorial input states.
/// Canonicalization may change the reconstructed state by these small amounts.
public struct OrbitalElementConversionOptions: Equatable, Sendable {
    public let circularEccentricityTolerance: Double
    public let equatorialSineTolerance: Double
    public init(circularEccentricityTolerance: Double = 1e-12,
                equatorialSineTolerance: Double = 1e-12) throws {
        guard circularEccentricityTolerance.isFinite, circularEccentricityTolerance >= 0,
              circularEccentricityTolerance < 1,
              equatorialSineTolerance.isFinite, equatorialSineTolerance >= 0,
              equatorialSineTolerance < 1 else { throw OrbitalElementsError.invalidOptions }
        self.circularEccentricityTolerance = circularEccentricityTolerance
        self.equatorialSineTolerance = equatorialSineTolerance
    }
    public static let `default` = try! OrbitalElementConversionOptions()
}

public enum OrbitalAngleConvention: Equatable, Sendable {
    case classical
    /// Argument of periapsis is zero; trueAnomaly is argument of latitude.
    case circularInclined
    /// Node is zero; periapsis is measured from +X in the direction of motion.
    case eccentricEquatorial
    /// Node and periapsis are zero; trueAnomaly is longitude in direction of motion.
    case circularEquatorial
}

public enum OrbitalConic: Equatable, Sendable { case elliptic, parabolic, hyperbolic }

/// Instantaneous two-body osculating elements, in right-handed inertial axes
/// centred on centralBody. Lengths are metres and all angles are radians.
/// These are independent of the legacy Meeus mean planetary elements.
public struct OsculatingOrbitalElements {
    public let epoch: Epoch
    public let frame: Frame
    public let centralBody: BodyIdentifier
    public let gravitationalParameter: GravitationalParameter
    public let semiLatusRectumMetres: Double
    public let eccentricity: Double
    public let inclination: Double
    public let longitudeOfAscendingNode: Double
    public let argumentOfPeriapsis: Double
    /// Signed principal angle [-pi, pi]; for unbound orbits this selects the
    /// physically reachable branch, subject to 1 + e cos(nu) > 0.
    public let trueAnomaly: Double
    public let angleConvention: OrbitalAngleConvention

    public var conic: OrbitalConic {
        eccentricity < 1 ? .elliptic : eccentricity == 1 ? .parabolic : .hyperbolic
    }
    /// Nil for a parabola; negative for a hyperbola. Near unity, p is better conditioned.
    public var semiMajorAxisMetres: Double? {
        guard eccentricity != 1 else { return nil }
        return semiLatusRectumMetres / ((1 - eccentricity) * (1 + eccentricity))
    }
    public var periapsisDistanceMetres: Double { semiLatusRectumMetres / (1 + eccentricity) }

    public init(epoch: Epoch, frame: Frame, centralBody: BodyIdentifier,
                gravitationalParameter: GravitationalParameter,
                semiLatusRectumMetres: Double, eccentricity: Double,
                inclination: Double, longitudeOfAscendingNode: Double,
                argumentOfPeriapsis: Double, trueAnomaly: Double) throws {
        guard frame.pseudoInertial else { throw OrbitalElementsError.nonInertialFrame }
        try centralBody.validate()
        guard semiLatusRectumMetres.isFinite, semiLatusRectumMetres > 0,
              eccentricity.isFinite, eccentricity >= 0,
              inclination.isFinite, inclination >= 0, inclination <= .pi,
              longitudeOfAscendingNode.isFinite, argumentOfPeriapsis.isFinite,
              trueAnomaly.isFinite else { throw OrbitalElementsError.invalidElements }
        let anomaly = trueAnomaly.remainder(dividingBy: 2 * .pi)
        guard 1 + eccentricity * cos(anomaly) > 0 else { throw OrbitalElementsError.anomalyOutsideConic }
        // Canonicalize exact singularities, preserving the combined phase.
        let equatorial = inclination == 0 || inclination == .pi
        let circular = eccentricity == 0
        let direction = inclination == .pi ? -1.0 : 1.0
        let inputNode = Self.positiveAngle(longitudeOfAscendingNode)
        let inputPeriapsis = Self.positiveAngle(argumentOfPeriapsis)
        let node = equatorial ? 0 : inputNode
        let periapsis = equatorial ? inputPeriapsis + direction * inputNode : inputPeriapsis
        self.epoch = epoch; self.frame = frame; self.centralBody = centralBody
        self.gravitationalParameter = gravitationalParameter
        self.semiLatusRectumMetres = semiLatusRectumMetres; self.eccentricity = eccentricity
        self.inclination = inclination
        self.longitudeOfAscendingNode = Self.positiveAngle(node)
        self.argumentOfPeriapsis = circular ? 0 : Self.positiveAngle(periapsis)
        self.trueAnomaly = (anomaly + (circular ? periapsis : 0)).remainder(dividingBy: 2 * .pi)
        self.angleConvention = Self.convention(circular: circular, equatorial: equatorial)
    }

    /// The state must already be relative to centralBody, including velocity,
    /// in nonrotating axes. No origin translation or frame conversion is inferred.
    public init(state: CartesianState, centralBody: BodyIdentifier,
                gravitationalParameter: GravitationalParameter,
                options: OrbitalElementConversionOptions = .default) throws {
        guard state.frame.pseudoInertial else { throw OrbitalElementsError.nonInertialFrame }
        let r = state.position, v = state.velocity
        let radius = Self.norm(r), speed = Self.norm(v)
        guard radius > 0 else { throw OrbitalElementsError.zeroPosition }
        guard radius.isFinite, speed.isFinite else { throw OrbitalElementsError.numericalRangeExceeded }
        guard speed > 0 else { throw OrbitalElementsError.radialTrajectory }
        let radial = r / radius
        let angular = Self.cross(radial, v / speed)
        let angularNorm = Self.norm(angular)
        guard angularNorm > 0 else { throw OrbitalElementsError.radialTrajectory }
        let normal = angular / angularNorm
        let mu = gravitationalParameter.cubicMetresPerSecondSquared
        let h = radius * (speed * angularNorm)
        let p = (h / mu) * h
        let eccentricityVector = Self.cross(v, normal) * (h / mu) - radial
        var e = Self.norm(eccentricityVector)
        guard p.isFinite, p > 0, e.isFinite else { throw OrbitalElementsError.numericalRangeExceeded }
        let nodeVector = SIMD3<Double>(-normal.y, normal.x, 0)
        let sineInclination = Self.norm(nodeVector)
        var inclination = atan2(sineInclination, normal.z)
        let equatorial = sineInclination <= options.equatorialSineTolerance
        let circular = e <= options.circularEccentricityTolerance
        if circular { e = 0 }
        if equatorial { inclination = normal.z >= 0 ? 0 : .pi }
        let nodeAxis = equatorial ? SIMD3<Double>(1, 0, 0) : nodeVector / sineInclination
        let periapsisAxis = circular ? nodeAxis : eccentricityVector / Self.norm(eccentricityVector)
        let node = equatorial ? 0 : atan2(nodeAxis.y, nodeAxis.x)
        let periapsis = circular ? 0 : Self.angle(from: nodeAxis, to: periapsisAxis, normal: normal)
        let anomaly = Self.angle(from: periapsisAxis, to: radial, normal: normal)
        try self.init(epoch: state.epoch, frame: state.frame, centralBody: centralBody,
                      gravitationalParameter: gravitationalParameter, semiLatusRectumMetres: p,
                      eccentricity: e, inclination: inclination, longitudeOfAscendingNode: node,
                      argumentOfPeriapsis: periapsis, trueAnomaly: anomaly)
    }

    /// Reconstructs position and velocity at the element epoch. Acceleration is
    /// unavailable: extracting elements does not assert a two-body force model.
    public func cartesianState() throws -> CartesianState {
        let p = semiLatusRectumMetres, e = eccentricity, nu = trueAnomaly
        let denominator = 1 + e * cos(nu)
        guard denominator > 0 else { throw OrbitalElementsError.anomalyOutsideConic }
        let radius = p / denominator
        let velocityScale = sqrt(gravitationalParameter.cubicMetresPerSecondSquared) / sqrt(p)
        let cO = cos(longitudeOfAscendingNode), sO = sin(longitudeOfAscendingNode)
        let co = cos(argumentOfPeriapsis), so = sin(argumentOfPeriapsis)
        let ci = cos(inclination), si = sin(inclination)
        let x = SIMD3<Double>(cO * co - sO * so * ci, sO * co + cO * so * ci, so * si)
        let y = SIMD3<Double>(-cO * so - sO * co * ci, -sO * so + cO * co * ci, co * si)
        let position = radius * (cos(nu) * x + sin(nu) * y)
        let velocity = velocityScale * (-sin(nu) * x + (e + cos(nu)) * y)
        guard [position.x, position.y, position.z, velocity.x, velocity.y, velocity.z].allSatisfy(\.isFinite) else {
            throw OrbitalElementsError.numericalRangeExceeded
        }
        return try CartesianState(epoch: epoch, frame: frame, position: position, velocity: velocity)
    }

    private static func norm(_ x: SIMD3<Double>) -> Double { hypot(hypot(x.x, x.y), x.z) }
    private static func cross(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> SIMD3<Double> {
        SIMD3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x)
    }
    private static func dot(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> Double { a.x*b.x+a.y*b.y+a.z*b.z }
    private static func angle(from a: SIMD3<Double>, to b: SIMD3<Double>, normal: SIMD3<Double>) -> Double {
        atan2(dot(cross(a, b), normal), dot(a, b))
    }
    private static func positiveAngle(_ x: Double) -> Double {
        let value = x.remainder(dividingBy: 2 * .pi)
        return value < 0 ? value + 2 * .pi : value
    }
    private static func convention(circular: Bool, equatorial: Bool) -> OrbitalAngleConvention {
        if circular { return equatorial ? .circularEquatorial : .circularInclined }
        return equatorial ? .eccentricEquatorial : .classical
    }
}
