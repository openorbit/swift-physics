// SPDX-License-Identifier: Apache-2.0
import Foundation

/// Dimensions for nondimensional force/moment coefficients in one stated set of axes.
/// Moment lengths are per axis (typically span, chord, span for aircraft body axes).
public struct AerodynamicReferenceGeometry: Sendable {
    public let areaSquareMetres: Double
    public let momentLengthsMetres: SIMD3<Double>
    public let momentReferencePointMetres: SIMD3<Double>
    public init(areaSquareMetres: Double, momentLengthsMetres: SIMD3<Double>,
                momentReferencePointMetres: SIMD3<Double> = .zero) throws {
        guard areaSquareMetres.isFinite, areaSquareMetres > 0,
              [momentLengthsMetres.x,momentLengthsMetres.y,momentLengthsMetres.z].allSatisfy({ $0.isFinite && $0 > 0 }),
              aeroFinite(momentReferencePointMetres) else { throw AerodynamicError.invalidGeometry }
        self.areaSquareMetres = areaSquareMetres; self.momentLengthsMetres = momentLengthsMetres
        self.momentReferencePointMetres = momentReferencePointMetres
    }
}

/// Six signed Cartesian coefficients, including all moments about the supplied
/// geometry's reference point. Axes must match geometry and the returned loads.
public struct AerodynamicCoefficients: Sendable {
    public let force: SIMD3<Double>
    public let moment: SIMD3<Double>
    public init(force: SIMD3<Double>, moment: SIMD3<Double>) throws {
        guard aeroFinite(force), aeroFinite(moment) else { throw AerodynamicError.invalidCoefficients }
        self.force = force; self.moment = moment
    }
    /// F=q S Cforce, M=q S L Cmoment. No separate centre-of-pressure moment is added.
    public func loads(dynamicPressurePascals q: Double, geometry: AerodynamicReferenceGeometry) throws -> AerodynamicLoads {
        guard q.isFinite, q >= 0 else { throw AerodynamicError.invalidFlow }
        let scale = q*geometry.areaSquareMetres
        guard scale.isFinite else { throw AerodynamicError.nonFiniteResult }
        return try AerodynamicLoads(forceNewtons: force*scale,
            momentNewtonMetres: (moment*scale)*geometry.momentLengthsMetres,
            referencePointMetres: geometry.momentReferencePointMetres)
    }
}

/// Physical force and moment about an explicit reference point, all in the same axes.
/// No body identity, spacecraft state, assembly or integrator is owned here.
public struct AerodynamicLoads: Sendable {
    public let forceNewtons: SIMD3<Double>
    public let momentNewtonMetres: SIMD3<Double>
    public let referencePointMetres: SIMD3<Double>
    public init(forceNewtons: SIMD3<Double>, momentNewtonMetres: SIMD3<Double>,
                referencePointMetres: SIMD3<Double> = .zero) throws {
        guard aeroFinite(forceNewtons), aeroFinite(momentNewtonMetres), aeroFinite(referencePointMetres) else {
            throw AerodynamicError.nonFiniteResult
        }
        self.forceNewtons = forceNewtons; self.momentNewtonMetres = momentNewtonMetres
        self.referencePointMetres = referencePointMetres
    }
    /// Transfer an already complete moment: Mnew = Mold + (oldPoint−newPoint)×F.
    public func referred(to pointMetres: SIMD3<Double>) throws -> Self {
        guard aeroFinite(pointMetres) else { throw AerodynamicError.invalidGeometry }
        return try Self(forceNewtons: forceNewtons,
            momentNewtonMetres: momentNewtonMetres+aeroCross(referencePointMetres-pointMetres,forceNewtons),
            referencePointMetres: pointMetres)
    }
    /// Sum loads after expressing both moments about this load's reference point.
    public func adding(_ other: Self) throws -> Self {
        let shifted = try other.referred(to: referencePointMetres)
        return try Self(forceNewtons: forceNewtons+shifted.forceNewtons,
            momentNewtonMetres: momentNewtonMetres+shifted.momentNewtonMetres,referencePointMetres: referencePointMetres)
    }
}

/// Surface-local heat-flux contract for future correlations and external thermal models.
/// Positive values represent heat into the surface. This type does not evolve temperature.
public struct AerodynamicHeatFlux: Sendable {
    public let wattsPerSquareMetre: Double
    public let surfacePositionMetres: SIMD3<Double>
    public init(wattsPerSquareMetre: Double, surfacePositionMetres: SIMD3<Double>) throws {
        guard wattsPerSquareMetre.isFinite, wattsPerSquareMetre >= 0, aeroFinite(surfacePositionMetres) else {
            throw AerodynamicError.invalidParameter
        }
        self.wattsPerSquareMetre = wattsPerSquareMetre; self.surfacePositionMetres = surfacePositionMetres
    }
}
