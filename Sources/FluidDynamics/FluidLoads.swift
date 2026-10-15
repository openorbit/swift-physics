//
// SPDX-License-Identifier: Apache-2.0
//
// Copyright 2026 Mattias Holm
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
import Foundation

/// Dimensions for nondimensional force/moment coefficients in one stated set of axes.
/// Moment lengths are per axis (typically span, chord, span for aircraft body axes).
public struct FluidReferenceGeometry: Sendable {
    public let areaSquareMetres: Double
    public let momentLengthsMetres: SIMD3<Double>
    public let momentReferencePointMetres: SIMD3<Double>
    public init(areaSquareMetres: Double, momentLengthsMetres: SIMD3<Double>,
                momentReferencePointMetres: SIMD3<Double> = .zero) throws {
        guard areaSquareMetres.isFinite, areaSquareMetres > 0,
              [momentLengthsMetres.x,momentLengthsMetres.y,momentLengthsMetres.z].allSatisfy({ $0.isFinite && $0 > 0 }),
              fluidFinite(momentReferencePointMetres) else { throw FluidError.invalidGeometry }
        self.areaSquareMetres = areaSquareMetres; self.momentLengthsMetres = momentLengthsMetres
        self.momentReferencePointMetres = momentReferencePointMetres
    }
}

/// Six signed Cartesian coefficients, including all moments about the supplied
/// geometry's reference point. Axes must match geometry and the returned loads.
public struct FluidCoefficients: Sendable, Codable {
    public let force: SIMD3<Double>
    public let moment: SIMD3<Double>
    public init(force: SIMD3<Double>, moment: SIMD3<Double>) throws {
        guard fluidFinite(force), fluidFinite(moment) else { throw FluidError.invalidCoefficients }
        self.force = force; self.moment = moment
    }
    private enum CodingKeys: String, CodingKey { case force, moment }
    public init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        try self.init(force: c.decode(SIMD3<Double>.self,forKey: .force),moment: c.decode(SIMD3<Double>.self,forKey: .moment))
    }
    /// F=q S Cforce, M=q S L Cmoment. No separate centre-of-pressure moment is added.
    public func loads(dynamicPressurePascals q: Double, geometry: FluidReferenceGeometry) throws -> FluidLoads {
        guard q.isFinite, q >= 0 else { throw FluidError.invalidFlow }
        let scale = q*geometry.areaSquareMetres
        guard scale.isFinite else { throw FluidError.nonFiniteResult }
        return try FluidLoads(forceNewtons: force*scale,
            momentNewtonMetres: (moment*scale)*geometry.momentLengthsMetres,
            referencePointMetres: geometry.momentReferencePointMetres)
    }
}

/// Physical force and moment about an explicit reference point, all in the same axes.
/// No body identity, spacecraft state, assembly or integrator is owned here.
public struct FluidLoads: Sendable {
    public let forceNewtons: SIMD3<Double>
    public let momentNewtonMetres: SIMD3<Double>
    public let referencePointMetres: SIMD3<Double>
    public init(forceNewtons: SIMD3<Double>, momentNewtonMetres: SIMD3<Double>,
                referencePointMetres: SIMD3<Double> = .zero) throws {
        guard fluidFinite(forceNewtons), fluidFinite(momentNewtonMetres), fluidFinite(referencePointMetres) else {
            throw FluidError.nonFiniteResult
        }
        self.forceNewtons = forceNewtons; self.momentNewtonMetres = momentNewtonMetres
        self.referencePointMetres = referencePointMetres
    }
    /// Transfer an already complete moment: Mnew = Mold + (oldPoint−newPoint)×F.
    public func referred(to pointMetres: SIMD3<Double>) throws -> Self {
        guard fluidFinite(pointMetres) else { throw FluidError.invalidGeometry }
        return try Self(forceNewtons: forceNewtons,
            momentNewtonMetres: momentNewtonMetres+fluidCross(referencePointMetres-pointMetres,forceNewtons),
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
public struct SurfaceHeatFlux: Sendable {
    public let wattsPerSquareMetre: Double
    public let surfacePositionMetres: SIMD3<Double>
    public init(wattsPerSquareMetre: Double, surfacePositionMetres: SIMD3<Double>) throws {
        guard wattsPerSquareMetre.isFinite, wattsPerSquareMetre >= 0, fluidFinite(surfacePositionMetres) else {
            throw FluidError.invalidParameter
        }
        self.wattsPerSquareMetre = wattsPerSquareMetre; self.surfacePositionMetres = surfacePositionMetres
    }
}
