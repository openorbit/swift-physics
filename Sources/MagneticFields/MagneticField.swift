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

public enum MagneticFieldError: Error, Equatable {
    case invalidParameters, invalidPosition, invalidDate, outsideCoverage, singularField, nonFiniteResult
}

/// Position relative to the body centre in the provider's stated body-fixed axes.
/// Decimal year is a caller-supplied calendar date, not elapsed Julian years from J2000.
public struct MagneticFieldQuery: Sendable {
    public let positionMetres: SIMD3<Double>
    public let decimalYear: Double
    public init(positionMetres: SIMD3<Double>, decimalYear: Double) throws {
        guard [positionMetres.x,positionMetres.y,positionMetres.z].allSatisfy(\.isFinite) else {
            throw MagneticFieldError.invalidPosition
        }
        guard decimalYear.isFinite else { throw MagneticFieldError.invalidDate }
        self.positionMetres = positionMetres; self.decimalYear = decimalYear
    }
}

public struct MagneticFluxDensity: Sendable {
    /// Magnetic flux density B in tesla, in the same axes as the query position.
    public let vectorTesla: SIMD3<Double>
    public var magnitudeTesla: Double { hypot(hypot(vectorTesla.x,vectorTesla.y),vectorTesla.z) }
    public init(vectorTesla: SIMD3<Double>) throws {
        guard [vectorTesla.x,vectorTesla.y,vectorTesla.z].allSatisfy(\.isFinite),
              hypot(hypot(vectorTesla.x,vectorTesla.y),vectorTesla.z).isFinite else {
            throw MagneticFieldError.nonFiniteResult
        }
        self.vectorTesla = vectorTesla
    }
}

/// Coverage is geocentric radius, not altitude. nil date coverage means time-independent.
public struct MagneticFieldMetadata: Sendable {
    public let name: String
    public let source: String
    public let version: String
    public let axesDescription: String
    public let assumptions: String
    public let radiusRangeMetres: ClosedRange<Double>
    public let decimalYearRange: ClosedRange<Double>?
    public init(name: String, source: String, version: String, axesDescription: String,
                assumptions: String, radiusRangeMetres: ClosedRange<Double>,
                decimalYearRange: ClosedRange<Double>? = nil) throws {
        guard [name,source,version,axesDescription].allSatisfy({ !$0.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty }),
              radiusRangeMetres.lowerBound.isFinite, radiusRangeMetres.upperBound.isFinite,
              radiusRangeMetres.lowerBound >= 0,
              decimalYearRange.map({ $0.lowerBound.isFinite && $0.upperBound.isFinite }) ?? true else {
            throw MagneticFieldError.invalidParameters
        }
        self.name = name; self.source = source; self.version = version
        self.axesDescription = axesDescription; self.assumptions = assumptions
        self.radiusRangeMetres = radiusRangeMetres; self.decimalYearRange = decimalYearRange
    }
    public func validate(_ query: MagneticFieldQuery) throws {
        let p = query.positionMetres
        let radius = hypot(hypot(p.x,p.y),p.z)
        guard radius.isFinite, radiusRangeMetres.contains(radius),
              decimalYearRange.map({ $0.contains(query.decimalYear) }) ?? true else {
            throw MagneticFieldError.outsideCoverage
        }
    }
}

/// Pluggable field evaluator: dipoles, prescribed fields or future WMM/IGRF harmonics.
/// Providers perform no origin translations, frame transforms or calendar conversions.
public protocol MagneticFieldModel: Sendable {
    var metadata: MagneticFieldMetadata { get }
    func field(at query: MagneticFieldQuery) throws -> MagneticFluxDensity
}
public extension MagneticFieldModel {
    /// Enforces metadata coverage even for a caller-supplied provider.
    func checkedField(at query: MagneticFieldQuery) throws -> MagneticFluxDensity {
        try metadata.validate(query)
        return try field(at: query)
    }
}

/// B = (a/r)^3 [3 n (b·n) − b], where b = μ₀ m/(4π a³).
/// b points along the signed dipole moment; at the magnetic equator B points opposite b.
/// Supports tilted and displaced dipoles. No magnetopause, solar wind, induction or multipoles.
public struct DipoleMagneticField: MagneticFieldModel {
    public let metadata: MagneticFieldMetadata
    public let referenceRadiusMetres: Double
    public let momentAxisFieldTesla: SIMD3<Double>
    public let offsetMetres: SIMD3<Double>
    public init(referenceRadiusMetres: Double, momentAxisFieldTesla: SIMD3<Double>,
                offsetMetres: SIMD3<Double> = .zero, metadata: MagneticFieldMetadata) throws {
        guard referenceRadiusMetres.isFinite, referenceRadiusMetres > 0,
              [momentAxisFieldTesla.x,momentAxisFieldTesla.y,momentAxisFieldTesla.z,
               offsetMetres.x,offsetMetres.y,offsetMetres.z].allSatisfy(\.isFinite),
              hypot(hypot(momentAxisFieldTesla.x,momentAxisFieldTesla.y),momentAxisFieldTesla.z).isFinite else {
            throw MagneticFieldError.invalidParameters
        }
        self.referenceRadiusMetres = referenceRadiusMetres; self.momentAxisFieldTesla = momentAxisFieldTesla
        self.offsetMetres = offsetMetres; self.metadata = metadata
    }
    public func field(at query: MagneticFieldQuery) throws -> MagneticFluxDensity {
        try metadata.validate(query)
        let p = query.positionMetres-offsetMetres
        let radius = hypot(hypot(p.x,p.y),p.z)
        guard radius.isFinite else { throw MagneticFieldError.nonFiniteResult }
        guard radius > 0 else { throw MagneticFieldError.singularField }
        let n = p/radius, ratio = referenceRadiusMetres/radius
        // Normalize the moment before dot products to avoid avoidable overflow.
        let magnitude = hypot(hypot(momentAxisFieldTesla.x,momentAxisFieldTesla.y),momentAxisFieldTesla.z)
        if magnitude == 0 { return try MagneticFluxDensity(vectorTesla: .zero) }
        let axis = momentAxisFieldTesla/magnitude
        let dot = (axis*n).sum()
        let scale = magnitude*ratio*ratio*ratio
        guard scale.isFinite else { throw MagneticFieldError.nonFiniteResult }
        return try MagneticFluxDensity(vectorTesla: (3*dot*n-axis)*scale)
    }
}

/// Explicit caller-prescribed field, useful for controlled experiments and external fields.
public struct ConstantMagneticField: MagneticFieldModel {
    public let metadata: MagneticFieldMetadata
    private let sample: MagneticFluxDensity
    public init(vectorTesla: SIMD3<Double>, metadata: MagneticFieldMetadata) throws {
        self.sample = try MagneticFluxDensity(vectorTesla: vectorTesla); self.metadata = metadata
    }
    public func field(at query: MagneticFieldQuery) throws -> MagneticFluxDensity {
        try metadata.validate(query); return sample
    }
}
