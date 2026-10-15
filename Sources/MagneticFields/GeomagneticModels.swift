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

public struct WMMMagneticField: MagneticFieldModel {
    private let evaluator: GeomagneticHarmonics
    public var metadata: MagneticFieldMetadata { evaluator.metadata }
    /// Standard WMM .COF text; validity is explicit, normally five years.
    public init(coefficientText: String, validityYears: Double = 5,
                source: String = "Caller-supplied WMM coefficient text") throws {
        evaluator = try GeomagneticHarmonics(coefficients: .wmm(coefficientText,validityYears: validityYears),source: source)
    }
    public init(coefficientFile: URL, validityYears: Double = 5) throws {
        try self.init(coefficientText: String(contentsOf: coefficientFile,encoding: .utf8),
                      validityYears: validityYears,source: coefficientFile.absoluteString)
    }
    public static func wmm2025() throws -> Self {
        try Self(coefficientText: String(contentsOf: Bundle.module.url(forResource: "WMM2025",withExtension: "COF")!,encoding: .utf8),
                 source: "https://www.ncei.noaa.gov/products/world-magnetic-model; WMM2025 2024-11-13")
    }
    public func field(at query: MagneticFieldQuery) throws -> MagneticFluxDensity { try evaluator.field(at: query) }
}

public struct IGRFMagneticField: MagneticFieldModel {
    private let evaluator: GeomagneticHarmonics
    public var metadata: MagneticFieldMetadata { evaluator.metadata }
    /// Traditional IGRF g/h table with epoch columns and a final SV interval column.
    public init(coefficientText: String, source: String = "Caller-supplied IGRF coefficient text") throws {
        evaluator = try GeomagneticHarmonics(coefficients: .igrf(coefficientText),source: source)
    }
    public init(coefficientFile: URL) throws {
        try self.init(coefficientText: String(contentsOf: coefficientFile,encoding: .utf8),source: coefficientFile.absoluteString)
    }
    public static func igrf14() throws -> Self {
        try Self(coefficientText: String(contentsOf: Bundle.module.url(forResource: "igrf14coeffs",withExtension: "txt")!,encoding: .utf8),
                 source: "https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt; IGRF-14 2024-11-18")
    }
    public func field(at query: MagneticFieldQuery) throws -> MagneticFluxDensity { try evaluator.field(at: query) }
}

/// Geodetic latitude/longitude in radians, height above WGS84 ellipsoid in metres.
/// Mean-sea-level/geoid heights must be converted by the caller.
public struct GeodeticMagneticQuery: Sendable {
    public let latitudeRadians: Double
    public let longitudeRadians: Double
    public let ellipsoidHeightMetres: Double
    public let decimalYear: Double
    public init(latitudeRadians: Double, longitudeRadians: Double, ellipsoidHeightMetres: Double,
                decimalYear: Double) throws {
        guard [latitudeRadians,longitudeRadians,ellipsoidHeightMetres].allSatisfy(\.isFinite),
              abs(latitudeRadians) <= .pi/2, (-1000.0...850000.0).contains(ellipsoidHeightMetres) else {
            throw MagneticFieldError.invalidPosition
        }
        guard decimalYear.isFinite else { throw MagneticFieldError.invalidDate }
        self.latitudeRadians = latitudeRadians; self.longitudeRadians = longitudeRadians
        self.ellipsoidHeightMetres = ellipsoidHeightMetres; self.decimalYear = decimalYear
    }
    public func cartesianQuery() throws -> MagneticFieldQuery {
        let s = sin(latitudeRadians), c = cos(latitudeRadians)
        let a = 6_378_137.0, inverseFlattening = 298.257223563
        let f = 1/inverseFlattening, e2 = f*(2-f)
        let n = a/sqrt(1-e2*s*s), h = ellipsoidHeightMetres
        return try MagneticFieldQuery(positionMetres: SIMD3((n+h)*c*cos(longitudeRadians),
            (n+h)*c*sin(longitudeRadians),(n*(1-e2)+h)*s),decimalYear: decimalYear)
    }
}

public struct GeomagneticElements: Sendable {
    /// Local north/east/down components in tesla.
    public let northEastDownTesla: SIMD3<Double>
    public var horizontalIntensityTesla: Double { hypot(northEastDownTesla.x,northEastDownTesla.y) }
    public var totalIntensityTesla: Double { hypot(horizontalIntensityTesla,northEastDownTesla.z) }
    /// East-positive declination; undefined for zero horizontal field.
    public var declinationRadians: Double? {
        horizontalIntensityTesla > 0 ? atan2(northEastDownTesla.y,northEastDownTesla.x) : nil
    }
    public var inclinationRadians: Double? {
        totalIntensityTesla > 0 ? atan2(northEastDownTesla.z,horizontalIntensityTesla) : nil
    }
}
public extension MagneticFieldModel {
    /// WGS84 geographic convenience evaluation; also usable with dipole/constant models
    /// whose axes match Earth-fixed XYZ. At geographic poles the supplied longitude
    /// defines the local north/east basis; compass declination has no unique geographic meaning.
    func elements(at query: GeodeticMagneticQuery) throws -> GeomagneticElements {
        let b = try checkedField(at: query.cartesianQuery()).vectorTesla
        let s = sin(query.latitudeRadians), c = cos(query.latitudeRadians)
        let sl = sin(query.longitudeRadians), cl = cos(query.longitudeRadians)
        let ned = SIMD3(-s*cl*b.x-s*sl*b.y+c*b.z,
                         -sl*b.x+cl*b.y,-c*cl*b.x-c*sl*b.y-s*b.z)
        _ = try MagneticFluxDensity(vectorTesla: ned)
        return GeomagneticElements(northEastDownTesla: ned)
    }
}
