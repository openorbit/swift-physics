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

public enum MagneticBody: Int32, CaseIterable, Sendable {
    case mercury = 199, earth = 399, jupiter = 599, saturn = 699, uranus = 799, neptune = 899
}

/// Deliberately aligned, centred ideal dipoles, using nominal NASA field strengths.
/// No measured pole longitude, tilt, offset or secular variation is implied.
public struct BodyMagneticField: MagneticFieldModel {
    public let body: MagneticBody
    public let dipole: DipoleMagneticField
    public var metadata: MagneticFieldMetadata { dipole.metadata }
    public func field(at query: MagneticFieldQuery) throws -> MagneticFluxDensity { try dipole.field(at: query) }
    public static func model(for body: MagneticBody) -> BodyMagneticField {
        let radius: Double, gauss: Double, name: String
        switch body {
        case .mercury: (radius,gauss,name) = (2_439_700,-0.002,"mercury")
        case .earth: (radius,gauss,name) = (6_378_000,-0.306,"earth")
        case .jupiter: (radius,gauss,name) = (71_398_000,4.30,"jupiter")
        case .saturn: (radius,gauss,name) = (60_330_000,0.215,"saturn")
        case .uranus: (radius,gauss,name) = (25_600_000,0.228,"uranus")
        case .neptune: (radius,gauss,name) = (24_765_000,0.142,"neptune")
        }
        let metadata = try! MagneticFieldMetadata(name: "Nominal aligned \(name) dipole",
            source: "https://nssdc.gsfc.nasa.gov/planetary/factsheet/\(name)fact.html",
            version: "Idealized catalogue v1; reference strengths reviewed 2026-10-14",
            axesDescription: "Body-fixed right-handed axes; +Z north rotation pole, X/Y equatorial",
            assumptions: "Centred aligned dipole; signed -Z moment for Earth/Mercury, +Z for giant presets. Giant polarities are prescribed simplifications. Ignores measured tilts/offsets (especially large for Uranus/Neptune), harmonics, secular variation and external fields. Radius coverage 1...10 reference radii is a chosen simulation bound, not magnetopause geometry.",
            radiusRangeMetres: radius...10*radius)
        return BodyMagneticField(body: body, dipole: try! DipoleMagneticField(referenceRadiusMetres: radius,
            momentAxisFieldTesla: SIMD3(0,0,gauss*1e-4), metadata: metadata))
    }
    /// nil means no supplied intrinsic model, not a zero total magnetic field.
    public static func model(naifBodyCode: Int32) -> BodyMagneticField? {
        MagneticBody(rawValue: naifBodyCode).map { model(for: $0) }
    }
}
