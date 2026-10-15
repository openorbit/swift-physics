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
import FluidDynamics
import Foundation

/// Published baseline and approximate composition scalings from Carroll–Brandis 2022.
/// Saturn/Uranus factors were derived with hot-wall CFD and also used by the authors
/// for cold-wall trajectory comparisons; they are not independent cold-wall fits.
public enum HydrogenHeliumHeatingMixture: Sendable {
    case neptuneReference, saturnApproximation, uranusApproximation
    public var gasIdentifier: String {
        switch self {
        case .neptuneReference: "h2-he-81.55-18.45-mole"
        case .saturnApproximation: "h2-he-saturn-11-percent-he-mole"
        case .uranusApproximation: "h2-he-uranus-15-percent-he-mole"
        }
    }
    public var heliumMoleFraction: Double {
        switch self { case .neptuneReference: 0.1845; case .saturnApproximation: 0.11; case .uranusApproximation: 0.15 }
    }
    var scaling: Double {
        switch self { case .neptuneReference: 1; case .saturnApproximation: 0.952; case .uranusApproximation: 0.980 }
    }
}

public struct HydrogenHeliumHeatingResult: Sendable {
    public let heatFlux: SurfaceHeatFlux
    public let fittedStagnationPressurePascals: Double
    public let fittedEnthalpyDifferenceJoulesPerKilogram: Double
    /// K₂, with source units kg/(s m^1.5 atm^0.5).
    public let enthalpyCoefficient: Double
    public let usesHighSpeedFit: Bool
}

/// Carroll–Brandis 2022 cold-wall nonlinear correlation, equations 4–6, table 3.
/// Fully catalytic, laminar, non-ablating spherical surface at 300 K; freestream
/// temperature 141.2 K. No radiation. There is no implicit Jupiter extrapolation.
public struct HydrogenHeliumColdWallHeating: ConvectiveHeatingModel {
    public let mixture: HydrogenHeliumHeatingMixture
    public let source = "https://ntrs.nasa.gov/api/citations/20220018610/downloads/AIAA_Aviation_2022_Paper-6.pdf"
    public let version = "Carroll–Brandis 2022 cold-wall fit"
    public init(mixture: HydrogenHeliumHeatingMixture) { self.mixture = mixture }

    /// Conservative piecewise log-density envelope through table-2 end points.
    /// Interior values interpolate coverage, not CFD data or a certified accuracy bound.
    public func supportedSpeedRange(densityKilogramsPerCubicMetre density: Double) throws -> ClosedRange<Double> {
        guard density.isFinite, (1e-5...0.1).contains(density) else { throw FluidError.outsideCoverage }
        let levels: [(Double, Double, Double)] = [(-5,15000,45000),(-4,7000,45000),
            (-3,7000,27000),(-2,7000,23000),(-1,7000,23000)]
        let logarithm = log10(density)
        for i in 0..<levels.count-1 {
            let a = levels[i], b = levels[i+1]
            if logarithm <= b.0 {
                let t = max(0,min(1,logarithm-a.0))
                return (a.1+(b.1-a.1)*t)...(a.2+(b.2-a.2)*t)
            }
        }
        return 7000...23000
    }
    public func heatFlux(for query: StagnationHeatingQuery) throws -> SurfaceHeatFlux {
        try evaluate(query).heatFlux
    }
    public func evaluate(_ query: StagnationHeatingQuery) throws -> HydrogenHeliumHeatingResult {
        guard query.gasIdentifier == mixture.gasIdentifier else { throw FluidError.unsupportedQuery }
        // Nil explicitly adopts the model's fixed 300 K wall. A specified wall must match.
        guard query.wall.map({ abs($0.temperatureKelvin-300) <= 1e-6 }) ?? true,
              query.gasState == nil else { throw FluidError.unsupportedQuery }
        let flow = query.flow, v = flow.speedMetresPerSecond, rho = flow.densityKilogramsPerCubicMetre
        guard (0.2...1).contains(query.noseRadiusMetres) else { throw FluidError.outsideCoverage }
        guard let temperature = flow.temperatureKelvin, let mach = flow.machNumber,
              let kn = try flow.knudsenNumber(referenceLengthMetres: query.noseRadiusMetres) else {
            throw FluidError.missingFlowProperty
        }
        guard abs(temperature-141.2) <= 1e-6 else { throw FluidError.unsupportedQuery }
        guard try supportedSpeedRange(densityKilogramsPerCubicMetre: rho).contains(v),
              mach >= 5, kn <= 0.003 else { throw FluidError.outsideCoverage }
        let u = v/45000, radius = query.noseRadiusMetres, densityTerm = log10(rho)/(-5)
        let high = v >= 29000
        let scaledFlux: Double
        if high {
            scaledFlux = 59.084*(u-0.2199)*pow(radius,-0.0429)*pow(densityTerm,-0.0287)
        } else {
            scaledFlux = 73.180*(u*u-0.0583*u)*pow(radius,-0.0184)*pow(densityTerm,-0.0254)
        }
        let pressure = 0.9474*rho*v*v // Figure 3a; Pa before converting to source atm.
        let difference = 0.5024*v*v-166.9*v
        let coefficient = scaledFlux/difference*1e6*mixture.scaling
        let flux = scaledFlux*1e6*sqrt((pressure/101325)/radius)*mixture.scaling
        guard flux.isFinite, coefficient.isFinite, flux >= 0 else { throw FluidError.nonFiniteResult }
        return HydrogenHeliumHeatingResult(heatFlux: try SurfaceHeatFlux(wattsPerSquareMetre: flux,
            surfacePositionMetres: query.surfacePositionMetres), fittedStagnationPressurePascals: pressure,
            fittedEnthalpyDifferenceJoulesPerKilogram: difference,
            enthalpyCoefficient: coefficient, usesHighSpeedFit: high)
    }
}
