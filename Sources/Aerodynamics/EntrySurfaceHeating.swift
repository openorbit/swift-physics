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

public enum AerothermalBoundaryLayer: String, Sendable, Codable { case laminar, turbulent }

/// Edge and wall enthalpies share a dataset-defined reference and composition.
/// Adiabatic-wall enthalpy includes recovery; no perfect-gas cp*T substitution is made.
public struct EnthalpySurfaceHeatingQuery: Sendable {
    public let edgeFlow: FluidFlow
    public let gasIdentifier: String
    public let distanceFromLeadingEdgeMetres: Double
    public let adiabaticWallEnthalpyJoulesPerKilogram: Double
    public let wallEnthalpyJoulesPerKilogram: Double
    public let surfacePositionMetres: SIMD3<Double>
    public init(edgeFlow: FluidFlow, gasIdentifier: String, distanceFromLeadingEdgeMetres: Double,
                adiabaticWallEnthalpyJoulesPerKilogram: Double, wallEnthalpyJoulesPerKilogram: Double,
                surfacePositionMetres: SIMD3<Double>) throws {
        guard !gasIdentifier.isEmpty, distanceFromLeadingEdgeMetres.isFinite, distanceFromLeadingEdgeMetres > 0,
              adiabaticWallEnthalpyJoulesPerKilogram.isFinite, adiabaticWallEnthalpyJoulesPerKilogram > 0,
              wallEnthalpyJoulesPerKilogram.isFinite, wallEnthalpyJoulesPerKilogram >= 0,
              fluidFinite(surfacePositionMetres) else { throw FluidError.invalidParameter }
        self.edgeFlow = edgeFlow; self.gasIdentifier = gasIdentifier
        self.distanceFromLeadingEdgeMetres = distanceFromLeadingEdgeMetres
        self.adiabaticWallEnthalpyJoulesPerKilogram = adiabaticWallEnthalpyJoulesPerKilogram
        self.wallEnthalpyJoulesPerKilogram = wallEnthalpyJoulesPerKilogram
        self.surfacePositionMetres = surfacePositionMetres
    }
}
public struct EnthalpySurfaceHeatingResult: Sendable {
    /// Signed heat into the surface, W/m².
    public let wattsPerSquareMetre: Double
    public let stantonNumber: Double
    public let surfacePositionMetres: SIMD3<Double>
}
public protocol EnthalpySurfaceHeatingModel: Sendable {
    func evaluate(_ query: EnthalpySurfaceHeatingQuery) throws -> EnthalpySurfaceHeatingResult
}


/// Additional thermodynamic/geometry coverage for a supplied Stanton dataset.
/// Bounds must enclose only states established by the dataset's assumptions.
public struct EnthalpySurfaceHeatingCoverage: Sendable {
    public let edgeTemperatureKelvin: ClosedRange<Double>
    public let edgePressurePascals: ClosedRange<Double>
    public let adiabaticWallEnthalpyJoulesPerKilogram: ClosedRange<Double>
    public let distanceFromLeadingEdgeMetres: ClosedRange<Double>
    public init(edgeTemperatureKelvin: ClosedRange<Double>, edgePressurePascals: ClosedRange<Double>,
                adiabaticWallEnthalpyJoulesPerKilogram: ClosedRange<Double>, distanceFromLeadingEdgeMetres: ClosedRange<Double>) throws {
        guard [edgeTemperatureKelvin,edgePressurePascals,adiabaticWallEnthalpyJoulesPerKilogram,distanceFromLeadingEdgeMetres]
            .allSatisfy({ $0.lowerBound.isFinite && $0.lowerBound > 0 && $0.upperBound.isFinite }) else {
            throw FluidError.invalidParameter
        }
        self.edgeTemperatureKelvin = edgeTemperatureKelvin; self.edgePressurePascals = edgePressurePascals
        self.adiabaticWallEnthalpyJoulesPerKilogram = adiabaticWallEnthalpyJoulesPerKilogram
        self.distanceFromLeadingEdgeMetres = distanceFromLeadingEdgeMetres
    }
}

/// Compressible/high-enthalpy local convection from supplied Stanton data, not a
/// universal gas correlation. Axes: Mach, log10(Re_x), wall/adiabatic-wall enthalpy.
/// q = rho_edge V_edge St (h_aw - h_wall). Recovery, chemistry, catalytic state,
/// pressure-gradient and geometry assumptions must be identified in the dataset.
public struct TabulatedEnthalpySurfaceHeating: EnthalpySurfaceHeatingModel {
    public let gasIdentifier: String
    public let boundaryLayer: AerothermalBoundaryLayer
    public let stanton: FluidScalarTable
    public let coverage: EnthalpySurfaceHeatingCoverage
    public let maximumKnudsenNumber: Double
    public init(gasIdentifier: String, boundaryLayer: AerothermalBoundaryLayer,
                stanton: FluidScalarTable, coverage: EnthalpySurfaceHeatingCoverage, maximumKnudsenNumber: Double = 0.01) throws {
        guard !gasIdentifier.isEmpty, stanton.valueUnits == "1", stanton.values.allSatisfy({ $0 >= 0 }),
              stanton.axes.map(\.name) == ["mach","log10Re","wallEnthalpyRatio"],
              stanton.axes.allSatisfy({ $0.units == "1" }),
              stanton.axes[0].coordinates.first! >= 0, stanton.axes[2].coordinates.first! >= 0,
              maximumKnudsenNumber.isFinite, maximumKnudsenNumber > 0, maximumKnudsenNumber <= 0.01 else {
            throw FluidError.malformedData
        }
        self.gasIdentifier = gasIdentifier; self.boundaryLayer = boundaryLayer
        self.stanton = stanton; self.coverage = coverage; self.maximumKnudsenNumber = maximumKnudsenNumber
    }
    public func evaluate(_ query: EnthalpySurfaceHeatingQuery) throws -> EnthalpySurfaceHeatingResult {
        guard query.gasIdentifier == gasIdentifier else { throw FluidError.unsupportedQuery }
        let flow = query.edgeFlow, x = query.distanceFromLeadingEdgeMetres
        guard let temperature = flow.temperatureKelvin, let pressure = flow.pressurePascals,
              let mach = flow.machNumber, let re = try flow.reynoldsNumber(referenceLengthMetres: x),
              let kn = try flow.knudsenNumber(referenceLengthMetres: x) else { throw FluidError.missingFlowProperty }
        guard re > 0, flow.speedMetresPerSecond > 0, kn <= maximumKnudsenNumber,
              coverage.edgeTemperatureKelvin.contains(temperature), coverage.edgePressurePascals.contains(pressure),
              coverage.adiabaticWallEnthalpyJoulesPerKilogram.contains(query.adiabaticWallEnthalpyJoulesPerKilogram),
              coverage.distanceFromLeadingEdgeMetres.contains(x) else { throw FluidError.outsideCoverage }
        let st = try stanton.value(at: [mach,log10(re),query.wallEnthalpyJoulesPerKilogram/query.adiabaticWallEnthalpyJoulesPerKilogram])
        let flux = flow.densityKilogramsPerCubicMetre*flow.speedMetresPerSecond*st
            * (query.adiabaticWallEnthalpyJoulesPerKilogram-query.wallEnthalpyJoulesPerKilogram)
        guard flux.isFinite else { throw FluidError.nonFiniteResult }
        return EnthalpySurfaceHeatingResult(wattsPerSquareMetre: flux,stantonNumber: st,surfacePositionMetres: query.surfacePositionMetres)
    }
}

/// Leading-edge stagnation state plus effective sweep relative to incident flow.
public struct LeadingEdgeHeatingQuery: Sendable {
    public let stagnation: StagnationHeatingQuery
    public let effectiveSweepRadians: Double
    public init(stagnation: StagnationHeatingQuery, effectiveSweepRadians: Double) throws {
        guard effectiveSweepRadians.isFinite, abs(effectiveSweepRadians) < .pi/2 else {
            throw FluidError.invalidGeometry
        }
        self.stagnation = stagnation; self.effectiveSweepRadians = effectiveSweepRadians
    }
}
public protocol LeadingEdgeConvectiveHeatingModel: Sendable {
    func heatFlux(for query: LeadingEdgeHeatingQuery) throws -> SurfaceHeatFlux
}

/// Idealized laminar swept-cylinder stagnation line from a supplied sphere reference.
/// NASA CR-181922 eq. 1: qSL=qSphere/sqrt(2)*sqrt(RSphere/RLE)*cos(sweep)^1.1.
/// Earth-air only. Only the stagnation line, not a circumferential distribution or wing acreage.
public struct SweptCylinderStagnationHeating: LeadingEdgeConvectiveHeatingModel {
    public let referenceSphereModel: any ConvectiveHeatingModel
    public let referenceSphereRadiusMetres: Double
    public let coverage: StagnationHeatingCoverage
    public let maximumSweepRadians: Double
    public var source: String { "https://ntrs.nasa.gov/api/citations/19900007739/downloads/19900007739.pdf" }
    public init(referenceSphereModel: any ConvectiveHeatingModel, referenceSphereRadiusMetres: Double,
                coverage: StagnationHeatingCoverage, maximumSweepRadians: Double = .pi/3) throws {
        guard referenceSphereRadiusMetres.isFinite, referenceSphereRadiusMetres > 0,
              maximumSweepRadians.isFinite, maximumSweepRadians >= 0, maximumSweepRadians <= .pi/3 else {
            throw FluidError.invalidParameter
        }
        self.referenceSphereModel = referenceSphereModel; self.referenceSphereRadiusMetres = referenceSphereRadiusMetres
        self.coverage = coverage; self.maximumSweepRadians = maximumSweepRadians
    }
    public func heatFlux(for query: LeadingEdgeHeatingQuery) throws -> SurfaceHeatFlux {
        try heatFlux(for: query.stagnation,effectiveSweepRadians: query.effectiveSweepRadians)
    }
    /// Query nose radius denotes the actual leading-edge radius. Effective sweep is
    /// relative to the incident flow, not automatically the planform sweep angle.
    public func heatFlux(for query: StagnationHeatingQuery, effectiveSweepRadians: Double) throws -> SurfaceHeatFlux {
        // The published geometry correction is an Earth-air approximation.
        guard query.gasIdentifier == "earth-air" else { throw FluidError.unsupportedQuery }
        let flow = query.flow
        guard effectiveSweepRadians.isFinite, abs(effectiveSweepRadians) <= maximumSweepRadians,
              coverage.noseRadiusMetres.contains(query.noseRadiusMetres),
              coverage.speedMetresPerSecond.contains(flow.speedMetresPerSecond),
              coverage.densityKilogramsPerCubicMetre.contains(flow.densityKilogramsPerCubicMetre) else {
            throw FluidError.outsideCoverage
        }
        guard let mach = flow.machNumber, let kn = try flow.knudsenNumber(referenceLengthMetres: query.noseRadiusMetres) else {
            throw FluidError.missingFlowProperty
        }
        guard coverage.machNumber.contains(mach), kn <= coverage.maximumKnudsenNumber else { throw FluidError.outsideCoverage }
        let sphereQuery = try StagnationHeatingQuery(flow: flow,gasIdentifier: query.gasIdentifier,
            noseRadiusMetres: referenceSphereRadiusMetres,surfacePositionMetres: query.surfacePositionMetres,
            wall: query.wall,gasState: query.gasState)
        let reference = try referenceSphereModel.heatFlux(for: sphereQuery)
        let flux = reference.wattsPerSquareMetre/sqrt(2)*sqrt(referenceSphereRadiusMetres/query.noseRadiusMetres)
            * pow(cos(effectiveSweepRadians),1.1)
        return try SurfaceHeatFlux(wattsPerSquareMetre: flux,surfacePositionMetres: query.surfacePositionMetres)
    }
}

/// Incoming shock-layer radiation from supplied gas-specific data; no radiation
/// transport solver or empirical universal planetary coefficient is implied.
/// Axes: speed (m/s), density (kg/m^3), sphere nose radius (m). Cold-wall queries only.
public struct TabulatedRadiativeHeating: RadiativeHeatingModel {
    public let gasIdentifier: String
    public let flux: FluidScalarTable
    public let freestreamTemperatureKelvin: ClosedRange<Double>
    public let coverage: StagnationHeatingCoverage
    public init(gasIdentifier: String, flux: FluidScalarTable, coverage: StagnationHeatingCoverage,
                freestreamTemperatureKelvin: ClosedRange<Double>) throws {
        guard freestreamTemperatureKelvin.lowerBound.isFinite, freestreamTemperatureKelvin.lowerBound > 0,
              freestreamTemperatureKelvin.upperBound.isFinite, !gasIdentifier.isEmpty, flux.valueUnits == "W/m^2", flux.values.allSatisfy({ $0 >= 0 }),
              flux.axes.map(\.name) == ["speed","density","radius"],
              flux.axes.map(\.units) == ["m/s","kg/m^3","m"],
              flux.axes.allSatisfy({ $0.coordinates.first! > 0 }) else { throw FluidError.malformedData }
        self.gasIdentifier = gasIdentifier; self.flux = flux; self.coverage = coverage
        self.freestreamTemperatureKelvin = freestreamTemperatureKelvin
    }
    public func heatFlux(for query: StagnationHeatingQuery) throws -> SurfaceHeatFlux {
        guard query.gasIdentifier == gasIdentifier, query.wall == nil, query.gasState == nil else {
            throw FluidError.unsupportedQuery
        }
        let flow = query.flow
        guard let temperature = flow.temperatureKelvin, let mach = flow.machNumber,
              let kn = try flow.knudsenNumber(referenceLengthMetres: query.noseRadiusMetres) else {
            throw FluidError.missingFlowProperty
        }
        guard freestreamTemperatureKelvin.contains(temperature), coverage.speedMetresPerSecond.contains(flow.speedMetresPerSecond),
              coverage.densityKilogramsPerCubicMetre.contains(flow.densityKilogramsPerCubicMetre),
              coverage.noseRadiusMetres.contains(query.noseRadiusMetres), coverage.machNumber.contains(mach),
              kn <= coverage.maximumKnudsenNumber else { throw FluidError.outsideCoverage }
        return try SurfaceHeatFlux(wattsPerSquareMetre: flux.value(at: [flow.speedMetresPerSecond,
            flow.densityKilogramsPerCubicMetre,query.noseRadiusMetres]),surfacePositionMetres: query.surfacePositionMetres)
    }
}
