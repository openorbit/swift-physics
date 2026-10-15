// SPDX-License-Identifier: Apache-2.0
import FluidDynamics
import Foundation

/// Molecular weights and Mγ transport products from NASA TR R-376, table I.
public enum SuttonGravesBaseGas: String, Sendable, CaseIterable {
    case nitrogen, oxygen, hydrogen, helium, neon, argon, carbonDioxide, ammonia, methane
    public var molecularWeight: Double {
        switch self {
        case .nitrogen: 28.014; case .oxygen: 32; case .hydrogen: 2.016
        case .helium: 4.003; case .neon: 20.183; case .argon: 39.948
        case .carbonDioxide: 44.011; case .ammonia: 17.031; case .methane: 16.043
        }
    }
    var molecularWeightTransportProduct: Double {
        switch self {
        case .nitrogen: 1.024; case .oxygen: 1.321; case .hydrogen: 0.137
        case .helium: 0.434; case .neon: 1.504; case .argon: 1.612
        case .carbonDioxide: 1.285; case .ammonia: 0.784; case .methane: 0.537
        }
    }
}

/// Equation 42 mixture K, distinct from the coefficient of sqrt(ρ/R) V³.
/// Uses wall Prandtl number 0.69 and the report's rounded table-I products.
public struct SuttonGravesMixtureParameters: Sendable {
    public let gasIdentifier: String
    public let massFractions: [SuttonGravesBaseGas: Double]
    /// Source units kg/(s m^1.5 atm^0.5). Evaluators convert supplied Pa to atm.
    public let enthalpyCoefficient: Double
    public let source = "https://ntrs.nasa.gov/citations/19720003329"
    public let version = "NASA TR R-376, equation 42 and table I"
    public init(gasIdentifier: String, massFractions: [SuttonGravesBaseGas: Double]) throws {
        guard !gasIdentifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              !massFractions.isEmpty, massFractions.values.allSatisfy({ $0.isFinite && $0 >= 0 }),
              abs(massFractions.values.reduce(0,+)-1) <= 1e-10 else { throw FluidError.invalidParameter }
        self.gasIdentifier = gasIdentifier; self.massFractions = massFractions
        let sum = SuttonGravesBaseGas.allCases.reduce(0.0) {
            $0 + (massFractions[$1] ?? 0)/$1.molecularWeightTransportProduct
        }
        enthalpyCoefficient = 0.1106/sqrt(sum)
    }
    public init(gasIdentifier: String, moleFractions: [SuttonGravesBaseGas: Double]) throws {
        guard !moleFractions.isEmpty, moleFractions.values.allSatisfy({ $0.isFinite && $0 >= 0 }),
              abs(moleFractions.values.reduce(0,+)-1) <= 1e-10 else { throw FluidError.invalidParameter }
        let mass = moleFractions.reduce(0.0) { $0+$1.value*$1.key.molecularWeight }
        try self.init(gasIdentifier: gasIdentifier, massFractions: moleFractions.reduce(into: [:]) {
            $0[$1.key] = $1.value*$1.key.molecularWeight/mass
        })
    }
    /// Idealized 96.5% CO₂ / 3.5% N₂ by mole, ignoring trace constituents.
    /// Composition source: NASA Venus fact sheet; this is a derived mixture, not
    /// a fitted or flight-validated Venus coefficient.
    public static var venusCarbonDioxideNitrogen: Self {
        try! Self(gasIdentifier: "venus-co2-n2-96.5-3.5-mole",
                  moleFractions: [.carbonDioxide: 0.965, .nitrogen: 0.035])
    }
}

/// Cold-wall enthalpy-form mixture approximation q=K sqrt(pStag[atm]/R) hEdge.
/// Requires caller-supplied shock-layer pressure/enthalpy; wall enthalpy is neglected.
public struct SuttonGravesMixtureHeating: ConvectiveHeatingModel {
    public let parameters: SuttonGravesMixtureParameters
    public let coverage: StagnationHeatingCoverage
    public init(parameters: SuttonGravesMixtureParameters, coverage: StagnationHeatingCoverage) {
        self.parameters = parameters; self.coverage = coverage
    }
    public func heatFlux(for query: StagnationHeatingQuery) throws -> SurfaceHeatFlux {
        guard query.gasIdentifier == parameters.gasIdentifier, query.wall == nil else { throw FluidError.unsupportedQuery }
        guard coverage.noseRadiusMetres.contains(query.noseRadiusMetres) else { throw FluidError.outsideCoverage }
        if query.flow.speedMetresPerSecond == 0 {
            return try SurfaceHeatFlux(wattsPerSquareMetre: 0, surfacePositionMetres: query.surfacePositionMetres)
        }
        guard let state = query.gasState, let mach = query.flow.machNumber,
              let kn = try query.flow.knudsenNumber(referenceLengthMetres: query.noseRadiusMetres) else {
            throw FluidError.missingFlowProperty
        }
        guard coverage.speedMetresPerSecond.contains(query.flow.speedMetresPerSecond),
              coverage.densityKilogramsPerCubicMetre.contains(query.flow.densityKilogramsPerCubicMetre),
              coverage.machNumber.contains(mach), kn <= coverage.maximumKnudsenNumber,
              (101.325...10_132_500).contains(state.pressurePascals),
              (2.3e6...116.2e6).contains(state.edgeSpecificEnthalpyJoulesPerKilogram) else { throw FluidError.outsideCoverage }
        let flux = parameters.enthalpyCoefficient*sqrt((state.pressurePascals/101325)/query.noseRadiusMetres)
            * state.edgeSpecificEnthalpyJoulesPerKilogram
        guard flux.isFinite else { throw FluidError.nonFiniteResult }
        return try SurfaceHeatFlux(wattsPerSquareMetre: flux, surfacePositionMetres: query.surfacePositionMetres)
    }
}
