// SPDX-License-Identifier: Apache-2.0
import FluidDynamics
import Foundation

/// Caller-supplied wall state for correlations which require it. No thermal evolution.
public struct AerothermalWallState: Sendable {
    public let temperatureKelvin: Double
    public init(temperatureKelvin: Double) throws {
        guard temperatureKelvin.isFinite, temperatureKelvin > 0 else { throw FluidError.invalidParameter }
        self.temperatureKelvin = temperatureKelvin
    }
}

/// Flow at a caller-identified stagnation point of an axisymmetric blunt body.
/// Gas identifier names the parameter dataset/composition, not an astronomical body.
/// Richer gas/surface state can be carried by a replacement model's own provider.
public struct StagnationHeatingQuery: Sendable {
    public let flow: FluidFlow
    public let gasIdentifier: String
    public let noseRadiusMetres: Double
    public let surfacePositionMetres: SIMD3<Double>
    public let wall: AerothermalWallState?
    public let gasState: StagnationGasState?
    public init(flow: FluidFlow, gasIdentifier: String, noseRadiusMetres: Double,
                surfacePositionMetres: SIMD3<Double>, wall: AerothermalWallState? = nil,
                gasState: StagnationGasState? = nil) throws {
        guard !gasIdentifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              noseRadiusMetres.isFinite, noseRadiusMetres > 0, fluidFinite(surfacePositionMetres) else {
            throw FluidError.invalidParameter
        }
        self.flow = flow; self.gasIdentifier = gasIdentifier; self.noseRadiusMetres = noseRadiusMetres
        self.surfacePositionMetres = surfacePositionMetres; self.wall = wall
        self.gasState = gasState
    }
}

/// Incoming convective flux only; wall response and emitted radiation are external.
public protocol ConvectiveHeatingModel: Sendable {
    func heatFlux(for query: StagnationHeatingQuery) throws -> SurfaceHeatFlux
}
/// Incoming shock-layer radiative flux only. Absence of a model means unknown,
/// not zero radiation; this interface does not provide a radiative correlation.
public protocol RadiativeHeatingModel: Sendable {
    func heatFlux(for query: StagnationHeatingQuery) throws -> SurfaceHeatFlux
}

/// Parameter k for q = k sqrt(ρ/Rn) V³ with kg/m³, metres, m/s and W/m².
/// k has SI dimensions kg^0.5/m. Gas and provenance are mandatory.
public struct SuttonGravesParameters: Sendable {
    public let gasIdentifier: String
    public let coefficientSI: Double
    public let source: String
    public let version: String
    public init(gasIdentifier: String, coefficientSI: Double, source: String, version: String) throws {
        guard [gasIdentifier,source,version].allSatisfy({ !$0.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty }),
              coefficientSI.isFinite, coefficientSI > 0 else { throw FluidError.invalidParameter }
        self.gasIdentifier = gasIdentifier; self.coefficientSI = coefficientSI; self.source = source; self.version = version
    }
    /// Explicit Earth-air cold-wall approximation from NASA's 2012 TFAWS course.
    /// Selection does not establish accuracy or coverage for a given trajectory.
    public static var earthAirColdWall: Self {
        try! Self(gasIdentifier: "earth-air", coefficientSI: 1.7415e-4,
            source: "https://tfaws.nasa.gov/TFAWS12/Proceedings/Aerothermodynamics%20Course.pdf",
            version: "TFAWS 2012")
    }
    /// Published cold-wall approximation for a 97% CO₂ / 3% N₂ Mars mixture by mass.
    public static var marsCarbonDioxideNitrogenColdWall: Self {
        try! Self(gasIdentifier: "mars-co2-n2-97-3-mass", coefficientSI: 1.83e-4,
            source: "https://ntrs.nasa.gov/api/citations/20200002354/downloads/20200002354.pdf",
            version: "West and Brandis 2020, equation 2")
    }
}

/// Supplied stagnation pressure and boundary-layer-edge enthalpy for mixture
/// correlations. Neither is inferred from an ideal-gas atmosphere profile.
public struct StagnationGasState: Sendable {
    public let pressurePascals: Double
    public let edgeSpecificEnthalpyJoulesPerKilogram: Double
    public init(pressurePascals: Double, edgeSpecificEnthalpyJoulesPerKilogram: Double) throws {
        guard pressurePascals.isFinite, pressurePascals > 0,
              edgeSpecificEnthalpyJoulesPerKilogram.isFinite, edgeSpecificEnthalpyJoulesPerKilogram > 0 else {
            throw FluidError.invalidParameter
        }
        self.pressurePascals = pressurePascals
        self.edgeSpecificEnthalpyJoulesPerKilogram = edgeSpecificEnthalpyJoulesPerKilogram
    }
}

/// Caller-declared applicability bounds, not a calibrated accuracy guarantee.
public struct StagnationHeatingCoverage: Sendable {
    public let speedMetresPerSecond: ClosedRange<Double>
    public let densityKilogramsPerCubicMetre: ClosedRange<Double>
    public let noseRadiusMetres: ClosedRange<Double>
    public let machNumber: ClosedRange<Double>
    public let maximumKnudsenNumber: Double
    public init(speedMetresPerSecond: ClosedRange<Double>, densityKilogramsPerCubicMetre: ClosedRange<Double>,
                noseRadiusMetres: ClosedRange<Double>, machNumber: ClosedRange<Double>,
                maximumKnudsenNumber: Double = 0.01) throws {
        guard [speedMetresPerSecond,densityKilogramsPerCubicMetre,noseRadiusMetres,machNumber]
            .allSatisfy({ $0.lowerBound.isFinite && $0.lowerBound > 0 && $0.upperBound.isFinite }),
              machNumber.lowerBound >= 5, maximumKnudsenNumber.isFinite,
              maximumKnudsenNumber > 0, maximumKnudsenNumber <= 0.01 else { throw FluidError.invalidParameter }
        self.speedMetresPerSecond = speedMetresPerSecond; self.densityKilogramsPerCubicMetre = densityKilogramsPerCubicMetre
        self.noseRadiusMetres = noseRadiusMetres; self.machNumber = machNumber; self.maximumKnudsenNumber = maximumKnudsenNumber
    }
}

/// Laminar cold-wall stagnation-point estimate, not a surface heating distribution.
/// No wall-temperature correction, catalytic-state correction or turbulent multiplier.
public struct SuttonGravesHeating: ConvectiveHeatingModel {
    public let parameters: SuttonGravesParameters
    public let coverage: StagnationHeatingCoverage
    public init(parameters: SuttonGravesParameters, coverage: StagnationHeatingCoverage) {
        self.parameters = parameters; self.coverage = coverage
    }
    public func heatFlux(for query: StagnationHeatingQuery) throws -> SurfaceHeatFlux {
        guard query.gasIdentifier == parameters.gasIdentifier else { throw FluidError.unsupportedQuery }
        // A supplied finite-temperature wall needs a different correlation.
        guard query.wall == nil else { throw FluidError.unsupportedQuery }
        let flow = query.flow
        guard coverage.noseRadiusMetres.contains(query.noseRadiusMetres) else { throw FluidError.outsideCoverage }
        if flow.speedMetresPerSecond == 0 {
            return try SurfaceHeatFlux(wattsPerSquareMetre: 0, surfacePositionMetres: query.surfacePositionMetres)
        }
        guard let mach = flow.machNumber,
              let kn = try flow.knudsenNumber(referenceLengthMetres: query.noseRadiusMetres) else {
            throw FluidError.missingFlowProperty
        }
        guard coverage.speedMetresPerSecond.contains(flow.speedMetresPerSecond),
              coverage.densityKilogramsPerCubicMetre.contains(flow.densityKilogramsPerCubicMetre),
              coverage.machNumber.contains(mach), kn <= coverage.maximumKnudsenNumber else { throw FluidError.outsideCoverage }
        let v = flow.speedMetresPerSecond
        let flux = ((parameters.coefficientSI*sqrt(flow.densityKilogramsPerCubicMetre)/sqrt(query.noseRadiusMetres))*v)*v*v
        guard flux.isFinite else { throw FluidError.nonFiniteResult }
        return try SurfaceHeatFlux(wattsPerSquareMetre: flux, surfacePositionMetres: query.surfacePositionMetres)
    }
}
