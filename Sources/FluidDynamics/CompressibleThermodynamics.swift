// SPDX-License-Identifier: Apache-2.0
import Foundation

public struct FluidThermodynamicState: Sendable {
    public let pressurePascals: Double
    public let temperatureKelvin: Double
    public let speedOfSoundMetresPerSecond: Double
    public init(pressurePascals: Double, temperatureKelvin: Double, speedOfSoundMetresPerSecond: Double) throws {
        guard [pressurePascals,temperatureKelvin,speedOfSoundMetresPerSecond].allSatisfy({ $0.isFinite && $0 > 0 }) else {
            throw FluidError.invalidFlow
        }
        self.pressurePascals = pressurePascals; self.temperatureKelvin = temperatureKelvin
        self.speedOfSoundMetresPerSecond = speedOfSoundMetresPerSecond
    }
}

/// Thermodynamic closure for conservative compressible flow. The energy reference
/// must yield positive internal energy in the supported domain; sound speed is the
/// thermodynamic acoustic speed. Chemistry/phase assumptions belong to the provider.
public protocol CompressibleEquationOfState: Sendable {
    var identifier: String { get }
    func thermodynamics(densityKilogramsPerCubicMetre: Double, specificInternalEnergyJoulesPerKilogram: Double) throws -> FluidThermodynamicState
    func specificInternalEnergyJoulesPerKilogram(densityKilogramsPerCubicMetre: Double, pressurePascals: Double) throws -> Double
}

/// Calorically perfect ideal gas with explicit temperature bounds and provenance.
/// p=(gamma-1) rho e; e=cv T. No dissociation, ionization or non-equilibrium chemistry.
public struct CaloricallyPerfectGasEOS: CompressibleEquationOfState {
    public let identifier: String
    public let specificGasConstantJoulesPerKilogramKelvin: Double
    public let specificHeatRatio: Double
    public let temperatureRangeKelvin: ClosedRange<Double>
    public let source: String
    public init(identifier: String, specificGasConstantJoulesPerKilogramKelvin: Double, specificHeatRatio: Double,
                temperatureRangeKelvin: ClosedRange<Double>, source: String) throws {
        guard !identifier.isEmpty, !source.isEmpty, specificGasConstantJoulesPerKilogramKelvin.isFinite,
              specificGasConstantJoulesPerKilogramKelvin > 0, specificHeatRatio.isFinite, specificHeatRatio > 1, specificHeatRatio <= 2,
              temperatureRangeKelvin.lowerBound.isFinite, temperatureRangeKelvin.lowerBound > 0, temperatureRangeKelvin.upperBound.isFinite else {
            throw FluidError.invalidParameter
        }
        self.identifier = identifier; self.specificGasConstantJoulesPerKilogramKelvin = specificGasConstantJoulesPerKilogramKelvin
        self.specificHeatRatio = specificHeatRatio; self.temperatureRangeKelvin = temperatureRangeKelvin; self.source = source
    }
    public func thermodynamics(densityKilogramsPerCubicMetre rho: Double, specificInternalEnergyJoulesPerKilogram e: Double) throws -> FluidThermodynamicState {
        guard rho.isFinite, rho > 0, e.isFinite, e > 0 else { throw FluidError.invalidFlow }
        let temperature = e*(specificHeatRatio-1)/specificGasConstantJoulesPerKilogramKelvin
        guard temperature.isFinite, temperatureRangeKelvin.contains(temperature) else { throw FluidError.outsideCoverage }
        return try FluidThermodynamicState(pressurePascals: (specificHeatRatio-1)*rho*e,temperatureKelvin: temperature,
            speedOfSoundMetresPerSecond: sqrt(specificHeatRatio*(specificHeatRatio-1)*e))
    }
    public func specificInternalEnergyJoulesPerKilogram(densityKilogramsPerCubicMetre rho: Double, pressurePascals p: Double) throws -> Double {
        guard rho.isFinite, rho > 0, p.isFinite, p > 0 else { throw FluidError.invalidFlow }
        let e = p/rho/(specificHeatRatio-1)
        _ = try thermodynamics(densityKilogramsPerCubicMetre: rho,specificInternalEnergyJoulesPerKilogram: e)
        return e
    }
}

public struct CompressiblePrimitiveState: Sendable {
    public let densityKilogramsPerCubicMetre: Double
    public let velocityMetresPerSecond: SIMD3<Double>
    public let thermodynamics: FluidThermodynamicState
}
public struct ConservedFluidState: Sendable, Equatable {
    public let densityKilogramsPerCubicMetre: Double
    public let momentumDensityKilogramsPerSquareMetreSecond: SIMD3<Double>
    public let totalEnergyDensityJoulesPerCubicMetre: Double
    public init(densityKilogramsPerCubicMetre: Double, momentumDensityKilogramsPerSquareMetreSecond: SIMD3<Double>,
                totalEnergyDensityJoulesPerCubicMetre: Double) throws {
        guard densityKilogramsPerCubicMetre.isFinite, densityKilogramsPerCubicMetre > 0,
              fluidFinite(momentumDensityKilogramsPerSquareMetreSecond), totalEnergyDensityJoulesPerCubicMetre.isFinite,
              totalEnergyDensityJoulesPerCubicMetre > 0 else { throw FluidError.invalidFlow }
        self.densityKilogramsPerCubicMetre = densityKilogramsPerCubicMetre
        self.momentumDensityKilogramsPerSquareMetreSecond = momentumDensityKilogramsPerSquareMetreSecond
        self.totalEnergyDensityJoulesPerCubicMetre = totalEnergyDensityJoulesPerCubicMetre
    }
    public init(densityKilogramsPerCubicMetre rho: Double, velocityMetresPerSecond v: SIMD3<Double>, pressurePascals: Double,
                equationOfState: any CompressibleEquationOfState) throws {
        guard fluidFinite(v) else { throw FluidError.invalidFlow }
        let e = try equationOfState.specificInternalEnergyJoulesPerKilogram(densityKilogramsPerCubicMetre: rho,pressurePascals: pressurePascals)
        try self.init(densityKilogramsPerCubicMetre: rho,momentumDensityKilogramsPerSquareMetreSecond: rho*v,
            totalEnergyDensityJoulesPerCubicMetre: rho*(e+0.5*(v*v).sum()))
        _ = try primitive(using: equationOfState)
    }
    public func primitive(using eos: any CompressibleEquationOfState) throws -> CompressiblePrimitiveState {
        let v = momentumDensityKilogramsPerSquareMetreSecond/densityKilogramsPerCubicMetre
        guard fluidFinite(v) else { throw FluidError.nonFiniteResult }
        let e = totalEnergyDensityJoulesPerCubicMetre/densityKilogramsPerCubicMetre-0.5*(v*v).sum()
        let thermo = try eos.thermodynamics(densityKilogramsPerCubicMetre: densityKilogramsPerCubicMetre,specificInternalEnergyJoulesPerKilogram: e)
        return CompressiblePrimitiveState(densityKilogramsPerCubicMetre: densityKilogramsPerCubicMetre,velocityMetresPerSecond: v,thermodynamics: thermo)
    }
    package var vector: SIMD8<Double> {
        SIMD8(densityKilogramsPerCubicMetre,momentumDensityKilogramsPerSquareMetreSecond.x,
              momentumDensityKilogramsPerSquareMetreSecond.y,momentumDensityKilogramsPerSquareMetreSecond.z,
              totalEnergyDensityJoulesPerCubicMetre,0,0,0)
    }
    package init(vector: SIMD8<Double>) throws {
        try self.init(densityKilogramsPerCubicMetre: vector[0],momentumDensityKilogramsPerSquareMetreSecond: SIMD3(vector[1],vector[2],vector[3]),
            totalEnergyDensityJoulesPerCubicMetre: vector[4])
    }
}

public extension CompressiblePrimitiveState {
    /// All velocities share axes/frame. Viscosity and molecular mean free path are
    /// transport inputs, not inferred from the pressure/energy closure.
    func relativeFlow(bodyVelocityMetresPerSecond: SIMD3<Double>, dynamicViscosityPascalSeconds: Double? = nil,
                      meanFreePathMetres: Double? = nil) throws -> FluidFlow {
        try FluidFlow(relativeVelocityMetresPerSecond: bodyVelocityMetresPerSecond-velocityMetresPerSecond,
            densityKilogramsPerCubicMetre: densityKilogramsPerCubicMetre,pressurePascals: thermodynamics.pressurePascals,
            temperatureKelvin: thermodynamics.temperatureKelvin,speedOfSoundMetresPerSecond: thermodynamics.speedOfSoundMetresPerSecond,
            dynamicViscosityPascalSeconds: dynamicViscosityPascalSeconds,meanFreePathMetres: meanFreePathMetres)
    }
}
