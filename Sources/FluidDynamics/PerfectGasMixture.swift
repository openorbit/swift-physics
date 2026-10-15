// SPDX-License-Identifier: Apache-2.0
import Foundation

/// Calorically perfect species over a caller-declared temperature interval.
/// No temperature-dependent cp, dissociation, ionization or transport mixing law.
public struct PerfectGasSpecies: Sendable {
    public let identifier: String
    public let molarMassKilogramsPerMole: Double
    public let specificHeatJoulesPerKilogramKelvin: Double
    public let temperatureRangeKelvin: ClosedRange<Double>
    public let source: String
    public var specificGasConstantJoulesPerKilogramKelvin: Double { 8.31446261815324/molarMassKilogramsPerMole }
    public init(identifier: String, molarMassKilogramsPerMole: Double,
                specificHeatJoulesPerKilogramKelvin: Double, temperatureRangeKelvin: ClosedRange<Double>, source: String) throws {
        guard !identifier.isEmpty, !source.isEmpty, molarMassKilogramsPerMole.isFinite, molarMassKilogramsPerMole > 0,
              specificHeatJoulesPerKilogramKelvin.isFinite,
              specificHeatJoulesPerKilogramKelvin > 8.31446261815324/molarMassKilogramsPerMole,
              temperatureRangeKelvin.lowerBound.isFinite, temperatureRangeKelvin.lowerBound > 0,
              temperatureRangeKelvin.upperBound.isFinite else { throw FluidError.invalidParameter }
        self.identifier = identifier; self.molarMassKilogramsPerMole = molarMassKilogramsPerMole
        self.specificHeatJoulesPerKilogramKelvin = specificHeatJoulesPerKilogramKelvin
        self.temperatureRangeKelvin = temperatureRangeKelvin; self.source = source
    }
}
public struct PerfectGasMassFraction: Sendable {
    public let species: PerfectGasSpecies
    public let massFraction: Double
    public init(species: PerfectGasSpecies, massFraction: Double) throws {
        guard massFraction.isFinite, massFraction > 0, massFraction <= 1 else { throw FluidError.invalidParameter }
        self.species = species; self.massFraction = massFraction
    }
}

/// Explicit composition, mass-weighted cp and R. Transport remains separately
/// supplied at the evaluation temperature; Earth-air defaults are never inferred.
public struct PerfectGasMixture: Sendable {
    public let identifier: String
    public let components: [PerfectGasMassFraction]
    public let source: String
    public let temperatureRangeKelvin: ClosedRange<Double>
    public let specificGasConstantJoulesPerKilogramKelvin: Double
    public let specificHeatJoulesPerKilogramKelvin: Double
    public var specificHeatRatio: Double {
        specificHeatJoulesPerKilogramKelvin/(specificHeatJoulesPerKilogramKelvin-specificGasConstantJoulesPerKilogramKelvin)
    }
    public init(identifier: String, components: [PerfectGasMassFraction], source: String) throws {
        let total = components.reduce(0) { $0+$1.massFraction }
        guard !identifier.isEmpty, !source.isEmpty, !components.isEmpty, abs(total-1) <= 1e-10,
              Set(components.map { $0.species.identifier }).count == components.count else { throw FluidError.invalidParameter }
        let lower = components.map { $0.species.temperatureRangeKelvin.lowerBound }.max()!
        let upper = components.map { $0.species.temperatureRangeKelvin.upperBound }.min()!
        guard lower <= upper else { throw FluidError.outsideCoverage }
        let r = components.reduce(0) { $0+$1.massFraction/total*$1.species.specificGasConstantJoulesPerKilogramKelvin }
        let cp = components.reduce(0) { $0+$1.massFraction/total*$1.species.specificHeatJoulesPerKilogramKelvin }
        guard r.isFinite, cp.isFinite, cp > r else { throw FluidError.nonFiniteResult }
        self.identifier = identifier; self.components = try components.map { try PerfectGasMassFraction(species: $0.species,massFraction: $0.massFraction/total) }
        self.source = source; self.temperatureRangeKelvin = lower...upper
        self.specificGasConstantJoulesPerKilogramKelvin = r; self.specificHeatJoulesPerKilogramKelvin = cp
    }
    public func densityKilogramsPerCubicMetre(temperatureKelvin: Double, pressurePascals: Double) throws -> Double {
        guard temperatureKelvin.isFinite, temperatureRangeKelvin.contains(temperatureKelvin) else { throw FluidError.outsideCoverage }
        guard pressurePascals.isFinite, pressurePascals >= 0 else { throw FluidError.invalidParameter }
        let density = pressurePascals/(specificGasConstantJoulesPerKilogramKelvin*temperatureKelvin)
        guard density.isFinite else { throw FluidError.nonFiniteResult }
        return density
    }
    public func speedOfSoundMetresPerSecond(temperatureKelvin: Double) throws -> Double {
        guard temperatureKelvin.isFinite, temperatureRangeKelvin.contains(temperatureKelvin) else { throw FluidError.outsideCoverage }
        let speed = sqrt(specificHeatRatio*specificGasConstantJoulesPerKilogramKelvin*temperatureKelvin)
        guard speed.isFinite else { throw FluidError.nonFiniteResult }
        return speed
    }
    /// All supplied transport values must describe this mixture at this temperature.
    public func transport(temperatureKelvin: Double, pressurePascals: Double,
                          dynamicViscosityPascalSeconds: Double, thermalConductivityWattsPerMetreKelvin: Double,
                          transportSource: String) throws -> ThermalTransportProperties {
        guard !transportSource.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else { throw FluidError.invalidParameter }
        return try ThermalTransportProperties(fluidIdentifier: identifier,source: source+"; transport: "+transportSource,
            temperatureKelvin: temperatureKelvin,
            densityKilogramsPerCubicMetre: densityKilogramsPerCubicMetre(temperatureKelvin: temperatureKelvin,pressurePascals: pressurePascals),
            dynamicViscosityPascalSeconds: dynamicViscosityPascalSeconds,
            thermalConductivityWattsPerMetreKelvin: thermalConductivityWattsPerMetreKelvin,
            specificHeatJoulesPerKilogramKelvin: specificHeatJoulesPerKilogramKelvin)
    }
}
