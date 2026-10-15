// SPDX-License-Identifier: Apache-2.0
import Foundation
/// Caller-supplied fluid properties at the film temperature (edge + wall)/2.
/// Composition and property evaluation remain with the gas-property provider.
public struct ThermalTransportProperties: Sendable {
    public let fluidIdentifier: String
    public let source: String
    public let temperatureKelvin: Double
    public let densityKilogramsPerCubicMetre: Double
    public let dynamicViscosityPascalSeconds: Double
    public let thermalConductivityWattsPerMetreKelvin: Double
    public let specificHeatJoulesPerKilogramKelvin: Double
    public var prandtlNumber: Double {
        (dynamicViscosityPascalSeconds/thermalConductivityWattsPerMetreKelvin)*specificHeatJoulesPerKilogramKelvin
    }
    public init(fluidIdentifier: String, source: String, temperatureKelvin: Double,
                densityKilogramsPerCubicMetre: Double, dynamicViscosityPascalSeconds: Double,
                thermalConductivityWattsPerMetreKelvin: Double, specificHeatJoulesPerKilogramKelvin: Double) throws {
        guard [fluidIdentifier,source].allSatisfy({ !$0.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty }),
              [temperatureKelvin,densityKilogramsPerCubicMetre,dynamicViscosityPascalSeconds,
               thermalConductivityWattsPerMetreKelvin,specificHeatJoulesPerKilogramKelvin]
                .allSatisfy({ $0.isFinite && $0 > 0 }) else { throw FluidError.invalidParameter }
        self.fluidIdentifier = fluidIdentifier; self.source = source; self.temperatureKelvin = temperatureKelvin
        self.densityKilogramsPerCubicMetre = densityKilogramsPerCubicMetre
        self.dynamicViscosityPascalSeconds = dynamicViscosityPascalSeconds
        self.thermalConductivityWattsPerMetreKelvin = thermalConductivityWattsPerMetreKelvin
        self.specificHeatJoulesPerKilogramKelvin = specificHeatJoulesPerKilogramKelvin
        guard prandtlNumber.isFinite, prandtlNumber > 0 else { throw FluidError.nonFiniteResult }
    }
}

