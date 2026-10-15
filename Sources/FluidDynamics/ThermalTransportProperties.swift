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

