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
import FluidDynamics

/// State-local material data with matched temperature/absolute-pressure axes.
/// Composition belongs to the supplied identifier and dataset assumptions. The
/// rectangular domain must be supported single-phase liquid; no extrapolation.
public struct TabulatedLiquidModel: LiquidPropertyModel {
    public let identifier: String
    public let density: FluidScalarTable
    public let viscosity: FluidScalarTable
    public let soundSpeed: FluidScalarTable
    public let vaporPressure: FluidScalarTable
    public init(identifier: String, density: FluidScalarTable, viscosity: FluidScalarTable,
                soundSpeed: FluidScalarTable, vaporPressure: FluidScalarTable) throws {
        guard !identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              density.axes.map(\.name) == ["temperature","pressure"], density.axes.map(\.units) == ["K","Pa"],
              density.axes[0].coordinates.first! > 0, density.axes[1].coordinates.first! > 0,
              viscosity.axes == density.axes, soundSpeed.axes == density.axes, vaporPressure.axes == density.axes,
              density.valueUnits == "kg/m^3", viscosity.valueUnits == "Pa*s", soundSpeed.valueUnits == "m/s", vaporPressure.valueUnits == "Pa",
              [density,viscosity,soundSpeed].allSatisfy({ $0.values.allSatisfy({ $0 > 0 }) }), vaporPressure.values.allSatisfy({ $0 >= 0 }) else {
            throw FluidError.malformedData
        }
        let pressures = density.axes[1].coordinates
        for index in vaporPressure.values.indices {
            guard pressures[index%pressures.count] > vaporPressure.values[index] else { throw FluidError.outsideCoverage }
        }
        self.identifier = identifier; self.density = density; self.viscosity = viscosity
        self.soundSpeed = soundSpeed; self.vaporPressure = vaporPressure
    }
    public func properties(temperatureKelvin: Double, absolutePressurePascals: Double) throws -> LiquidProperties {
        let query = [temperatureKelvin,absolutePressurePascals]
        let vapor = try vaporPressure.value(at: query)
        guard absolutePressurePascals > vapor else { throw FluidError.outsideCoverage }
        let provenance = [density,viscosity,soundSpeed,vaporPressure].map { $0.source+" ("+$0.version+")" }.joined(separator: "; ")
        return try LiquidProperties(identifier: identifier,source: provenance,
            densityKilogramsPerCubicMetre: density.value(at: query),dynamicViscosityPascalSeconds: viscosity.value(at: query),
            speedOfSoundMetresPerSecond: soundSpeed.value(at: query),vaporPressurePascals: vapor)
    }
}
