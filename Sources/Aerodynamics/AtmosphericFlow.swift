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
import Atmospheres

public extension FluidFlow {
    init(relativeVelocityMetresPerSecond: SIMD3<Double>, properties: AtmosphericProperties,
                meanFreePathMetres: Double? = nil) throws {
        try self.init(relativeVelocityMetresPerSecond: relativeVelocityMetresPerSecond,
            densityKilogramsPerCubicMetre: properties.densityKilogramsPerCubicMetre,
            pressurePascals: properties.pressurePascals,temperatureKelvin: properties.temperatureKelvin,
            speedOfSoundMetresPerSecond: properties.speedOfSoundMetresPerSecond,
            dynamicViscosityPascalSeconds: properties.dynamicViscosityPascalSeconds,meanFreePathMetres: meanFreePathMetres)
    }
}
