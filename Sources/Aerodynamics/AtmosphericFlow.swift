// SPDX-License-Identifier: Apache-2.0
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
