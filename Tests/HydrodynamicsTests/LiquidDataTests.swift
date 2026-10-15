// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import FluidDynamics
import Hydrodynamics

@Suite struct LiquidDataTests {
    func table(units: String, value: (Double,Double) -> Double) throws -> FluidScalarTable {
        let temperatures = [290.0,300], pressures = [1e5,1e6]
        return try FluidScalarTable(axes: [FluidTableAxis(name: "temperature",units: "K",coordinates: temperatures),
            FluidTableAxis(name: "pressure",units: "Pa",coordinates: pressures)],
            values: temperatures.flatMap { t in pressures.map { value(t,$0) } },valueUnits: units,
            source: "Synthetic affine liquid validation data",version: "1",assumptions: "Single synthetic composition; interpolation test only")
    }
    @Test func pressureTemperatureInterpolationAndPhaseCoverage() throws {
        let density = try table(units: "kg/m^3") { 1000-0.1*($0-290)+1e-8*($1-1e5) }
        let viscosity = try table(units: "Pa*s") { t,_ in 0.001*(1-0.01*(t-290)) }
        let sound = try table(units: "m/s") { t,_ in 1500+0.01*(t-290) }
        let vapor = try table(units: "Pa") { t,_ in 1000+10*(t-290) }
        let model: any LiquidPropertyModel = try TabulatedLiquidModel(identifier: "Synthetic liquid",density: density,
            viscosity: viscosity,soundSpeed: sound,vaporPressure: vapor)
        let p = try model.properties(temperatureKelvin: 295,absolutePressurePascals: 550000)
        #expect(abs(p.densityKilogramsPerCubicMetre-999.5045) < 1e-10)
        #expect(abs(p.dynamicViscosityPascalSeconds-0.00095) < 1e-15)
        #expect(abs(p.speedOfSoundMetresPerSecond-1500.05) < 1e-10)
        #expect(p.vaporPressurePascals == 1050)
        #expect(throws: FluidError.outsideCoverage) { try model.properties(temperatureKelvin: 94,absolutePressurePascals: 150000) }
        #expect(throws: FluidError.outsideCoverage) {
            try TabulatedLiquidModel(identifier: "bad-phase",density: density,viscosity: viscosity,soundSpeed: sound,
                vaporPressure: table(units: "Pa") { _,_ in 2e5 })
        }
        #expect(throws: FluidError.malformedData) {
            try TabulatedLiquidModel(identifier: "bad-units",density: table(units: "g/cm^3") { _,_ in 1 },
                viscosity: viscosity,soundSpeed: sound,vaporPressure: vapor)
        }
    }
}
