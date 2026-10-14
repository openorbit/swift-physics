import Foundation
import Testing
@testable import Atmospheres
@Suite struct ProfileTests {
    @Test func isaReferenceAndBoundaries() throws {
        let model = ISAAtmosphere()
        let sea = try model.properties(at: 0)
        #expect(abs(sea.pressurePascals-101325) < 1e-8)
        #expect(abs(sea.densityKilogramsPerCubicMetre-1.225) < 1e-6)
        #expect(abs(sea.speedOfSoundMetresPerSecond!-340.294) < 0.001)
        for (h,t,p) in [(11000.0,216.65,22632.06),(20000,216.65,5474.889),
                        (32000,228.65,868.0187),(47000,270.65,110.9063),
                        (51000,270.65,66.93887),(71000,214.65,3.956420),(84852,186.946,0.3733836)] {
            let geometric = 6356766*h/(6356766-h)
            let sample = try model.properties(at: geometric)
            #expect(abs(sample.temperatureKelvin-t) < 1e-8)
            // Published tables use slightly different rounded gas constants.
            #expect(abs(sample.pressurePascals/p-1) < 1e-5)
            if h < 84852 {
                let below = try model.properties(at: geometric-0.001)
                let above = try model.properties(at: geometric+0.001)
                #expect(abs(above.pressurePascals/below.pressurePascals-1) < 1e-6)
            }
        }
        #expect(throws: AtmosphericModelError.altitudeOutsideCoverage) { try model.properties(at: .nan) }
        #expect(throws: AtmosphericModelError.altitudeOutsideCoverage) { try model.properties(at: 90000) }
    }
    @Test func catalogueAndGiantDescent() throws {
        for body in AtmosphericBody.allCases {
            let model = BodyAtmosphere.model(for: body)
            let zero = try model.properties(at: 0)
            let high = try model.properties(at: model.altitudeRangeMetres.upperBound)
            #expect(zero.pressurePascals > high.pressurePascals)
            #expect(high.densityKilogramsPerCubicMetre > 0)
            #expect(!model.source.isEmpty)
            if [.jupiter,.saturn,.uranus,.neptune].contains(body) {
                let deep = try model.properties(at: model.altitudeRangeMetres.lowerBound)
                #expect(abs(deep.pressurePascals/1e7-1) < 1e-12)
                #expect(deep.temperatureKelvin > zero.temperatureKelvin)
                #expect(deep.densityKilogramsPerCubicMetre > zero.densityKilogramsPerCubicMetre)
            }
        }
        #expect(BodyAtmosphere.model(naifBodyCode: 5) == nil)
        #expect(try BodyAtmosphere.model(for: .callisto).properties(at: 0).speedOfSoundMetresPerSecond == nil)
        #expect(try BodyAtmosphere.model(for: .mars).properties(at: 0).dynamicViscosityPascalSeconds == nil)
    }
    @Test func isothermalAnalyticAndValidation() throws {
        let model = try HydrostaticAtmosphere(levels: [.init(altitudeMetres: 0,temperatureKelvin: 300),
            .init(altitudeMetres: 1000,temperatureKelvin: 300)], referencePressurePascals: 100000,
            specificGasConstant: 300, gravityMetresPerSecondSquared: 10,
            altitudeDatum: "Test", source: "Synthetic", assumptions: "Test")
        #expect(abs(try model.properties(at: 1000).pressurePascals-100000*exp(-1.0/9)) < 1e-8)
        #expect(throws: AtmosphericModelError.invalidParameters) {
            try HydrostaticAtmosphere(levels: [], referencePressurePascals: -1,
                specificGasConstant: 300, gravityMetresPerSecondSquared: 10,
                altitudeDatum: "Test", source: "Test", assumptions: "")
        }
    }
}
