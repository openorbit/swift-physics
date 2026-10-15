// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Aerodynamics

@Suite struct CoverageDataTests {
    func axis(_ name: String, _ coordinates: [Double], units: String = "1") throws -> AerodynamicTableAxis {
        try AerodynamicTableAxis(name: name,units: units,coordinates: coordinates)
    }
    func table(_ axes: [AerodynamicTableAxis], units: String = "1", function: ([Double]) -> Double) throws -> AerodynamicScalarTable {
        var points: [[Double]] = [[]]
        for axis in axes { points = points.flatMap { p in axis.coordinates.map { p+[$0] } } }
        return try AerodynamicScalarTable(axes: axes,values: points.map(function),valueUnits: units,
            source: "Synthetic affine validation fixture",version: "1",assumptions: "Tests interpolation only; no physical calibration")
    }
    @Test func checkedGridInterpolationAndDecodedValidation() throws {
        let axes = try [axis("a",[1,2]),axis("b",[3,5,8]),axis("c",[7])]
        let t = try table(axes) { 2*$0[0]-3*$0[1]+$0[2] }
        #expect(abs(try t.value(at: [1.4,6,7])+8.2) < 1e-12)
        #expect(try t.value(at: [2,8,7]) == -13)
        #expect(throws: AerodynamicError.outsideCoverage) { try t.value(at: [0,6,7]) }
        #expect(throws: AerodynamicError.invalidParameter) { try t.value(at: [.nan,6,7]) }
        #expect(throws: AerodynamicError.self) { try axis("bad",[2,1]) }
        #expect(throws: AerodynamicError.self) { try axis("bad",[1,1]) }
        #expect(throws: AerodynamicError.self) { try axis("overflow",[-1e308,1e308]) }
        let decoded = try JSONDecoder().decode(AerodynamicScalarTable.self,from: JSONEncoder().encode(t))
        #expect(try decoded.value(at: [1.4,6,7]) == t.value(at: [1.4,6,7]))
        var object = try #require(JSONSerialization.jsonObject(with: JSONEncoder().encode(t)) as? [String:Any])
        object["values"] = [1.0]
        let invalid = try JSONSerialization.data(withJSONObject: object)
        #expect(throws: AerodynamicError.incompleteGrid) { try JSONDecoder().decode(AerodynamicScalarTable.self,from: invalid) }
    }
    @Test func transitionalTractionIsContinuousWithinSuppliedDataset() throws {
        let axes = try [axis("log10Kn",[-3,0,2]),axis("incidenceCosine",[-1,1]),
                        axis("speedRatio",[0,3]),axis("wallTemperatureRatio",[1,3])]
        let p = try table(axes) { 1+0.1*$0[0]+0.2*$0[1]+0.3*$0[2]+0.4*$0[3] }
        let shear = try table(axes) { _ in 0.5 }
        let gas = try FreeMolecularGas(identifier: "synthetic",specificGasConstantJoulesPerKilogramKelvin: 4,source: "test")
        let model = try TabulatedRarefiedTraction(gas: gas,wallTemperatureKelvin: 8,pressure: p,shear: shear)
        let panel = try AerodynamicSurfacePanel(areaSquareMetres: 2,centreMetres: SIMD3(0,0,1),outwardNormal: SIMD3(1,0,0))
        for logKn in [-3.0,-1,-1e-8,0,1e-8,1,2] {
            let flow = try AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(3,4,0),densityKilogramsPerCubicMetre: 2,
                temperatureKelvin: 4,meanFreePathMetres: pow(10,logKn))
            let t = try model.traction(panel: panel,flow: flow,referenceLengthMetres: 1)
            #expect(abs(t.gaugePressurePascals-32*(2.295+0.1*logKn)) < 1e-11)
            #expect(abs(t.shearPascals.y+16) < 1e-12)
            #expect(t.shearPascals.x == 0 && t.shearPascals.z == 0)
            let load = try SurfaceTractionAerodynamics(panels: [panel],referenceLengthMetres: 1,model: model).evaluate(flow: flow).loads
            #expect(abs(load.momentNewtonMetres.y+2*t.gaugePressurePascals) < 1e-10)
        }
        let outside = try AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(3,4,0),densityKilogramsPerCubicMetre: 2,
            temperatureKelvin: 4,meanFreePathMetres: 1000)
        #expect(throws: AerodynamicError.outsideCoverage) { try model.traction(panel: panel,flow: outside,referenceLengthMetres: 1) }
        #expect(throws: AerodynamicError.malformedData) {
            try TabulatedRarefiedTraction(gas: gas,wallTemperatureKelvin: 8,pressure: p,shear: table(axes) { _ in -1 })
        }
    }
    @Test func explicitMachDispatchAcrossAllFlightRegimes() throws {
        let bounds = [0.0...0.8,0.8...1.2,1.2...5,5...25]
        let models = try bounds.enumerated().map { index, bound in
            try CoveredBodyCoefficients(model: ConstantBodyCoefficients(value: AerodynamicCoefficients(
                force: SIMD3(-Double(index+1),0,0),moment: .zero),source: "Synthetic regime test",version: "1"),
                coverage: AerodynamicFlowCoverage(machNumber: bound,reynoldsNumber: 1...1e10,
                    angleOfAttackRadians: -.pi...(.pi),sideslipRadians: -1...1,knudsenNumber: 0...0.01),
                reynoldsReferenceLengthMetres: 1,momentLengthsMetres: SIMD3(repeating: 1))
        }
        let model = try MachRegimeBodyCoefficients(models: models,machThresholds: [0.8,1.2,5],source: "test dispatch",version: "1")
        for (mach, expected) in [(0.3,1.0),(0.8,2),(1.2,3),(5,4),(20,4)] {
            let flow = try AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(mach*100,0,0),densityKilogramsPerCubicMetre: 1,
                speedOfSoundMetresPerSecond: 100,dynamicViscosityPascalSeconds: 1e-5,meanFreePathMetres: 1e-7)
            let query = try BodyAerodynamicQuery(flow: flow,reynoldsReferenceLengthMetres: 1,momentLengthsMetres: SIMD3(repeating: 1))
            #expect(try model.coefficients(for: query).force.x == -expected)
            #expect(throws: AerodynamicError.unsupportedQuery) {
                try model.coefficients(for: BodyAerodynamicQuery(flow: flow,reynoldsReferenceLengthMetres: 2,momentLengthsMetres: SIMD3(repeating: 1)))
            }
        }
    }
    @Test func explicitCompositionAndTransport() throws {
        let a = try PerfectGasSpecies(identifier: "a",molarMassKilogramsPerMole: 8.31446261815324/200,
            specificHeatJoulesPerKilogramKelvin: 1000,temperatureRangeKelvin: 200...500,source: "synthetic")
        let b = try PerfectGasSpecies(identifier: "b",molarMassKilogramsPerMole: 8.31446261815324/400,
            specificHeatJoulesPerKilogramKelvin: 1500,temperatureRangeKelvin: 250...600,source: "synthetic")
        let components = try [PerfectGasMassFraction(species: a,massFraction: 0.25),PerfectGasMassFraction(species: b,massFraction: 0.75)]
        let mixture = try PerfectGasMixture(identifier: "mix",components: components,source: "synthetic")
        #expect(abs(mixture.specificGasConstantJoulesPerKilogramKelvin-350) < 1e-12)
        #expect(mixture.specificHeatJoulesPerKilogramKelvin == 1375)
        #expect(mixture.temperatureRangeKelvin == 250...500)
        let density = try mixture.densityKilogramsPerCubicMetre(temperatureKelvin: 300,pressurePascals: 105000)
        #expect(abs(density-1) < 1e-12)
        #expect(abs(try mixture.speedOfSoundMetresPerSecond(temperatureKelvin: 300)-sqrt(1375.0/1025*350*300)) < 1e-10)
        let transport = try mixture.transport(temperatureKelvin: 300,pressurePascals: 105000,
            dynamicViscosityPascalSeconds: 2e-5,thermalConductivityWattsPerMetreKelvin: 0.03,transportSource: "supplied")
        #expect(transport.gasIdentifier == "mix")
        #expect(abs(transport.prandtlNumber-1375*2e-5/0.03) < 1e-12)
        #expect(throws: AerodynamicError.outsideCoverage) { try mixture.speedOfSoundMetresPerSecond(temperatureKelvin: 800) }
        #expect(throws: AerodynamicError.invalidParameter) { try PerfectGasMixture(identifier: "bad",components: [components[0]],source: "test") }
    }
}
