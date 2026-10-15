// SPDX-License-Identifier: Apache-2.0
import FluidDynamics
import Foundation
import Testing
import Aerodynamics

@Suite struct EntrySurfaceHeatingTests {
    let data = CoverageDataTests()
    func flow(speed: Double = 2000, density: Double = 0.001, path: Double = 1e-6) throws -> FluidFlow {
        try FluidFlow(relativeVelocityMetresPerSecond: SIMD3(speed,0,0),densityKilogramsPerCubicMetre: density,
            pressurePascals: 100,temperatureKelvin: 400,speedOfSoundMetresPerSecond: 400,dynamicViscosityPascalSeconds: 2e-5,meanFreePathMetres: path)
    }
    func coverage() throws -> StagnationHeatingCoverage {
        try StagnationHeatingCoverage(speedMetresPerSecond: 1000...5000,densityKilogramsPerCubicMetre: 1e-4...0.01,
            noseRadiusMetres: 0.01...1,machNumber: 5...20)
    }
    @Test func suppliedEnthalpyDataPreservesRecoveryAndWallCooling() throws {
        let axes = try [data.axis("mach",[5,10]),data.axis("log10Re",[4,7]),data.axis("wallEnthalpyRatio",[0,1,2])]
        let table = try data.table(axes) { 0.001+0.0001*$0[0]+0.00001*$0[1]+0.0002*$0[2] }
        let model = try TabulatedEnthalpySurfaceHeating(gasIdentifier: "test",boundaryLayer: .turbulent,stanton: table,coverage: EnthalpySurfaceHeatingCoverage(
            edgeTemperatureKelvin: 200...1000,edgePressurePascals: 1...1000,
            adiabaticWallEnthalpyJoulesPerKilogram: 1e6...3e6,distanceFromLeadingEdgeMetres: 0.1...10))
        for wall in [0.0,1e6,2e6,3e6] {
            let result = try model.evaluate(EnthalpySurfaceHeatingQuery(edgeFlow: flow(),gasIdentifier: "test",
                distanceFromLeadingEdgeMetres: 1,adiabaticWallEnthalpyJoulesPerKilogram: 2e6,
                wallEnthalpyJoulesPerKilogram: wall,surfacePositionMetres: SIMD3(1,2,3)))
            let st = 0.00155+0.0002*wall/2e6
            #expect(abs(result.stantonNumber-st) < 1e-14)
            #expect(abs(result.wattsPerSquareMetre-2*st*(2e6-wall)) < 1e-8)
        }
        #expect(throws: FluidError.unsupportedQuery) {
            try model.evaluate(EnthalpySurfaceHeatingQuery(edgeFlow: flow(),gasIdentifier: "other",distanceFromLeadingEdgeMetres: 1,
                adiabaticWallEnthalpyJoulesPerKilogram: 2e6,wallEnthalpyJoulesPerKilogram: 1e6,surfacePositionMetres: .zero))
        }
        #expect(throws: FluidError.outsideCoverage) {
            try model.evaluate(EnthalpySurfaceHeatingQuery(edgeFlow: flow(),gasIdentifier: "test",distanceFromLeadingEdgeMetres: 1,
                adiabaticWallEnthalpyJoulesPerKilogram: 1e8,wallEnthalpyJoulesPerKilogram: 1e6,surfacePositionMetres: .zero))
        }
        #expect(throws: FluidError.outsideCoverage) {
            try model.evaluate(EnthalpySurfaceHeatingQuery(edgeFlow: flow(path: 1),gasIdentifier: "test",distanceFromLeadingEdgeMetres: 1,
                adiabaticWallEnthalpyJoulesPerKilogram: 2e6,wallEnthalpyJoulesPerKilogram: 1e6,surfacePositionMetres: .zero))
        }
    }
    struct SphereReference: ConvectiveHeatingModel {
        func heatFlux(for query: StagnationHeatingQuery) throws -> SurfaceHeatFlux {
            guard query.noseRadiusMetres == 0.3048, query.gasIdentifier == "earth-air", query.wall == nil else { throw FluidError.unsupportedQuery }
            return try SurfaceHeatFlux(wattsPerSquareMetre: 1000,surfacePositionMetres: query.surfacePositionMetres)
        }
    }
    @Test func sweptLeadingEdgeReferenceRadiusAndSweep() throws {
        let model = try SweptCylinderStagnationHeating(referenceSphereModel: SphereReference(),referenceSphereRadiusMetres: 0.3048,coverage: coverage())
        let query = try StagnationHeatingQuery(flow: flow(),gasIdentifier: "earth-air",noseRadiusMetres: 0.0762,surfacePositionMetres: SIMD3(1,2,3))
        let pluggable: any LeadingEdgeConvectiveHeatingModel = model
        let unswept = try pluggable.heatFlux(for: LeadingEdgeHeatingQuery(stagnation: query,effectiveSweepRadians: 0))
        #expect(abs(unswept.wattsPerSquareMetre-1414.213562373095) < 1e-10)
        let swept = try model.heatFlux(for: query,effectiveSweepRadians: .pi/3)
        #expect(abs(swept.wattsPerSquareMetre-659.753955386447) < 1e-9)
        #expect(try model.heatFlux(for: query,effectiveSweepRadians: -.pi/3).wattsPerSquareMetre == swept.wattsPerSquareMetre)
        #expect(swept.surfacePositionMetres == SIMD3(1,2,3))
        #expect(throws: FluidError.outsideCoverage) { try model.heatFlux(for: query,effectiveSweepRadians: 1.1) }
        #expect(throws: FluidError.invalidGeometry) { try LeadingEdgeHeatingQuery(stagnation: query,effectiveSweepRadians: .pi/2) }
        let wrongGas = try StagnationHeatingQuery(flow: flow(),gasIdentifier: "other",noseRadiusMetres: 0.0762,surfacePositionMetres: .zero)
        #expect(throws: FluidError.unsupportedQuery) { try model.heatFlux(for: wrongGas,effectiveSweepRadians: 0) }
        let rarefied = try StagnationHeatingQuery(flow: flow(path: 1),gasIdentifier: "earth-air",noseRadiusMetres: 0.0762,surfacePositionMetres: .zero)
        #expect(throws: FluidError.outsideCoverage) { try model.heatFlux(for: rarefied,effectiveSweepRadians: 0) }
    }
    @Test func suppliedRadiationIsIndependentIncomingFlux() throws {
        let axes = try [data.axis("speed",[1000,5000],units: "m/s"),data.axis("density",[1e-4,0.01],units: "kg/m^3"),
                        data.axis("radius",[0.01,1],units: "m")]
        let table = try data.table(axes,units: "W/m^2") { 100*$0[0]+1000*$0[1]+10*$0[2] }
        let model = try TabulatedRadiativeHeating(gasIdentifier: "test",flux: table,coverage: coverage(),freestreamTemperatureKelvin: 300...500)
        let query = try StagnationHeatingQuery(flow: flow(),gasIdentifier: "test",noseRadiusMetres: 0.5,surfacePositionMetres: .zero)
        #expect(abs(try model.heatFlux(for: query).wattsPerSquareMetre-200006) < 1e-7)
        let hotFlow = try FluidFlow(relativeVelocityMetresPerSecond: SIMD3(2000,0,0),densityKilogramsPerCubicMetre: 0.001,
            temperatureKelvin: 600,speedOfSoundMetresPerSecond: 400,meanFreePathMetres: 1e-6)
        let hotQuery = try StagnationHeatingQuery(flow: hotFlow,gasIdentifier: "test",noseRadiusMetres: 0.5,surfacePositionMetres: .zero)
        #expect(throws: FluidError.outsideCoverage) { try model.heatFlux(for: hotQuery) }
        let wallQuery = try StagnationHeatingQuery(flow: flow(),gasIdentifier: "test",noseRadiusMetres: 0.5,
            surfacePositionMetres: .zero,wall: AerothermalWallState(temperatureKelvin: 300))
        #expect(throws: FluidError.unsupportedQuery) { try model.heatFlux(for: wallQuery) }
    }
    @Test func fullyTurbulentFlatPlateAndIntegratedDistribution() throws {
        let model = TurbulentFlatPlateHeating()
        let fixture = SurfaceHeatingTests()
        let end = try model.evaluate(fixture.query(x: 2))
        // At Pr=1, Re=2e6. Numerical constant reference from local turbulent correlation.
        #expect(abs(end.localNusseltNumber-3251.739208186111) < 1e-6)
        #expect(end.wattsPerSquareMetre > 0)
        var lastError = Double.infinity
        let exactIntegral = end.wattsPerSquareMetre*2/pow(2,0.8)*(pow(2,0.8)-1)/0.8
        for count in [16,64] {
            var integral = 0.0
            for i in 0..<count {
                let x = 1+(Double(i)+0.5)/Double(count)
                integral += try model.evaluate(fixture.query(x: x)).wattsPerSquareMetre/Double(count)
            }
            let error = abs(integral-exactIntegral)
            #expect(error < lastError); lastError = error
        }
        #expect(lastError/exactIntegral < 1e-5)

        let start = try model.evaluate(fixture.query(x: 1))
        #expect(abs(start.wattsPerSquareMetre/end.wattsPerSquareMetre-pow(2,0.2)) < 1e-12)
        #expect(try model.evaluate(fixture.query(x: 2,wall: 330,edge: 300)).wattsPerSquareMetre < 0)
        #expect(throws: FluidError.outsideCoverage) { try model.evaluate(fixture.query(x: 0.1)) }
    }
}
