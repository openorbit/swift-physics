// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Aerodynamics

@Suite struct PlanetaryHeatingTests {
    func coverage() throws -> StagnationHeatingCoverage {
        try StagnationHeatingCoverage(speedMetresPerSecond: 1000...10000,
            densityKilogramsPerCubicMetre: 1e-6...0.1, noseRadiusMetres: 0.1...5, machNumber: 5...50)
    }
    @Test func marsAndVenusMixtures() throws {
        let mars = SuttonGravesHeating(parameters: .marsCarbonDioxideNitrogenColdWall, coverage: try coverage())
        let flow = try AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(3000,0,0),
            densityKilogramsPerCubicMetre: 0.01, speedOfSoundMetresPerSecond: 300, meanFreePathMetres: 1e-4)
        let query = try StagnationHeatingQuery(flow: flow, gasIdentifier: mars.parameters.gasIdentifier,
            noseRadiusMetres: 1, surfacePositionMetres: .zero)
        #expect(abs((try mars.heatFlux(for: query)).wattsPerSquareMetre-494100) < 1e-8)
        let venus = SuttonGravesMixtureParameters.venusCarbonDioxideNitrogen
        #expect(abs(venus.massFractions[.carbonDioxide]!-0.977434008541) < 1e-6)
        // Compare Eq42's transport approximation with table II's independent
        // boundary-layer computer fit for air; agreement is approximate, not exact.
        let air = try SuttonGravesMixtureParameters(gasIdentifier: "air", massFractions: [.oxygen: 0.232,.nitrogen: 0.768])
        #expect(abs(air.enthalpyCoefficient/0.1113-1) < 0.035)
        let heating = SuttonGravesMixtureHeating(parameters: venus, coverage: try coverage())
        let state = try StagnationGasState(pressurePascals: 101325, edgeSpecificEnthalpyJoulesPerKilogram: 10e6)
        let venusQuery = try StagnationHeatingQuery(flow: flow, gasIdentifier: venus.gasIdentifier,
            noseRadiusMetres: 1, surfacePositionMetres: SIMD3(1,0,0), gasState: state)
        #expect(abs((try heating.heatFlux(for: venusQuery)).wattsPerSquareMetre-venus.enthalpyCoefficient*10e6) < 1e-8)
        #expect(throws: AerodynamicError.missingFlowProperty) {
            try heating.heatFlux(for: StagnationHeatingQuery(flow: flow, gasIdentifier: venus.gasIdentifier,
                noseRadiusMetres: 1, surfacePositionMetres: .zero))
        }
        #expect(throws: AerodynamicError.invalidParameter) {
            try SuttonGravesMixtureParameters(gasIdentifier: "invalid", massFractions: [.nitrogen: 0.5])
        }
    }
    func query(v: Double, rho: Double, radius: Double, mixture: HydrogenHeliumHeatingMixture = .neptuneReference,
               temperature: Double = 141.2, wall: Double = 300, path: Double? = 1e-5) throws -> StagnationHeatingQuery {
        try StagnationHeatingQuery(flow: AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(v,0,0),
            densityKilogramsPerCubicMetre: rho, temperatureKelvin: temperature,
            speedOfSoundMetresPerSecond: 1000, meanFreePathMetres: path), gasIdentifier: mixture.gasIdentifier,
            noseRadiusMetres: radius, surfacePositionMetres: SIMD3(1,0,0), wall: AerothermalWallState(temperatureKelvin: wall))
    }
    @Test func hydrogenHeliumPublishedColdWallReferences() throws {
        struct References: Decodable { let rows: [[Double]] }
        let url = try #require(Bundle.module.url(forResource: "H2HeColdWallReferences", withExtension: "json"))
        let rows = try JSONDecoder().decode(References.self,from: Data(contentsOf: url)).rows
        #expect(rows.count == 180)
        let model = HydrogenHeliumColdWallHeating(mixture: .neptuneReference)
        for row in rows {
            let result = try model.evaluate(query(v: row[0],rho: row[1],radius: row[2]))
            // Published fitted correlation is approximate; compare with independent CFD data.
            #expect(abs(result.heatFlux.wattsPerSquareMetre/row[3]-1) < 0.20)
            #expect(abs(result.enthalpyCoefficient/row[4]-1) < 0.20)
        }
    }
    @Test func hydrogenHeliumCoverageAndCompositionScaling() throws {
        let model = HydrogenHeliumColdWallHeating(mixture: .neptuneReference)
        let reference = try model.evaluate(query(v: 19000,rho: 1e-4,radius: 0.4))
        for (mixture,scale) in [(HydrogenHeliumHeatingMixture.saturnApproximation,0.952),(.uranusApproximation,0.980)] {
            let result = try HydrogenHeliumColdWallHeating(mixture: mixture).evaluate(query(v: 19000,rho: 1e-4,radius: 0.4,mixture: mixture))
            #expect(abs(result.heatFlux.wattsPerSquareMetre/reference.heatFlux.wattsPerSquareMetre-scale) < 1e-12)
        }
        #expect(try model.evaluate(query(v: 28999,rho: 1e-4,radius: 0.4)).usesHighSpeedFit == false)
        #expect(try model.evaluate(query(v: 29000,rho: 1e-4,radius: 0.4)).usesHighSpeedFit)
        #expect(throws: AerodynamicError.outsideCoverage) { try model.evaluate(query(v: 45000,rho: 0.1,radius: 0.4)) }
        #expect(throws: AerodynamicError.outsideCoverage) { try model.evaluate(query(v: 7000,rho: 1e-5,radius: 0.4)) }
        #expect(throws: AerodynamicError.unsupportedQuery) { try model.evaluate(query(v: 19000,rho: 1e-4,radius: 0.4,temperature: 150)) }
        #expect(throws: AerodynamicError.unsupportedQuery) { try model.evaluate(query(v: 19000,rho: 1e-4,radius: 0.4,wall: 1000)) }
        #expect(throws: AerodynamicError.missingFlowProperty) { try model.evaluate(query(v: 19000,rho: 1e-4,radius: 0.4,path: nil)) }
    }
}
