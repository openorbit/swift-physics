// SPDX-License-Identifier: Apache-2.0
import Testing
import Aerodynamics

@Suite struct AerothermalHeatingTests {
    func model(parameters: SuttonGravesParameters = .earthAirColdWall) throws -> SuttonGravesHeating {
        SuttonGravesHeating(parameters: parameters, coverage: try StagnationHeatingCoverage(
            speedMetresPerSecond: 1000...10000, densityKilogramsPerCubicMetre: 1e-6...0.1,
            noseRadiusMetres: 0.1...10, machNumber: 5...50))
    }
    func query(speed: Double = 3535, density: Double = 3.1459e-4, radius: Double = 1,
               path: Double? = 1e-4, gas: String = "earth-air",
               wall: AerothermalWallState? = nil) throws -> StagnationHeatingQuery {
        try StagnationHeatingQuery(flow: AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(speed,0,0),
            densityKilogramsPerCubicMetre: density, speedOfSoundMetresPerSecond: 300, meanFreePathMetres: path),
            gasIdentifier: gas, noseRadiusMetres: radius, surfacePositionMetres: SIMD3(2,0,0), wall: wall)
    }
    @Test func nasaShuttleReferenceAndScaling() throws {
        let heating = try model()
        let base = try heating.heatFlux(for: query())
        // NASA TFAWS 2012 course: 13.6 W/cm², rounded. SI output is W/m².
        #expect(abs(base.wattsPerSquareMetre/10000-13.6) < 0.05)
        #expect(base.surfacePositionMetres == SIMD3(2,0,0))
        #expect(abs((try heating.heatFlux(for: query(speed: 7070))).wattsPerSquareMetre/base.wattsPerSquareMetre-8) < 1e-12)
        #expect(abs((try heating.heatFlux(for: query(density: 4*3.1459e-4))).wattsPerSquareMetre/base.wattsPerSquareMetre-2) < 1e-12)
        #expect(abs((try heating.heatFlux(for: query(radius: 4))).wattsPerSquareMetre/base.wattsPerSquareMetre-0.5) < 1e-12)
    }
    @Test func explicitGasSelectionAndCustomParameters() throws {
        let earth = try model()
        #expect(throws: AerodynamicError.unsupportedQuery) { try earth.heatFlux(for: query(gas: "titan-mixture")) }
        let parameters = try SuttonGravesParameters(gasIdentifier: "synthetic", coefficientSI: 2e-4,
            source: "Synthetic algebraic reference", version: "1")
        let custom = try model(parameters: parameters)
        // sqrt(0.01/1)=0.1; V³=27×10⁹; k=0.0002 -> 540000 W/m².
        #expect(abs((try custom.heatFlux(for: query(speed: 3000,density: 0.01,gas: "synthetic"))).wattsPerSquareMetre-540000) < 1e-8)
    }
    @Test func missingPropertiesCoverageAndWallState() throws {
        let heating = try model()
        #expect(throws: AerodynamicError.missingFlowProperty) { try heating.heatFlux(for: query(path: nil)) }
        #expect(throws: AerodynamicError.outsideCoverage) { try heating.heatFlux(for: query(path: 0.02)) }
        #expect(throws: AerodynamicError.outsideCoverage) { try heating.heatFlux(for: query(speed: 900)) }
        #expect(throws: AerodynamicError.outsideCoverage) { try heating.heatFlux(for: query(density: 0)) }
        #expect(throws: AerodynamicError.outsideCoverage) { try heating.heatFlux(for: query(radius: 20)) }
        #expect(throws: AerodynamicError.unsupportedQuery) {
            try heating.heatFlux(for: query(wall: AerothermalWallState(temperatureKelvin: 1000)))
        }
        #expect(try heating.heatFlux(for: query(speed: 0,path: nil)).wattsPerSquareMetre == 0)
        #expect(throws: AerodynamicError.invalidParameter) { try query(radius: 0) }
    }
}
