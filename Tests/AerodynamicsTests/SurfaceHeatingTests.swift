// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Aerodynamics

@Suite struct SurfaceHeatingTests {
    func query(x: Double = 0.1, speed: Double = 10, wall: Double = 300, edge: Double = 330,
               pr: Double = 1, sound: Double? = 340, path: Double? = 1e-7,
               film: Double? = nil) throws -> SurfaceHeatingQuery {
        try SurfaceHeatingQuery(edgeFlow: AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(speed,0,0),
            densityKilogramsPerCubicMetre: 2,temperatureKelvin: edge,
            speedOfSoundMetresPerSecond: sound,meanFreePathMetres: path),
            distanceFromLeadingEdgeMetres: x,surfacePositionMetres: SIMD3(x,2,3),
            wall: AerothermalWallState(temperatureKelvin: wall),
            transport: AerothermalTransportProperties(gasIdentifier: "synthetic",source: "test constant properties",
                temperatureKelvin: film ?? (wall+edge)/2,densityKilogramsPerCubicMetre: 1,
                dynamicViscosityPascalSeconds: 1e-5,thermalConductivityWattsPerMetreKelvin: 0.01,
                specificHeatJoulesPerKilogramKelvin: pr*1000))
    }
    @Test func similarityReferenceTable() throws {
        // Independent numerical similarity solutions: NPTEL lecture 8, table 2.2.
        for (pr, coefficient) in [(0.6,0.276),(0.8,0.307),(1,0.332),(7,0.645),(15,0.835)] {
            let result = try LaminarFlatPlateHeating().evaluate(query(pr: pr))
            #expect(abs(result.localNusseltNumber/sqrt(result.localReynoldsNumber)/coefficient-1) < 0.02)
            // Reynolds number uses film density, not the different edge density.
            #expect(abs(result.localReynoldsNumber-100000) < 1e-8)
            #expect(result.surfacePositionMetres == SIMD3(0.1,2,3))
        }
    }
    @Test func localDistributionIntegratesToMeanPlateCorrelation() throws {
        let model = LaminarFlatPlateHeating(), length = 0.4
        let end = try model.evaluate(query(x: length))
        // Integral from the leading edge: mean Nu_L = 0.664 sqrt(Re_L) at Pr=1.
        // Use x=L*t² to integrate the integrable leading-edge singularity.
        var integral = 0.0
        let count = 16
        for i in 0..<count {
            let t = (Double(i)+0.5)/Double(count)
            // Integrate t=0.1...1 to keep every station inside Reynolds coverage.
            let boundedT = 0.1+0.9*t
            let result = try model.evaluate(query(x: length*boundedT*boundedT))
            integral += result.wattsPerSquareMetre*2*length*boundedT*0.9/Double(count)
        }
        let expected = 2*end.wattsPerSquareMetre*length*0.9
        #expect(abs(integral/expected-1) < 1e-12)
        let upstream = try model.evaluate(query(x: 0.1))
        #expect(abs(upstream.wattsPerSquareMetre/end.wattsPerSquareMetre-2) < 1e-12)
    }
    @Test func coolingAndEquilibriumPreserveSign() throws {
        let model = LaminarFlatPlateHeating()
        let heating = try model.evaluate(query())
        let cooling = try model.evaluate(query(wall: 330,edge: 300))
        #expect(heating.wattsPerSquareMetre > 0)
        #expect(abs(heating.wattsPerSquareMetre+cooling.wattsPerSquareMetre) < 1e-12)
        #expect(try model.evaluate(query(wall: 315,edge: 315)).wattsPerSquareMetre == 0)
    }
    @Test func coverageAndRequiredProperties() throws {
        let model = LaminarFlatPlateHeating()
        for q in [try query(speed: 0),try query(speed: 110),try query(x: 1),try query(x: 0.0001),
                  try query(pr: 0.5),try query(pr: 16),try query(path: 0.002),try query(wall: 100)] {
            #expect(throws: AerodynamicError.outsideCoverage) { try model.evaluate(q) }
        }
        #expect(throws: AerodynamicError.missingFlowProperty) { try model.evaluate(query(sound: nil)) }
        #expect(throws: AerodynamicError.missingFlowProperty) { try model.evaluate(query(path: nil)) }
        #expect(throws: AerodynamicError.unsupportedQuery) { try model.evaluate(query(film: 300)) }
        #expect(throws: AerodynamicError.self) { try query(x: 0) }
        #expect(throws: AerodynamicError.self) { try query(pr: .nan) }
    }
}
