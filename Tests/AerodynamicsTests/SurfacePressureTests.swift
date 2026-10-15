// SPDX-License-Identifier: Apache-2.0
import FluidDynamics
import Foundation
import Testing
import Aerodynamics

@Suite struct SurfacePressureTests {
    func flow(_ velocity: SIMD3<Double> = SIMD3(1000,0,0), path: Double? = 0.001) throws -> FluidFlow {
        try FluidFlow(relativeVelocityMetresPerSecond: velocity, densityKilogramsPerCubicMetre: 0.02,
                            speedOfSoundMetresPerSecond: 100, meanFreePathMetres: path)
    }
    @Test func normalShockReferenceAndInclinedPanel() throws {
        let model = try ModifiedNewtonianPressure(specificHeatRatio: 1.4)
        // Independently evaluated normal-shock + downstream stagnation reference, M=10.
        let maximum = try model.stagnationPressureCoefficient(machNumber: 10)
        #expect(abs(maximum-1.831670977387537) < 1e-12)
        // γ=2, M=5 has an exact rational result from the normal-shock relations.
        let otherGas = try ModifiedNewtonianPressure(specificHeatRatio: 2)
        #expect(abs((try otherGas.stagnationPressureCoefficient(machNumber: 5))-20141.0/12100) < 1e-12)
        let panel = try FluidSurfacePanel(areaSquareMetres: 2, centreMetres: SIMD3(0,1,0),
                                               outwardNormal: SIMD3(1,0,1))
        let surface = try SurfacePressureIntegration(panels: [panel], referenceLengthMetres: 1, model: model)
        let result = try surface.evaluate(flow: flow())
        #expect(abs(result.pressureCoefficients[0]-maximum/2) < 1e-12)
        let component = -10000*maximum/sqrt(2)
        #expect(abs(result.loads.forceNewtons.x-component) < 1e-8)
        #expect(abs(result.loads.forceNewtons.z-component) < 1e-8)
        #expect(abs(result.loads.momentNewtonMetres.x-component) < 1e-8)
        #expect(abs(result.loads.momentNewtonMetres.z+component) < 1e-8)
    }
    @Test func sphereMatchesAnalyticDragWithConvergence() throws {
        let model = try ModifiedNewtonianPressure(specificHeatRatio: 1.4)
        let maximum = try model.stagnationPressureCoefficient(machNumber: 10)
        var previous = Double.infinity
        for count in [8,16,32] {
            var panels = [FluidSurfacePanel]()
            for i in 0..<count {
                let mu = (Double(i)+0.5)/Double(count)
                for j in 0..<16 {
                    let phi = 2*Double.pi*(Double(j)+0.5)/16
                    let normal = SIMD3(mu,sqrt(1-mu*mu)*cos(phi),sqrt(1-mu*mu)*sin(phi))
                    panels.append(try FluidSurfacePanel(areaSquareMetres: 2*Double.pi/Double(count*16),
                        centreMetres: normal, outwardNormal: normal))
                }
            }
            let surface = try SurfacePressureIntegration(panels: panels, referenceLengthMetres: 1, model: model)
            let result = try surface.evaluate(flow: flow())
            // Unit sphere: CD=2 CpMax ∫₀¹ μ³ dμ = CpMax/2, based on πR².
            let expected = -10000*Double.pi*maximum/2
            let error = abs(result.loads.forceNewtons.x-expected)
            #expect(error < previous/3); previous = error
            #expect(abs(result.loads.forceNewtons.y)+abs(result.loads.forceNewtons.z) < 1e-8)
            #expect(abs(result.loads.momentNewtonMetres.x)+abs(result.loads.momentNewtonMetres.y)
                + abs(result.loads.momentNewtonMetres.z) < 1e-8)
        }
    }
    @Test func exposureRestAndCoverage() throws {
        let model = try ModifiedNewtonianPressure(specificHeatRatio: 1.4)
        for normal in [SIMD3<Double>(-1,0,0), SIMD3(0,1,0)] {
            let panel = try FluidSurfacePanel(areaSquareMetres: 1, centreMetres: .zero, outwardNormal: normal)
            #expect(try model.pressureCoefficient(panel: panel, flow: flow(), referenceLengthMetres: 1) == 0)
        }
        let panel = try FluidSurfacePanel(areaSquareMetres: 1, centreMetres: .zero,
            outwardNormal: SIMD3(1,0,0), exposedToFlow: false)
        #expect(try model.pressureCoefficient(panel: panel, flow: flow(), referenceLengthMetres: 1) == 0)
        #expect(throws: FluidError.outsideCoverage) {
            try model.pressureCoefficient(panel: panel, flow: flow(path: 0.02), referenceLengthMetres: 1)
        }
        #expect(throws: FluidError.missingFlowProperty) {
            try model.pressureCoefficient(panel: panel, flow: flow(path: nil), referenceLengthMetres: 1)
        }
        #expect(throws: FluidError.outsideCoverage) { try model.stagnationPressureCoefficient(machNumber: 4.9) }
        #expect(throws: FluidError.invalidGeometry) {
            try FluidSurfacePanel(areaSquareMetres: 1, centreMetres: .zero, outwardNormal: .zero)
        }
        let surface = try SurfacePressureIntegration(panels: [panel], referenceLengthMetres: 1, model: model)
        #expect(try surface.evaluate(flow: flow(.zero,path: nil)).loads.forceNewtons == .zero)
    }
}
