// SPDX-License-Identifier: Apache-2.0
import FluidDynamics
import Foundation
import Testing
import Aerodynamics

@Suite struct FreeMolecularTests {
    func model(specular: Double = 0, wall: Double = 2) throws -> MaxwellianFreeMolecularTraction {
        try MaxwellianFreeMolecularTraction(gas: FreeMolecularGas(identifier: "Synthetic single species",
            specificGasConstantJoulesPerKilogramKelvin: 2,source: "Analytic kinetic test"),
            surface: GasSurfaceReflection(specularFraction: specular,wallTemperatureKelvin: wall))
    }
    func flow(_ v: SIMD3<Double>, path: Double? = 20, temperature: Double? = 2) throws -> FluidFlow {
        try FluidFlow(relativeVelocityMetresPerSecond: v,densityKilogramsPerCubicMetre: 1,
                            temperatureKelvin: temperature,meanFreePathMetres: path)
    }
    func panel() throws -> FluidSurfacePanel {
        try FluidSurfacePanel(areaSquareMetres: 2,centreMetres: SIMD3(0,1,0),outwardNormal: SIMD3(1,0,0))
    }
    @Test func gaussianMomentumMatchesIndependentQuadrature() throws {
        for normalSpeed in [-25.0,-3.0,0,3,25] {
            // Independent Simpson integration over incoming normal molecular velocities.
            let sigma = 2.0, lower = -max(20*sigma,normalSpeed+12*sigma), count = 4096
            let step = -lower/Double(count)
            var mass = 0.0, pressure = 0.0
            for i in 0...count {
                let c = lower+Double(i)*step
                let density = exp(-pow((c+normalSpeed)/sigma,2)/2)/(sigma*sqrt(2*Double.pi))
                let weight = i == 0 || i == count ? 1.0 : (i%2 == 0 ? 2.0 : 4.0)
                mass += weight*(-c)*density; pressure += weight*c*c*density
            }
            mass *= step/3; pressure *= step/3
            for specular in [0.0,0.35,1] {
                let value = try model(specular: specular,wall: 8).traction(panel: panel(),
                    flow: flow(SIMD3(normalSpeed,4,0)),referenceLengthMetres: 1)
                let expectedPressure = (1+specular)*pressure+(1-specular)*mass*sqrt(2*8*Double.pi/2)-4
                #expect(abs(value.gaugePressurePascals-expectedPressure) < 1e-7)
                #expect(abs(value.shearPascals.y+(1-specular)*mass*4) < 1e-7)
                #expect(value.shearPascals.x == 0)
            }
        }
    }
    @Test func thermalRestAndLeewardImpacts() throws {
        let rest = try model().traction(panel: panel(),flow: flow(.zero),referenceLengthMetres: 1)
        #expect(abs(rest.gaugePressurePascals) < 1e-12)
        #expect(rest.shearPascals == .zero)
        let hotWall = try model(wall: 8).traction(panel: panel(),flow: flow(.zero),referenceLengthMetres: 1)
        #expect(abs(hotWall.gaugePressurePascals-2) < 1e-12)
        let specular = try model(specular: 1,wall: 8).traction(panel: panel(),flow: flow(.zero),referenceLengthMetres: 1)
        #expect(abs(specular.gaugePressurePascals) < 1e-12)
        let leeward = try model().traction(panel: panel(),flow: flow(SIMD3(-3,4,0)),referenceLengthMetres: 1)
        #expect(leeward.gaugePressurePascals > -4)
        #expect(leeward.shearPascals.y < 0)
    }
    @Test func specularBeamLimitAndMomentTransfer() throws {
        let normal = SIMD3<Double>(1/sqrt(2),1/sqrt(2),0)
        let surfacePanel = try FluidSurfacePanel(areaSquareMetres: 2,centreMetres: SIMD3(0,0,1),outwardNormal: normal)
        let evaluator = try SurfaceTractionIntegration(panels: [surfacePanel],referenceLengthMetres: 1,model: model(specular: 1))
        let result = try evaluator.evaluate(flow: flow(SIMD3(1000,0,0)))
        // 2 rho (V·n)^2 + rho RT, after subtracting ambient rho RT.
        #expect(abs(result.tractions[0].gaugePressurePascals-(1e6+4)) < 1e-8)
        #expect(result.tractions[0].shearPascals == .zero)
        #expect(abs(result.loads.forceNewtons.x+sqrt(2)*(1e6+4)) < 1e-8)
        #expect(abs(result.loads.momentNewtonMetres.y+sqrt(2)*(1e6+4)) < 1e-8)
    }
    @Test func continuumAdapterAndCoverage() throws {
        let continuum = try ModifiedNewtonianPressure(specificHeatRatio: 1.4)
        let local = try FluidFlow(relativeVelocityMetresPerSecond: SIMD3(1000,0,0),
            densityKilogramsPerCubicMetre: 0.02,speedOfSoundMetresPerSecond: 100,meanFreePathMetres: 0.001)
        let original = try SurfacePressureIntegration(panels: [panel()],referenceLengthMetres: 1,model: continuum).evaluate(flow: local)
        let adapted = try SurfaceTractionIntegration(panels: [panel()],referenceLengthMetres: 1,
            model: PressureOnlyTraction(continuum)).evaluate(flow: local)
        #expect(original.loads.forceNewtons == adapted.loads.forceNewtons)
        #expect(original.loads.momentNewtonMetres == adapted.loads.momentNewtonMetres)
        #expect(throws: FluidError.outsideCoverage) { try model().traction(panel: panel(),flow: flow(.zero,path: 1),referenceLengthMetres: 1) }
        #expect(throws: FluidError.missingFlowProperty) { try model().traction(panel: panel(),flow: flow(.zero,path: nil),referenceLengthMetres: 1) }
        #expect(throws: FluidError.missingFlowProperty) { try model().traction(panel: panel(),flow: flow(.zero,temperature: nil),referenceLengthMetres: 1) }
        let oblique = try FluidSurfacePanel(areaSquareMetres: 1,centreMetres: .zero,
                                                  outwardNormal: SIMD3(1,2,3))
        let evaluator = try SurfaceTractionIntegration(panels: [oblique],referenceLengthMetres: 1,model: model())
        let parallel = try evaluator.evaluate(flow: flow(oblique.outwardNormal*1000))
        #expect(abs(parallel.tractions[0].shearPascals.x)+abs(parallel.tractions[0].shearPascals.y)
            + abs(parallel.tractions[0].shearPascals.z) < 1e-8)
    }
}
