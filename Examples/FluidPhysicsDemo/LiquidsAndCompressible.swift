// SPDX-License-Identifier: Apache-2.0
import Foundation
import FluidDynamics
import Hydrodynamics

// Illustrative constant properties, not a temperature/composition-resolved property database.
func runLiquidExamples() throws {
    let hull = try ConvexDisplacementHull.box(halfExtentsMetres: SIMD3(2,1,0.5))
    for (name,density,gravity,viscosity,temperature,pressure) in [
        ("Earth water illustration",1000.0,9.80665,0.001,293.15,101325.0),
        ("Titan hydrocarbon illustration",660.0,1.352,0.0003,94.0,150000.0)
    ] {
        let liquid = try LiquidProperties(identifier: name,source: "Illustrative constant-property example; not a lake composition/EOS preset",
            densityKilogramsPerCubicMetre: density,dynamicViscosityPascalSeconds: viscosity,speedOfSoundMetresPerSecond: 1500,vaporPressurePascals: 2000)
        let environment = try HydrostaticLiquid(properties: liquid,up: SIMD3(0,0,-1),gravityMetresPerSecondSquared: SIMD3(0,0,gravity),
            surfaceAbsolutePressurePascals: pressure,maximumDepthMetres: 100)
        let immersed = try hull.displacement(in: environment), loads = try immersed.buoyancy(in: environment)
        let heel = 0.1, up = SIMD3(0,sin(heel),-cos(heel))
        let heeled = try HydrostaticLiquid(properties: liquid,up: up,gravityMetresPerSecondSquared: -up*gravity,
            surfaceAbsolutePressurePascals: pressure,maximumDepthMetres: 100)
        let restoring = try hull.displacement(in: heeled).buoyancy(in: heeled,momentReferencePointMetres: SIMD3(0,0,0.25))
        print("\(name): displaced \(immersed.volumeCubicMetres) m³, buoyancy \(loads.forceNewtons) N; heeled buoyancy moment \(restoring.momentNewtonMetres) N m")
        let geometry = try CartesianFluidGrid(dimensions: SIMD3(repeating: 4),minimumCornerMetres: .zero,cellSizeMetres: 1)
        var submergedGrid = try IncompressibleFlowGrid(configuration: IncompressibleFlowConfiguration(grid: geometry,liquid: liquid))
        let report = try submergedGrid.advance(timeStepSeconds: 0.01,accelerationMetresPerSecondSquared: { _ in SIMD3(0,0,gravity) })
        let reference = try LiquidCFDPressureReference(positionMetres: SIMD3(0.5,0.5,0.5),absolutePressurePascals: pressure)
        print("  Fully filled liquid CFD: gauge-calibrated pressure at +3 m depth \(try reference.absolutePressurePascals(in: submergedGrid,at: SIMD3(0.5,0.5,3.5))) Pa; divergence \(report.maximumDivergenceAfterPerSecond) 1/s. T=\(temperature) K is an illustrative material condition.")
    }
}

func runCompressibleExample() throws {
    let eos = try CaloricallyPerfectGasEOS(identifier: "Synthetic shock-tube gas",specificGasConstantJoulesPerKilogramKelvin: 1,
        specificHeatRatio: 1.4,temperatureRangeKelvin: 1e-4...1e4,source: "Standard ideal-gas Sod numerical benchmark")
    let geometry = try CartesianFluidGrid(dimensions: SIMD3(160,1,1),minimumCornerMetres: .zero,cellSizeMetres: 1.0/160)
    var grid = try CompressibleFlowGrid(configuration: CompressibleFlowConfiguration(grid: geometry,equationOfState: eos),initialState: { p in
        try ConservedFluidState(densityKilogramsPerCubicMetre: p.x < 0.5 ? 1 : 0.125,velocityMetresPerSecond: .zero,
            pressurePascals: p.x < 0.5 ? 1 : 0.1,equationOfState: eos)
    })
    let report = try grid.advance(timeStepSeconds: 0.2,maximumSubsteps: 10000)
    let star = try grid.state(at: SIMD3(0.6,0.002,0.002))
    print("Compressible Sod tube: \(report.substeps) acoustic-CFL substeps, pressure near left star region \(star.thermodynamics.pressurePascals) Pa, velocity \(star.velocityMetresPerSecond.x) m/s. Inviscid shock benchmark; wall heat flux unavailable.")
}
