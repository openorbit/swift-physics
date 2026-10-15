// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import FluidDynamics
import Hydrodynamics

@Suite struct SubmersibleTests {
    @Test func passiveDampingAddedInertiaAndCombinedLoads() throws {
        let damping = try HydrodynamicDamping(linearForce: SIMD3(1,2,3),quadraticForce: SIMD3(2,3,4),
            linearMoment: SIMD3(4,5,6),quadraticMoment: SIMD3(1,2,3),source: "Synthetic passive coefficients")
        let v = SIMD3(2.0,-3,1), w = SIMD3(-1.0,2,-0.5)
        let loads = try damping.evaluate(relativeVelocityMetresPerSecond: v,relativeAngularVelocityRadiansPerSecond: w)
        #expect((loads.forceNewtons*v).sum()+(loads.momentNewtonMetres*w).sum() < 0)
        let reverse = try damping.evaluate(relativeVelocityMetresPerSecond: -v,relativeAngularVelocityRadiansPerSecond: -w)
        #expect(reverse.forceNewtons == -loads.forceNewtons && reverse.momentNewtonMetres == -loads.momentNewtonMetres)
        let inertia = try HydrodynamicAddedInertia(massKilograms: SIMD3(2,3,4),momentKilogramSquareMetres: SIMD3(5,6,7))
        #expect(try inertia.accelerationReaction(relativeLinearAccelerationMetresPerSecondSquared: v,
            relativeAngularAccelerationRadiansPerSecondSquared: w).forceNewtons == SIMD3(-4,9,-4))
        let p = try HydrostaticsTests().properties()
        let body = SubmergedBodyHydrodynamics(buoyancy: try Buoyancy(displacedVolumeCubicMetres: 2),damping: damping)
        let flow = try FluidFlow(relativeVelocityMetresPerSecond: v,liquid: p,temperatureKelvin: 293,absolutePressurePascals: 100000)
        let combined = try body.evaluate(liquid: p,flow: flow,relativeAngularVelocityRadiansPerSecond: w,gravityMetresPerSecondSquared: SIMD3(0,0,10))
        #expect(combined.forceNewtons == loads.forceNewtons+SIMD3(0,0,-20000))
    }
    @Test func liquidPropertiesAndAbsolutePressureCFDReference() throws {
        let liquid = try HydrostaticsTests().properties()
        let model = try ConstantLiquidModel(value: liquid,temperatureRangeKelvin: 290...300,pressureRangePascals: 5000...1e7)
        #expect(try model.properties(temperatureKelvin: 293,absolutePressurePascals: 100000).densityKilogramsPerCubicMetre == 1000)
        #expect(throws: FluidError.outsideCoverage) { try model.properties(temperatureKelvin: 94,absolutePressurePascals: 100000) }
        #expect(abs(try liquid.cavitationNumber(absolutePressurePascals: 102000,speedMetresPerSecond: 10)-2) < 1e-12)
        let geometry = try CartesianFluidGrid(dimensions: SIMD3(repeating: 4),minimumCornerMetres: .zero,cellSizeMetres: 1)
        var grid = try IncompressibleFlowGrid(configuration: IncompressibleFlowConfiguration(grid: geometry,liquid: liquid))
        _ = try grid.advance(timeStepSeconds: 0.01,accelerationMetresPerSecondSquared: { _ in SIMD3(0,0,10) })
        let reference = try LiquidCFDPressureReference(positionMetres: SIMD3(0.5,0.5,0.5),absolutePressurePascals: 100000)
        #expect(abs(try reference.absolutePressurePascals(in: grid,at: SIMD3(0.5,0.5,3.5))-130000) < 1e-6)
        #expect(throws: FluidError.outsideCoverage) { try FluidFlow(relativeVelocityMetresPerSecond: .zero,liquid: liquid,temperatureKelvin: 293,absolutePressurePascals: 1000) }
    }
}
