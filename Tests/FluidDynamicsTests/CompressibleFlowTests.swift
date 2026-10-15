// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import FluidDynamics

@Suite struct CompressibleFlowTests {
    func gas() throws -> CaloricallyPerfectGasEOS {
        try CaloricallyPerfectGasEOS(identifier: "Synthetic gamma=1.4 gas",specificGasConstantJoulesPerKilogramKelvin: 1,
            specificHeatRatio: 1.4,temperatureRangeKelvin: 1e-4...1e4,source: "Analytic Euler references")
    }
    func state(rho: Double = 1, velocity: SIMD3<Double> = .zero, pressure: Double = 1) throws -> ConservedFluidState {
        try ConservedFluidState(densityKilogramsPerCubicMetre: rho,velocityMetresPerSecond: velocity,pressurePascals: pressure,equationOfState: gas())
    }
    @Test func uniformSupersonicStateAndEOSCoverage() throws {
        let eos = try gas(), value = try state(velocity: SIMD3(5,1,-2))
        let geometry = try CartesianFluidGrid(dimensions: SIMD3(8,4,4),minimumCornerMetres: .zero,cellSizeMetres: 0.25)
        var grid = try CompressibleFlowGrid(configuration: CompressibleFlowConfiguration(grid: geometry,equationOfState: eos),initialState: { _ in value })
        let before = try grid.conservedTotals(), report = try grid.advance(timeStepSeconds: 0.01)
        #expect(report.maximumCourantNumber <= 0.4+1e-15)
        for cell in grid.cells { #expect(abs(cell.densityKilogramsPerCubicMetre-1) < 1e-14) }
        let after = try grid.conservedTotals()
        #expect(abs(after.energyJoules-before.energyJoules) < 1e-12)
        #expect(throws: FluidError.outsideCoverage) { try eos.thermodynamics(densityKilogramsPerCubicMetre: 1,specificInternalEnergyJoulesPerKilogram: 1e9) }
    }
    @Test func periodicConservationAndSolidActionReaction() throws {
        let eos = try gas(), geometry = try CartesianFluidGrid(dimensions: SIMD3(repeating: 8),minimumCornerMetres: .zero,cellSizeMetres: 1)
        let config = try CompressibleFlowConfiguration(grid: geometry,equationOfState: eos,boundaries: Array(repeating: .periodic,count: 6))
        var grid = try CompressibleFlowGrid(configuration: config,solidAtCellCentre: { p in
            p.x > 3 && p.x < 5 && p.y > 3 && p.y < 5 && p.z > 3 && p.z < 5
        },initialState: { _ in try state(velocity: SIMD3(1,0,0)) })
        let before = try grid.conservedTotals(), report = try grid.advance(timeStepSeconds: 0.05)
        let after = try grid.conservedTotals(), impulse = report.averageLoadsOnSolids.forceNewtons*0.05
        #expect(abs(after.massKilograms-before.massKilograms) < 1e-9)
        #expect(abs(after.energyJoules-before.energyJoules) < 1e-9)
        #expect(abs(after.momentumKilogramMetresPerSecond.x-before.momentumKilogramMetresPerSecond.x+impulse.x
                    + report.boundaryMomentumOutflowKilogramMetresPerSecond.x) < 1e-9)
        #expect(impulse.x > 0)
    }
    func sod(count: Int) throws -> CompressibleFlowGrid {
        let geometry = try CartesianFluidGrid(dimensions: SIMD3(count,1,1),minimumCornerMetres: .zero,cellSizeMetres: 1/Double(count))
        var grid = try CompressibleFlowGrid(configuration: CompressibleFlowConfiguration(grid: geometry,equationOfState: gas()),initialState: { p in
            try state(rho: p.x < 0.5 ? 1 : 0.125,pressure: p.x < 0.5 ? 1 : 0.1)
        })
        _ = try grid.advance(timeStepSeconds: 0.2,maximumSubsteps: 10000)
        return grid
    }
    // Independently specified exact Sod solution: standard star state and self-similar fan.
    func exactSodDensity(x: Double) -> Double {
        let xi = (x-0.5)/0.2, gamma = 1.4, c = sqrt(gamma), u = 0.92745262004895, p = 0.303130178050647
        let starLeft = pow(p,1/gamma), starRight = 0.125*((p/0.1)+(gamma-1)/(gamma+1))/((gamma-1)/(gamma+1)*(p/0.1)+1)
        let tail = u-c*pow(p,(gamma-1)/(2*gamma)), shock = sqrt((gamma+1)/(2*gamma)*(p/0.1)+(gamma-1)/(2*gamma))*sqrt(gamma*0.1/0.125)
        if xi < -c { return 1 }
        if xi < tail { return pow(2/(gamma+1)-(gamma-1)/(gamma+1)*xi/c,2/(gamma-1)) }
        if xi < u { return starLeft }
        if xi < shock { return starRight }
        return 0.125
    }
    @Test func sodShockTubeConvergesToIndependentExactReference() throws {
        var previous = Double.infinity
        for count in [80,160] {
            let result = try sod(count: count)
            let error = result.cells.enumerated().reduce(0.0) { value, item in
                value+abs(item.element.densityKilogramsPerCubicMetre-exactSodDensity(x: (Double(item.offset)+0.5)/Double(count)))/Double(count)
            }
            #expect(error < previous); previous = error
            for cell in result.cells { #expect(try cell.primitive(using: gas()).thermodynamics.pressurePascals > 0) }
        }
        #expect(previous < 0.04)
        let reference = try sod(count: 160).state(at: SIMD3(0.6,0.002,0.002))
        #expect(abs(reference.thermodynamics.pressurePascals-0.30313) < 0.025)
        #expect(abs(reference.velocityMetresPerSecond.x-0.92745) < 0.05)
    }
    @Test func failedBudgetAndInvalidEOSStateAreAtomic() throws {
        let geometry = try CartesianFluidGrid(dimensions: SIMD3(8,1,1),minimumCornerMetres: .zero,cellSizeMetres: 0.1)
        var grid = try CompressibleFlowGrid(configuration: CompressibleFlowConfiguration(grid: geometry,equationOfState: gas()),initialState: { _ in try state() })
        let before = grid.cells
        #expect(throws: FluidError.workBudgetExceeded) { try grid.advance(timeStepSeconds: 1,maximumSubsteps: 1) }
        #expect(grid.cells == before)
        #expect(throws: FluidError.invalidFlow) { try ConservedFluidState(densityKilogramsPerCubicMetre: 0,velocityMetresPerSecond: .zero,pressurePascals: 1,equationOfState: gas()) }
        var boundaries = Array(repeating: CompressibleBoundary.transmissive,count: 6); boundaries[0] = .periodic
        #expect(throws: FluidError.invalidParameter) { try CompressibleFlowConfiguration(grid: geometry,equationOfState: gas(),boundaries: boundaries) }
    }
}
