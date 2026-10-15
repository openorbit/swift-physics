// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Aerodynamics

@Suite struct LocalCFDTests {
    func configuration(count: Int = 8, viscosity: Double = 0, boundary: SIMD3<Double> = .zero) throws -> LocalCFDConfiguration {
        try LocalCFDConfiguration(dimensions: SIMD3(repeating: count),minimumCornerMetres: .zero,
            cellSizeMetres: 1,densityKilogramsPerCubicMetre: 2,speedOfSoundMetresPerSecond: 340,meanFreePathMetres: 1e-7,
            kinematicViscositySquareMetresPerSecond: viscosity,boundaryVelocityMetresPerSecond: boundary)
    }
    @Test func uniformFlowRemainsUniformOnRectangularGrid() throws {
        let velocity = SIMD3(2.0,1,-0.5)
        let config = try LocalCFDConfiguration(dimensions: SIMD3(4,6,8),minimumCornerMetres: SIMD3(-2,-3,-4),
            cellSizeMetres: 1,densityKilogramsPerCubicMetre: 1,speedOfSoundMetresPerSecond: 340,meanFreePathMetres: 1e-7,
            kinematicViscositySquareMetresPerSecond: 0.02,boundaryVelocityMetresPerSecond: velocity)
        var grid = try LocalCFDGrid(configuration: config)
        let report = try grid.advance(timeStepSeconds: 0.1)
        #expect(report.substeps == 1)
        #expect(report.maximumDivergenceAfterPerSecond < 1e-12)
        let sampled = try grid.velocityMetresPerSecond(at: .zero)
        #expect(abs(sampled.x-2)+abs(sampled.y-1)+abs(sampled.z+0.5) < 1e-12)
        #expect(try grid.pressureLoadsOnSolids().forceNewtons == .zero)
    }
    @Test func manufacturedGradientIsRemovedAndPressureMatchesPoissonSolution() throws {
        let wave = Double.pi/8, amplitude = 2*sin(Double.pi/16)
        var grid = try LocalCFDGrid(configuration: configuration()) { _ in false } initialVelocity: { p in
            SIMD3(-amplitude*sin(wave*p.x)*cos(wave*p.y)*cos(wave*p.z),
                  -amplitude*cos(wave*p.x)*sin(wave*p.y)*cos(wave*p.z),
                  -amplitude*cos(wave*p.x)*cos(wave*p.y)*sin(wave*p.z))
        }
        let report = try grid.project(timeStepSeconds: 0.1)
        #expect(report.maximumDivergenceBeforePerSecond > 0.1)
        #expect(report.maximumDivergenceAfterPerSecond < 1e-8)
        #expect(report.relativePressureResidual < 1e-8)
        #expect(grid.maximumSpeedBoundMetresPerSecond < 1e-8)
        for k in 0..<8 { for j in 0..<8 { for i in 0..<8 {
            let expected = 20*cos(wave*(Double(i)+0.5))*cos(wave*(Double(j)+0.5))*cos(wave*(Double(k)+0.5))
            #expect(abs(grid.pressurePascals[i+8*(j+8*k)]-expected) < 1e-7)
        } } }
    }
    @Test func solenoidalFieldSurvivesProjectionAndViscosityDissipates() throws {
        let wave = Double.pi/8
        var grid = try LocalCFDGrid(configuration: configuration(viscosity: 0.1),initialVelocity: { p in
            SIMD3(sin(wave*p.x)*cos(wave*p.y),-cos(wave*p.x)*sin(wave*p.y),0)
        })
        let before = grid.maximumSpeedBoundMetresPerSecond
        #expect(grid.divergencePerSecond().map(abs).max()! < 1e-14)
        _ = try grid.project(timeStepSeconds: 0.1)
        #expect(abs(grid.maximumSpeedBoundMetresPerSecond-before) < 1e-12)
        let point = SIMD3(3.5,2.5,4.5), original = try grid.velocityMetresPerSecond(at: point)
        var inviscid = try LocalCFDGrid(configuration: configuration(),initialVelocity: { p in
            SIMD3(sin(wave*p.x)*cos(wave*p.y),-cos(wave*p.x)*sin(wave*p.y),0)
        })
        _ = try inviscid.advance(timeStepSeconds: 0.1)
        let report = try grid.advance(timeStepSeconds: 0.1)
        #expect(report.maximumDivergenceAfterPerSecond < 1e-7)
        #expect(grid.maximumSpeedBoundMetresPerSecond < before)
        let change = try grid.velocityMetresPerSecond(at: point)-inviscid.velocityMetresPerSecond(at: point)
        let expected = -original*(0.1*0.1*8*pow(sin(wave/2),2))
        #expect(abs(change.x-expected.x)+abs(change.y-expected.y)+abs(change.z-expected.z) < 1e-8)
    }
    @Test func obstaclesProjectionPressureLoadsAndReferencePoint() throws {
        var grid = try LocalCFDGrid(configuration: configuration(boundary: SIMD3(2,0,0)),solidAtCellCentre: { p in
            p.x > 3 && p.x < 5 && p.y > 3 && p.y < 5 && p.z > 3 && p.z < 5
        })
        let report = try grid.project(timeStepSeconds: 0.1)
        #expect(report.maximumDivergenceBeforePerSecond > 1)
        #expect(report.maximumDivergenceAfterPerSecond < 1e-7)
        #expect(throws: AerodynamicError.unsupportedQuery) { try grid.velocityMetresPerSecond(at: SIMD3(4,4,4)) }
        let loads = try grid.pressureLoadsOnSolids(momentReferencePointMetres: SIMD3(4,4,4))
        #expect(abs(loads.forceNewtons.y)+abs(loads.forceNewtons.z) < 1e-8)
        #expect(loads.forceNewtons.x > 0)
        let shifted = try grid.pressureLoadsOnSolids(momentReferencePointMetres: .zero)
        let expected = try loads.referred(to: .zero)
        #expect(abs(shifted.momentNewtonMetres.y-expected.momentNewtonMetres.y) < 1e-8)
        _ = try grid.advance(timeStepSeconds: 0.02)
    }
    @Test func hydrostaticPressureGivesDisplacedFluidBuoyancy() throws {
        var grid = try LocalCFDGrid(configuration: configuration(),solidAtCellCentre: { p in
            p.x > 3 && p.x < 5 && p.y > 3 && p.y < 5 && p.z > 3 && p.z < 5
        })
        _ = try grid.advance(timeStepSeconds: 0.1,accelerationMetresPerSecondSquared: { _ in SIMD3(0,0,3) })
        #expect(grid.maximumSpeedBoundMetresPerSecond < 1e-8)
        let loads = try grid.pressureLoadsOnSolids(momentReferencePointMetres: SIMD3(4,4,4))
        // Eight cubic metres displaced, rho=2 kg/m³ and gravity=3 m/s².
        #expect(abs(loads.forceNewtons.z+48) < 1e-7)
        #expect(abs(loads.forceNewtons.x)+abs(loads.forceNewtons.y) < 1e-7)
    }
    @Test func rotorForcingProducesDownstreamFlowAndSwirl() throws {
        var grid = try LocalCFDGrid(configuration: configuration(count: 12))
        let disk = try ActuatorDiskPropwash(centreMetres: SIMD3(6,6,6),downstreamAxis: SIMD3(1,0,0),radiusMetres: 2,
            thrustNewtons: 50,airTorqueNewtonMetres: 5,densityKilogramsPerCubicMetre: 2,upstreamAxialSpeedMetresPerSecond: 0)
        let source = try ActuatorDiskAccelerationField(disk: disk,thicknessMetres: 2)
        for _ in 0..<3 { _ = try grid.advance(timeStepSeconds: 0.05,accelerationMetresPerSecondSquared: source.accelerationMetresPerSecondSquared) }
        let downstream = try grid.velocityMetresPerSecond(at: SIMD3(7.5,6.5,6))
        #expect(downstream.x > 0)
        #expect(downstream.z > 0)
        #expect(grid.divergencePerSecond().map(abs).max()! < 1e-7)
    }
    @Test func limitsAndFailuresDoNotMutateGrid() throws {
        var grid = try LocalCFDGrid(configuration: configuration(),initialVelocity: { p in SIMD3(p.x/8,p.y/8,p.z/8) })
        let before = try grid.velocityMetresPerSecond(at: SIMD3(4,4,4))
        #expect(throws: AerodynamicError.nonConvergence) { try grid.project(timeStepSeconds: 0.1,maximumIterations: 1) }
        #expect(try grid.velocityMetresPerSecond(at: SIMD3(4,4,4)) == before)
        #expect(throws: AerodynamicError.outsideCoverage) { try grid.advance(timeStepSeconds: 10,maximumSubsteps: 1) }
        #expect(try grid.velocityMetresPerSecond(at: SIMD3(4,4,4)) == before)
        #expect(throws: AerodynamicError.invalidFlow) {
            try grid.advance(timeStepSeconds: 0.1,accelerationMetresPerSecondSquared: { _ in SIMD3(.nan,0,0) })
        }
        #expect(throws: AerodynamicError.nonFiniteResult) { try grid.project(timeStepSeconds: .leastNonzeroMagnitude) }
        #expect(try grid.velocityMetresPerSecond(at: SIMD3(4,4,4)) == before)
        #expect(throws: AerodynamicError.outsideCoverage) { try grid.velocityMetresPerSecond(at: SIMD3(-1,4,4)) }
        #expect(throws: AerodynamicError.outsideCoverage) { try configuration(boundary: SIMD3(200,0,0)) }
        #expect(throws: AerodynamicError.invalidGeometry) { try LocalCFDGrid(configuration: configuration(),solidAtCellCentre: { _ in true }) }
        #expect(throws: AerodynamicError.invalidParameter) { try configuration(count: 3) }
        #expect(throws: AerodynamicError.outsideCoverage) {
            try LocalCFDConfiguration(minimumCornerMetres: .zero,cellSizeMetres: 1,
                densityKilogramsPerCubicMetre: 1,speedOfSoundMetresPerSecond: 340,meanFreePathMetres: 1)
        }
    }
}
