// SPDX-License-Identifier: Apache-2.0
import FluidDynamics
import Foundation
import Testing
import Aerodynamics

@Suite struct FiniteWingTests {
    private struct LinearPolar: AirfoilCoefficientModel {
        var offset = 0.0
        func coefficients(angleOfAttackRadians: Double,reynoldsNumber: Double,machNumber: Double) throws -> AirfoilCoefficients {
            AirfoilCoefficients(lift: offset+2 * .pi*angleOfAttackRadians,drag:0.02,moment:-0.03)
        }
    }
    private func model(aspectRatio:Double = 8, offset:Double = 0, iterations:Int = 64) throws -> FiniteWingAirfoil {
        try FiniteWingAirfoil(airfoil:LinearPolar(offset:offset),aspectRatio:aspectRatio,spanEfficiency:0.8,
            effectiveAngleRangeRadians:-0.3...0.3,operatingAngleRangeRadians:-0.2...0.2,maximumIterations:iterations)
    }
    @Test(arguments:[-0.15,0.0,0.15]) func matchesClassicalFiniteWingSlope(alpha:Double) throws {
        let finite = try model()
        let result = try finite.solve(angleOfAttackRadians:alpha,reynoldsNumber:1e6,machNumber:0.2)
        let a0 = 2*Double.pi
        let expectedLift = a0*alpha/(1+a0/(Double.pi*0.8*8))
        #expect(abs(result.coefficients.lift-expectedLift) < 1e-9)
        #expect(abs(result.inducedDragCoefficient-expectedLift*expectedLift/(Double.pi*0.8*8)) < 1e-9)
        #expect(abs(result.coefficients.drag-0.02-result.inducedDragCoefficient) < 1e-14)
        #expect(result.coefficients.moment == -0.03)
        #expect(abs(result.angleResidualRadians) <= finite.angleToleranceRadians)
        #expect(result.iterations <= 64)
    }
    @Test func camberAndAspectRatioLimits() throws {
        let cambered = try model(offset:0.4)
        let result = try cambered.solve(angleOfAttackRadians:0,reynoldsNumber:1e6,machNumber:0.1)
        #expect(abs(result.coefficients.lift-0.4/(1+2/(0.8*8))) < 1e-9)
        #expect(result.effectiveAngleOfAttackRadians < 0)
        let low = try model(aspectRatio:4).solve(angleOfAttackRadians:0.1,reynoldsNumber:1e6,machNumber:0.1)
        let high = try model(aspectRatio:10000).solve(angleOfAttackRadians:0.1,reynoldsNumber:1e6,machNumber:0.1)
        #expect(low.coefficients.lift < high.coefficients.lift)
        #expect(low.inducedDragCoefficient > high.inducedDragCoefficient)
        #expect(abs(high.coefficients.lift-2 * .pi*0.1) < 0.001)
    }
    @Test func boundsAndWorkLimits() throws {
        let finite = try model()
        #expect(throws:FluidError.outsideCoverage) { try finite.solve(angleOfAttackRadians:0.3,reynoldsNumber:1e6,machNumber:0.1) }
        #expect(throws:FluidError.outsideCoverage) { try finite.solve(angleOfAttackRadians:0.1,reynoldsNumber:1e6,machNumber:0.4) }
        #expect(throws:FiniteWingError.iterationLimit) { try model(iterations:1).solve(angleOfAttackRadians:0.1,reynoldsNumber:1e6,machNumber:0.1) }
        let narrow = try FiniteWingAirfoil(airfoil:LinearPolar(),aspectRatio:8,spanEfficiency:1,
            effectiveAngleRangeRadians:-0.01...0.01,operatingAngleRangeRadians:-0.2...0.2)
        #expect(throws:FiniteWingError.unbracketedSolution) { try narrow.solve(angleOfAttackRadians:0.1,reynoldsNumber:1e6,machNumber:0.1) }
        #expect(throws:FluidError.invalidParameter) {
            try FiniteWingAirfoil(airfoil:LinearPolar(),aspectRatio:0,spanEfficiency:1,
                effectiveAngleRangeRadians:-0.2...0.2,operatingAngleRangeRadians:-0.2...0.2)
        }
    }
    @Test func uniformStripIntegrationAvoidsResolutionDependentCorrection() throws {
        let finite = try model(), section = AirfoilWingSectionModel(airfoil:finite)
        let angle = 0.1
        let flow = try FluidFlow(relativeVelocityMetresPerSecond:SIMD3(10*cos(angle),0,10*sin(angle)),
            densityKilogramsPerCubicMetre:2,speedOfSoundMetresPerSecond:100,dynamicViscosityPascalSeconds:0.001)
        let expected = try finite.coefficients(angleOfAttackRadians:angle,reynoldsNumber:20000,machNumber:0.1)
        for count in [1,8,32] {
            let wing = try SegmentedWing.rectangular(spanMetres:16,chordMetres:2,segmentCount:count,model:section)
            let loads = try wing.evaluate(flow:flow).loads
            let direction = flow.relativeVelocityMetresPerSecond/flow.speedMetresPerSecond
            let lift = SIMD3(sin(angle),0,-cos(angle))
            #expect(abs((loads.forceNewtons*lift).sum()-3200*expected.lift) < 1e-7)
            #expect(abs(-(loads.forceNewtons*direction).sum()-3200*expected.drag) < 1e-7)
        }
    }
}
