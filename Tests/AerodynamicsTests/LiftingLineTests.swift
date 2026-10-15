// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Aerodynamics

@Suite struct LiftingLineTests {
    func flow(alpha: Double = 0.1, speed: Double = 50, sideways: Double = 0) throws -> AerodynamicFlow {
        try AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(speed*cos(alpha),sideways,speed*sin(alpha)),
            densityKilogramsPerCubicMetre: 1,speedOfSoundMetresPerSecond: 300,dynamicViscosityPascalSeconds: 0.001)
    }
    func wing(count: Int = 8, ellipse: Bool = true, asymmetricTwist: Double = 0,
              centre: SIMD3<Double> = .zero, profileDrag: Double = 0) throws -> LiftingLineWing {
        try LiftingLineWing(spanMetres: 8,referenceAreaSquareMetres: ellipse ? 4*Double.pi : 16,
            stationCount: count,momentReferenceLineCentreMetres: centre,source: "Analytic reference geometry",version: "1") { fraction in
                try LiftingLineSection(chordMetres: ellipse ? 2*sqrt(1-fraction*fraction) : 2,
                    liftSlopePerRadian: 2*Double.pi,geometricIncidenceRadians: -asymmetricTwist*fraction,
                    profileDragCoefficient: profileDrag,effectiveAngleRangeRadians: -0.2...0.2,reynoldsRange: 100...1e8)
            }
    }
    @Test(arguments: [4,8,16,32]) func ellipticalWingMatchesClassicalReference(count: Int) throws {
        let wing = try wing(count: count)
        let result = try wing.evaluate(flow: flow())
        let slope = 2*Double.pi/(1+2/wing.aspectRatio)
        let expectedCL = slope*0.1
        #expect(abs(result.liftCoefficient-expectedCL) < 1e-12)
        #expect(abs(result.inducedDragCoefficient-expectedCL*expectedCL/(Double.pi*wing.aspectRatio)) < 1e-12)
        #expect(abs(try #require(result.spanEfficiency)-1) < 1e-12)
        #expect(result.maximumEquationResidualRadians < 1e-12)
        for station in result.stations {
            #expect(abs(station.sectionLiftCoefficient-2*Double.pi*station.effectiveAngleRadians) < 1e-12)
        }
        #expect(abs(result.loads.momentNewtonMetres.x)+abs(result.loads.momentNewtonMetres.z) < 1e-9)
    }
    @Test func asymmetricTwistProducesAnalyticSecondHarmonicAndRoll() throws {
        let wing = try wing(asymmetricTwist: 0.04)
        let result = try wing.evaluate(flow: flow(alpha: 0))
        let lambda = 4*8/(2*Double.pi*2)
        let expectedA2 = 0.04/(2*(lambda+2))
        #expect(abs(result.fourierCoefficients[1]-expectedA2) < 1e-12)
        #expect(abs(result.liftCoefficient) < 1e-12)
        #expect(result.spanEfficiency == nil)
        #expect(result.inducedDragCoefficient > 0)
        let expectedRoll = Double.pi*1250*pow(8,3)*expectedA2/4
        #expect(abs(result.loads.momentNewtonMetres.x-expectedRoll) < 1e-9)
    }
    @Test func rectangularWingConvergenceAndProfileDrag() throws {
        let reference = try wing(count: 64,ellipse: false,profileDrag: 0.02).evaluate(flow: flow())
        var lastError = Double.infinity
        for count in [8,16,32] {
            let result = try wing(count: count,ellipse: false,profileDrag: 0.02).evaluate(flow: flow())
            let error = abs(result.liftCoefficient-reference.liftCoefficient)
            #expect(error < lastError); lastError = error
            #expect(abs(result.profileDragCoefficient-0.02) < 1e-12)
            #expect((try #require(result.spanEfficiency)) < 1)
            #expect((try #require(result.spanEfficiency)) > 0.8)
        }
        let atLine = try wing(centre: SIMD3(2,0,0)).evaluate(flow: flow(),momentReferencePointMetres: SIMD3(2,0,0))
        let atOrigin = try wing(centre: SIMD3(2,0,0)).evaluate(flow: flow())
        #expect(abs(atOrigin.loads.momentNewtonMetres.y-atLine.loads.momentNewtonMetres.y+2*atOrigin.loads.forceNewtons.z) < 1e-9)
    }
    @Test func coverageAndMissingProperties() throws {
        let wing = try wing()
        #expect(throws: AerodynamicError.outsideCoverage) { try wing.evaluate(flow: flow(speed: 100)) }
        #expect(throws: AerodynamicError.outsideCoverage) { try wing.evaluate(flow: flow(sideways: 1)) }
        #expect(throws: AerodynamicError.outsideCoverage) { try wing.evaluate(flow: flow(alpha: 0.34)) }
        #expect(throws: AerodynamicError.missingFlowProperty) {
            try wing.evaluate(flow: AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(50,0,0),densityKilogramsPerCubicMetre: 1))
        }
    }
}
