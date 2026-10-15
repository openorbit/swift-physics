// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Aerodynamics

@Suite struct NonlinearLiftingLineTests {
    struct Polar: WingSectionCoefficientModel {
        var nonlinear = false
        func coefficients(for query: WingSectionQuery) throws -> AirfoilCoefficients {
            guard let alpha = query.angleOfAttackRadians else { throw AerodynamicError.missingFlowProperty }
            return AirfoilCoefficients(lift: nonlinear ? 3*sin(2*alpha) : 2*Double.pi*alpha,drag: 0.01,moment: -0.02)
        }
    }
    func flow(alpha: Double = 0.08, speed: Double = 30) throws -> AerodynamicFlow {
        try AerodynamicFlow(relativeVelocityMetresPerSecond: SIMD3(speed*cos(alpha),0,speed*sin(alpha)),
            densityKilogramsPerCubicMetre: 1,speedOfSoundMetresPerSecond: 340,dynamicViscosityPascalSeconds: 1.8e-5)
    }
    func segments(count: Int = 8, span: Double = 8, centre: SIMD3<Double> = .zero, nonlinear: Bool = false) throws -> [WingSegment] {
        try SegmentedWing.rectangular(spanMetres: span,chordMetres: 1,segmentCount: count,
            centreMetres: centre,model: Polar(nonlinear: nonlinear)).segments
    }
    func solver(_ segments: [WingSegment], length: Double = 800, core: Double = 0.0001,
                iterations: Int = 1000) throws -> NonlinearLiftingLine {
        try NonlinearLiftingLine(segments: segments,wakeLengthMetres: length,vortexCoreRadiusMetres: core,
            relaxation: 0.05,tolerance: 1e-9,maximumIterations: iterations)
    }
    @Test func singleHorseshoeMatchesAnalyticSemiInfiniteWake() throws {
        let wing = try solver(segments(count: 1),length: 80000,core: 1e-6)
        let result = try wing.evaluate { _ in try flow() }
        let gamma = result.circulationSquareMetresPerSecond[0]
        // Two semi-infinite trailing legs, each distance b/2: total w=Gamma/(pi*b).
        #expect(abs(result.inducedAirVelocityMetresPerSecond[0].z/(gamma/(Double.pi*8))-1) < 1e-8)
        #expect(abs(result.inducedAirVelocityMetresPerSecond[0].x) < 1e-12)
        // Independent scalar bisection of nonlinear section/circulation compatibility.
        var lower = 0.0, upper = 30.0
        for _ in 0..<80 {
            let g = (lower+upper)/2, vx = 30*cos(0.08), vz = 30*sin(0.08)-g/(Double.pi*8)
            let target = 0.5*hypot(vx,vz)*2*Double.pi*atan2(vz,vx)
            if g > target { upper = g } else { lower = g }
        }
        #expect(abs(gamma-(lower+upper)/2) < 1e-7)
        #expect(result.maximumNormalizedCirculationResidual <= 1e-9)
    }
    @Test func rectangularConvergenceAgainstFourierLiftingLine() throws {
        let referenceWing = try LiftingLineWing(spanMetres: 8,referenceAreaSquareMetres: 8,stationCount: 64,
            source: "linear rectangular reference",version: "1") { _ in
                try LiftingLineSection(chordMetres: 1,liftSlopePerRadian: 2*Double.pi,
                    effectiveAngleRangeRadians: -0.2...0.2,reynoldsRange: 1...1e9)
            }
        let alpha = 0.01
        let reference = try referenceWing.evaluate(flow: flow(alpha: alpha))
        var previous = Double.infinity
        for count in [8,16,32,64] {
            let result = try solver(segments(count: count)).evaluate { _ in try flow(alpha: alpha) }
            let cl = -result.loads.forceNewtons.z/(450*8)
            let error = abs(cl-reference.liftCoefficient)
            #expect(error < previous); previous = error
            #expect(abs(result.loads.momentNewtonMetres.x) < 1e-7)
            #expect(abs(result.loads.momentNewtonMetres.z) < 1e-7)
        }
        #expect(previous/reference.liftCoefficient < 0.02)
    }
    @Test func coupledTailNonuniformFlowAndExternalWake() throws {
        let main = try segments(nonlinear: true)
        let tail = try segments(count: 4,span: 3,centre: SIMD3(-4,0,0.5),nonlinear: true)
        let isolated = try solver(tail).evaluate { _ in try flow() }
        let coupled = try solver(main+tail).evaluate { _ in try flow() }
        let tailLift = coupled.segments.suffix(tail.count).reduce(0.0) { $0-$1.loads.forceNewtons.z }
        #expect(tailLift < -isolated.loads.forceNewtons.z)
        let mainSolver = try solver(main)
        let mainResult = try mainSolver.evaluate { _ in try flow() }
        let wake = try mainSolver.inducedAirVelocity(at: SIMD3(-4,0,0.5),circulation: mainResult.circulationSquareMetresPerSecond)
        #expect(wake.z > 0)
        let nonuniform = try mainSolver.evaluate { segment in
            try flow(alpha: 0.08+0.005*segment.positionMetres.y,speed: 30+segment.positionMetres.y)
        }
        #expect(abs(nonuniform.loads.momentNewtonMetres.x) > 1)
        let shifted = try mainSolver.evaluate(momentReferencePointMetres: SIMD3(1,2,3)) { _ in try flow() }
        let expected = try mainResult.loads.referred(to: SIMD3(1,2,3))
        #expect(hypot(hypot((shifted.loads.momentNewtonMetres-expected.momentNewtonMetres).x,
                           (shifted.loads.momentNewtonMetres-expected.momentNewtonMetres).y),
                          (shifted.loads.momentNewtonMetres-expected.momentNewtonMetres).z) < 1e-8)
    }
    @Test func wakeLengthAndCoreSensitivity() throws {
        let s = try segments(count: 16)
        let reference = try solver(s,length: 8000,core: 1e-6).evaluate { _ in try flow() }
        var previous = Double.infinity
        for length in [80.0,800] {
            let result = try solver(s,length: length,core: 1e-6).evaluate { _ in try flow() }
            let error = abs(result.loads.forceNewtons.z-reference.loads.forceNewtons.z)
            #expect(error < previous); previous = error
        }
        let finiteCore = try solver(s,length: 8000,core: 0.001).evaluate { _ in try flow() }
        #expect(abs((finiteCore.loads.forceNewtons.z-reference.loads.forceNewtons.z)/reference.loads.forceNewtons.z) < 1e-4)
    }
    @Test func failureIsExplicitAndNoInvalidFlowFallback() throws {
        let s = try segments()
        #expect(throws: AerodynamicError.nonConvergence) {
            try solver(s,iterations: 1).evaluate { _ in try flow() }
        }
        #expect(throws: AerodynamicError.outsideCoverage) { try solver(s).evaluate { _ in try flow(speed: 200) } }
        #expect(throws: AerodynamicError.invalidGeometry) { try solver([s[0],s[0]]) }
        #expect(throws: AerodynamicError.invalidGeometry) { try solver(s,core: 1) }
        let wing = try solver(s)
        #expect(throws: AerodynamicError.invalidParameter) { try wing.inducedAirVelocity(at: .zero,circulation: []) }
    }
}
