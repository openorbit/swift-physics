// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Math
@testable import CelMek

@Suite struct LambertSolverTests {
    let frame = Frame(name: "Lambert", pseudoInertial: true)
    func gm(_ value: Double = 1) throws -> GravitationalParameter {
        try GravitationalParameter(cubicMetresPerSecondSquared: value, convention: .newtonian,
            provenance: PhysicalDataProvenance(source: "Analytic test", version: "1"))
    }
    func solve(_ r1: SIMD3<Double>, _ r2: SIMD3<Double>, _ time: Double,
        direction: LambertDirection = .positive, options: LambertOptions = try! LambertOptions(absoluteTimeToleranceSeconds: 1e-13)) throws -> LambertSolution {
        try LambertSolver.solve(departurePosition: r1, arrivalPosition: r2, departureEpoch: .j2000,
            arrivalEpoch: Epoch.j2000.advanced(by: CelMek.Duration(seconds: time)), frame: frame,
            centralBody: .sun, gravitationalParameter: gm(), referencePole: SIMD3(0,0,1), direction: direction, options: options)
    }
    func check(_ solution: LambertSolution, _ v1: SIMD3<Double>, _ v2: SIMD3<Double>, tolerance: Double = 1e-9) throws {
        #expect(routeNorm(solution.departureState.velocity-v1) < tolerance)
        #expect(routeNorm(solution.arrivalState.velocity-v2) < tolerance)
        let endpoint = try solution.trajectory().state(at: solution.arrivalState.epoch)
        #expect(routeNorm(endpoint.position-solution.arrivalState.position) < tolerance * 10)
        #expect(routeNorm(endpoint.velocity-v2) < tolerance * 10)
        #expect(abs(solution.timeResidualSeconds) < 2e-9)
    }
    @Test func circularBothDirectionsAndLongWay() throws {
        let positive = try solve(SIMD3(1,0,0), SIMD3(0,1,0), .pi/2)
        try check(positive, SIMD3(0,1,0), SIMD3(-1,0,0))
        #expect(positive.branch == .shortWay)
        let negative = try solve(SIMD3(1,0,0), SIMD3(0,1,0), 3 * .pi/2, direction: .negative)
        try check(negative, SIMD3(0,-1,0), SIMD3(1,0,0))
        #expect(negative.branch == .longWay)
        #expect(negative.orbitalNormal.z == -1)
    }
    @Test func analyticConics() throws {
        // Independent Kepler eccentric-anomaly equations, a=1, mu=1.
        let e = 0.4, b = sqrt(1-e*e), e1 = 0.2, e2 = 2.1
        func ellipse(_ anomaly: Double) -> (SIMD3<Double>, SIMD3<Double>) {
            (SIMD3(cos(anomaly)-e,b*sin(anomaly),0), SIMD3(-sin(anomaly),b*cos(anomaly),0)/(1-e*cos(anomaly)))
        }
        let first = ellipse(e1), second = ellipse(e2)
        try check(solve(first.0, second.0, e2-e*sin(e2)-e1+e*sin(e1)), first.1, second.1)
        // Hyperbolic anomaly equations, |a|=1, e=1.5, mu=1.
        let eh = 1.5, bh = sqrt(eh*eh-1), h1 = -0.3, h2 = 0.6
        func hyperbola(_ anomaly: Double) -> (SIMD3<Double>, SIMD3<Double>) {
            (SIMD3(eh-cosh(anomaly),bh*sinh(anomaly),0), SIMD3(-sinh(anomaly),bh*cosh(anomaly),0)/(eh*cosh(anomaly)-1))
        }
        let hf = hyperbola(h1), hs = hyperbola(h2)
        try check(solve(hf.0, hs.0, eh*sinh(h2)-h2-eh*sinh(h1)+h1), hf.1, hs.1)
        // Barker's equation for a parabola, p=1 and mu=1.
        let f1 = -0.4, f2 = 0.7, d1 = tan(f1/2), d2 = tan(f2/2)
        func parabola(_ f: Double) -> (SIMD3<Double>, SIMD3<Double>) {
            (SIMD3(cos(f),sin(f),0)/(1+cos(f)), SIMD3(-sin(f),1+cos(f),0))
        }
        let pf = parabola(f1), ps = parabola(f2)
        try check(solve(pf.0, ps.0, (d2+d2*d2*d2/3-d1-d1*d1*d1/3)/2), pf.1, ps.1)
    }
    @Test func inclinedNearCollinearAndTimeReversal() throws {
        let r1 = SIMD3<Double>(1,0,0), r2 = SIMD3<Double>(0,0,1)
        let end = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: .pi/2))
        let inclined = try LambertSolver.solve(departurePosition: r1, arrivalPosition: r2,
            departureEpoch: .j2000, arrivalEpoch: end, frame: frame, centralBody: .sun,
            gravitationalParameter: gm(), referencePole: SIMD3(0,-1,0), direction: .positive)
        try check(inclined, SIMD3(0,0,1), SIMD3(-1,0,0))
        let angle = 1e-4
        let near = try solve(SIMD3(1,0,0), SIMD3(cos(angle),sin(angle),0), angle,
            options: LambertOptions(absoluteTimeToleranceSeconds: 1e-14))
        try check(near, SIMD3(0,1,0), SIMD3(-sin(angle),cos(angle),0), tolerance: 1e-7)
        let nearPi = Double.pi - 1e-4
        let almostOpposite = try solve(SIMD3(1,0,0), SIMD3(cos(nearPi),sin(nearPi),0), nearPi)
        try check(almostOpposite, SIMD3(0,1,0), SIMD3(-sin(nearPi),cos(nearPi),0), tolerance: 1e-6)
        let reversed = try solve(SIMD3(0,1,0), SIMD3(1,0,0), .pi/2, direction: .negative)
        try check(reversed, SIMD3(1,0,0), SIMD3(0,-1,0))
    }
    @Test func longWayConicsAndSIUnits() throws {
        let e = 0.3, b = sqrt(1-e*e), first = 0.1, last = 5.2
        func ellipse(_ anomaly: Double) -> (SIMD3<Double>, SIMD3<Double>) {
            (SIMD3(cos(anomaly)-e,b*sin(anomaly),0), SIMD3(-sin(anomaly),b*cos(anomaly),0)/(1-e*cos(anomaly)))
        }
        let ef = ellipse(first), es = ellipse(last)
        let elliptic = try solve(ef.0, es.0, last-e*sin(last)-first+e*sin(first))
        #expect(elliptic.branch == .longWay)
        try check(elliptic, ef.1, es.1)
        let f1 = -2.0, f2 = 2.0, d1 = tan(f1/2), d2 = tan(f2/2)
        func parabola(_ f: Double) -> (SIMD3<Double>, SIMD3<Double>) {
            (SIMD3(cos(f),sin(f),0)/(1+cos(f)), SIMD3(-sin(f),1+cos(f),0))
        }
        let pf = parabola(f1), ps = parabola(f2)
        let parabolic = try solve(pf.0, ps.0, (d2+d2*d2*d2/3-d1-d1*d1*d1/3)/2)
        #expect(parabolic.branch == .longWay)
        try check(parabolic, pf.1, ps.1)
        let radius = 7e6, mu = 3.986004418e14, speed = sqrt(mu/radius)
        let seconds = Double.pi/2 * radius/speed
        let arrival = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: seconds))
        let si = try LambertSolver.solve(departurePosition: SIMD3(radius,0,0), arrivalPosition: SIMD3(0,radius,0),
            departureEpoch: .j2000, arrivalEpoch: arrival, frame: frame, centralBody: .earth,
            gravitationalParameter: gm(mu), referencePole: SIMD3(0,0,1), direction: .positive)
        #expect(routeNorm(si.departureState.velocity-SIMD3(0,speed,0)) < 1e-7)
        #expect(routeNorm(si.arrivalState.velocity-SIMD3(-speed,0,0)) < 1e-7)
        let actual = try si.trajectory().state(at: arrival)
        #expect(routeNorm(actual.position-si.arrivalState.position) < 1e-3)
    }
    @Test func explicitFailuresAndBudgets() throws {
        #expect(throws: LambertError.degenerateGeometry) { try solve(SIMD3(1,0,0), SIMD3(-1,0,0), 1) }
        #expect(throws: LambertError.degenerateGeometry) { try solve(SIMD3(1,0,0), SIMD3(2,0,0), 1) }
        #expect(throws: LambertError.invalidInput) { try solve(.zero, SIMD3(1,0,0), 1) }
        #expect(throws: LambertError.invalidInput) { try solve(SIMD3(1,0,0), SIMD3(0,1,0), -1) }
        #expect(throws: LambertError.invalidOptions) { try LambertOptions(maximumIterations: 0) }
        do {
            _ = try solve(SIMD3(1,0,0), SIMD3(0,1,0), .pi/2, options: LambertOptions(maximumIterations: 1))
            Issue.record("Expected iteration limit")
        } catch LambertError.convergenceFailure(let iterations, _) { #expect(iterations == 1) }
        do {
            _ = try solve(SIMD3(1,0,0), SIMD3(0,1,0), 1e4, options: LambertOptions(maximumBracketExpansions: 1))
            Issue.record("Expected bracket limit")
        } catch LambertError.bracketFailure(let evaluations) { #expect(evaluations == 2) }
        let end = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: 1))
        #expect(throws: LambertError.unsupportedRevolutions) {
            try LambertSolver.solve(departurePosition: SIMD3(1,0,0), arrivalPosition: SIMD3(0,1,0),
                departureEpoch: .j2000, arrivalEpoch: end, frame: frame, centralBody: .sun,
                gravitationalParameter: gm(), referencePole: SIMD3(0,0,1), direction: .positive, revolutions: 33)
        }
        #expect(throws: LambertError.ambiguousDirection) {
            try LambertSolver.solve(departurePosition: SIMD3(1,0,0), arrivalPosition: SIMD3(0,1,0),
                departureEpoch: .j2000, arrivalEpoch: end, frame: frame, centralBody: .sun,
                gravitationalParameter: gm(), referencePole: SIMD3(1,0,0), direction: .positive)
        }
    }
}
