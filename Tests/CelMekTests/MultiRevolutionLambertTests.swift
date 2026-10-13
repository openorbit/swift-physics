// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import Math
@testable import CelMek

@Suite struct MultiRevolutionLambertTests {
    let frame = Frame(name: "Multi-revolution Lambert", pseudoInertial: true)
    func end(_ seconds: Double) throws -> Epoch { try Epoch.j2000.advanced(by: CelMek.Duration(seconds: seconds)) }
    func gm() throws -> GravitationalParameter {
        try GravitationalParameter(cubicMetresPerSecondSquared: 1, convention: .newtonian,
            provenance: PhysicalDataProvenance(source: "Analytic conics", version: "1"))
    }
    func branches(_ r1: SIMD3<Double>, _ r2: SIMD3<Double>, _ time: Double, _ count: Int,
        direction: LambertDirection = .positive, options: LambertOptions = .default) throws -> [LambertSolution] {
        try LambertSolver.solveBranches(departurePosition: r1, arrivalPosition: r2,
            departureEpoch: .j2000, arrivalEpoch: end(time), frame: frame, centralBody: .sun,
            gravitationalParameter: gm(), referencePole: SIMD3(0,0,1), direction: direction,
            revolutions: count, options: options)
    }
    func verify(_ solution: LambertSolution, _ time: Double) throws {
        let trajectory = try solution.trajectory()
        let endpoint = try trajectory.state(at: solution.arrivalState.epoch)
        #expect(routeNorm(endpoint.position-solution.arrivalState.position) < 1e-7)
        #expect(routeNorm(endpoint.velocity-solution.arrivalState.velocity) < 1e-7)
        #expect(abs(solution.timeResidualSeconds) < 2e-9)
        let samples = 128*(solution.revolutions+1)
        var previous = atan2(solution.departureState.position.y, solution.departureState.position.x), sweep = 0.0
        for index in 1...samples {
            let p = try trajectory.state(at: end(time*Double(index)/Double(samples))).position
            let angle = atan2(p.y,p.x)
            sweep += (angle-previous).remainder(dividingBy: 2 * .pi)
            previous = angle
        }
        let signedSweep = solution.direction == .positive ? sweep : -sweep
        #expect(abs(signedSweep-(2 * .pi * Double(solution.revolutions)+solution.transferAngleRadians)) < 1e-7)
    }
    @Test func bothRootsDirectionsAndRevolutionCounts() throws {
        for count in [1,2,5,32] {
            for direction in [LambertDirection.positive, .negative] {
                let partial = direction == .positive ? Double.pi/2 : 3 * Double.pi/2
                let time = 2 * Double.pi * Double(count) + partial
                let solutions = try branches(SIMD3(1,0,0), SIMD3(0,1,0), time, count, direction: direction)
                #expect(solutions.count == 2)
                #expect(solutions[0].multiRevolutionBranch == .lowerUniversalParameter)
                #expect(solutions[1].multiRevolutionBranch == .upperUniversalParameter)
                #expect(solutions[0].universalParameter < solutions[1].universalParameter)
                #expect(routeNorm(solutions[0].departureState.velocity-solutions[1].departureState.velocity) > 1e-3)
                let expected = SIMD3<Double>(0,direction == .positive ? 1 : -1,0)
                #expect(solutions.map { routeNorm($0.departureState.velocity-expected) }.min()! < 1e-8)
                for solution in solutions { try verify(solution,time) }
            }
        }
    }
    @Test func independentEccentricOrbitWithTwoRevolutions() throws {
        let e = 0.4, b = sqrt(1-e*e), first = 0.2, last = 2.1
        let r1 = SIMD3<Double>(cos(first)-e,b*sin(first),0), r2 = SIMD3<Double>(cos(last)-e,b*sin(last),0)
        let v1 = SIMD3<Double>(-sin(first),b*cos(first),0)/(1-e*cos(first))
        let v2 = SIMD3<Double>(-sin(last),b*cos(last),0)/(1-e*cos(last))
        let time = last-e*sin(last)-first+e*sin(first)+4 * Double.pi
        let solutions = try branches(r1,r2,time,2)
        let known = solutions.min { routeNorm($0.departureState.velocity-v1) < routeNorm($1.departureState.velocity-v1) }!
        #expect(routeNorm(known.departureState.velocity-v1) < 1e-8)
        #expect(routeNorm(known.arrivalState.velocity-v2) < 1e-8)
        for solution in solutions { try verify(solution,time) }
    }
    @Test func minimumCoalescenceAndInfeasibleTime() throws {
        let supplied = try branches(SIMD3(1,0,0),SIMD3(0,1,0),10,1)
        let minimum = try #require(supplied[0].minimumFlightTimeSeconds)
        #expect(minimum > 0 && minimum < 10)
        let merged = try branches(SIMD3(1,0,0),SIMD3(0,1,0),minimum,1)
        #expect(merged.count == 1)
        #expect(merged[0].branchesCoalesced)
        try verify(merged[0],minimum)
        let near = try branches(SIMD3(1,0,0),SIMD3(0,1,0),minimum+1e-5,1)
        #expect(near.count == 2)
        #expect(near.allSatisfy { !$0.branchesCoalesced })
        for solution in near { try verify(solution,minimum+1e-5) }
        do {
            _ = try branches(SIMD3(1,0,0),SIMD3(0,1,0),minimum-0.01,1)
            Issue.record("Expected infeasible flight time")
        } catch LambertError.belowMinimumFlightTime(let reported) { #expect(abs(reported-minimum) < 1e-10) }
    }
    @Test func explicitBranchAndWorkLimits() throws {
        let arrival = try end(10)
        #expect(throws: LambertError.branchRequired) {
            try LambertSolver.solve(departurePosition: SIMD3(1,0,0), arrivalPosition: SIMD3(0,1,0),
                departureEpoch: .j2000, arrivalEpoch: arrival, frame: frame, centralBody: .sun,
                gravitationalParameter: gm(), referencePole: SIMD3(0,0,1), direction: .positive, revolutions: 1)
        }
        #expect(throws: LambertError.invalidBranch) {
            try LambertSolver.solve(departurePosition: SIMD3(1,0,0), arrivalPosition: SIMD3(0,1,0),
                departureEpoch: .j2000, arrivalEpoch: arrival, frame: frame, centralBody: .sun,
                gravitationalParameter: gm(), referencePole: SIMD3(0,0,1), direction: .positive,
                multiRevolutionBranch: .lowerUniversalParameter)
        }
        #expect(throws: LambertError.minimumSearchFailure(iterations: 1)) {
            try branches(SIMD3(1,0,0),SIMD3(0,1,0),10,1,options: LambertOptions(maximumIterations: 1))
        }
        do {
            _ = try branches(SIMD3(1,0,0),SIMD3(0,1,0),1e6,1,options: LambertOptions(maximumBracketExpansions: 1))
            Issue.record("Expected bracket failure")
        } catch LambertError.bracketFailure(let evaluations) { #expect(evaluations == 1) }
    }
}
