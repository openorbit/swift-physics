import Foundation
import Testing
@testable import CelMek

@Test func ellipticKeplerKnownAnomalies() throws {
    for e in [0.0, 0.1, 0.9, 0.999999999999] {
        for expected in [-3.0, -1.0, -0.001, 0, 0.001, 1, 3] {
            let m = expected - e * sin(expected)
            let result = try KeplerSolver.elliptic(meanAnomaly: m, eccentricity: e)
            #expect(abs(result.anomaly - expected) < 2e-10)
            #expect(result.iterations <= 128)
        }
    }
    let m = 20 * Double.pi + 0.8
    let result = try KeplerSolver.elliptic(meanAnomaly: m, eccentricity: 0.7)
    #expect(abs(result.anomaly - 0.7 * sin(result.anomaly) - m) < 1e-12)
    #expect(try KeplerSolver.elliptic(meanAnomaly: 0, eccentricity: 0.999).anomaly == 0)
}

@Test func hyperbolicKeplerKnownAnomalies() throws {
    for e in [1.000000000001, 1.1, 2.0, 100.0] {
        for expected in [-10.0, -1, -0.001, 0, 0.001, 1, 10] {
            let m = e * sinh(expected) - expected
            let result = try KeplerSolver.hyperbolic(meanAnomaly: m, eccentricity: e)
            #expect(abs(result.anomaly - expected) < 2e-10)
        }
    }
}

@Test func keplerSolverRejectsInvalidInputsAndReportsExhaustion() throws {
    #expect(throws: KeplerSolverError.invalidEccentricity) {
        try KeplerSolver.elliptic(meanAnomaly: 1, eccentricity: 1)
    }
    #expect(throws: KeplerSolverError.invalidEccentricity) {
        try KeplerSolver.hyperbolic(meanAnomaly: 1, eccentricity: 0.5)
    }
    #expect(throws: KeplerSolverError.nonFiniteMeanAnomaly) {
        try KeplerSolver.elliptic(meanAnomaly: .nan, eccentricity: 0.5)
    }
    #expect(throws: KeplerSolverError.invalidOptions) { try KeplerSolverOptions(maximumIterations: 0) }
    #expect(throws: KeplerSolverError.invalidOptions) { try KeplerSolverOptions(absoluteTolerance: .infinity) }
    let options = try KeplerSolverOptions(maximumIterations: 1)
    do {
        _ = try KeplerSolver.elliptic(meanAnomaly: 0.1, eccentricity: 0.99, options: options)
        Issue.record("Expected iteration exhaustion")
    } catch let KeplerSolverError.convergenceFailure(iterations, low, high) {
        #expect(iterations == 1)
        #expect(low < high)
    }
    #expect(throws: KeplerSolverError.numericalRangeExceeded) {
        try KeplerSolver.hyperbolic(meanAnomaly: .greatestFiniteMagnitude, eccentricity: 2)
    }
}

@Test func keplerHighPrecisionReferenceCases() throws {
    // Generated independently by Tools/generate_kepler_fixtures.py using
    // 80-digit Decimal arithmetic and bisection of the defining equations.
    let cases: [(Bool, Double, Double, Double)] = [
        (false, 0.9, 0.1, 0.6308435275631535),
        (false, 0.9999999999999999, 1e-18, 1.8169983969771785e-06),
        (false, 0.999999, 1e-08, 0.003407264597719929),
        (false, 0.5, 2.0, 2.3542427582227807),
        (true, 1.0000000000000002, 1e-18, 1.8168762011230227e-06),
        (true, 1.000001, 1e-08, 0.0034072615353025817),
        (true, 2.0, 10.0, 2.5348145176603545),
        (true, 1.1, 1000.0, 7.513077572718448),
    ]
    for (hyperbolic, e, m, expected) in cases {
        for sign in [-1.0, 1.0] {
            let result = try hyperbolic
                ? KeplerSolver.hyperbolic(meanAnomaly: sign * m, eccentricity: e)
                : KeplerSolver.elliptic(meanAnomaly: sign * m, eccentricity: e)
            #expect(abs(result.anomaly - sign * expected) <= 2e-14 * (1 + abs(expected)))
        }
    }
}
