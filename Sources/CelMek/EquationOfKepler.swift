//
// SPDX-License-Identifier: Apache-2.0
//
// Copyright 2026 Mattias Holm
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
import Foundation

public enum KeplerSolverError: Error, Equatable {
    case invalidEccentricity
    case nonFiniteMeanAnomaly
    case invalidOptions
    case numericalRangeExceeded
    case convergenceFailure(iterations: Int, lowerBound: Double, upperBound: Double)
}

/// Tolerances apply to the solved anomaly, in radians. The relative tolerance
/// scales with its magnitude; machine representability also limits resolution.
public struct KeplerSolverOptions: Equatable, Sendable {
    public let absoluteTolerance: Double
    public let relativeTolerance: Double
    public let maximumIterations: Int

    public init(absoluteTolerance: Double = 1e-14,
                relativeTolerance: Double = 1e-14,
                maximumIterations: Int = 128) throws {
        guard absoluteTolerance.isFinite, absoluteTolerance > 0,
              relativeTolerance.isFinite, relativeTolerance >= 0,
              maximumIterations > 0 else { throw KeplerSolverError.invalidOptions }
        self.absoluteTolerance = absoluteTolerance
        self.relativeTolerance = relativeTolerance
        self.maximumIterations = maximumIterations
    }

    public static let `default` = try! KeplerSolverOptions()
}

public struct KeplerSolution: Equatable, Sendable {
    /// Elliptic E or hyperbolic H, in radians.
    public let anomaly: Double
    public let iterations: Int
}

public enum KeplerSolver {
    /// Solves M = E - e sin(E), for 0 <= e < 1. The result retains the
    /// revolution offset of M. Very large M loses phase information in Double.
    public static func elliptic(meanAnomaly: Double, eccentricity: Double,
                                options: KeplerSolverOptions = .default) throws -> KeplerSolution {
        guard eccentricity.isFinite, eccentricity >= 0, eccentricity < 1 else {
            throw KeplerSolverError.invalidEccentricity
        }
        guard meanAnomaly.isFinite else { throw KeplerSolverError.nonFiniteMeanAnomaly }
        if eccentricity == 0 { return KeplerSolution(anomaly: meanAnomaly, iterations: 0) }
        let reduced = meanAnomaly.remainder(dividingBy: 2 * .pi)
        let magnitude = abs(reduced)
        if magnitude == 0 { return KeplerSolution(anomaly: meanAnomaly, iterations: 0) }
        let root = try solve(lower: magnitude, upper: .pi, options: options) { x in
            ((1 - eccentricity) * x + eccentricity * sineDifference(x) - magnitude,
             (1 - eccentricity) + 2 * eccentricity * pow(sin(x / 2), 2))
        }
        let signed = reduced < 0 ? -root.anomaly : root.anomaly
        // Add the small correction to M to preserve the revolution offset.
        return KeplerSolution(anomaly: meanAnomaly + (signed - reduced), iterations: root.iterations)
    }

    /// Solves M = e sinh(H) - H, for e > 1. M is dimensionless and is
    /// not periodic. Parabolic Barker/universal-variable equations are separate.
    public static func hyperbolic(meanAnomaly: Double, eccentricity: Double,
                                  options: KeplerSolverOptions = .default) throws -> KeplerSolution {
        guard eccentricity.isFinite, eccentricity > 1 else { throw KeplerSolverError.invalidEccentricity }
        guard meanAnomaly.isFinite else { throw KeplerSolverError.nonFiniteMeanAnomaly }
        if meanAnomaly == 0 { return KeplerSolution(anomaly: meanAnomaly, iterations: 0) }
        let magnitude = abs(meanAnomaly)
        // For extreme inputs, fail explicitly rather than overflow an intermediate.
        let upper = asinh(magnitude / eccentricity) + 2
        let root = try solve(lower: 0, upper: upper, options: options) { x in
            let value = (eccentricity - 1) * x + eccentricity * sinhDifference(x)
            let derivative = (eccentricity - 1) + 2 * eccentricity * pow(sinh(x / 2), 2)
            return (value - magnitude, derivative)
        }
        return KeplerSolution(anomaly: meanAnomaly < 0 ? -root.anomaly : root.anomaly,
                              iterations: root.iterations)
    }

    private static func solve(lower: Double, upper: Double, options: KeplerSolverOptions,
                              evaluate: (Double) -> (Double, Double)) throws -> KeplerSolution {
        var low = lower
        var high = upper
        let endpoints = (evaluate(low).0, evaluate(high).0)
        guard endpoints.0.isFinite, endpoints.1.isFinite else { throw KeplerSolverError.numericalRangeExceeded }
        if endpoints.0 == 0 { return KeplerSolution(anomaly: low, iterations: 0) }
        if endpoints.1 == 0 { return KeplerSolution(anomaly: high, iterations: 0) }
        guard endpoints.0 < 0, endpoints.1 > 0 else { throw KeplerSolverError.numericalRangeExceeded }
        var x = low + (high - low) / 2
        for iteration in 1...options.maximumIterations {
            let (residual, derivative) = evaluate(x)
            guard residual.isFinite, derivative.isFinite, derivative > 0 else {
                throw KeplerSolverError.numericalRangeExceeded
            }
            if residual == 0 { return KeplerSolution(anomaly: x, iterations: iteration) }
            if residual < 0 { low = x } else { high = x }
            let midpoint = low + (high - low) / 2
            let tolerance = options.absoluteTolerance + options.relativeTolerance * abs(midpoint)
            if (high - low) / 2 <= tolerance || midpoint == low || midpoint == high {
                return KeplerSolution(anomaly: midpoint, iterations: iteration)
            }
            let newton = x - residual / derivative
            // Require progress within the middle half of the bracket. Bisection
            // otherwise guarantees contraction even close to a parabolic orbit.
            let quarter = (high - low) / 4
            x = newton > low + quarter && newton < high - quarter ? newton : midpoint
        }
        throw KeplerSolverError.convergenceFailure(iterations: options.maximumIterations,
                                                   lowerBound: low, upperBound: high)
    }

    static func sineDifference(_ x: Double) -> Double {
        if abs(x) >= 0.5 { return x - sin(x) }
        var term = x * x * x / 6
        var sum = term
        for n in 2...10 {
            term *= -x * x / Double((2 * n) * (2 * n + 1))
            sum += term
        }
        return sum
    }

    static func sinhDifference(_ x: Double) -> Double {
        if abs(x) >= 0.5 { return sinh(x) - x }
        var term = x * x * x / 6
        var sum = term
        for n in 2...10 {
            term *= x * x / Double((2 * n) * (2 * n + 1))
            sum += term
        }
        return sum
    }
}
