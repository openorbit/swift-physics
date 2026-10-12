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

/// Distinct dynamical phase parameters. Barker's D = tan(nu/2) is
/// dimensionless, whereas elliptic E and hyperbolic H are in radians.
public enum ConicAnomaly: Equatable, Sendable {
    case eccentric(Double)
    case hyperbolic(Double)
    case barker(Double)

    /// Principal phase; elliptic revolution offsets are not retained here.
    public static func from(trueAnomaly: Double, eccentricity e: Double) throws -> ConicAnomaly {
        try Self.validate(e)
        guard trueAnomaly.isFinite else { throw OrbitalElementsError.invalidElements }
        let nu = trueAnomaly.remainder(dividingBy: 2 * .pi)
        let denominator = try OsculatingOrbitalElements.radialDenominator(eccentricity: e, trueAnomaly: nu)
        if e < 1 {
            return .eccentric(2 * atan2(sqrt(1 - e) * sin(nu / 2), sqrt(1 + e) * cos(nu / 2)))
        }
        if e == 1 { return .barker(tan(nu / 2)) }
        let sinhH = (sqrt(e - 1) * sqrt(e + 1)) * (sin(nu) / denominator)
        guard sinhH.isFinite else { throw OrbitalElementsError.numericalRangeExceeded }
        return .hyperbolic(asinh(sinhH))
    }

    /// Solves the appropriate equation. For e==1 the input is Barker's
    /// B = D + D³/3, not an elliptic/hyperbolic mean anomaly.
    public static func from(meanAnomaly: Double, eccentricity e: Double,
                            options: KeplerSolverOptions = .default) throws -> ConicAnomaly {
        try Self.validate(e)
        guard meanAnomaly.isFinite else { throw KeplerSolverError.nonFiniteMeanAnomaly }
        if e < 1 {
            return .eccentric(try KeplerSolver.elliptic(meanAnomaly: meanAnomaly, eccentricity: e,
                                                       options: options).anomaly)
        }
        if e > 1 {
            return .hyperbolic(try KeplerSolver.hyperbolic(meanAnomaly: meanAnomaly, eccentricity: e,
                                                         options: options).anomaly)
        }
        // D=2*sinh(asinh(3B/2)/3), evaluated without overflowing 3B/2.
        let magnitude = abs(meanAnomaly)
        if magnitude < 1e-8 { return .barker(meanAnomaly) }
        let argument = magnitude > Double.greatestFiniteMagnitude / 1.5
            ? log(magnitude) + log(3.0) : asinh(1.5 * magnitude)
        let d = 2 * sinh(argument / 3)
        return .barker(meanAnomaly < 0 ? -d : d)
    }

    public func trueAnomaly(eccentricity e: Double) throws -> Double {
        try Self.validate(e)
        let result: Double
        switch self {
        case .eccentric(let anomaly):
            guard e < 1, anomaly.isFinite else { throw OrbitalElementsError.invalidElements }
            let x = anomaly.remainder(dividingBy: 2 * .pi)
            result = 2 * atan2(sqrt(1 + e) * sin(x / 2), sqrt(1 - e) * cos(x / 2))
        case .hyperbolic(let anomaly):
            guard e > 1, anomaly.isFinite else { throw OrbitalElementsError.invalidElements }
            result = 2 * atan2(sqrt(e + 1) * tanh(anomaly / 2), sqrt(e - 1))
        case .barker(let d):
            guard e == 1, d.isFinite else { throw OrbitalElementsError.invalidElements }
            result = 2 * atan(d)
        }
        // For extremely large unbound parameters, the angle can round to the
        // asymptote: reject rather than silently return an unreachable state.
        _ = try OsculatingOrbitalElements.radialDenominator(eccentricity: e, trueAnomaly: result)
        return result
    }

    /// Elliptic/hyperbolic M, or Barker B=D+D³/3 for a parabola.
    public func meanAnomaly(eccentricity e: Double) throws -> Double {
        try Self.validate(e)
        let result: Double
        switch self {
        case .eccentric(let anomaly):
            guard e < 1, anomaly.isFinite else { throw OrbitalElementsError.invalidElements }
            let phase = anomaly.remainder(dividingBy: 2 * .pi)
            let reduced = (1 - e) * phase + e * KeplerSolver.sineDifference(phase)
            result = (anomaly - phase) + reduced
        case .hyperbolic(let anomaly):
            guard e > 1, anomaly.isFinite else { throw OrbitalElementsError.invalidElements }
            result = (e - 1) * anomaly + e * KeplerSolver.sinhDifference(anomaly)
        case .barker(let d):
            guard e == 1, d.isFinite else { throw OrbitalElementsError.invalidElements }
            result = d + (d / 3) * d * d
        }
        guard result.isFinite else { throw OrbitalElementsError.numericalRangeExceeded }
        return result
    }

    private static func validate(_ e: Double) throws {
        guard e.isFinite, e >= 0 else { throw OrbitalElementsError.invalidElements }
    }
}
