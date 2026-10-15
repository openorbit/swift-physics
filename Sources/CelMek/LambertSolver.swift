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
import Math

public enum LambertError: Error, Equatable {
    case invalidInput, invalidOptions, unsupportedRevolutions, degenerateGeometry
    case ambiguousDirection, numericalRangeExceeded, branchRequired, invalidBranch
    case belowMinimumFlightTime(minimumSeconds: Double)
    case minimumSearchFailure(iterations: Int)
    case bracketFailure(evaluations: Int)
    case convergenceFailure(iterations: Int, timeResidualSeconds: Double)
}
/// Sign of angular momentum projected on the caller's reference pole.
public enum LambertDirection: Equatable, Sendable { case positive, negative }
public enum LambertBranch: Equatable, Sendable { case shortWay, longWay }
/// Root on the decreasing/increasing side of the universal time curve's minimum.
/// These names do not imply an energy, focus-position or delta-v ranking.
public enum LambertMultiRevolutionBranch: Equatable, Sendable {
    case lowerUniversalParameter, upperUniversalParameter
}
public struct LambertOptions: Equatable, Sendable {
    public let absoluteTimeToleranceSeconds: Double
    public let relativeTimeTolerance: Double
    /// Separate cap for root bisection and (when needed) minimum-time search.
    public let maximumIterations: Int
    public let maximumBracketExpansions: Int
    public init(absoluteTimeToleranceSeconds: Double = 1e-9, relativeTimeTolerance: Double = 1e-13,
        maximumIterations: Int = 128, maximumBracketExpansions: Int = 64) throws {
        guard absoluteTimeToleranceSeconds.isFinite, absoluteTimeToleranceSeconds > 0,
              relativeTimeTolerance.isFinite, relativeTimeTolerance >= 0,
              (1...10_000).contains(maximumIterations), (1...1_000).contains(maximumBracketExpansions) else {
            throw LambertError.invalidOptions
        }
        self.absoluteTimeToleranceSeconds = absoluteTimeToleranceSeconds
        self.relativeTimeTolerance = relativeTimeTolerance
        self.maximumIterations = maximumIterations
        self.maximumBracketExpansions = maximumBracketExpansions
    }
    public static let `default` = try! LambertOptions()
}
public struct LambertSolution {
    public let departureState: CartesianState
    public let arrivalState: CartesianState
    public let centralBody: BodyIdentifier
    public let gravitationalParameter: GravitationalParameter
    public let branch: LambertBranch
    public let direction: LambertDirection
    public let transferAngleRadians: Double
    public let orbitalNormal: SIMD3<Double>
    public let iterations: Int
    public let bracketEvaluations: Int
    /// Signed universal time-equation residual, not a position-error bound.
    public let timeResidualSeconds: Double
    public let revolutions: Int
    public let multiRevolutionBranch: LambertMultiRevolutionBranch?
    public let minimumFlightTimeSeconds: Double?
    public let minimumSearchIterations: Int
    public let branchesCoalesced: Bool
    public let universalParameter: Double
    public func trajectory(options: TwoBodyPropagationOptions = .default) throws -> TwoBodyPropagator {
        try TwoBodyPropagator(initialState: departureState, centralBody: centralBody,
            gravitationalParameter: gravitationalParameter, options: options)
    }
}

/// Scaled universal-variable Lambert boundary-value solution.
/// Positions are central-body-relative metres in one explicitly supplied inertial frame.
public enum LambertSolver {
    public static func solve(departurePosition: SIMD3<Double>, arrivalPosition: SIMD3<Double>,
        departureEpoch: Epoch, arrivalEpoch: Epoch, frame: Frame, centralBody: BodyIdentifier,
        gravitationalParameter gm: GravitationalParameter, referencePole: SIMD3<Double>,
        direction: LambertDirection, revolutions: Int = 0,
        multiRevolutionBranch: LambertMultiRevolutionBranch? = nil, options: LambertOptions = .default) throws -> LambertSolution {
        guard (0...32).contains(revolutions) else { throw LambertError.unsupportedRevolutions }
        if revolutions > 0 && multiRevolutionBranch == nil { throw LambertError.branchRequired }
        if revolutions == 0 && multiRevolutionBranch != nil { throw LambertError.invalidBranch }
        guard frame.pseudoInertial else { throw PropagationError.nonInertialFrame }
        guard gm.convention == .newtonian || gm.convention == .ttCompatible else { throw PropagationError.incompatibleGMConvention }
        try centralBody.validate()
        guard arrivalEpoch > departureEpoch,
              [departurePosition, arrivalPosition, referencePole].allSatisfy({ [$0.x,$0.y,$0.z].allSatisfy(\.isFinite) }) else {
            throw LambertError.invalidInput
        }
        let radius1 = routeNorm(departurePosition), radius2 = routeNorm(arrivalPosition), poleNorm = routeNorm(referencePole)
        guard radius1 > 0, radius2 > 0, poleNorm > 0,
              radius1.isFinite, radius2.isFinite, poleNorm.isFinite else { throw LambertError.invalidInput }
        let u = departurePosition / radius1, v = arrivalPosition / radius2
        let cross = SIMD3(u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x)
        let sine = routeNorm(cross), cosine = max(-1.0, min(1.0, (u*v).sum()))
        // The plane and f/g reconstruction become ill-conditioned at 0 and pi.
        guard sine > 1e-8 else { throw LambertError.degenerateGeometry }
        let projected = (cross / sine * (referencePole / poleNorm)).sum()
        guard abs(projected) > 1e-12 else { throw LambertError.ambiguousDirection }
        let short = (projected > 0) == (direction == .positive)
        let angle = short ? atan2(sine, cosine) : 2 * Double.pi - atan2(sine, cosine)
        let scale = max(radius1, radius2), r1 = radius1 / scale, r2 = radius2 / scale
        let speedScale = sqrt(gm.cubicMetresPerSecondSquared) / sqrt(scale), timeScale = scale / speedScale
        let seconds = try arrivalEpoch.duration(since: departureEpoch).seconds, tau = seconds / timeScale
        guard speedScale.isFinite, speedScale > 0, timeScale.isFinite, timeScale > 0, tau.isFinite, tau > 0 else {
            throw LambertError.numericalRangeExceeded
        }
        // Equivalent to sin(theta)*sqrt(r1*r2/(1-cos(theta))), without small-angle division.
        let a = (short ? 1.0 : -1.0) * sqrt(r1 * r2) * (cosine < 0 ? sine / sqrt(1 - cosine) : sqrt(1 + cosine))
        let cosHalf = cosine < 0 ? sine / sqrt(2 * (1 - cosine)) : sqrt((1 + cosine)/2)
        let sinHalf = cosine > 0 ? sine / sqrt(2 * (1 + cosine)) : sqrt((1 - cosine)/2)
        let yAtZero = short ? pow(sqrt(r1)-sqrt(r2),2) +
            2 * sqrt(r1*r2) * sinHalf*sinHalf / (1 + cosHalf) : r1+r2-sqrt(2)*a
        func value(_ z: Double) throws -> (time: Double, y: Double)? {
            let c: Double, s: Double
            if abs(z) < 0.5 {
                var cc = 0.5, ss = 1.0 / 6, tc = cc, ts = ss
                for n in 1...16 {
                    tc *= -z / Double((2*n+1)*(2*n+2)); ts *= -z / Double((2*n+2)*(2*n+3))
                    cc += tc; ss += ts
                }
                c = cc; s = ss
            } else {
                let x = sqrt(abs(z))
                if z > 0 {
                    let local = x - 2 * Double.pi * Double(revolutions)
                    c = 2 * pow(sin(local/2),2) / z; s = (x-sin(local))/(x*x*x)
                }
                else { c = 2 * pow(sinh(x/2),2) / -z; s = (sinh(x)-x)/(x*x*x) }
            }
            guard c.isFinite, s.isFinite, c > 0 else { throw LambertError.numericalRangeExceeded }
            // Use the local anomaly within this revolution interval. The sign
            // from sqrt(C) cancels the completed-turn parity. Half-angle differences
            // avoid subtracting almost equal radii on short, small-angle arcs.
            let quarter = (sqrt(abs(z)) - (z > 0 ? 2 * Double.pi * Double(revolutions) : 0))/4
            let correction = z >= 0 ? 2 * pow(sin(quarter),2) : -2 * pow(sinh(quarter),2)
            let y = yAtZero + sqrt(2)*a*correction
            if y <= 0 { return nil } // Below the physical time-of-flight domain.
            let time = pow(y/c,1.5)*s + a*sqrt(y)
            guard y.isFinite, time.isFinite, time >= 0 else { throw LambertError.numericalRangeExceeded }
            return (time,y)
        }
        let tolerance = options.absoluteTimeToleranceSeconds / timeScale + options.relativeTimeTolerance * tau
        guard tolerance.isFinite, tolerance > 0 else { throw LambertError.numericalRangeExceeded }
        var lower = 0.0, upper = 0.0, bracketEvaluations = 0
        var minimumTime: Double?, minimumIterations = 0, coalesced = false
        var startingResidual = 0.0
        if revolutions == 0 {
            let limit = 4 * Double.pi * Double.pi
            bracketEvaluations = 1
            guard let zero = try value(0) else { throw LambertError.numericalRangeExceeded }
            startingResidual = zero.time - tau
            if zero.time < tau {
                upper = limit / 2
                var found = false
                for _ in 0..<options.maximumBracketExpansions {
                    bracketEvaluations += 1
                    if let test = try value(upper), test.time >= tau { found = true; break }
                    let next = (upper + limit)/2
                    if next == upper || next >= limit { break }
                    upper = next
                }
                guard found else { throw LambertError.bracketFailure(evaluations: bracketEvaluations) }
            } else {
                lower = -4
                var found = false
                for _ in 0..<options.maximumBracketExpansions {
                    bracketEvaluations += 1
                    let test = try value(lower)
                    if test == nil || test!.time <= tau { found = true; break }
                    // Bound hyperbolic cancellation/overflow; extreme flights are explicit failures.
                    if lower <= -4096 { break }
                    lower *= 2
                }
                guard found else { throw LambertError.bracketFailure(evaluations: bracketEvaluations) }
            }
        } else {
            let edge1 = pow(2 * Double.pi * Double(revolutions),2)
            let edge2 = pow(2 * Double.pi * Double(revolutions+1),2)
            var left = edge1, right = edge2
            let ratio = (sqrt(5.0)-1)/2
            var x1 = right-ratio*(right-left), x2 = left+ratio*(right-left)
            guard var t1 = try value(x1), var t2 = try value(x2) else { throw LambertError.numericalRangeExceeded }
            var minimumFound = false
            for iteration in 1...options.maximumIterations {
                minimumIterations = iteration
                if right-left <= 1e-12 * max(1, abs((left+right)/2)) { minimumFound = true; break }
                if t1.time < t2.time {
                    right = x2; x2 = x1; t2 = t1
                    x1 = right-ratio*(right-left)
                    guard let test = try value(x1) else { throw LambertError.numericalRangeExceeded }; t1 = test
                } else {
                    left = x1; x1 = x2; t1 = t2
                    x2 = left+ratio*(right-left)
                    guard let test = try value(x2) else { throw LambertError.numericalRangeExceeded }; t2 = test
                }
            }
            guard minimumFound else { throw LambertError.minimumSearchFailure(iterations: minimumIterations) }
            let minimumZ = (left+right)/2
            guard let minimum = try value(minimumZ) else { throw LambertError.numericalRangeExceeded }
            minimumTime = minimum.time * timeScale
            startingResidual = minimum.time-tau
            guard tau >= minimum.time-tolerance else {
                throw LambertError.belowMinimumFlightTime(minimumSeconds: minimum.time*timeScale)
            }
            coalesced = abs(minimum.time-tau) <= tolerance
            if coalesced { lower = minimumZ; upper = minimumZ }
            else {
                let lowerBranch = multiRevolutionBranch == .lowerUniversalParameter
                let edge = lowerBranch ? edge1 : edge2
                var outer = (edge+minimumZ)/2, found = false
                for _ in 0..<options.maximumBracketExpansions {
                    bracketEvaluations += 1
                    if let test = try value(outer), test.time >= tau { found = true; break }
                    let next = (edge+outer)/2
                    if next == outer || next == edge { break }
                    outer = next
                }
                guard found else { throw LambertError.bracketFailure(evaluations: bracketEvaluations) }
                lower = lowerBranch ? outer : minimumZ
                upper = lowerBranch ? minimumZ : outer
            }
        }
        var performed = 0
        var residual = startingResidual
        for iteration in 1...options.maximumIterations {
            performed = iteration
            let z = (lower+upper)/2
            guard let test = try value(z) else { lower = z; continue }
            residual = test.time - tau
            if abs(residual) <= tolerance {
                let f = 1 - test.y/r1, g = a*sqrt(test.y), gDot = 1 - test.y/r2
                guard g.isFinite, g != 0 else { throw LambertError.numericalRangeExceeded }
                let p1 = departurePosition / scale, p2 = arrivalPosition / scale
                let v1 = (p2-f*p1)/g * speedScale, v2 = (gDot*p2-p1)/g * speedScale
                guard [v1,v2].allSatisfy({ [$0.x,$0.y,$0.z].allSatisfy(\.isFinite) }) else {
                    throw LambertError.numericalRangeExceeded
                }
                return LambertSolution(departureState: try CartesianState(epoch: departureEpoch, frame: frame,
                    position: departurePosition, velocity: v1), arrivalState: try CartesianState(epoch: arrivalEpoch,
                    frame: frame, position: arrivalPosition, velocity: v2), centralBody: centralBody,
                    gravitationalParameter: gm, branch: short ? .shortWay : .longWay, direction: direction,
                    transferAngleRadians: angle, orbitalNormal: cross / sine * (short ? 1 : -1),
                    iterations: iteration, bracketEvaluations: bracketEvaluations, timeResidualSeconds: residual * timeScale,
                    revolutions: revolutions, multiRevolutionBranch: multiRevolutionBranch,
                    minimumFlightTimeSeconds: minimumTime, minimumSearchIterations: minimumIterations,
                    branchesCoalesced: coalesced, universalParameter: z)
            }
            let increasing = revolutions == 0 || multiRevolutionBranch == .upperUniversalParameter
            if (residual < 0) == increasing { lower = z } else { upper = z }
            if lower == upper || (lower+upper)/2 == z { break }
        }
        throw LambertError.convergenceFailure(iterations: performed, timeResidualSeconds: residual * timeScale)
    }

    /// Enumerates both roots for a requested positive revolution count; one result
    /// when they coalesce within the requested time tolerance. Budgets apply per solve.
    public static func solveBranches(departurePosition: SIMD3<Double>, arrivalPosition: SIMD3<Double>,
        departureEpoch: Epoch, arrivalEpoch: Epoch, frame: Frame, centralBody: BodyIdentifier,
        gravitationalParameter: GravitationalParameter, referencePole: SIMD3<Double>,
        direction: LambertDirection, revolutions: Int, options: LambertOptions = .default) throws -> [LambertSolution] {
        guard revolutions > 0 else { throw LambertError.unsupportedRevolutions }
        let lower = try solve(departurePosition: departurePosition, arrivalPosition: arrivalPosition,
            departureEpoch: departureEpoch, arrivalEpoch: arrivalEpoch, frame: frame, centralBody: centralBody,
            gravitationalParameter: gravitationalParameter, referencePole: referencePole, direction: direction,
            revolutions: revolutions, multiRevolutionBranch: .lowerUniversalParameter, options: options)
        if lower.branchesCoalesced { return [lower] }
        let upper = try solve(departurePosition: departurePosition, arrivalPosition: arrivalPosition,
            departureEpoch: departureEpoch, arrivalEpoch: arrivalEpoch, frame: frame, centralBody: centralBody,
            gravitationalParameter: gravitationalParameter, referencePole: referencePole, direction: direction,
            revolutions: revolutions, multiRevolutionBranch: .upperUniversalParameter, options: options)
        return [lower,upper]
    }

}
