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

public enum PropagationError: Error, Equatable {
    case nonInertialFrame, incompatibleGMConvention, invalidOptions
    case zeroPosition, radialTrajectory, numericalRangeExceeded
    case bracketFailure
    case convergenceFailure(iterations: Int)
}

/// A trajectory queried at an absolute epoch, independently of previous queries.
public protocol Propagator: StateProvider {
    var referenceEpoch: Epoch { get }
    var centralBody: BodyIdentifier { get }
}

public struct TwoBodyPropagationOptions: Equatable, Sendable {
    /// Universal time-equation residual tolerances, not physical position bounds.
    public let absoluteTimeToleranceSeconds: Double
    public let relativeTimeTolerance: Double
    public let maximumIterations: Int
    public let maximumBracketExpansions: Int
    public init(absoluteTimeToleranceSeconds: Double = 1e-9,
                relativeTimeTolerance: Double = 1e-13,
                maximumIterations: Int = 128, maximumBracketExpansions: Int = 128) throws {
        guard absoluteTimeToleranceSeconds.isFinite, absoluteTimeToleranceSeconds > 0,
              relativeTimeTolerance.isFinite, relativeTimeTolerance >= 0,
              maximumIterations > 0, maximumBracketExpansions > 0 else {
            throw PropagationError.invalidOptions
        }
        self.absoluteTimeToleranceSeconds = absoluteTimeToleranceSeconds
        self.relativeTimeTolerance = relativeTimeTolerance
        self.maximumIterations = maximumIterations
        self.maximumBracketExpansions = maximumBracketExpansions
    }
    public static let `default` = try! TwoBodyPropagationOptions()
}

public struct TwoBodyPropagationResult {
    public let state: CartesianState
    public let iterations: Int
    public let bracketExpansions: Int
}

/// Newtonian point-mass motion in central-body-relative, nonrotating axes.
/// Uses scaled universal variables without a near-parabolic conic switch.
public struct TwoBodyPropagator: Propagator {
    public let initialState: CartesianState
    public let centralBody: BodyIdentifier
    public let gravitationalParameter: GravitationalParameter
    public let options: TwoBodyPropagationOptions
    public var frame: Frame { initialState.frame }
    public var referenceEpoch: Epoch { initialState.epoch }
    private let radiusScale: Double
    private let velocityScale: Double
    private let timeScale: Double
    private let radial: SIMD3<Double>
    private let scaledVelocity: SIMD3<Double>
    private let alpha: Double
    private let radialVelocity: Double

    public init(initialState: CartesianState, centralBody: BodyIdentifier,
                gravitationalParameter: GravitationalParameter,
                options: TwoBodyPropagationOptions = .default) throws {
        guard initialState.frame.pseudoInertial else { throw PropagationError.nonInertialFrame }
        // Epoch differences and state derivatives use uniform atomic/TT seconds.
        guard gravitationalParameter.convention == .newtonian || gravitationalParameter.convention == .ttCompatible else {
            throw PropagationError.incompatibleGMConvention
        }
        try centralBody.validate()
        let r = Self.norm(initialState.position)
        guard r > 0 else { throw PropagationError.zeroPosition }
        let scale = sqrt(gravitationalParameter.cubicMetresPerSecondSquared) / sqrt(r)
        let t = r / scale
        let u = initialState.position / r, v = initialState.velocity / scale
        let a = 2 - (v * v).sum(), b = (u * v).sum()
        guard r.isFinite, scale.isFinite, scale > 0, t.isFinite, t > 0, a.isFinite, b.isFinite else {
            throw PropagationError.numericalRangeExceeded
        }
        guard Self.norm(Self.cross(u, v)) > 0 else { throw PropagationError.radialTrajectory }
        self.initialState = initialState; self.centralBody = centralBody
        self.gravitationalParameter = gravitationalParameter; self.options = options
        radiusScale = r; velocityScale = scale; timeScale = t
        radial = u; scaledVelocity = v; alpha = a; radialVelocity = b
    }

    public func state(at epoch: Epoch) throws -> CartesianState { try propagation(at: epoch).state }

    public func propagation(at epoch: Epoch) throws -> TwoBodyPropagationResult {
        var tau = try epoch.duration(since: referenceEpoch).seconds / timeScale
        guard tau.isFinite else { throw PropagationError.numericalRangeExceeded }
        // Bound trigonometric arguments and f/g cancellation on multi-orbit queries.
        if alpha > 0 {
            let period = (2 * Double.pi / alpha) / sqrt(alpha)
            if period.isFinite { tau = tau.remainder(dividingBy: period) }
        }
        if tau == 0 {
            return TwoBodyPropagationResult(state: try makeState(epoch: epoch, position: initialState.position,
                                                                 velocity: initialState.velocity),
                                             iterations: 0, bracketExpansions: 0)
        }
        let sign = tau < 0 ? -1.0 : 1.0
        let magnitude = abs(tau)
        func evaluate(_ positiveAnomaly: Double) -> (residual: Double, derivative: Double)? {
            let x = sign * positiveAnomaly
            guard let functions = Self.stumpff(alpha * x * x) else { return nil }
            let x2c = x * x * functions.c, x3s = x * x * x * functions.s
            let value = radialVelocity * x2c + (1 - alpha) * x3s + x
            let derivative = radialVelocity * x * (1 - alpha * x * x * functions.s) + (1 - alpha) * x2c + 1
            guard value.isFinite, derivative.isFinite else { return nil }
            return (sign * value - magnitude, derivative)
        }
        var low = 0.0, high = 1.0, expansions = 0
        while let upper = evaluate(high), upper.residual < 0 {
            guard expansions < options.maximumBracketExpansions else { throw PropagationError.bracketFailure }
            high *= 2; expansions += 1
            guard high.isFinite else { throw PropagationError.numericalRangeExceeded }
        }
        let tolerance = options.absoluteTimeToleranceSeconds / timeScale + options.relativeTimeTolerance * magnitude
        var x = min(magnitude, high / 2)
        var root: Double?
        var iterations = 0
        for iteration in 1...options.maximumIterations {
            iterations = iteration
            guard let value = evaluate(x) else {
                high = x; x = low + (high - low) / 2
                continue
            }
            guard value.derivative > 0 else { throw PropagationError.numericalRangeExceeded }
            if abs(value.residual) <= tolerance { root = sign * x; break }
            if value.residual < 0 { low = x } else { high = x }
            let midpoint = low + (high - low) / 2
            if midpoint == low || midpoint == high { root = sign * midpoint; break }
            let newton = x - value.residual / value.derivative
            let quarter = (high - low) / 4
            x = newton > low + quarter && newton < high - quarter ? newton : midpoint
        }
        guard let chi = root else { throw PropagationError.convergenceFailure(iterations: iterations) }
        guard let functions = Self.stumpff(alpha * chi * chi) else { throw PropagationError.numericalRangeExceeded }
        let x2c = chi * chi * functions.c, x3s = chi * chi * chi * functions.s
        let f = 1 - x2c, g = tau - x3s
        let position = f * radial + g * scaledVelocity
        let radius = Self.norm(position)
        guard radius.isFinite, radius > 0 else { throw PropagationError.numericalRangeExceeded }
        let fdot = (alpha * x3s - chi) / radius, gdot = 1 - x2c / radius
        let velocity = fdot * radial + gdot * scaledVelocity
        let state = try makeState(epoch: epoch, position: radiusScale * position, velocity: velocityScale * velocity)
        return TwoBodyPropagationResult(state: state, iterations: iterations, bracketExpansions: expansions)
    }

    private func makeState(epoch: Epoch, position: SIMD3<Double>, velocity: SIMD3<Double>) throws -> CartesianState {
        let r = Self.norm(position)
        let magnitude = (gravitationalParameter.cubicMetresPerSecondSquared / r) / r
        let acceleration = -magnitude * (position / r)
        guard [position, velocity, acceleration].allSatisfy({ [$0.x, $0.y, $0.z].allSatisfy(\.isFinite) }) else {
            throw PropagationError.numericalRangeExceeded
        }
        return try CartesianState(epoch: epoch, frame: frame, position: position, velocity: velocity, acceleration: acceleration)
    }

    // Series near zero; half-angle C avoids cancellation elsewhere. Large
    // hyperbolic arguments return nil so root bracketing can bisect inward.
    private static func stumpff(_ z: Double) -> (c: Double, s: Double)? {
        guard z.isFinite else { return nil }
        if abs(z) < 0.5 {
            var c = 0.5, s = 1.0 / 6, tc = c, ts = s
            for n in 1...16 {
                tc *= -z / Double((2*n+1)*(2*n+2))
                ts *= -z / Double((2*n+2)*(2*n+3))
                c += tc; s += ts
            }
            return (c, s)
        }
        let x = sqrt(abs(z))
        let c: Double, s: Double
        if z > 0 {
            c = 2 * pow(sin(x / 2), 2) / z
            s = (x - sin(x)) / (x * x * x)
        } else {
            c = 2 * pow(sinh(x / 2), 2) / (-z)
            s = (sinh(x) - x) / (x * x * x)
        }
        return c.isFinite && s.isFinite ? (c, s) : nil
    }
    private static func norm(_ x: SIMD3<Double>) -> Double { hypot(hypot(x.x, x.y), x.z) }
    private static func cross(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> SIMD3<Double> {
        SIMD3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x)
    }
}
