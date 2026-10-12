// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

/// A prescribed thrust cone in native route axes. The axis provider is pure,
/// evaluated at RK stages; this does not integrate attitude or certify slew rates.
public struct WaypointThrustCone {
    public let frame: Frame
    public let halfAngleRadians: Double
    private let axis: (Epoch) throws -> SIMD3<Double>
    public init(frame: Frame, halfAngleRadians: Double,
                axis: @escaping (Epoch) throws -> SIMD3<Double>) throws {
        guard halfAngleRadians.isFinite, halfAngleRadians >= 0, halfAngleRadians <= .pi else {
            throw WaypointRoutePropagationError.invalidLimits
        }
        self.frame = frame; self.halfAngleRadians = halfAngleRadians; self.axis = axis
    }
    func check(_ acceleration: SIMD3<Double>, at epoch: Epoch, leg: Int) throws {
        let length = routeNorm(acceleration)
        guard length > 0 else { return } // No attitude requirement while coasting.
        let direction = try axis(epoch)
        let scale = max(abs(direction.x), max(abs(direction.y), abs(direction.z)))
        guard [direction.x, direction.y, direction.z].allSatisfy(\.isFinite), scale.isFinite, scale > 0 else { throw WaypointRoutePropagationError.invalidLimits }
        let scaled = direction/scale
        let unit = scaled/routeNorm(scaled)
        let cosine = (unit * (acceleration/length)).sum()
        guard cosine + 64*Double.ulpOfOne >= cos(halfAngleRadians) else {
            throw WaypointRoutePropagationError.attitudeLimitExceeded(leg: leg)
        }
    }
}

/// Conservative stored-heat envelope: dH/dt = gain*thrust - coolingRate*H.
/// Heat is carried across spatial/temporal boundaries using each interval's
/// engine ceiling, rather than inferring a thermal model from acceleration limits.
public struct WaypointThermalEnvelope: Equatable, Sendable {
    public let initialHeatJoules, maximumHeatJoules: Double
    public let heatingWattsPerNewton, coolingRatePerSecond: Double
    public init(initialHeatJoules: Double = 0, maximumHeatJoules: Double,
                heatingWattsPerNewton: Double, coolingRatePerSecond: Double) throws {
        guard [initialHeatJoules, maximumHeatJoules, heatingWattsPerNewton, coolingRatePerSecond]
            .allSatisfy(\.isFinite),
            initialHeatJoules >= 0, maximumHeatJoules > 0, initialHeatJoules <= maximumHeatJoules,
            heatingWattsPerNewton >= 0, coolingRatePerSecond >= 0 else {
            throw WaypointRoutePropagationError.invalidLimits
        }
        self.initialHeatJoules = initialHeatJoules; self.maximumHeatJoules = maximumHeatJoules
        self.heatingWattsPerNewton = heatingWattsPerNewton; self.coolingRatePerSecond = coolingRatePerSecond
    }
    func advance(heat: Double, duration: Double, thrustCeiling: Double, leg: Int) throws -> Double {
        let power = heatingWattsPerNewton * thrustCeiling
        let x = coolingRatePerSecond * duration
        let integral: Double
        if x == 0 { integral = duration }
        else if x.isInfinite { integral = 1/coolingRatePerSecond }
        else { integral = duration * (-expm1(-x)/x) }
        let added = power * integral
        var end = heat * exp(-x) + added
        // An upper envelope, including a small allowance for floating arithmetic.
        if duration > 0 && (power > 0 || coolingRatePerSecond > 0) {
            end += 64*Double.ulpOfOne*max(heat, added)
        }
        guard end.isFinite else { throw SpacecraftDynamicsError.numericalRangeExceeded }
        guard max(heat, end) <= maximumHeatJoules else {
            throw WaypointRoutePropagationError.thermalBoundExceeded(leg: leg)
        }
        return end
    }
}

public struct WaypointRoutePhysicalConstraints {
    public let thrustCone: WaypointThrustCone?
    public let thermalEnvelope: WaypointThermalEnvelope?
    public static let unconstrained = WaypointRoutePhysicalConstraints()
    public init(thrustCone: WaypointThrustCone? = nil, thermalEnvelope: WaypointThermalEnvelope? = nil) {
        self.thrustCone = thrustCone; self.thermalEnvelope = thermalEnvelope
    }
    func validate(frame: Frame) throws {
        guard thrustCone == nil || thrustCone!.frame === frame else { throw SpacecraftDynamicsError.frameMismatch }
    }
}
