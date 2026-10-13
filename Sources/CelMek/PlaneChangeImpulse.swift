// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum PlaneChangeError: Error, Equatable {
    case invalidSpeed, invalidAngle, invalidDirections, numericalRangeExceeded
}
public enum PlaneChangeSense: Equatable, Sendable { case positive, negative }

/// Simultaneous tangential speed and orbital-plane change at a common node.
public struct PlaneChangeImpulse {
    public let beforeSpeedMetresPerSecond, afterSpeedMetresPerSecond, angleRadians: Double
    public let deltaVelocityMetresPerSecond: Double
    public init(beforeSpeedMetresPerSecond: Double, afterSpeedMetresPerSecond: Double, angleRadians: Double) throws {
        guard [beforeSpeedMetresPerSecond, afterSpeedMetresPerSecond].allSatisfy({ $0.isFinite && $0 >= 0 }) else {
            throw PlaneChangeError.invalidSpeed
        }
        guard angleRadians.isFinite, angleRadians >= 0, angleRadians <= .pi else { throw PlaneChangeError.invalidAngle }
        let lateral = (sqrt(beforeSpeedMetresPerSecond)*sin(angleRadians/2))*(2*sqrt(afterSpeedMetresPerSecond))
        let delta = hypot(afterSpeedMetresPerSecond-beforeSpeedMetresPerSecond, lateral)
        guard delta.isFinite else { throw PlaneChangeError.numericalRangeExceeded }
        self.beforeSpeedMetresPerSecond = beforeSpeedMetresPerSecond
        self.afterSpeedMetresPerSecond = afterSpeedMetresPerSecond
        self.angleRadians = angleRadians; deltaVelocityMetresPerSecond = delta
    }
    /// Positive sense rotates velocity by the right-hand rule about the supplied
    /// outward radius. Both velocities are tangential; radial changes are excluded.
    public func maneuver(epoch: Epoch, frame: Frame, radialDirection: SIMD3<Double>,
                         progradeDirection: SIMD3<Double>, sense: PlaneChangeSense = .positive) throws -> ImpulsiveManeuver {
        guard frame.pseudoInertial else { throw PropagationError.nonInertialFrame }
        let radial = try planeUnit(radialDirection), supplied = try planeUnit(progradeDirection)
        let dot = (radial*supplied).sum()
        guard abs(dot) <= 1e-12 else { throw PlaneChangeError.invalidDirections }
        let tangent = try planeUnit(supplied-radial*dot)
        let rotation = Quaternion(angle: sense == .positive ? angleRadians : -angleRadians, axis: radial)
        let final = rotation.matrix*tangent
        return try ImpulsiveManeuver(epoch: epoch, frame: frame,
            deltaVelocityMetresPerSecond: final*afterSpeedMetresPerSecond-tangent*beforeSpeedMetresPerSecond)
    }
}
private func planeUnit(_ v: SIMD3<Double>) throws -> SIMD3<Double> {
    let n = hypot(hypot(v.x,v.y),v.z)
    guard n.isFinite, n > 0 else { throw PlaneChangeError.invalidDirections }
    return v/n
}
