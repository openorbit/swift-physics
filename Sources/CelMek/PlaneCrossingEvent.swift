// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Mattias Holm
import Foundation
import Math

/// Signed distance to an origin-centred fixed plane in an explicit inertial frame.
/// Use EventDetectedPropagation for bounded window searches over moving sites.
public struct PlaneCrossingEvent: PropagationEventDetector {
    public let identifier: String
    public let frame: Frame
    public let normal: SIMD3<Double>
    public let direction: EventCrossingDirection
    public let action: PropagationEventAction
    public init(identifier: String, frame: Frame, normal: SIMD3<Double>, direction: EventCrossingDirection = .any,
                action: PropagationEventAction = .record) throws {
        let length = hypot(hypot(normal.x,normal.y),normal.z)
        guard !identifier.isEmpty, length.isFinite, length > 0 else { throw PropagationEventError.invalidDetector }
        self.identifier = identifier; self.frame = frame; self.normal = normal/length
        self.direction = direction; self.action = action
    }
    public func value(at state: SpacecraftState) throws -> Double {
        guard state.frame === frame else { throw SpacecraftDynamicsError.frameMismatch }
        return (state.cartesian.position*normal).sum()
    }
}
