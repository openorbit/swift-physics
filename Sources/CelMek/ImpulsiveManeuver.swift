// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum ManeuverError: Error, Equatable {
    case nonFiniteDeltaVelocity, beforeReferenceEpoch, duplicateEpoch, frameTreeRequired
}

/// An instantaneous free-vector velocity change. Units are metres/second in
/// the explicitly supplied axes, evaluated at the burn epoch.
public struct ImpulsiveManeuver {
    public let epoch: Epoch
    public let frame: Frame
    public let deltaVelocityMetresPerSecond: SIMD3<Double>
    public init(epoch: Epoch, frame: Frame, deltaVelocityMetresPerSecond: SIMD3<Double>) throws {
        guard [deltaVelocityMetresPerSecond.x, deltaVelocityMetresPerSecond.y,
               deltaVelocityMetresPerSecond.z].allSatisfy(\.isFinite) else {
            throw ManeuverError.nonFiniteDeltaVelocity
        }
        self.epoch = epoch; self.frame = frame
        self.deltaVelocityMetresPerSecond = deltaVelocityMetresPerSecond
    }
}

public enum ManeuverEpochSide: Sendable { case before, after }

public struct ResolvedImpulsiveManeuver {
    public let maneuver: ImpulsiveManeuver
    /// Delta-v rotated into the propagation axes, without translation or frame-rate terms.
    public let deltaVelocityInPropagationFrame: SIMD3<Double>
    public let preBurnState: CartesianState
    public let postBurnState: CartesianState
}

/// Immutable piecewise two-body trajectory. Burns are resolved once, in time
/// order during construction. Query order cannot change any segment.
public struct ManeuveredTwoBodyPropagator: Propagator {
    public let base: TwoBodyPropagator
    public let resolvedManeuvers: [ResolvedImpulsiveManeuver]
    public var referenceEpoch: Epoch { base.referenceEpoch }
    public var frame: Frame { base.frame }
    public var centralBody: BodyIdentifier { base.centralBody }
    private let segments: [TwoBodyPropagator]

    public init(base: TwoBodyPropagator, maneuvers: [ImpulsiveManeuver], using tree: FrameTree? = nil) throws {
        let ordered = maneuvers.sorted { $0.epoch < $1.epoch }
        var records: [ResolvedImpulsiveManeuver] = [], segments = [base]
        var previousEpoch: Epoch?
        for maneuver in ordered {
            guard maneuver.epoch >= base.referenceEpoch else { throw ManeuverError.beforeReferenceEpoch }
            guard maneuver.epoch != previousEpoch else { throw ManeuverError.duplicateEpoch }
            let before = try segments.last!.state(at: maneuver.epoch)
            let deltaVelocity: SIMD3<Double>
            if maneuver.frame === base.frame { deltaVelocity = maneuver.deltaVelocityMetresPerSecond }
            else {
                guard let tree else { throw ManeuverError.frameTreeRequired }
                let transform = try tree.transform(from: maneuver.frame, to: base.frame, at: maneuver.epoch)
                deltaVelocity = transform.rotation.matrix * maneuver.deltaVelocityMetresPerSecond
            }
            let seed = try CartesianState(epoch: maneuver.epoch, frame: base.frame, position: before.position,
                                          velocity: before.velocity + deltaVelocity)
            let segment = try TwoBodyPropagator(initialState: seed, centralBody: base.centralBody,
                gravitationalParameter: base.gravitationalParameter, options: base.options)
            let after = try segment.state(at: maneuver.epoch)
            records.append(ResolvedImpulsiveManeuver(maneuver: maneuver,
                deltaVelocityInPropagationFrame: deltaVelocity, preBurnState: before, postBurnState: after))
            segments.append(segment); previousEpoch = maneuver.epoch
        }
        self.base = base; self.resolvedManeuvers = records; self.segments = segments
    }

    /// Right-continuous: queries at a burn epoch return the post-burn state.
    public func state(at epoch: Epoch) throws -> CartesianState { try state(at: epoch, side: .after) }

    public func state(at epoch: Epoch, side: ManeuverEpochSide) throws -> CartesianState {
        // Upper/lower-bound search chooses a fixed segment in O(log burns).
        var low = 0, high = resolvedManeuvers.count
        while low < high {
            let middle = low + (high - low) / 2
            let burnEpoch = resolvedManeuvers[middle].maneuver.epoch
            let included = burnEpoch < epoch || (burnEpoch == epoch && side == .after)
            if included { low = middle + 1 } else { high = middle }
        }
        return try segments[low].state(at: epoch)
    }
}
