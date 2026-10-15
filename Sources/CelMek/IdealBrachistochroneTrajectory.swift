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

public enum BrachistochroneError: Error, Equatable {
    case invalidAcceleration, nonFiniteEndpoint, nonInertialFrame
    case numericalRangeExceeded, outsideCoverage
}

public enum BrachistochronePhase: Equatable, Sendable { case accelerating, braking, arrived }

/// Minimum-time, straight-line, rest-to-rest motion with bounded acceleration
/// magnitude in a Newtonian inertial frame. No gravity, mass or fuel model.
/// Immutable and independent of query order; coverage is departure through arrival.
public struct IdealBrachistochroneTrajectory: StateProvider {
    public let frame: Frame
    public let departureEpoch: Epoch
    public let switchEpoch: Epoch
    public let arrivalEpoch: Epoch
    public let departurePositionMetres: SIMD3<Double>
    public let arrivalPositionMetres: SIMD3<Double>
    public let maximumAccelerationMetresPerSecondSquared: Double
    public let distanceMetres: Double
    public let duration: Duration
    public let peakSpeedMetresPerSecond: Double
    /// Integral of acceleration magnitude, not a propellant estimate.
    public let totalDeltaVelocityMetresPerSecond: Double
    public let accelerationDirection: SIMD3<Double>
    public var brakingDirection: SIMD3<Double> { -accelerationDirection }
    public var coverage: ClosedRange<Epoch> { departureEpoch...arrivalEpoch }
    private let halfDuration: Duration

    public init(departureEpoch: Epoch, frame: Frame,
                departurePositionMetres start: SIMD3<Double>, arrivalPositionMetres end: SIMD3<Double>,
                maximumAccelerationMetresPerSecondSquared acceleration: Double) throws {
        guard acceleration.isFinite, acceleration > 0 else { throw BrachistochroneError.invalidAcceleration }
        guard frame.pseudoInertial else { throw BrachistochroneError.nonInertialFrame }
        guard [start.x,start.y,start.z,end.x,end.y,end.z].allSatisfy(\.isFinite) else {
            throw BrachistochroneError.nonFiniteEndpoint
        }
        let delta = end-start, distance = hypot(hypot(delta.x,delta.y),delta.z)
        guard distance.isFinite else { throw BrachistochroneError.numericalRangeExceeded }
        // Split square roots avoid overflow/underflow in distance / acceleration.
        let halfSeconds = distance == 0 ? 0 : sqrt(distance)/sqrt(acceleration)
        let peak = distance == 0 ? 0 : sqrt(distance)*sqrt(acceleration)
        guard halfSeconds.isFinite, peak.isFinite, (2*peak).isFinite,
              distance == 0 || (halfSeconds > 0 && peak > 0) else {
            throw BrachistochroneError.numericalRangeExceeded
        }
        let half = try Duration(seconds: halfSeconds), total = try half.adding(half)
        let change = try departureEpoch.advanced(by: half), arrival = try change.advanced(by: half)
        guard distance == 0 || (change > departureEpoch && arrival > change) else {
            throw BrachistochroneError.numericalRangeExceeded
        }
        self.frame = frame; self.departureEpoch = departureEpoch; switchEpoch = change; arrivalEpoch = arrival
        departurePositionMetres = start; arrivalPositionMetres = end
        maximumAccelerationMetresPerSecondSquared = acceleration; distanceMetres = distance
        halfDuration = half; duration = total; peakSpeedMetresPerSecond = peak
        totalDeltaVelocityMetresPerSecond = 2*peak
        accelerationDirection = distance == 0 ? .zero : delta/distance
    }

    /// Right-continuous at the switch (braking) and arrival (zero acceleration).
    public func phase(at epoch: Epoch) throws -> BrachistochronePhase {
        guard coverage.contains(epoch) else { throw BrachistochroneError.outsideCoverage }
        if epoch == arrivalEpoch { return .arrived }
        return epoch < switchEpoch ? .accelerating : .braking
    }

    public func state(at epoch: Epoch) throws -> CartesianState {
        let phase = try phase(at: epoch)
        if phase == .arrived {
            return try CartesianState(epoch: epoch,frame: frame,position: arrivalPositionMetres,
                velocity: .zero,acceleration: .zero)
        }
        let accelerating = phase == .accelerating
        // Evaluate braking backwards from the endpoint, preserving small remaining
        // distances and speeds rather than subtracting two large elapsed times.
        let interval = try accelerating ? epoch.duration(since: departureEpoch) : arrivalEpoch.duration(since: epoch)
        let fraction = interval.seconds/halfDuration.seconds
        let displacement = 0.5*fraction*fraction
        let endWeight = accelerating ? displacement : 1-displacement
        let startWeight = accelerating ? 1-displacement : displacement
        let position = departurePositionMetres*startWeight + arrivalPositionMetres*endWeight
        let velocity = accelerationDirection*(peakSpeedMetresPerSecond*fraction)
        let acceleration = accelerationDirection*(accelerating ? maximumAccelerationMetresPerSecondSquared : -maximumAccelerationMetresPerSecondSquared)
        return try CartesianState(epoch: epoch,frame: frame,position: position,velocity: velocity,acceleration: acceleration)
    }
}
