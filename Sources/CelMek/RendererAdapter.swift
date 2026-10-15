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

public enum RendererAdapterError: Error, Equatable, Sendable {
    case invalidScale, epochMismatch, frameMismatch, coordinateOverflow
}

/// Render-only coordinates in the source frame's right-handed axes.
/// Derivatives are per simulation second; unavailable acceleration stays nil.
public struct RenderState: Sendable, Equatable {
    public let position: SIMD3<Float>
    public let velocity: SIMD3<Float>
    public let acceleration: SIMD3<Float>?
}

/// A snapshot of a moving render origin at one epoch. Convert states into the
/// origin's frame with FrameTree before adapting them. This translates axes;
/// it does not rotate them or define a new physical reference frame.
public struct RendererAdapter {
    public let origin: CartesianState
    public let metresPerRenderUnit: Double

    public init(origin: CartesianState, metresPerRenderUnit: Double = 1) throws {
        guard metresPerRenderUnit.isFinite, metresPerRenderUnit > 0 else {
            throw RendererAdapterError.invalidScale
        }
        self.origin = origin
        self.metresPerRenderUnit = metresPerRenderUnit
    }

    /// Subtracts the origin in Double before scaling and converting to Float.
    /// Rebuild the adapter when the origin epoch changes. Float rounding is
    /// intentional; nonfinite intermediate or output coordinates are rejected.
    public func adapt(_ state: CartesianState) throws -> RenderState {
        guard state.epoch == origin.epoch else { throw RendererAdapterError.epochMismatch }
        guard state.frame === origin.frame else { throw RendererAdapterError.frameMismatch }
        let acceleration: SIMD3<Float>?
        if let value = state.acceleration, let reference = origin.acceleration {
            acceleration = try convert(value - reference)
        } else { acceleration = nil }
        return try RenderState(position: convert(state.position - origin.position),
                               velocity: convert(state.velocity - origin.velocity),
                               acceleration: acceleration)
    }

    private func convert(_ relative: SIMD3<Double>) throws -> SIMD3<Float> {
        let scaled = relative / metresPerRenderUnit
        let result = SIMD3<Float>(Float(scaled.x), Float(scaled.y), Float(scaled.z))
        guard scaled.x.isFinite, scaled.y.isFinite, scaled.z.isFinite,
              result.x.isFinite, result.y.isFinite, result.z.isFinite else {
            throw RendererAdapterError.coordinateOverflow
        }
        return result
    }
}
