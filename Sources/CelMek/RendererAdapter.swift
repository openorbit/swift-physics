// SPDX-License-Identifier: Apache-2.0

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
