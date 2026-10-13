// SPDX-License-Identifier: Apache-2.0
import Testing
import CelMek
import Math

@Suite struct RendererAdapterTests {
    let frame = Frame(name: "Rendering", pseudoInertial: true)
    func state(_ position: SIMD3<Double>, velocity: SIMD3<Double> = .zero,
               acceleration: SIMD3<Double>? = nil) throws -> CartesianState {
        try CartesianState(epoch: .j2000, frame: frame, position: position,
                           velocity: velocity, acceleration: acceleration)
    }

    @Test func preservesNearbyDetailAndMovingOriginDerivatives() throws {
        let base = SIMD3<Double>(1e12, -2e12, 3e12)
        let origin = try state(base, velocity: SIMD3(10, 20, 30), acceleration: SIMD3(1, 2, 3))
        let target = try state(base + SIMD3(1, -2, 3), velocity: SIMD3(14, 15, 36), acceleration: SIMD3(3, 1, 7))
        let render = try RendererAdapter(origin: origin, metresPerRenderUnit: 2).adapt(target)
        #expect(render.position == SIMD3<Float>(0.5, -1, 1.5))
        #expect(render.velocity == SIMD3<Float>(2, -2.5, 3))
        #expect(render.acceleration == SIMD3<Float>(1, -0.5, 2))
        #expect(Float(base.x + 1) == Float(base.x)) // Premature conversion loses the metre.
        #expect(try RendererAdapter(origin: origin).adapt(origin).position == .zero)
        #expect(target.position == base + SIMD3(1, -2, 3))
    }

    @Test func missingAccelerationAndSnapshotContracts() throws {
        let origin = try state(.zero)
        let adapter = try RendererAdapter(origin: origin)
        #expect(try adapter.adapt(state(.zero, acceleration: .zero)).acceleration == nil)
        let acceleratingOrigin = try RendererAdapter(origin: state(.zero, acceleration: .zero))
        #expect(try acceleratingOrigin.adapt(origin).acceleration == nil)
        let later = try CartesianState(epoch: Epoch.j2000.advanced(by: CelMek.Duration(seconds: 1)),
                                       frame: frame, position: .zero, velocity: .zero)
        #expect(throws: RendererAdapterError.epochMismatch) { try adapter.adapt(later) }
        let other = try CartesianState(epoch: .j2000, frame: Frame(name: frame.name, pseudoInertial: true),
                                       position: .zero, velocity: .zero)
        #expect(throws: RendererAdapterError.frameMismatch) { try adapter.adapt(other) }
    }

    @Test func rejectsInvalidScaleAndOverflow() throws {
        for scale in [0, -1, Double.nan, .infinity] {
            #expect(throws: RendererAdapterError.invalidScale) {
                try RendererAdapter(origin: state(.zero), metresPerRenderUnit: scale)
            }
        }
        let adapter = try RendererAdapter(origin: state(.zero))
        #expect(throws: RendererAdapterError.coordinateOverflow) {
            try adapter.adapt(state(SIMD3(repeating: Double.greatestFiniteMagnitude)))
        }
        let remote = try state(SIMD3(repeating: -Double.greatestFiniteMagnitude))
        #expect(throws: RendererAdapterError.coordinateOverflow) {
            try RendererAdapter(origin: remote).adapt(state(SIMD3(repeating: Double.greatestFiniteMagnitude)))
        }
        #expect(throws: RendererAdapterError.coordinateOverflow) {
            try RendererAdapter(origin: state(.zero), metresPerRenderUnit: Double.leastNonzeroMagnitude)
                .adapt(state(SIMD3(1, 0, 0)))
        }
    }
}
