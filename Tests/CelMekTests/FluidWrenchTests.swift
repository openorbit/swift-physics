// SPDX-License-Identifier: Apache-2.0
import FluidDynamics
import Atmospheres
import Math
import Testing
@testable import CelMek

@Suite struct FluidWrenchTests {
    @Test func rotatesForceTransfersMomentAndSubtractsAirVelocity() throws {
        let frame = Frame(name: "Aero inertial", pseudoInertial: true)
        let body = Frame(parent: frame, name: "Aero body", pseudoInertial: false)
        let spacecraft = try SpacecraftState(cartesian: CartesianState(epoch: .j2000, frame: frame,
            position: SIMD3(6371000, 0, 0), velocity: SIMD3(0, 737.1, 0)),
            mass: SpacecraftMassProperties(dryMassKilograms: 100))
        let attitude = try AttitudeState(epoch: .j2000, referenceFrame: frame, bodyFrame: body,
            bodyToReference: Quaternion(angle: .pi/2, axis: SIMD3(0, 0, 1)))
        let state = try CoupledSpacecraftState(spacecraft: spacecraft, attitude: attitude)
        let context = try RotationalEvaluationContext(inertia: RigidBodyInertia(
            principalMomentsKilogramSquareMetres: SIMD3(2, 3, 4)), bodyFrame: body,
            referenceContext: ForceEvaluationContext(frame: frame, origin: .earth))
        let ambient = try ISAAtmosphere().properties(at: 0)
        let flow = try state.fluidFlow(mediumVelocityReferenceMetresPerSecond: SIMD3(0, 637.1, 0),
                                            properties: ambient)
        #expect(abs(flow.relativeVelocityMetresPerSecond.x-100) < 1e-10)
        #expect(abs(flow.relativeVelocityMetresPerSecond.y) < 1e-10)
        let sample = try AtmosphereSample(epoch: .j2000, frame: frame, origin: .earth,
            densityKilogramsPerCubicMetre: ambient.densityKilogramsPerCubicMetre,
            velocityMetresPerSecond: SIMD3(0, 637.1, 0))
        let checked = try state.fluidFlow(atmosphere: sample, origin: .earth, properties: ambient)
        #expect(checked.relativeVelocityMetresPerSecond == flow.relativeVelocityMetresPerSecond)
        #expect(throws: EnvironmentalForceError.originMismatch) {
            try state.fluidFlow(atmosphere: sample, origin: .mars, properties: ambient)
        }
        let adapter = try FluidWrench(bodyFrame: body, centreOfMassBodyMetres: SIMD3(1, 0, 0)) { _, _ in
            try FluidLoads(forceNewtons: SIMD3(-10, 0, -20),
                momentNewtonMetres: SIMD3(0, 3, 0), referencePointMetres: SIMD3(2, 0, 0))
        }
        let wrench = try adapter.contribution(for: state, in: context)
        #expect(abs(wrench.forceReferenceNewtons.x) < 1e-12)
        #expect(abs(wrench.forceReferenceNewtons.y+10) < 1e-12)
        #expect(wrench.forceReferenceNewtons.z == -20)
        #expect(wrench.torqueBodyNewtonMetres == SIMD3(0, 23, 0))
        #expect(wrench.propellantMassRateKilogramsPerSecond == 0)
        let derivative = try CoupledSpacecraftDynamics(context: context, wrenchModels: [adapter]).derivative(at: state)
        #expect(abs(derivative.translation.velocityRateMetresPerSecondSquared.y+0.1) < 1e-12)
        #expect(derivative.translation.velocityRateMetresPerSecondSquared.z == -0.2)
        #expect(abs(derivative.rotation.angularAccelerationBodyRadiansPerSecondSquared.y-23.0/3) < 1e-12)
        let wrongBody = Frame(name: "Other body", pseudoInertial: false)
        let wrong = try FluidWrench(bodyFrame: wrongBody) { _, _ in
            Issue.record("Frame mismatch must fail before evaluation")
            return try FluidLoads(forceNewtons: .zero, momentNewtonMetres: .zero)
        }
        #expect(throws: AttitudeError.frameMismatch) { try wrong.contribution(for: state, in: context) }
    }
}
