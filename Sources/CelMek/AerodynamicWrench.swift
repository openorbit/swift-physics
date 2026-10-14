// SPDX-License-Identifier: Apache-2.0
import Aerodynamics
import Atmospheres
import Math

/// Adapts a pure aerodynamic evaluation into one coupled force/torque contribution.
/// Returned loads must use the attitude's body axes; reference points and the centre
/// of mass use the same body-coordinate origin. The evaluator is called at every
/// dynamics stage and must not advance subsystem state or a simulation clock.
public struct AerodynamicWrench: SpacecraftWrenchModel {
    public let bodyFrame: Frame
    public let centreOfMassBodyMetres: SIMD3<Double>
    private let evaluator: (CoupledSpacecraftState, RotationalEvaluationContext) throws -> AerodynamicLoads

    public init(bodyFrame: Frame, centreOfMassBodyMetres: SIMD3<Double> = .zero,
                evaluate: @escaping (CoupledSpacecraftState, RotationalEvaluationContext) throws -> AerodynamicLoads) throws {
        guard attitudeFinite(centreOfMassBodyMetres) else { throw AttitudeError.invalidVector }
        self.bodyFrame = bodyFrame
        self.centreOfMassBodyMetres = centreOfMassBodyMetres
        evaluator = evaluate
    }

    public func contribution(for state: CoupledSpacecraftState,
                             in context: RotationalEvaluationContext) throws -> SpacecraftWrenchContribution {
        try validateRotationalState(state.attitude, spacecraft: state.spacecraft, context: context)
        guard bodyFrame === context.bodyFrame else { throw AttitudeError.frameMismatch }
        let loads = try evaluator(state, context).referred(to: centreOfMassBodyMetres)
        return try SpacecraftWrenchContribution(epoch: state.epoch,
            referenceFrame: context.referenceFrame, bodyFrame: bodyFrame,
            forceReferenceNewtons: state.attitude.referenceVector(fromBody: loads.forceNewtons),
            torqueBodyNewtonMetres: loads.momentNewtonMetres)
    }
}

public extension CoupledSpacecraftState {
    /// Checked atmosphere sample variant, requiring a matching epoch, reference
    /// frame and explicitly supplied atmosphere origin.
    func aerodynamicFlow(atmosphere sample: AtmosphereSample, origin: BodyIdentifier,
                         properties: AtmosphericProperties) throws -> AerodynamicFlow {
        guard sample.epoch == epoch else { throw AttitudeError.epochMismatch }
        guard sample.frame === spacecraft.frame else { throw AttitudeError.frameMismatch }
        guard sample.origin == origin else { throw EnvironmentalForceError.originMismatch }
        return try aerodynamicFlow(airVelocityReferenceMetresPerSecond: sample.velocityMetresPerSecond,
                                   properties: properties)
    }
    /// Construct centre-of-mass flow from vehicle and air velocities in the same
    /// inertial reference frame. Air velocity must already include co-rotation;
    /// thermodynamic properties are sampled separately by the caller.
    func aerodynamicFlow(airVelocityReferenceMetresPerSecond: SIMD3<Double>,
                         properties: AtmosphericProperties) throws -> AerodynamicFlow {
        try AerodynamicFlow(relativeVelocityMetresPerSecond:
            attitude.bodyVector(fromReference: spacecraft.cartesian.velocity-airVelocityReferenceMetresPerSecond),
            properties: properties)
    }
}
