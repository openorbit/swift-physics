// SPDX-License-Identifier: Apache-2.0
import FluidDynamics
import Atmospheres
import Math

/// Adapts a pure aerodynamic evaluation into one coupled force/torque contribution.
/// Returned loads must use the attitude's body axes; reference points and the centre
/// of mass use the same body-coordinate origin. The evaluator is called at every
/// dynamics stage and must not advance subsystem state or a simulation clock.
public struct FluidWrench: SpacecraftWrenchModel {
    public let bodyFrame: Frame
    public let centreOfMassBodyMetres: SIMD3<Double>
    private let evaluator: (CoupledSpacecraftState, RotationalEvaluationContext) throws -> FluidLoads

    public init(bodyFrame: Frame, centreOfMassBodyMetres: SIMD3<Double> = .zero,
                evaluate: @escaping (CoupledSpacecraftState, RotationalEvaluationContext) throws -> FluidLoads) throws {
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
    func fluidFlow(atmosphere sample: AtmosphereSample, origin: BodyIdentifier,
                         properties: AtmosphericProperties) throws -> FluidFlow {
        guard sample.epoch == epoch else { throw AttitudeError.epochMismatch }
        guard sample.frame === spacecraft.frame else { throw AttitudeError.frameMismatch }
        guard sample.origin == origin else { throw EnvironmentalForceError.originMismatch }
        return try fluidFlow(mediumVelocityReferenceMetresPerSecond: sample.velocityMetresPerSecond,
                                   densityKilogramsPerCubicMetre: properties.densityKilogramsPerCubicMetre,
            pressurePascals: properties.pressurePascals,temperatureKelvin: properties.temperatureKelvin,
            speedOfSoundMetresPerSecond: properties.speedOfSoundMetresPerSecond,
            dynamicViscosityPascalSeconds: properties.dynamicViscosityPascalSeconds)
    }
    /// Atmosphere convenience; general fluid sampling uses the scalar-properties overload.
    func fluidFlow(mediumVelocityReferenceMetresPerSecond: SIMD3<Double>,
                   properties: AtmosphericProperties) throws -> FluidFlow {
        try fluidFlow(mediumVelocityReferenceMetresPerSecond: mediumVelocityReferenceMetresPerSecond,
            densityKilogramsPerCubicMetre: properties.densityKilogramsPerCubicMetre,
            pressurePascals: properties.pressurePascals,temperatureKelvin: properties.temperatureKelvin,
            speedOfSoundMetresPerSecond: properties.speedOfSoundMetresPerSecond,
            dynamicViscosityPascalSeconds: properties.dynamicViscosityPascalSeconds)
    }
    /// Medium and vehicle velocities share reference frame and epoch; rotate their
    /// difference to body axes. Liquid currents and atmospheric rotation are explicit.
    func fluidFlow(mediumVelocityReferenceMetresPerSecond: SIMD3<Double>,
                   densityKilogramsPerCubicMetre: Double, pressurePascals: Double? = nil,
                   temperatureKelvin: Double? = nil, speedOfSoundMetresPerSecond: Double? = nil,
                   dynamicViscosityPascalSeconds: Double? = nil) throws -> FluidFlow {
        try FluidFlow(relativeVelocityMetresPerSecond:
            attitude.bodyVector(fromReference: spacecraft.cartesian.velocity-mediumVelocityReferenceMetresPerSecond),
            densityKilogramsPerCubicMetre: densityKilogramsPerCubicMetre,
            pressurePascals: pressurePascals,temperatureKelvin: temperatureKelvin,
            speedOfSoundMetresPerSecond: speedOfSoundMetresPerSecond,
            dynamicViscosityPascalSeconds: dynamicViscosityPascalSeconds)
    }
}
