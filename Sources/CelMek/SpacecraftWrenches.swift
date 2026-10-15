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

public struct CoupledSpacecraftState {
    public let spacecraft: SpacecraftState
    public let attitude: AttitudeState
    public var epoch: Epoch { spacecraft.epoch }
    public init(spacecraft: SpacecraftState, attitude: AttitudeState) throws {
        guard spacecraft.epoch == attitude.epoch else { throw AttitudeError.epochMismatch }
        guard spacecraft.frame === attitude.referenceFrame else { throw AttitudeError.frameMismatch }
        self.spacecraft = spacecraft; self.attitude = attitude
    }
}
public struct SpacecraftWrenchContribution {
    public let epoch: Epoch
    public let referenceFrame: Frame
    public let bodyFrame: Frame
    public let forceReferenceNewtons: SIMD3<Double>
    public let torqueBodyNewtonMetres: SIMD3<Double>
    public let propellantMassRateKilogramsPerSecond: Double
    public init(epoch: Epoch, referenceFrame: Frame, bodyFrame: Frame, forceReferenceNewtons: SIMD3<Double>,
        torqueBodyNewtonMetres: SIMD3<Double>, propellantMassRateKilogramsPerSecond: Double = 0) throws {
        guard attitudeFinite(forceReferenceNewtons), attitudeFinite(torqueBodyNewtonMetres), propellantMassRateKilogramsPerSecond.isFinite else { throw AttitudeError.invalidVector }
        self.epoch = epoch; self.referenceFrame = referenceFrame; self.bodyFrame = bodyFrame
        self.forceReferenceNewtons = forceReferenceNewtons; self.torqueBodyNewtonMetres = torqueBodyNewtonMetres
        self.propellantMassRateKilogramsPerSecond = propellantMassRateKilogramsPerSecond
    }
}
/// Force and torque from one physical model, counted once in coupled dynamics.
public protocol SpacecraftWrenchModel {
    func contribution(for state: CoupledSpacecraftState, in context: RotationalEvaluationContext) throws -> SpacecraftWrenchContribution
}
/// Continuous body-mounted engine. Position is measured from the fixed centre of
/// mass; nozzle direction is force on the spacecraft, not the exhaust direction.
public struct BodyMountedThruster: SpacecraftWrenchModel {
    public let bodyFrame: Frame
    public let thrustNewtons: Double
    public let directionBody: SIMD3<Double>
    public let positionBodyMetres: SIMD3<Double>
    public let effectiveExhaustVelocityMetresPerSecond: Double
    public let propellantConsumptionKilogramsPerSecond: Double
    public init(bodyFrame: Frame, thrustNewtons: Double, directionBody: SIMD3<Double>, positionBodyMetres: SIMD3<Double> = .zero,
        effectiveExhaustVelocityMetresPerSecond: Double) throws {
        guard thrustNewtons.isFinite, thrustNewtons >= 0 else { throw FiniteThrustError.invalidThrust }
        guard effectiveExhaustVelocityMetresPerSecond.isFinite, effectiveExhaustVelocityMetresPerSecond > 0 else { throw FiniteThrustError.invalidExhaustVelocity }
        guard attitudeFinite(positionBodyMetres) else { throw AttitudeError.invalidVector }
        let flow = thrustNewtons/effectiveExhaustVelocityMetresPerSecond
        guard flow.isFinite, thrustNewtons == 0 || flow > 0 else { throw FiniteThrustError.invalidExhaustVelocity }
        self.bodyFrame = bodyFrame; self.thrustNewtons = thrustNewtons; self.directionBody = try attitudeUnit(directionBody)
        self.positionBodyMetres = positionBodyMetres; self.effectiveExhaustVelocityMetresPerSecond = effectiveExhaustVelocityMetresPerSecond
        propellantConsumptionKilogramsPerSecond = flow
    }
    public func contribution(for state: CoupledSpacecraftState, in context: RotationalEvaluationContext) throws -> SpacecraftWrenchContribution {
        guard bodyFrame === state.attitude.bodyFrame, bodyFrame === context.bodyFrame,
              state.spacecraft.frame === context.referenceFrame else { throw AttitudeError.frameMismatch }
        guard thrustNewtons == 0 || state.spacecraft.mass.propellantMassKilograms > 0 else { throw NumericalPropagationError.propellantDepleted }
        let force = directionBody*thrustNewtons
        return try SpacecraftWrenchContribution(epoch: state.epoch,referenceFrame: context.referenceFrame,bodyFrame: bodyFrame,
            forceReferenceNewtons: state.attitude.referenceVector(fromBody: force),torqueBodyNewtonMetres: attitudeCross(positionBodyMetres,force),
            propellantMassRateKilogramsPerSecond: -propellantConsumptionKilogramsPerSecond)
    }
    /// Translation driven by a prescribed/replayed attitude. Lever-arm torque is
    /// available through coupled dynamics; this adapter integrates translation only.
    public func forceModel(attitude: any AttitudeProvider) throws -> BodyMountedThrustForceModel {
        try BodyMountedThrustForceModel(thruster: self,attitude: attitude)
    }
}
public struct BodyMountedThrustForceModel: FuelConsumingForceModel {
    public let thruster: BodyMountedThruster
    public let attitude: any AttitudeProvider
    public init(thruster: BodyMountedThruster, attitude: any AttitudeProvider) throws {
        guard thruster.bodyFrame === attitude.bodyFrame, attitude.referenceFrame.pseudoInertial else { throw AttitudeError.frameMismatch }
        self.thruster = thruster; self.attitude = attitude
    }
    public func contribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
        guard thruster.thrustNewtons == 0 || state.mass.propellantMassKilograms > 0 else { throw NumericalPropagationError.propellantDepleted }
        return try depletionBoundaryContribution(for: state,in: context)
    }
    public func depletionBoundaryContribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
        guard state.frame === attitude.referenceFrame, context.frame === attitude.referenceFrame else { throw AttitudeError.frameMismatch }
        let axes = try attitude.checkedAttitude(at: state.epoch)
        return try ForceContribution(epoch: state.epoch,frame: state.frame,
            accelerationMetresPerSecondSquared: axes.referenceVector(fromBody: thruster.directionBody)*(thruster.thrustNewtons/state.mass.totalMassKilograms),
            propellantMassRateKilogramsPerSecond: -thruster.propellantConsumptionKilogramsPerSecond)
    }
}
public struct CoupledSpacecraftDerivative {
    public let translation: SpacecraftDerivative
    public let rotation: RigidBodyDerivative
}
/// Fixed mass geometry/inertia; propellant depletion changes translation mass,
/// but never silently changes the tensor or centre of mass.
public struct CoupledSpacecraftDynamics {
    public let rotational: RigidBodyDynamics
    public let translational: SpacecraftDynamics
    public let wrenchModels: [any SpacecraftWrenchModel]
    public var context: RotationalEvaluationContext { rotational.context }
    public init(context: RotationalEvaluationContext, forceModels: [any ForceModel] = [], torqueModels: [any TorqueModel] = [],
        wrenchModels: [any SpacecraftWrenchModel] = []) throws {
        rotational = RigidBodyDynamics(context: context,models: torqueModels)
        translational = try SpacecraftDynamics(models: forceModels,context: context.referenceContext)
        self.wrenchModels = wrenchModels
    }
    public func derivative(at state: CoupledSpacecraftState) throws -> CoupledSpacecraftDerivative {
        guard state.attitude.bodyFrame === context.bodyFrame, state.attitude.referenceFrame === context.referenceFrame else { throw AttitudeError.frameMismatch }
        let base = try translational.derivative(at: state.spacecraft)
        var acceleration = base.velocityRateMetresPerSecondSquared, flow = base.propellantMassRateKilogramsPerSecond, torque = SIMD3<Double>.zero
        for model in wrenchModels {
            let wrench = try model.contribution(for: state,in: context)
            guard wrench.epoch == state.epoch else { throw AttitudeError.epochMismatch }
            guard wrench.referenceFrame === context.referenceFrame, wrench.bodyFrame === context.bodyFrame else { throw AttitudeError.frameMismatch }
            acceleration += wrench.forceReferenceNewtons/state.spacecraft.mass.totalMassKilograms
            flow += wrench.propellantMassRateKilogramsPerSecond; torque += wrench.torqueBodyNewtonMetres
        }
        guard attitudeFinite(acceleration), flow.isFinite, attitudeFinite(torque) else { throw AttitudeError.numericalRangeExceeded }
        let rotation = try rotational.derivative(at: state.attitude,spacecraft: state.spacecraft,additionalTorqueBodyNewtonMetres: torque)
        return CoupledSpacecraftDerivative(translation: SpacecraftDerivative(epoch: state.epoch,frame: context.referenceFrame,
            positionRateMetresPerSecond: state.spacecraft.cartesian.velocity,velocityRateMetresPerSecondSquared: acceleration,
            propellantMassRateKilogramsPerSecond: flow),rotation: rotation)
    }
}
public struct MagneticFieldSample {
    public let epoch: Epoch
    public let frame: Frame
    public let origin: BodyIdentifier
    public let fieldTeslas: SIMD3<Double>
    public init(epoch: Epoch, frame: Frame, origin: BodyIdentifier, fieldTeslas: SIMD3<Double>) throws {
        try origin.validate()
        guard attitudeFinite(fieldTeslas) else { throw AttitudeError.invalidVector }
        self.epoch = epoch; self.frame = frame; self.origin = origin; self.fieldTeslas = fieldTeslas
    }
}
public protocol MagneticFieldProvider {
    func field(at state: SpacecraftState) throws -> MagneticFieldSample
}
public struct UniformMagneticField: MagneticFieldProvider {
    public let frame: Frame
    public let origin: BodyIdentifier
    public let fieldTeslas: SIMD3<Double>
    public init(frame: Frame, origin: BodyIdentifier, fieldTeslas: SIMD3<Double>) throws {
        let sample = try MagneticFieldSample(epoch: .j2000,frame: frame,origin: origin,fieldTeslas: fieldTeslas)
        self.frame = frame; self.origin = origin; self.fieldTeslas = sample.fieldTeslas
    }
    public func field(at state: SpacecraftState) throws -> MagneticFieldSample {
        guard state.frame === frame else { throw AttitudeError.frameMismatch }
        return try MagneticFieldSample(epoch: state.epoch,frame: frame,origin: origin,fieldTeslas: fieldTeslas)
    }
}
/// Bounded commanded dipole (A m²), τBody = dipoleBody × fieldBody. Also
/// represents a fixed residual dipole with a constant command. No field model,
/// electronics, hysteresis or power consumption is inferred.
public struct MagnetorquerTorque: TorqueModel {
    public let fieldProvider: any MagneticFieldProvider
    public let maximumDipoleAmpereSquareMetres: Double
    private let command: (AttitudeState,SpacecraftState) throws -> SIMD3<Double>
    public init(fieldProvider: any MagneticFieldProvider, maximumDipoleAmpereSquareMetres: Double,
        command: @escaping (AttitudeState,SpacecraftState) throws -> SIMD3<Double>) throws {
        guard maximumDipoleAmpereSquareMetres.isFinite, maximumDipoleAmpereSquareMetres > 0 else { throw AttitudeError.invalidConfiguration }
        self.fieldProvider = fieldProvider; self.maximumDipoleAmpereSquareMetres = maximumDipoleAmpereSquareMetres; self.command = command
    }
    public func contribution(for attitude: AttitudeState, spacecraft: SpacecraftState?, in context: RotationalEvaluationContext) throws -> TorqueContribution {
        try validateRotationalState(attitude,spacecraft: spacecraft,context: context)
        guard let spacecraft else { throw AttitudeError.invalidConfiguration }
        let field = try fieldProvider.field(at: spacecraft)
        guard field.epoch == attitude.epoch else { throw AttitudeError.epochMismatch }
        guard field.frame === context.referenceFrame else { throw AttitudeError.frameMismatch }
        guard field.origin == context.referenceContext.origin else { throw GravityModelError.originMismatch }
        var dipole = try command(attitude,spacecraft)
        guard attitudeFinite(dipole) else { throw AttitudeError.invalidVector }
        let scale = max(abs(dipole.x),abs(dipole.y),abs(dipole.z))
        if scale > 0 {
            let scaled = dipole/scale, length = routeNorm(scaled)
            if scale > maximumDipoleAmpereSquareMetres/length { dipole = scaled*(maximumDipoleAmpereSquareMetres/length) }
        }
        let torque = attitudeCross(dipole,try attitude.bodyVector(fromReference: field.fieldTeslas))
        return try TorqueContribution(epoch: attitude.epoch,bodyFrame: context.bodyFrame,torqueBodyNewtonMetres: torque)
    }
}
