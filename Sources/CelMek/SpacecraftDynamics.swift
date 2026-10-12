// SPDX-License-Identifier: Apache-2.0
import Math

public enum SpacecraftDynamicsError: Error, Equatable {
    case invalidMass, invalidMassRate, nonFiniteAcceleration
    case frameMismatch, epochMismatch, nonInertialFrame, numericalRangeExceeded
}

/// Translational mass properties. Dry mass includes everything not represented
/// by this propellant reservoir; inertia tensors belong to rotational dynamics.
public struct SpacecraftMassProperties: Equatable, Sendable {
    public let dryMassKilograms: Double
    public let propellantMassKilograms: Double
    public var totalMassKilograms: Double { dryMassKilograms + propellantMassKilograms }
    public init(dryMassKilograms: Double, propellantMassKilograms: Double = 0) throws {
        guard dryMassKilograms.isFinite, dryMassKilograms > 0,
              propellantMassKilograms.isFinite, propellantMassKilograms >= 0,
              (dryMassKilograms+propellantMassKilograms).isFinite else {
            throw SpacecraftDynamicsError.invalidMass
        }
        self.dryMassKilograms = dryMassKilograms; self.propellantMassKilograms = propellantMassKilograms
    }
}

/// Position, velocity, epoch and frame plus the instantaneous mass inventory.
/// Stored Cartesian acceleration is optional observational data; force evaluation
/// computes its own derivative and does not add that acceleration again.
public struct SpacecraftState: Equatable {
    public let cartesian: CartesianState
    public let mass: SpacecraftMassProperties
    public var epoch: Epoch { cartesian.epoch }
    public var frame: Frame { cartesian.frame }
    public init(cartesian: CartesianState, mass: SpacecraftMassProperties) {
        self.cartesian = cartesian; self.mass = mass
    }
}

/// Origin and environmental providers for one dynamics evaluation. Model/time
/// coverage remains the responsibility of the supplied ephemeris providers.
/// Models must explicitly convert queried ephemerides to the state's coordinates.
public struct ForceEvaluationContext {
    public let frame: Frame
    public let origin: BodyIdentifier
    public let ephemerides: EphemerisCollection?
    public init(frame: Frame, origin: BodyIdentifier, ephemerides: EphemerisCollection? = nil) throws {
        try origin.validate()
        self.frame = frame; self.origin = origin; self.ephemerides = ephemerides
    }
}

/// One force model's contribution in the state's native axes. Derivatives use
/// atomic SI seconds: m/s² and kg/s. Negative mass rate consumes propellant;
/// positive rates describe an explicitly modelled inflow. No clipping is applied.
public struct ForceContribution {
    public let epoch: Epoch
    public let frame: Frame
    public let accelerationMetresPerSecondSquared: SIMD3<Double>
    public let propellantMassRateKilogramsPerSecond: Double
    public init(epoch: Epoch, frame: Frame, accelerationMetresPerSecondSquared: SIMD3<Double>,
                propellantMassRateKilogramsPerSecond: Double = 0) throws {
        guard [accelerationMetresPerSecondSquared.x,accelerationMetresPerSecondSquared.y,
               accelerationMetresPerSecondSquared.z].allSatisfy(\.isFinite) else {
            throw SpacecraftDynamicsError.nonFiniteAcceleration
        }
        guard propellantMassRateKilogramsPerSecond.isFinite else { throw SpacecraftDynamicsError.invalidMassRate }
        self.epoch = epoch; self.frame = frame
        self.accelerationMetresPerSecondSquared = accelerationMetresPerSecondSquared
        self.propellantMassRateKilogramsPerSecond = propellantMassRateKilogramsPerSecond
    }
}

/// Pure evaluation at an explicit state. Numerical solvers may repeat evaluations
/// or visit epochs out of order; models must not consume fuel or advance hidden time.
/// A model may use prescribed gravity acceleration or divide thrust force by mass.
public protocol ForceModel {
    func contribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution
}

/// Right-hand side for integration of Cartesian position, velocity and
/// propellant mass. Dry mass is constant in this initial translational model.
public struct SpacecraftDerivative {
    public let epoch: Epoch
    public let frame: Frame
    public let positionRateMetresPerSecond: SIMD3<Double>
    public let velocityRateMetresPerSecondSquared: SIMD3<Double>
    public let propellantMassRateKilogramsPerSecond: Double
}

/// Ordered, immutable force composition for inertial translational dynamics.
/// Empty models mean unforced motion, not implicit central gravity. Provider
/// thread safety remains each model's responsibility; no global environment/cache.
public struct SpacecraftDynamics {
    public let context: ForceEvaluationContext
    private let models: [any ForceModel]
    public init(models: [any ForceModel], context: ForceEvaluationContext) throws {
        guard context.frame.pseudoInertial else { throw SpacecraftDynamicsError.nonInertialFrame }
        self.models = models; self.context = context
    }
    public func derivative(at state: SpacecraftState) throws -> SpacecraftDerivative {
        guard state.frame === context.frame else { throw SpacecraftDynamicsError.frameMismatch }
        var acceleration = SIMD3<Double>.zero, massRate = 0.0
        for model in models {
            let contribution = try model.contribution(for: state,in: context)
            guard contribution.epoch == state.epoch else { throw SpacecraftDynamicsError.epochMismatch }
            guard contribution.frame === state.frame else { throw SpacecraftDynamicsError.frameMismatch }
            acceleration += contribution.accelerationMetresPerSecondSquared
            massRate += contribution.propellantMassRateKilogramsPerSecond
            guard [acceleration.x,acceleration.y,acceleration.z,massRate].allSatisfy(\.isFinite) else {
                throw SpacecraftDynamicsError.numericalRangeExceeded
            }
        }
        return SpacecraftDerivative(epoch: state.epoch,frame: state.frame,
            positionRateMetresPerSecond: state.cartesian.velocity,
            velocityRateMetresPerSecondSquared: acceleration,propellantMassRateKilogramsPerSecond: massRate)
    }
}


/// Consuming forces expose a continuous dry-mass limit for root localization.
/// This must equal normal evaluation for positive fuel and have nonpositive mass
/// rate. It does not authorize actual firing from an empty reservoir.
public protocol FuelConsumingForceModel: ForceModel {
    func depletionBoundaryContribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution
}
public enum FuelBoundaryError: Error, Equatable {
    case unsupportedMassFlow, positiveConsumptionRate
}
extension SpacecraftDynamics {
    func withoutFuelConsumers() throws -> SpacecraftDynamics {
        try SpacecraftDynamics(models: models.filter { !($0 is any FuelConsumingForceModel) }
            .map { FuelBoundaryModel(model: $0) as any ForceModel }, context: context)
    }
    func fuelBoundaryDerivative(at state: SpacecraftState) throws -> SpacecraftDerivative {
        let wrapped = models.map { FuelBoundaryModel(model: $0) as any ForceModel }
        return try SpacecraftDynamics(models: wrapped, context: context).derivative(at: state)
    }
}
private struct FuelBoundaryModel: ForceModel {
    let model: any ForceModel
    func contribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
        if let consumer = model as? any FuelConsumingForceModel {
            let result = try consumer.depletionBoundaryContribution(for: state, in: context)
            guard result.propellantMassRateKilogramsPerSecond <= 0 else { throw FuelBoundaryError.positiveConsumptionRate }
            return result
        }
        let result = try model.contribution(for: state, in: context)
        guard result.propellantMassRateKilogramsPerSecond == 0 else { throw FuelBoundaryError.unsupportedMassFlow }
        return result
    }
}
