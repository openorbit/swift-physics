// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum FiniteThrustError: Error, Equatable {
    case invalidThrust, invalidExhaustVelocity, invalidDirection, invalidInterval, tooManyBurns
}

/// Constant prescribed thrust in inertial axes, drawing from the single reservoir.
/// Effective exhaust velocity includes nozzle pressure effects: mass flow = thrust / ve.
public struct FiniteThrust: ForceModel {
    public let frame: Frame
    public let thrustNewtons: Double
    public let effectiveExhaustVelocityMetresPerSecond: Double
    public let direction: SIMD3<Double>
    public let propellantConsumptionKilogramsPerSecond: Double

    public init(frame: Frame, thrustNewtons: Double, effectiveExhaustVelocityMetresPerSecond: Double,
                direction: SIMD3<Double>) throws {
        guard frame.pseudoInertial else { throw SpacecraftDynamicsError.nonInertialFrame }
        guard thrustNewtons.isFinite, thrustNewtons >= 0 else { throw FiniteThrustError.invalidThrust }
        guard effectiveExhaustVelocityMetresPerSecond.isFinite, effectiveExhaustVelocityMetresPerSecond > 0 else {
            throw FiniteThrustError.invalidExhaustVelocity
        }
        guard [direction.x,direction.y,direction.z].allSatisfy(\.isFinite) else { throw FiniteThrustError.invalidDirection }
        let scale = max(abs(direction.x),abs(direction.y),abs(direction.z))
        guard scale > 0 else { throw FiniteThrustError.invalidDirection }
        let scaled = direction / scale
        let flow = thrustNewtons / effectiveExhaustVelocityMetresPerSecond
        guard flow.isFinite, thrustNewtons == 0 || flow > 0 else { throw FiniteThrustError.invalidExhaustVelocity }
        self.frame = frame
        self.thrustNewtons = thrustNewtons
        self.effectiveExhaustVelocityMetresPerSecond = effectiveExhaustVelocityMetresPerSecond
        self.direction = scaled / sqrt((scaled*scaled).sum())
        propellantConsumptionKilogramsPerSecond = flow
    }

    public func contribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
        guard state.frame === frame, context.frame === frame else { throw SpacecraftDynamicsError.frameMismatch }
        guard thrustNewtons == 0 || state.mass.propellantMassKilograms > 0 else {
            throw NumericalPropagationError.propellantDepleted
        }
        return try ForceContribution(epoch: state.epoch, frame: frame,
            accelerationMetresPerSecondSquared: direction * (thrustNewtons / state.mass.totalMassKilograms),
            propellantMassRateKilogramsPerSecond: -propellantConsumptionKilogramsPerSecond)
    }
}

/// Active on [start, end); independent engines may have overlapping burns.
public struct FiniteBurn {
    public let start: Epoch
    public let end: Epoch
    public let thrust: FiniteThrust
    public init(start: Epoch, end: Epoch, thrust: FiniteThrust) throws {
        guard end > start else { throw FiniteThrustError.invalidInterval }
        self.start = start; self.end = end; self.thrust = thrust
    }
}

/// Stateless absolute-epoch schedule. Use ScheduledThrustPropagator to ensure
/// numerical steps land on its discontinuities; a generic integrator can miss burns.
public struct ScheduledThrust: ForceModel {
    public let frame: Frame
    public let burns: [FiniteBurn]
    public init(frame: Frame, burns: [FiniteBurn]) throws {
        guard frame.pseudoInertial else { throw SpacecraftDynamicsError.nonInertialFrame }
        guard burns.count <= 10_000 else { throw FiniteThrustError.tooManyBurns }
        guard burns.allSatisfy({ $0.thrust.frame === frame }) else { throw SpacecraftDynamicsError.frameMismatch }
        self.frame = frame; self.burns = burns.sorted { $0.start < $1.start }
    }
    public func contribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
        guard state.frame === frame, context.frame === frame else { throw SpacecraftDynamicsError.frameMismatch }
        var acceleration = SIMD3<Double>.zero, flow = 0.0
        for burn in burns where state.epoch >= burn.start && state.epoch < burn.end {
            let contribution = try burn.thrust.contribution(for: state, in: context)
            acceleration += contribution.accelerationMetresPerSecondSquared
            flow += contribution.propellantMassRateKilogramsPerSecond
        }
        return try ForceContribution(epoch: state.epoch, frame: frame,
            accelerationMetresPerSecondSquared: acceleration, propellantMassRateKilogramsPerSecond: flow)
    }
}
