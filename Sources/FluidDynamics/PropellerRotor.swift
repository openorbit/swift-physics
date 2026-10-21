// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Mattias Holm
import Foundation

/// Quasi-steady coefficient rotor coupled to the existing actuator disk. Coefficients
/// are caller-supplied engineering assumptions, not a blade-element solution.
/// CT uses rho*A*(omega*R)^2; CP0 uses rho*A*(omega*R)^3.
public struct PropellerRotor:Sendable {
    public let radiusMetres,maximumThrustCoefficient,profilePowerCoefficient,inertiaKilogramMetresSquared:Double
    public init(radiusMetres:Double,maximumThrustCoefficient:Double,profilePowerCoefficient:Double,inertiaKilogramMetresSquared:Double) throws {
        guard [radiusMetres,maximumThrustCoefficient,profilePowerCoefficient,inertiaKilogramMetresSquared].allSatisfy({$0.isFinite && $0>0}) else {throw FluidError.invalidParameter}
        self.radiusMetres=radiusMetres;self.maximumThrustCoefficient=maximumThrustCoefficient
        self.profilePowerCoefficient=profilePowerCoefficient;self.inertiaKilogramMetresSquared=inertiaKilogramMetresSquared
    }
    /// Positive axial speed is climbing into still air. Descending/vortex-ring and
    /// windmilling regimes are deliberately not resolved. Coaxial interference is
    /// supplied as nonnegative added inflow, without adding its thrust twice.
    public func sample(angularSpeed:Double,collective:Double,density:Double,axialSpeed:Double=0,addedInflow:Double=0) throws -> PropellerRotorSample {
        guard [angularSpeed,collective,density,axialSpeed,addedInflow].allSatisfy(\.isFinite),angularSpeed>=0,(0...1).contains(collective),density>0,addedInflow>=0 else {throw FluidError.invalidParameter}
        let area=Double.pi*radiusMetres*radiusMetres,tip=angularSpeed*radiusMetres
        let thrust=density*area*tip*tip*maximumThrustCoefficient*collective
        let axial=max(0,axialSpeed)+addedInflow
        let disk=try ActuatorDiskSlipstream(centreMetres:.zero,downstreamAxis:SIMD3(0,0,-1),radiusMetres:radiusMetres,thrustNewtons:thrust,densityKilogramsPerCubicMetre:density,upstreamAxialSpeedMetresPerSecond:axial)
        let power=disk.idealAxialPowerWatts+profilePowerCoefficient*density*area*tip*tip*tip
        let torque=angularSpeed>0 ? power/angularSpeed : 0
        guard thrust.isFinite,power.isFinite,torque.isFinite else {throw FluidError.nonFiniteResult}
        return PropellerRotorSample(thrustNewtons:thrust,loadTorqueNewtonMetres:torque,shaftPowerWatts:power,inducedSpeed:disk.inducedSpeedAtDiskMetresPerSecond)
    }
}
public struct PropellerRotorSample:Sendable {
    public let thrustNewtons,loadTorqueNewtonMetres,shaftPowerWatts,inducedSpeed:Double
}
