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

/// Physical fluid velocity in the caller's common axes and reference frame.
/// Providers document whether they return total velocity or an induced increment.
public protocol FluidVelocityField: Sendable {
    func velocityMetresPerSecond(at positionMetres: SIMD3<Double>) throws -> SIMD3<Double>
}

/// A velocity increment, suitable for addition to one supplied background field.
public protocol InducedFluidVelocityField: FluidVelocityField {}

/// Prescribed-thrust axial actuator disk with a uniform contracting slipstream and
/// solid-body swirl carrying the supplied torque. Returns induced fluid velocity only.
/// No engine/RPM dynamics, blade loading, crossflow, windmilling or obstacle deflection.
public struct ActuatorDiskSlipstream: InducedFluidVelocityField {
    public let centreMetres: SIMD3<Double>
    /// Points downstream, along force exerted on the fluid (opposite aircraft thrust).
    public let downstreamAxis: SIMD3<Double>
    public let radiusMetres: Double
    public let thrustNewtons: Double
    /// Signed torque exerted on the fluid about downstreamAxis.
    public let fluidTorqueNewtonMetres: Double
    public let densityKilogramsPerCubicMetre: Double
    public let upstreamAxialSpeedMetresPerSecond: Double
    public let inducedSpeedAtDiskMetresPerSecond: Double
    public let massFlowKilogramsPerSecond: Double
    public var idealAxialPowerWatts: Double {
        thrustNewtons*(upstreamAxialSpeedMetresPerSecond+inducedSpeedAtDiskMetresPerSecond)
    }
    public init(centreMetres: SIMD3<Double>, downstreamAxis: SIMD3<Double>, radiusMetres: Double,
                thrustNewtons: Double, fluidTorqueNewtonMetres: Double = 0,
                densityKilogramsPerCubicMetre: Double, upstreamAxialSpeedMetresPerSecond: Double) throws {
        let length = fluidNorm(downstreamAxis)
        guard fluidFinite(centreMetres), fluidFinite(downstreamAxis), length.isFinite, length > 0,
              radiusMetres.isFinite, radiusMetres > 0, thrustNewtons.isFinite, thrustNewtons >= 0,
              fluidTorqueNewtonMetres.isFinite, densityKilogramsPerCubicMetre.isFinite, densityKilogramsPerCubicMetre > 0,
              upstreamAxialSpeedMetresPerSecond.isFinite, upstreamAxialSpeedMetresPerSecond >= 0,
              thrustNewtons > 0 || fluidTorqueNewtonMetres == 0 else { throw FluidError.invalidParameter }
        let area = Double.pi*radiusMetres*radiusMetres, v = upstreamAxialSpeedMetresPerSecond
        let loading = thrustNewtons/(densityKilogramsPerCubicMetre*area)
        let vi = thrustNewtons == 0 ? 0 : loading/(hypot(v,sqrt(2*loading))+v)
        let massFlow = densityKilogramsPerCubicMetre*area*(v+vi)
        guard area.isFinite, area > 0, loading.isFinite, vi.isFinite, massFlow.isFinite,
              thrustNewtons == 0 || (loading > 0 && vi > 0 && massFlow > 0),
              (thrustNewtons*(v+vi)).isFinite else { throw FluidError.nonFiniteResult }
        self.centreMetres = centreMetres; self.downstreamAxis = downstreamAxis/length; self.radiusMetres = radiusMetres
        self.thrustNewtons = thrustNewtons; self.fluidTorqueNewtonMetres = fluidTorqueNewtonMetres
        self.densityKilogramsPerCubicMetre = densityKilogramsPerCubicMetre
        self.upstreamAxialSpeedMetresPerSecond = v; self.inducedSpeedAtDiskMetresPerSecond = vi
        self.massFlowKilogramsPerSecond = massFlow
    }
    /// Chosen downstream interpolation: vi(z)=vi_disk*(1+z/sqrt(z²+R²)).
    /// Only disk/far-wake endpoints are fixed by momentum theory; intermediate
    /// development and top-hat cross-section are explicit idealizations.
    public func wakeRadiusMetres(downstreamDistanceMetres z: Double) throws -> Double {
        guard z.isFinite, z >= 0 else { throw FluidError.outsideCoverage }
        let vi = inducedSpeedAtDiskMetresPerSecond, v = upstreamAxialSpeedMetresPerSecond
        if thrustNewtons == 0 { return radiusMetres }
        let speed = v+vi*(1+z/hypot(z,radiusMetres))
        return radiusMetres*sqrt((v+vi)/speed)
    }
    public func velocityMetresPerSecond(at positionMetres: SIMD3<Double>) throws -> SIMD3<Double> {
        guard fluidFinite(positionMetres) else { throw FluidError.invalidParameter }
        let offset = positionMetres-centreMetres, z = (offset*downstreamAxis).sum()
        guard fluidFinite(offset), z.isFinite else { throw FluidError.nonFiniteResult }
        // Deliberately downstream-only. Upstream induction and the sharp wake edge
        // are not resolved; this is a prescribed slipstream, not a divergence-free CFD solution.
        if z < 0 || thrustNewtons == 0 { return .zero }
        let radial = offset-downstreamAxis*z, radius = try wakeRadiusMetres(downstreamDistanceMetres: z)
        if fluidNorm(radial) > radius { return .zero }
        let axial = inducedSpeedAtDiskMetresPerSecond*(1+z/hypot(z,radiusMetres))
        let omega = 2*fluidTorqueNewtonMetres/(massFlowKilogramsPerSecond*radius*radius)
        let result = downstreamAxis*axial+fluidCross(downstreamAxis,radial)*omega
        guard fluidFinite(result) else { throw FluidError.nonFiniteResult }
        return result
    }
    /// Equal/opposite prescribed propeller force and torque, including the lever arm.
    /// Use once; this is not an additional force to add to an existing propulsion model.
    public func reactionLoads(momentReferencePointMetres: SIMD3<Double> = .zero) throws -> FluidLoads {
        try FluidLoads(forceNewtons: -downstreamAxis*thrustNewtons,
            momentNewtonMetres: -downstreamAxis*fluidTorqueNewtonMetres,referencePointMetres: centreMetres)
            .referred(to: momentReferencePointMetres)
    }
}

/// Sum prescribed induced increments and a background fluid velocity. A CFD field
/// already contains its background; do not add the same background a second time.
public struct SuperposedFluidVelocityField: FluidVelocityField {
    public let backgroundMetresPerSecond: SIMD3<Double>
    public let increments: [any InducedFluidVelocityField]
    public init(backgroundMetresPerSecond: SIMD3<Double> = .zero, increments: [any InducedFluidVelocityField]) throws {
        guard fluidFinite(backgroundMetresPerSecond) else { throw FluidError.invalidFlow }
        self.backgroundMetresPerSecond = backgroundMetresPerSecond; self.increments = increments
    }
    public func velocityMetresPerSecond(at positionMetres: SIMD3<Double>) throws -> SIMD3<Double> {
        guard fluidFinite(positionMetres) else { throw FluidError.invalidParameter }
        var velocity = backgroundMetresPerSecond
        for field in increments { velocity += try field.velocityMetresPerSecond(at: positionMetres) }
        guard fluidFinite(velocity) else { throw FluidError.nonFiniteResult }
        return velocity
    }
}

/// Alternative to prescribing the analytic slipstream: apply the rotor's force
/// and torque to CFD fluid over a finite cylinder centred at the disk. The continuous
/// volume integral is thrust/torque; face-sampled grid quadrature has resolution
/// error and must be checked. Do not also add the analytic propwash velocity.
public struct ActuatorDiskAccelerationField: Sendable {
    public let disk: ActuatorDiskSlipstream
    public let thicknessMetres: Double
    public init(disk: ActuatorDiskSlipstream, thicknessMetres: Double) throws {
        guard thicknessMetres.isFinite, thicknessMetres > 0,
              (Double.pi*disk.radiusMetres*disk.radiusMetres*thicknessMetres).isFinite,
              Double.pi*disk.radiusMetres*disk.radiusMetres*thicknessMetres > 0 else {
            throw FluidError.invalidParameter
        }
        let fluidMass = disk.densityKilogramsPerCubicMetre*Double.pi*disk.radiusMetres*disk.radiusMetres*thicknessMetres
        guard fluidMass.isFinite, fluidMass > 0 else { throw FluidError.nonFiniteResult }
        self.disk = disk; self.thicknessMetres = thicknessMetres
    }
    public func accelerationMetresPerSecondSquared(at positionMetres: SIMD3<Double>) throws -> SIMD3<Double> {
        guard fluidFinite(positionMetres) else { throw FluidError.invalidParameter }
        let offset = positionMetres-disk.centreMetres, axial = (offset*disk.downstreamAxis).sum()
        guard fluidFinite(offset), axial.isFinite else { throw FluidError.nonFiniteResult }
        let radial = offset-disk.downstreamAxis*axial
        if abs(axial) > thicknessMetres/2 || fluidNorm(radial) > disk.radiusMetres { return .zero }
        let fluidMass = disk.densityKilogramsPerCubicMetre*Double.pi*disk.radiusMetres*disk.radiusMetres*thicknessMetres
        let axialAcceleration = disk.thrustNewtons/fluidMass
        let angularAcceleration = 2*disk.fluidTorqueNewtonMetres/(fluidMass*disk.radiusMetres*disk.radiusMetres)
        let result = disk.downstreamAxis*axialAcceleration+fluidCross(disk.downstreamAxis,radial)*angularAcceleration
        guard fluidFinite(result) else { throw FluidError.nonFiniteResult }
        return result
    }
}
