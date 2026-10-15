// SPDX-License-Identifier: Apache-2.0
import Foundation

/// Physical air velocity in the caller's common axes and reference frame.
/// Providers document whether they return total velocity or an induced increment.
public protocol AirVelocityField: Sendable {
    func velocityMetresPerSecond(at positionMetres: SIMD3<Double>) throws -> SIMD3<Double>
}

/// A velocity increment, suitable for addition to one supplied background field.
public protocol InducedAirVelocityField: AirVelocityField {}

/// Prescribed-thrust axial actuator disk with a uniform contracting slipstream and
/// solid-body swirl carrying the supplied torque. Returns induced air velocity only.
/// No engine/RPM dynamics, blade loading, crossflow, windmilling or obstacle deflection.
public struct ActuatorDiskPropwash: InducedAirVelocityField {
    public let centreMetres: SIMD3<Double>
    /// Points downstream, along force exerted on the air (opposite aircraft thrust).
    public let downstreamAxis: SIMD3<Double>
    public let radiusMetres: Double
    public let thrustNewtons: Double
    /// Signed torque exerted on the air about downstreamAxis.
    public let airTorqueNewtonMetres: Double
    public let densityKilogramsPerCubicMetre: Double
    public let upstreamAxialSpeedMetresPerSecond: Double
    public let inducedSpeedAtDiskMetresPerSecond: Double
    public let massFlowKilogramsPerSecond: Double
    public var idealAxialPowerWatts: Double {
        thrustNewtons*(upstreamAxialSpeedMetresPerSecond+inducedSpeedAtDiskMetresPerSecond)
    }
    public init(centreMetres: SIMD3<Double>, downstreamAxis: SIMD3<Double>, radiusMetres: Double,
                thrustNewtons: Double, airTorqueNewtonMetres: Double = 0,
                densityKilogramsPerCubicMetre: Double, upstreamAxialSpeedMetresPerSecond: Double) throws {
        let length = aeroNorm(downstreamAxis)
        guard aeroFinite(centreMetres), aeroFinite(downstreamAxis), length.isFinite, length > 0,
              radiusMetres.isFinite, radiusMetres > 0, thrustNewtons.isFinite, thrustNewtons >= 0,
              airTorqueNewtonMetres.isFinite, densityKilogramsPerCubicMetre.isFinite, densityKilogramsPerCubicMetre > 0,
              upstreamAxialSpeedMetresPerSecond.isFinite, upstreamAxialSpeedMetresPerSecond >= 0,
              thrustNewtons > 0 || airTorqueNewtonMetres == 0 else { throw AerodynamicError.invalidParameter }
        let area = Double.pi*radiusMetres*radiusMetres, v = upstreamAxialSpeedMetresPerSecond
        let loading = thrustNewtons/(densityKilogramsPerCubicMetre*area)
        let vi = thrustNewtons == 0 ? 0 : loading/(hypot(v,sqrt(2*loading))+v)
        let massFlow = densityKilogramsPerCubicMetre*area*(v+vi)
        guard area.isFinite, area > 0, loading.isFinite, vi.isFinite, massFlow.isFinite,
              thrustNewtons == 0 || (loading > 0 && vi > 0 && massFlow > 0),
              (thrustNewtons*(v+vi)).isFinite else { throw AerodynamicError.nonFiniteResult }
        self.centreMetres = centreMetres; self.downstreamAxis = downstreamAxis/length; self.radiusMetres = radiusMetres
        self.thrustNewtons = thrustNewtons; self.airTorqueNewtonMetres = airTorqueNewtonMetres
        self.densityKilogramsPerCubicMetre = densityKilogramsPerCubicMetre
        self.upstreamAxialSpeedMetresPerSecond = v; self.inducedSpeedAtDiskMetresPerSecond = vi
        self.massFlowKilogramsPerSecond = massFlow
    }
    /// Chosen downstream interpolation: vi(z)=vi_disk*(1+z/sqrt(z²+R²)).
    /// Only disk/far-wake endpoints are fixed by momentum theory; intermediate
    /// development and top-hat cross-section are explicit idealizations.
    public func wakeRadiusMetres(downstreamDistanceMetres z: Double) throws -> Double {
        guard z.isFinite, z >= 0 else { throw AerodynamicError.outsideCoverage }
        let vi = inducedSpeedAtDiskMetresPerSecond, v = upstreamAxialSpeedMetresPerSecond
        if thrustNewtons == 0 { return radiusMetres }
        let speed = v+vi*(1+z/hypot(z,radiusMetres))
        return radiusMetres*sqrt((v+vi)/speed)
    }
    public func velocityMetresPerSecond(at positionMetres: SIMD3<Double>) throws -> SIMD3<Double> {
        guard aeroFinite(positionMetres) else { throw AerodynamicError.invalidParameter }
        let offset = positionMetres-centreMetres, z = (offset*downstreamAxis).sum()
        guard aeroFinite(offset), z.isFinite else { throw AerodynamicError.nonFiniteResult }
        // Deliberately downstream-only. Upstream induction and the sharp wake edge
        // are not resolved; this is a prescribed slipstream, not a divergence-free CFD solution.
        if z < 0 || thrustNewtons == 0 { return .zero }
        let radial = offset-downstreamAxis*z, radius = try wakeRadiusMetres(downstreamDistanceMetres: z)
        if aeroNorm(radial) > radius { return .zero }
        let axial = inducedSpeedAtDiskMetresPerSecond*(1+z/hypot(z,radiusMetres))
        let omega = 2*airTorqueNewtonMetres/(massFlowKilogramsPerSecond*radius*radius)
        let result = downstreamAxis*axial+aeroCross(downstreamAxis,radial)*omega
        guard aeroFinite(result) else { throw AerodynamicError.nonFiniteResult }
        return result
    }
    /// Equal/opposite prescribed propeller force and torque, including the lever arm.
    /// Use once; this is not an additional force to add to an existing propulsion model.
    public func reactionLoads(momentReferencePointMetres: SIMD3<Double> = .zero) throws -> AerodynamicLoads {
        try AerodynamicLoads(forceNewtons: -downstreamAxis*thrustNewtons,
            momentNewtonMetres: -downstreamAxis*airTorqueNewtonMetres,referencePointMetres: centreMetres)
            .referred(to: momentReferencePointMetres)
    }
}

/// Sum prescribed induced increments and a background air velocity. A CFD field
/// already contains its background; do not add the same background a second time.
public struct SuperposedAirVelocityField: AirVelocityField {
    public let backgroundMetresPerSecond: SIMD3<Double>
    public let increments: [any InducedAirVelocityField]
    public init(backgroundMetresPerSecond: SIMD3<Double> = .zero, increments: [any InducedAirVelocityField]) throws {
        guard aeroFinite(backgroundMetresPerSecond) else { throw AerodynamicError.invalidFlow }
        self.backgroundMetresPerSecond = backgroundMetresPerSecond; self.increments = increments
    }
    public func velocityMetresPerSecond(at positionMetres: SIMD3<Double>) throws -> SIMD3<Double> {
        guard aeroFinite(positionMetres) else { throw AerodynamicError.invalidParameter }
        var velocity = backgroundMetresPerSecond
        for field in increments { velocity += try field.velocityMetresPerSecond(at: positionMetres) }
        guard aeroFinite(velocity) else { throw AerodynamicError.nonFiniteResult }
        return velocity
    }
}

/// Alternative to prescribing the analytic slipstream: apply the rotor's force
/// and torque to CFD air over a finite cylinder centred at the disk. The continuous
/// volume integral is thrust/torque; face-sampled grid quadrature has resolution
/// error and must be checked. Do not also add the analytic propwash velocity.
public struct ActuatorDiskAccelerationField: Sendable {
    public let disk: ActuatorDiskPropwash
    public let thicknessMetres: Double
    public init(disk: ActuatorDiskPropwash, thicknessMetres: Double) throws {
        guard thicknessMetres.isFinite, thicknessMetres > 0,
              (Double.pi*disk.radiusMetres*disk.radiusMetres*thicknessMetres).isFinite,
              Double.pi*disk.radiusMetres*disk.radiusMetres*thicknessMetres > 0 else {
            throw AerodynamicError.invalidParameter
        }
        let fluidMass = disk.densityKilogramsPerCubicMetre*Double.pi*disk.radiusMetres*disk.radiusMetres*thicknessMetres
        guard fluidMass.isFinite, fluidMass > 0 else { throw AerodynamicError.nonFiniteResult }
        self.disk = disk; self.thicknessMetres = thicknessMetres
    }
    public func accelerationMetresPerSecondSquared(at positionMetres: SIMD3<Double>) throws -> SIMD3<Double> {
        guard aeroFinite(positionMetres) else { throw AerodynamicError.invalidParameter }
        let offset = positionMetres-disk.centreMetres, axial = (offset*disk.downstreamAxis).sum()
        guard aeroFinite(offset), axial.isFinite else { throw AerodynamicError.nonFiniteResult }
        let radial = offset-disk.downstreamAxis*axial
        if abs(axial) > thicknessMetres/2 || aeroNorm(radial) > disk.radiusMetres { return .zero }
        let fluidMass = disk.densityKilogramsPerCubicMetre*Double.pi*disk.radiusMetres*disk.radiusMetres*thicknessMetres
        let axialAcceleration = disk.thrustNewtons/fluidMass
        let angularAcceleration = 2*disk.airTorqueNewtonMetres/(fluidMass*disk.radiusMetres*disk.radiusMetres)
        let result = disk.downstreamAxis*axialAcceleration+aeroCross(disk.downstreamAxis,radial)*angularAcceleration
        guard aeroFinite(result) else { throw AerodynamicError.nonFiniteResult }
        return result
    }
}
