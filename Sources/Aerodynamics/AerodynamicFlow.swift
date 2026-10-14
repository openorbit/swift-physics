// SPDX-License-Identifier: Apache-2.0
import Foundation
import Atmospheres

/// Pure local flow in caller-selected right-handed body/section axes:
/// +X forward, +Y right, +Z down. Velocity is body motion relative to the medium,
/// not the opposite air-on-body velocity. All vectors in a query must use the same axes.
public struct AerodynamicFlow: Sendable {
    public let relativeVelocityMetresPerSecond: SIMD3<Double>
    public let densityKilogramsPerCubicMetre: Double
    public let pressurePascals: Double?
    public let temperatureKelvin: Double?
    public let speedOfSoundMetresPerSecond: Double?
    public let dynamicViscosityPascalSeconds: Double?
    public let meanFreePathMetres: Double?
    public var speedMetresPerSecond: Double { aeroNorm(relativeVelocityMetresPerSecond) }
    public var dynamicPressurePascals: Double { (0.5*densityKilogramsPerCubicMetre*speedMetresPerSecond)*speedMetresPerSecond }
    public var machNumber: Double? { speedOfSoundMetresPerSecond.map { speedMetresPerSecond/$0 } }
    /// Positive incidence has a positive downward velocity component in body axes.
    public var angleOfAttackRadians: Double? {
        let v = relativeVelocityMetresPerSecond
        return v.x == 0 && v.z == 0 ? nil : atan2(v.z,v.x)
    }
    public var sideslipRadians: Double? {
        let v = relativeVelocityMetresPerSecond
        return speedMetresPerSecond == 0 ? nil : atan2(v.y,hypot(v.x,v.z))
    }
    public init(relativeVelocityMetresPerSecond: SIMD3<Double>, densityKilogramsPerCubicMetre: Double,
                pressurePascals: Double? = nil, temperatureKelvin: Double? = nil,
                speedOfSoundMetresPerSecond: Double? = nil, dynamicViscosityPascalSeconds: Double? = nil,
                meanFreePathMetres: Double? = nil) throws {
        guard aeroFinite(relativeVelocityMetresPerSecond), densityKilogramsPerCubicMetre.isFinite,
              densityKilogramsPerCubicMetre >= 0,
              pressurePascals.map({ $0.isFinite && $0 >= 0 }) ?? true,
              [temperatureKelvin,speedOfSoundMetresPerSecond,dynamicViscosityPascalSeconds,meanFreePathMetres]
                .compactMap({ $0 }).allSatisfy({ $0.isFinite && $0 > 0 }) else { throw AerodynamicError.invalidFlow }
        self.relativeVelocityMetresPerSecond = relativeVelocityMetresPerSecond
        self.densityKilogramsPerCubicMetre = densityKilogramsPerCubicMetre
        self.pressurePascals = pressurePascals; self.temperatureKelvin = temperatureKelvin
        self.speedOfSoundMetresPerSecond = speedOfSoundMetresPerSecond
        self.dynamicViscosityPascalSeconds = dynamicViscosityPascalSeconds; self.meanFreePathMetres = meanFreePathMetres
        guard speedMetresPerSecond.isFinite, dynamicPressurePascals.isFinite,
              machNumber.map(\.isFinite) ?? true else { throw AerodynamicError.nonFiniteResult }
    }
    public init(relativeVelocityMetresPerSecond: SIMD3<Double>, properties: AtmosphericProperties,
                meanFreePathMetres: Double? = nil) throws {
        try self.init(relativeVelocityMetresPerSecond: relativeVelocityMetresPerSecond,
            densityKilogramsPerCubicMetre: properties.densityKilogramsPerCubicMetre,
            pressurePascals: properties.pressurePascals,temperatureKelvin: properties.temperatureKelvin,
            speedOfSoundMetresPerSecond: properties.speedOfSoundMetresPerSecond,
            dynamicViscosityPascalSeconds: properties.dynamicViscosityPascalSeconds,meanFreePathMetres: meanFreePathMetres)
    }
    public func reynoldsNumber(referenceLengthMetres length: Double) throws -> Double? {
        try aeroPositiveLength(length)
        guard let viscosity = dynamicViscosityPascalSeconds else { return nil }
        let result = (densityKilogramsPerCubicMetre*speedMetresPerSecond)*length/viscosity
        guard result.isFinite else { throw AerodynamicError.nonFiniteResult }
        return result
    }
    public func knudsenNumber(referenceLengthMetres length: Double) throws -> Double? {
        try aeroPositiveLength(length)
        guard let path = meanFreePathMetres else { return nil }
        let result = path/length
        guard result.isFinite else { throw AerodynamicError.nonFiniteResult }
        return result
    }
    /// Local body velocity v + ω×r in a spatially uniform medium. Offset is from
    /// this flow's velocity reference point; atmosphere gradients/winds are sampled by callers.
    public func atLocalPoint(offsetMetres: SIMD3<Double>, angularVelocityRadiansPerSecond: SIMD3<Double>) throws -> Self {
        guard aeroFinite(offsetMetres), aeroFinite(angularVelocityRadiansPerSecond) else { throw AerodynamicError.invalidFlow }
        return try Self(relativeVelocityMetresPerSecond: relativeVelocityMetresPerSecond
            + aeroCross(angularVelocityRadiansPerSecond,offsetMetres),
            densityKilogramsPerCubicMetre: densityKilogramsPerCubicMetre,pressurePascals: pressurePascals,
            temperatureKelvin: temperatureKelvin,speedOfSoundMetresPerSecond: speedOfSoundMetresPerSecond,
            dynamicViscosityPascalSeconds: dynamicViscosityPascalSeconds,meanFreePathMetres: meanFreePathMetres)
    }
}

internal func aeroNorm(_ v: SIMD3<Double>) -> Double { hypot(hypot(v.x,v.y),v.z) }
internal func aeroFinite(_ v: SIMD3<Double>) -> Bool { [v.x,v.y,v.z].allSatisfy(\.isFinite) }
internal func aeroCross(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> SIMD3<Double> {
    SIMD3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x)
}
private func aeroPositiveLength(_ length: Double) throws {
    guard length.isFinite, length > 0 else { throw AerodynamicError.invalidGeometry }
}
