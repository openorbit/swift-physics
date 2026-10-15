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

/// Pure local flow in caller-selected right-handed body/section axes:
/// +X forward, +Y right, +Z down. Velocity is body motion relative to the medium,
/// not the opposite fluid-on-body velocity. All vectors in a query must use the same axes.
/// Subtract atmosphere velocity from body velocity in a common frame before rotating
/// the difference into these axes. A co-rotating atmosphere is at rest in the planet's
/// rotating frame, but has velocity Ω×r in a planet-centred inertial frame.
public struct FluidFlow: Sendable {
    public let relativeVelocityMetresPerSecond: SIMD3<Double>
    public let densityKilogramsPerCubicMetre: Double
    public let pressurePascals: Double?
    public let temperatureKelvin: Double?
    public let speedOfSoundMetresPerSecond: Double?
    public let dynamicViscosityPascalSeconds: Double?
    public let meanFreePathMetres: Double?
    public var speedMetresPerSecond: Double { fluidNorm(relativeVelocityMetresPerSecond) }
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
        guard fluidFinite(relativeVelocityMetresPerSecond), densityKilogramsPerCubicMetre.isFinite,
              densityKilogramsPerCubicMetre >= 0,
              pressurePascals.map({ $0.isFinite && $0 >= 0 }) ?? true,
              [temperatureKelvin,speedOfSoundMetresPerSecond,dynamicViscosityPascalSeconds,meanFreePathMetres]
                .compactMap({ $0 }).allSatisfy({ $0.isFinite && $0 > 0 }) else { throw FluidError.invalidFlow }
        self.relativeVelocityMetresPerSecond = relativeVelocityMetresPerSecond
        self.densityKilogramsPerCubicMetre = densityKilogramsPerCubicMetre
        self.pressurePascals = pressurePascals; self.temperatureKelvin = temperatureKelvin
        self.speedOfSoundMetresPerSecond = speedOfSoundMetresPerSecond
        self.dynamicViscosityPascalSeconds = dynamicViscosityPascalSeconds; self.meanFreePathMetres = meanFreePathMetres
        guard speedMetresPerSecond.isFinite, dynamicPressurePascals.isFinite,
              machNumber.map(\.isFinite) ?? true else { throw FluidError.nonFiniteResult }
    }
    public func reynoldsNumber(referenceLengthMetres length: Double) throws -> Double? {
        try fluidPositiveLength(length)
        guard let viscosity = dynamicViscosityPascalSeconds else { return nil }
        let result = (densityKilogramsPerCubicMetre*speedMetresPerSecond)*length/viscosity
        guard result.isFinite else { throw FluidError.nonFiniteResult }
        return result
    }
    public func knudsenNumber(referenceLengthMetres length: Double) throws -> Double? {
        try fluidPositiveLength(length)
        guard let path = meanFreePathMetres else { return nil }
        let result = path/length
        guard result.isFinite else { throw FluidError.nonFiniteResult }
        return result
    }
    /// Local relative velocity v + (ωBody − ωMedium)×r for a rigidly rotating
    /// medium. Both angular velocities must be measured relative to the same frame
    /// and expressed in these axes. The default medium rate describes uniform air.
    /// Offset is from this flow's velocity reference point. Callers sample additional
    /// wind gradients separately. In planet-fixed coordinates, co-rotating air has
    /// zero medium rate; in inertial coordinates, supply the planet's rotation rate.
    public func atLocalPoint(offsetMetres: SIMD3<Double>, angularVelocityRadiansPerSecond: SIMD3<Double>,
                             mediumAngularVelocityRadiansPerSecond: SIMD3<Double> = .zero) throws -> Self {
        guard fluidFinite(offsetMetres), fluidFinite(angularVelocityRadiansPerSecond),
              fluidFinite(mediumAngularVelocityRadiansPerSecond) else { throw FluidError.invalidFlow }
        return try Self(relativeVelocityMetresPerSecond: relativeVelocityMetresPerSecond
            + fluidCross(angularVelocityRadiansPerSecond-mediumAngularVelocityRadiansPerSecond,offsetMetres),
            densityKilogramsPerCubicMetre: densityKilogramsPerCubicMetre,pressurePascals: pressurePascals,
            temperatureKelvin: temperatureKelvin,speedOfSoundMetresPerSecond: speedOfSoundMetresPerSecond,
            dynamicViscosityPascalSeconds: dynamicViscosityPascalSeconds,meanFreePathMetres: meanFreePathMetres)
    }
}

package func fluidNorm(_ v: SIMD3<Double>) -> Double { hypot(hypot(v.x,v.y),v.z) }
package func fluidFinite(_ v: SIMD3<Double>) -> Bool { [v.x,v.y,v.z].allSatisfy(\.isFinite) }
package func fluidCross(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> SIMD3<Double> {
    SIMD3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x)
}
private func fluidPositiveLength(_ length: Double) throws {
    guard length.isFinite, length > 0 else { throw FluidError.invalidGeometry }
}

public extension FluidFlow {
    /// Preserve thermodynamic properties while applying caller-supplied wind/wake
    /// or other physical relative-velocity corrections in the same axes and frame.
    func replacingRelativeVelocity(_ velocityMetresPerSecond: SIMD3<Double>) throws -> Self {
        try Self(relativeVelocityMetresPerSecond: velocityMetresPerSecond,
            densityKilogramsPerCubicMetre: densityKilogramsPerCubicMetre,pressurePascals: pressurePascals,
            temperatureKelvin: temperatureKelvin,speedOfSoundMetresPerSecond: speedOfSoundMetresPerSecond,
            dynamicViscosityPascalSeconds: dynamicViscosityPascalSeconds,meanFreePathMetres: meanFreePathMetres)
    }
}
