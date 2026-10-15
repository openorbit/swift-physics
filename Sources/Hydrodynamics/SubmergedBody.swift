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
import FluidDynamics

/// Passive body-axis damping about one supplied velocity/moment reference point.
/// Coefficients are dimensional supplied data, not universal hull/fin coefficients.
public struct HydrodynamicDamping: Sendable {
    public let linearForce: SIMD3<Double>
    public let quadraticForce: SIMD3<Double>
    public let linearMoment: SIMD3<Double>
    public let quadraticMoment: SIMD3<Double>
    public let referencePointMetres: SIMD3<Double>
    public let source: String
    public init(linearForce: SIMD3<Double> = .zero, quadraticForce: SIMD3<Double> = .zero,
                linearMoment: SIMD3<Double> = .zero, quadraticMoment: SIMD3<Double> = .zero,
                referencePointMetres: SIMD3<Double> = .zero, source: String) throws {
        guard [linearForce,quadraticForce,linearMoment,quadraticMoment].allSatisfy({ fluidFinite($0) && $0.x >= 0 && $0.y >= 0 && $0.z >= 0 }),
              fluidFinite(referencePointMetres), !source.isEmpty else { throw FluidError.invalidParameter }
        self.linearForce = linearForce; self.quadraticForce = quadraticForce
        self.linearMoment = linearMoment; self.quadraticMoment = quadraticMoment
        self.referencePointMetres = referencePointMetres; self.source = source
    }
    /// Relative rates use the same axes and frame convention as the supplied data.
    public func evaluate(relativeVelocityMetresPerSecond v: SIMD3<Double>, relativeAngularVelocityRadiansPerSecond w: SIMD3<Double>,
                         momentReferencePointMetres: SIMD3<Double> = .zero) throws -> FluidLoads {
        guard fluidFinite(v), fluidFinite(w) else { throw FluidError.invalidFlow }
        let av = SIMD3(abs(v.x),abs(v.y),abs(v.z)), aw = SIMD3(abs(w.x),abs(w.y),abs(w.z))
        return try FluidLoads(forceNewtons: -(linearForce+quadraticForce*av)*v,
            momentNewtonMetres: -(linearMoment+quadraticMoment*aw)*w,referencePointMetres: referencePointMetres)
            .referred(to: momentReferencePointMetres)
    }
}

/// Diagonal added inertia for a fixed-axis small-motion approximation. General
/// rotating-frame/translation-rotation coupling requires a full marine dynamics model.
/// Incorporate in effective inertia OR use the acceleration reaction, never both.
public struct HydrodynamicAddedInertia: Sendable {
    public let massKilograms: SIMD3<Double>
    public let momentKilogramSquareMetres: SIMD3<Double>
    public init(massKilograms: SIMD3<Double>, momentKilogramSquareMetres: SIMD3<Double>) throws {
        guard [massKilograms,momentKilogramSquareMetres].allSatisfy({ fluidFinite($0) && $0.x >= 0 && $0.y >= 0 && $0.z >= 0 }) else {
            throw FluidError.invalidParameter
        }
        self.massKilograms = massKilograms; self.momentKilogramSquareMetres = momentKilogramSquareMetres
    }
    public func accelerationReaction(relativeLinearAccelerationMetresPerSecondSquared a: SIMD3<Double>,
                                     relativeAngularAccelerationRadiansPerSecondSquared alpha: SIMD3<Double>,
                                     referencePointMetres: SIMD3<Double> = .zero) throws -> FluidLoads {
        guard fluidFinite(a), fluidFinite(alpha) else { throw FluidError.invalidFlow }
        return try FluidLoads(forceNewtons: -massKilograms*a,momentNewtonMetres: -momentKilogramSquareMetres*alpha,
            referencePointMetres: referencePointMetres)
    }
}

/// Fully wetted supplied-state body loads. Buoyancy and damping are evaluated once;
/// current velocity, displaced geometry and rigid-body motion remain caller-owned.
public struct SubmergedBodyHydrodynamics: Sendable {
    public let buoyancy: Buoyancy
    public let damping: HydrodynamicDamping
    public init(buoyancy: Buoyancy, damping: HydrodynamicDamping) { self.buoyancy = buoyancy; self.damping = damping }
    public func evaluate(liquid: LiquidProperties, flow: FluidFlow, relativeAngularVelocityRadiansPerSecond: SIMD3<Double> = .zero,
                         gravityMetresPerSecondSquared: SIMD3<Double>, momentReferencePointMetres: SIMD3<Double> = .zero) throws -> FluidLoads {
        guard flow.densityKilogramsPerCubicMetre == liquid.densityKilogramsPerCubicMetre else { throw FluidError.unsupportedQuery }
        let pressure = try buoyancy.evaluate(densityKilogramsPerCubicMetre: liquid.densityKilogramsPerCubicMetre,
            field: .hydrostaticGravityMetresPerSecondSquared(gravityMetresPerSecondSquared),momentReferencePointMetres: momentReferencePointMetres)
        let resistance = try damping.evaluate(relativeVelocityMetresPerSecond: flow.relativeVelocityMetresPerSecond,
            relativeAngularVelocityRadiansPerSecond: relativeAngularVelocityRadiansPerSecond,momentReferencePointMetres: momentReferencePointMetres)
        return try pressure.adding(resistance)
    }
}
