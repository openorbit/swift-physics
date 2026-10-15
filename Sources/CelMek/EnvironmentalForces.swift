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

public enum EnvironmentalForceError: Error, Equatable {
    case invalidParameter, invalidDensity, invalidVelocity, invalidIllumination
    case originMismatch, missingEphemerides, singularSource, numericalRangeExceeded
}

/// Local medium properties in the same axes and physical origin as the state.
public struct AtmosphereSample {
    public let epoch: Epoch
    public let frame: Frame
    public let origin: BodyIdentifier
    public let densityKilogramsPerCubicMetre: Double
    public let velocityMetresPerSecond: SIMD3<Double>
    public init(epoch: Epoch, frame: Frame, origin: BodyIdentifier,
                densityKilogramsPerCubicMetre: Double, velocityMetresPerSecond: SIMD3<Double>) throws {
        try origin.validate()
        guard densityKilogramsPerCubicMetre.isFinite, densityKilogramsPerCubicMetre >= 0 else {
            throw EnvironmentalForceError.invalidDensity
        }
        guard [velocityMetresPerSecond.x,velocityMetresPerSecond.y,velocityMetresPerSecond.z].allSatisfy(\.isFinite) else {
            throw EnvironmentalForceError.invalidVelocity
        }
        self.epoch = epoch; self.frame = frame; self.origin = origin
        self.densityKilogramsPerCubicMetre = densityKilogramsPerCubicMetre
        self.velocityMetresPerSecond = velocityMetresPerSecond
    }
}

/// Pure local sampling; implementations supply winds/rotation, coverage and density.
public protocol AtmosphereProvider {
    func sample(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> AtmosphereSample
}

/// Explicit uniform medium for controlled scenarios; not a planetary atmosphere model.
public struct UniformAtmosphere: AtmosphereProvider {
    private let reference: AtmosphereSample
    public init(frame: Frame, origin: BodyIdentifier, densityKilogramsPerCubicMetre: Double,
                velocityMetresPerSecond: SIMD3<Double> = .zero) throws {
        reference = try AtmosphereSample(epoch: .j2000, frame: frame, origin: origin,
            densityKilogramsPerCubicMetre: densityKilogramsPerCubicMetre, velocityMetresPerSecond: velocityMetresPerSecond)
    }
    public func sample(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> AtmosphereSample {
        try AtmosphereSample(epoch: state.epoch, frame: reference.frame, origin: reference.origin,
            densityKilogramsPerCubicMetre: reference.densityKilogramsPerCubicMetre,
            velocityMetresPerSecond: reference.velocityMetresPerSecond)
    }
}

/// Isotropic constant-area drag: -rho Cd A |v-relative| v-relative / (2 mass).
public struct AtmosphericDrag: ForceModel {
    public let referenceAreaSquareMetres: Double
    public let dragCoefficient: Double
    public let atmosphere: any AtmosphereProvider
    public init(referenceAreaSquareMetres: Double, dragCoefficient: Double, atmosphere: any AtmosphereProvider) throws {
        try validateEnvironmentParameters([referenceAreaSquareMetres,dragCoefficient])
        self.referenceAreaSquareMetres = referenceAreaSquareMetres
        self.dragCoefficient = dragCoefficient; self.atmosphere = atmosphere
    }
    public func contribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
        try validateEnvironmentState(state, context)
        let medium = try atmosphere.sample(for: state, in: context)
        guard medium.frame === state.frame else { throw SpacecraftDynamicsError.frameMismatch }
        guard medium.epoch == state.epoch else { throw SpacecraftDynamicsError.epochMismatch }
        guard medium.origin == context.origin else { throw EnvironmentalForceError.originMismatch }
        let relative = state.cartesian.velocity - medium.velocityMetresPerSecond
        let speed = hypot(hypot(relative.x,relative.y),relative.z)
        guard speed.isFinite else { throw EnvironmentalForceError.numericalRangeExceeded }
        var acceleration = SIMD3<Double>.zero
        if speed > 0 && medium.densityKilogramsPerCubicMetre > 0 && dragCoefficient > 0 && referenceAreaSquareMetres > 0 {
            let magnitude = 0.5 * medium.densityKilogramsPerCubicMetre * dragCoefficient
                * referenceAreaSquareMetres / state.mass.totalMassKilograms * speed * speed
            acceleration = -(relative / speed) * magnitude
        }
        return try environmentalContribution(acceleration, state)
    }
}

/// Fraction of unobscured incident radiation, in 0...1. Caller models shadows.
public protocol RadiationIlluminationProvider {
    func fraction(for state: SpacecraftState, sourcePosition: SIMD3<Double>, in context: ForceEvaluationContext) throws -> Double
}
public struct ConstantRadiationIllumination: RadiationIlluminationProvider {
    public let illuminatedFraction: Double
    public init(illuminatedFraction: Double) throws {
        guard illuminatedFraction.isFinite, (0...1).contains(illuminatedFraction) else {
            throw EnvironmentalForceError.invalidIllumination
        }
        self.illuminatedFraction = illuminatedFraction
    }
    public func fraction(for state: SpacecraftState, sourcePosition: SIMD3<Double>, in context: ForceEvaluationContext) throws -> Double {
        illuminatedFraction
    }
}

/// Isotropic constant-area radiation pressure directed away from a point source.
/// Pressure scales by inverse-square source distance. No implicit eclipse model.
public struct SolarRadiationPressure: ForceModel {
    public let source: BodyIdentifier
    public let referencePressurePascals: SourcedPhysicalValue
    public let referenceDistanceMetres: SourcedPhysicalValue
    public let referenceAreaSquareMetres: Double
    public let radiationCoefficient: Double
    public let illumination: any RadiationIlluminationProvider
    public init(source: BodyIdentifier = .sun, referencePressurePascals: SourcedPhysicalValue,
                referenceDistanceMetres: SourcedPhysicalValue, referenceAreaSquareMetres: Double,
                radiationCoefficient: Double, illumination: any RadiationIlluminationProvider) throws {
        try source.validate()
        try validateEnvironmentParameters([referenceAreaSquareMetres,radiationCoefficient])
        self.source = source; self.referencePressurePascals = referencePressurePascals
        self.referenceDistanceMetres = referenceDistanceMetres
        self.referenceAreaSquareMetres = referenceAreaSquareMetres; self.radiationCoefficient = radiationCoefficient
        self.illumination = illumination
    }
    public func contribution(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> ForceContribution {
        try validateEnvironmentState(state, context)
        let sourcePosition: SIMD3<Double>
        if source == context.origin {
            sourcePosition = .zero
        } else {
            guard let collection = context.ephemerides else { throw EnvironmentalForceError.missingEphemerides }
            let provider = try collection.provider(for: source)
            guard provider.metadata.origin == context.origin else { throw EnvironmentalForceError.originMismatch }
            guard provider.frame === state.frame else { throw SpacecraftDynamicsError.frameMismatch }
            sourcePosition = try provider.checkedState(at: state.epoch).position
        }
        let displacement = state.cartesian.position - sourcePosition
        let distance = hypot(hypot(displacement.x,displacement.y),displacement.z)
        guard distance.isFinite else { throw EnvironmentalForceError.numericalRangeExceeded }
        guard distance > 0 else { throw EnvironmentalForceError.singularSource }
        let fraction = try illumination.fraction(for: state, sourcePosition: sourcePosition, in: context)
        guard fraction.isFinite, (0...1).contains(fraction) else { throw EnvironmentalForceError.invalidIllumination }
        var acceleration = SIMD3<Double>.zero
        if fraction > 0 && radiationCoefficient > 0 && referenceAreaSquareMetres > 0 {
            let ratio = referenceDistanceMetres.value / distance
            let magnitude = referencePressurePascals.value * fraction * radiationCoefficient
                * referenceAreaSquareMetres / state.mass.totalMassKilograms * ratio * ratio
            acceleration = displacement / distance * magnitude
        }
        return try environmentalContribution(acceleration, state)
    }
}

private func validateEnvironmentParameters(_ values: [Double]) throws {
    guard values.allSatisfy({ $0.isFinite && $0 >= 0 }) else { throw EnvironmentalForceError.invalidParameter }
}
private func validateEnvironmentState(_ state: SpacecraftState, _ context: ForceEvaluationContext) throws {
    guard state.frame === context.frame else { throw SpacecraftDynamicsError.frameMismatch }
    guard state.frame.pseudoInertial else { throw SpacecraftDynamicsError.nonInertialFrame }
}
private func environmentalContribution(_ acceleration: SIMD3<Double>, _ state: SpacecraftState) throws -> ForceContribution {
    guard [acceleration.x,acceleration.y,acceleration.z].allSatisfy(\.isFinite) else {
        throw EnvironmentalForceError.numericalRangeExceeded
    }
    return try ForceContribution(epoch: state.epoch, frame: state.frame, accelerationMetresPerSecondSquared: acceleration)
}
