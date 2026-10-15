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

/// State-local liquid properties, not a claim that one composition describes every sea.
/// Vapor pressure is absolute and must use the same temperature/composition as the sample.
public struct LiquidProperties: Sendable {
    public let identifier: String
    public let source: String
    public let densityKilogramsPerCubicMetre: Double
    public let dynamicViscosityPascalSeconds: Double
    public let speedOfSoundMetresPerSecond: Double
    public let vaporPressurePascals: Double
    public init(identifier: String, source: String, densityKilogramsPerCubicMetre: Double,
                dynamicViscosityPascalSeconds: Double, speedOfSoundMetresPerSecond: Double,
                vaporPressurePascals: Double) throws {
        guard !identifier.isEmpty, !source.isEmpty,
              [densityKilogramsPerCubicMetre,dynamicViscosityPascalSeconds,speedOfSoundMetresPerSecond]
                .allSatisfy({ $0.isFinite && $0 > 0 }), vaporPressurePascals.isFinite, vaporPressurePascals >= 0,
              (dynamicViscosityPascalSeconds/densityKilogramsPerCubicMetre).isFinite,
              dynamicViscosityPascalSeconds/densityKilogramsPerCubicMetre > 0 else {
            throw FluidError.invalidParameter
        }
        self.identifier = identifier; self.source = source; self.densityKilogramsPerCubicMetre = densityKilogramsPerCubicMetre
        self.dynamicViscosityPascalSeconds = dynamicViscosityPascalSeconds
        self.speedOfSoundMetresPerSecond = speedOfSoundMetresPerSecond; self.vaporPressurePascals = vaporPressurePascals
    }
    public var kinematicViscositySquareMetresPerSecond: Double { dynamicViscosityPascalSeconds/densityKilogramsPerCubicMetre }
    /// Diagnostic only: positive sigma does not establish a propeller cavitation margin.
    public func cavitationNumber(absolutePressurePascals: Double, speedMetresPerSecond: Double) throws -> Double {
        guard absolutePressurePascals.isFinite, absolutePressurePascals >= 0,
              speedMetresPerSecond.isFinite, speedMetresPerSecond > 0 else { throw FluidError.invalidParameter }
        let result = (absolutePressurePascals-vaporPressurePascals)/(0.5*densityKilogramsPerCubicMetre*speedMetresPerSecond*speedMetresPerSecond)
        guard result.isFinite else { throw FluidError.nonFiniteResult }; return result
    }
}
public protocol LiquidPropertyModel: Sendable {
    func properties(temperatureKelvin: Double, absolutePressurePascals: Double) throws -> LiquidProperties
}

/// Explicit constant-property approximation within declared temperature/pressure bounds.
/// No ocean composition, stratification, phase equilibrium or pressure EOS is inferred.
public struct ConstantLiquidModel: LiquidPropertyModel {
    public let value: LiquidProperties
    public let temperatureRangeKelvin: ClosedRange<Double>
    public let pressureRangePascals: ClosedRange<Double>
    public init(value: LiquidProperties, temperatureRangeKelvin: ClosedRange<Double>, pressureRangePascals: ClosedRange<Double>) throws {
        guard temperatureRangeKelvin.lowerBound.isFinite, temperatureRangeKelvin.lowerBound > 0,
              temperatureRangeKelvin.upperBound.isFinite, pressureRangePascals.lowerBound.isFinite,
              pressureRangePascals.lowerBound > value.vaporPressurePascals, pressureRangePascals.upperBound.isFinite else {
            throw FluidError.invalidParameter
        }
        self.value = value; self.temperatureRangeKelvin = temperatureRangeKelvin; self.pressureRangePascals = pressureRangePascals
    }
    public func properties(temperatureKelvin: Double, absolutePressurePascals: Double) throws -> LiquidProperties {
        guard temperatureKelvin.isFinite, absolutePressurePascals.isFinite,
              temperatureRangeKelvin.contains(temperatureKelvin), pressureRangePascals.contains(absolutePressurePascals) else {
            throw FluidError.outsideCoverage
        }
        return value
    }
}
public extension FluidFlow {
    init(relativeVelocityMetresPerSecond: SIMD3<Double>, liquid: LiquidProperties,
         temperatureKelvin: Double, absolutePressurePascals: Double) throws {
        guard absolutePressurePascals > liquid.vaporPressurePascals else { throw FluidError.outsideCoverage }
        try self.init(relativeVelocityMetresPerSecond: relativeVelocityMetresPerSecond,
            densityKilogramsPerCubicMetre: liquid.densityKilogramsPerCubicMetre,pressurePascals: absolutePressurePascals,
            temperatureKelvin: temperatureKelvin,speedOfSoundMetresPerSecond: liquid.speedOfSoundMetresPerSecond,
            dynamicViscosityPascalSeconds: liquid.dynamicViscosityPascalSeconds)
    }
}
