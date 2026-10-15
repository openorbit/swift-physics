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
import FluidDynamics
import Foundation
/// Idealized calorically perfect gas, windward modified-Newtonian pressure only.
/// Cp = CpMax max(0, normal·velocityDirection)². No shear, base suction, heating,
/// shock interaction or real-gas chemistry. Gas gamma is supplied, not Earth-specific.
public struct ModifiedNewtonianPressure: SurfacePressureCoefficientModel {
    public let specificHeatRatio: Double
    public let machRange: ClosedRange<Double>
    public let maximumKnudsenNumber: Double
    public init(specificHeatRatio: Double, machRange: ClosedRange<Double> = 5...25,
                maximumKnudsenNumber: Double = 0.01) throws {
        guard specificHeatRatio.isFinite, specificHeatRatio > 1, specificHeatRatio <= 2,
              machRange.lowerBound.isFinite, machRange.upperBound.isFinite,
              machRange.lowerBound >= 5, machRange.upperBound <= 50,
              maximumKnudsenNumber.isFinite, maximumKnudsenNumber > 0,
              maximumKnudsenNumber <= 0.01 else { throw FluidError.invalidParameter }
        self.specificHeatRatio = specificHeatRatio; self.machRange = machRange
        self.maximumKnudsenNumber = maximumKnudsenNumber
    }
    /// Normal shock followed by isentropic deceleration to stagnation.
    public func stagnationPressureCoefficient(machNumber m: Double) throws -> Double {
        guard m.isFinite, machRange.contains(m) else { throw FluidError.outsideCoverage }
        let g = specificHeatRatio, m2 = m*m
        let pressureRatio = (2*g*m2-(g-1))/(g+1)
        let downstreamMachSquared = ((g-1)*m2+2)/(2*g*m2-(g-1))
        let stagnationRatio = pressureRatio*exp((g/(g-1))*log1p((g-1)*downstreamMachSquared/2))
        let result = (stagnationRatio-1)/(0.5*g*m2)
        guard result.isFinite, result > 0 else { throw FluidError.nonFiniteResult }
        return result
    }
    public func pressureCoefficient(panel: FluidSurfacePanel, flow: FluidFlow,
                                    referenceLengthMetres: Double) throws -> Double {
        // Validate length even at rest, but rest needs no Mach/rarefaction inference.
        let kn = try flow.knudsenNumber(referenceLengthMetres: referenceLengthMetres)
        if flow.speedMetresPerSecond == 0 { return 0 }
        guard let mach = flow.machNumber, let kn else { throw FluidError.missingFlowProperty }
        guard kn <= maximumKnudsenNumber else { throw FluidError.outsideCoverage }
        let maximum = try stagnationPressureCoefficient(machNumber: mach)
        let incidence = max(0, min(1, (panel.outwardNormal*flow.relativeVelocityMetresPerSecond).sum()/flow.speedMetresPerSecond))
        return panel.exposedToFlow ? maximum*incidence*incidence : 0
    }
}
