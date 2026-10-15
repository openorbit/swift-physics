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

/// Two alternative descriptions of the same ambient pressure force. All vectors
/// use the load axes. Gravity must be physical gravity, not vehicle acceleration.
public enum BuoyancyField: Sendable {
    case hydrostaticGravityMetresPerSecondSquared(SIMD3<Double>)
    case pressureGradientPascalsPerMetre(SIMD3<Double>)
}

/// Local uniform-density/gradient displaced-volume pressure force. Volume and centre
/// of buoyancy are supplied current geometry; no gas mass or envelope evolution.
public struct Buoyancy: Sendable {
    public let displacedVolumeCubicMetres: Double
    public let centreOfBuoyancyMetres: SIMD3<Double>
    public init(displacedVolumeCubicMetres: Double, centreOfBuoyancyMetres: SIMD3<Double> = .zero) throws {
        guard displacedVolumeCubicMetres.isFinite, displacedVolumeCubicMetres >= 0,
              fluidFinite(centreOfBuoyancyMetres) else { throw FluidError.invalidGeometry }
        self.displacedVolumeCubicMetres = displacedVolumeCubicMetres; self.centreOfBuoyancyMetres = centreOfBuoyancyMetres
    }
    public func evaluate(densityKilogramsPerCubicMetre density: Double, field: BuoyancyField,
                         momentReferencePointMetres: SIMD3<Double> = .zero) throws -> FluidLoads {
        guard density.isFinite, density >= 0 else { throw FluidError.invalidFlow }
        let force: SIMD3<Double>
        switch field {
        case .hydrostaticGravityMetresPerSecondSquared(let gravity):
            guard fluidFinite(gravity) else { throw FluidError.invalidFlow }
            force = -gravity*(density*displacedVolumeCubicMetres)
        case .pressureGradientPascalsPerMetre(let gradient):
            guard fluidFinite(gradient) else { throw FluidError.invalidFlow }
            force = -gradient*displacedVolumeCubicMetres
        }
        return try FluidLoads(forceNewtons: force, momentNewtonMetres: .zero,
            referencePointMetres: centreOfBuoyancyMetres).referred(to: momentReferencePointMetres)
    }
}

