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

public extension IncompressibleFlowConfiguration {
    /// Fully filled single-phase liquid domain. No free-surface tracking is implied.
    init(grid: CartesianFluidGrid, liquid: LiquidProperties, boundaryVelocityMetresPerSecond: SIMD3<Double> = .zero) throws {
        try self.init(dimensions: grid.dimensions,minimumCornerMetres: grid.minimumCornerMetres,cellSizeMetres: grid.cellSizeMetres,
            densityKilogramsPerCubicMetre: liquid.densityKilogramsPerCubicMetre,
            speedOfSoundMetresPerSecond: liquid.speedOfSoundMetresPerSecond,
            kinematicViscositySquareMetresPerSecond: liquid.kinematicViscositySquareMetresPerSecond,
            boundaryVelocityMetresPerSecond: boundaryVelocityMetresPerSecond)
    }
}
/// Calibrates the arbitrary incompressible pressure gauge at a known fluid point.
/// Include hydrostatic pressure in that reference consistently with applied gravity.
public struct LiquidCFDPressureReference: Sendable {
    public let positionMetres: SIMD3<Double>
    public let absolutePressurePascals: Double
    public init(positionMetres: SIMD3<Double>, absolutePressurePascals: Double) throws {
        guard fluidFinite(positionMetres), absolutePressurePascals.isFinite, absolutePressurePascals >= 0 else { throw FluidError.invalidParameter }
        self.positionMetres = positionMetres; self.absolutePressurePascals = absolutePressurePascals
    }
    public func absolutePressurePascals(in grid: IncompressibleFlowGrid, at position: SIMD3<Double>) throws -> Double {
        let pressure = try absolutePressurePascals+grid.gaugePressurePascals(at: position)-grid.gaugePressurePascals(at: positionMetres)
        guard pressure.isFinite, pressure >= 0 else { throw FluidError.outsideCoverage }; return pressure
    }
}
