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
/// Combine displaced-volume buoyancy and optional supplied-shape drag exactly once.
/// Weight of envelope, lifting gas and payload is handled by external mass/gravity.
public struct BalloonAerodynamics: Sendable {
    public let buoyancy: Buoyancy
    public let drag: DragAreaForces?
    public init(buoyancy: Buoyancy, drag: DragAreaForces? = nil) {
        self.buoyancy = buoyancy; self.drag = drag
    }
    public func evaluate(flow: FluidFlow, field: BuoyancyField,
                         momentReferencePointMetres: SIMD3<Double> = .zero) throws -> FluidLoads {
        let pressureLoads = try buoyancy.evaluate(densityKilogramsPerCubicMetre: flow.densityKilogramsPerCubicMetre,
            field: field, momentReferencePointMetres: momentReferencePointMetres)
        guard let drag else { return pressureLoads }
        return try pressureLoads.adding(drag.evaluate(flow: flow, momentReferencePointMetres: momentReferencePointMetres))
    }
}
