// SPDX-License-Identifier: Apache-2.0
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
