// SPDX-License-Identifier: Apache-2.0
import FluidDynamics
import Foundation
/// Current-canopy steady drag. No opening loads, tether transmission or deployment.
public struct ParachuteAerodynamics: Sendable {
    public let drag: DragAreaForces
    public init(drag: DragAreaForces) { self.drag = drag }
    public func evaluate(flow: FluidFlow, momentReferencePointMetres: SIMD3<Double> = .zero) throws -> FluidLoads {
        try drag.evaluate(flow: flow, momentReferencePointMetres: momentReferencePointMetres)
    }
}
