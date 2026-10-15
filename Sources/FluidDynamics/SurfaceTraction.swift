// SPDX-License-Identifier: Apache-2.0
import Foundation

/// Signed pressure above ambient and tangential stress in the panel's body axes, Pa.
/// Force per unit area is −pressure*outwardNormal + shear. Negative pressure is suction.
public struct SurfaceTraction: Sendable {
    public let gaugePressurePascals: Double
    public let shearPascals: SIMD3<Double>
    public init(gaugePressurePascals: Double, shearPascals: SIMD3<Double>) throws {
        guard gaugePressurePascals.isFinite, fluidFinite(shearPascals) else { throw FluidError.nonFiniteResult }
        self.gaugePressurePascals = gaugePressurePascals; self.shearPascals = shearPascals
    }
}
public protocol SurfaceTractionModel: Sendable {
    func traction(panel: FluidSurfacePanel, flow: FluidFlow,
                  referenceLengthMetres: Double) throws -> SurfaceTraction
}
public struct SurfaceTractionResult: Sendable {
    public let tractions: [SurfaceTraction]
    public let loads: FluidLoads
}

/// Pressure and shear quadrature; no additional pressure-only load is implied.
public struct SurfaceTractionIntegration: Sendable {
    public let panels: [FluidSurfacePanel]
    public let referenceLengthMetres: Double
    public let model: any SurfaceTractionModel
    public init(panels: [FluidSurfacePanel], referenceLengthMetres: Double,
                model: any SurfaceTractionModel) throws {
        guard !panels.isEmpty, referenceLengthMetres.isFinite, referenceLengthMetres > 0 else { throw FluidError.invalidGeometry }
        self.panels = panels; self.referenceLengthMetres = referenceLengthMetres; self.model = model
    }
    public func evaluate(flow: FluidFlow, momentReferencePointMetres: SIMD3<Double> = .zero) throws -> SurfaceTractionResult {
        try evaluate(momentReferencePointMetres: momentReferencePointMetres) { _ in flow }
    }
    public func evaluate(momentReferencePointMetres: SIMD3<Double> = .zero,
                         flowAtPanel: (FluidSurfacePanel) throws -> FluidFlow) throws -> SurfaceTractionResult {
        var total = try FluidLoads(forceNewtons: .zero,momentNewtonMetres: .zero,referencePointMetres: momentReferencePointMetres)
        var tractions = [SurfaceTraction]()
        for panel in panels {
            let traction = try model.traction(panel: panel,flow: flowAtPanel(panel),referenceLengthMetres: referenceLengthMetres)
            let scale = max(abs(traction.shearPascals.x),abs(traction.shearPascals.y),abs(traction.shearPascals.z))
            if scale > 0 {
                guard abs(((traction.shearPascals/scale)*panel.outwardNormal).sum()) <= 1e-10 else {
                    throw FluidError.invalidCoefficients
                }
            }
            let force = (-panel.outwardNormal*traction.gaugePressurePascals+traction.shearPascals)*panel.areaSquareMetres
            total = try total.adding(FluidLoads(forceNewtons: force,momentNewtonMetres: .zero,
                referencePointMetres: panel.centreMetres))
            tractions.append(traction)
        }
        return SurfaceTractionResult(tractions: tractions,loads: total)
    }
}

/// Adapt a continuum pressure model with explicitly zero shear. This is an alternative
/// load evaluation, not a second force to add to the original pressure evaluator.
public struct PressureOnlyTraction: SurfaceTractionModel {
    public let pressureModel: any SurfacePressureCoefficientModel
    public init(_ pressureModel: any SurfacePressureCoefficientModel) { self.pressureModel = pressureModel }
    public func traction(panel: FluidSurfacePanel, flow: FluidFlow,
                         referenceLengthMetres: Double) throws -> SurfaceTraction {
        let coefficient = try pressureModel.pressureCoefficient(panel: panel,flow: flow,referenceLengthMetres: referenceLengthMetres)
        return try SurfaceTraction(gaugePressurePascals: coefficient*flow.dynamicPressurePascals,shearPascals: .zero)
    }
}
