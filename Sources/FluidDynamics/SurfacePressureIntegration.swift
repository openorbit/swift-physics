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

/// Flat surface quadrature element. Normal points outward; centre is the force
/// application point. All elements and flow vectors use the same body axes.
public struct FluidSurfacePanel: Sendable {
    public let areaSquareMetres: Double
    public let centreMetres: SIMD3<Double>
    public let outwardNormal: SIMD3<Double>
    /// Caller-supplied exposure; this library does not trace geometry occlusion.
    public let exposedToFlow: Bool
    public init(areaSquareMetres: Double, centreMetres: SIMD3<Double>,
                outwardNormal: SIMD3<Double>, exposedToFlow: Bool = true) throws {
        let length = fluidNorm(outwardNormal)
        guard areaSquareMetres.isFinite, areaSquareMetres > 0, fluidFinite(centreMetres),
              fluidFinite(outwardNormal), length.isFinite, length > 0 else { throw FluidError.invalidGeometry }
        self.areaSquareMetres = areaSquareMetres; self.centreMetres = centreMetres
        self.outwardNormal = outwardNormal/length; self.exposedToFlow = exposedToFlow
    }
}

/// Gauge pressure coefficient relative to ambient, not absolute pressure.
public protocol SurfacePressureCoefficientModel: Sendable {
    func pressureCoefficient(panel: FluidSurfacePanel, flow: FluidFlow,
                             referenceLengthMetres: Double) throws -> Double
}

public struct SurfacePressureResult: Sendable {
    public let pressureCoefficients: [Double]
    public let gaugePressuresPascals: [Double]
    public let loads: FluidLoads
}

/// O(N) independent-panel pressure integration. Ambient pressure is subtracted;
/// an open panel set represents aerodynamic gauge loads, not total structural pressure.
public struct SurfacePressureIntegration: Sendable {
    public let panels: [FluidSurfacePanel]
    public let referenceLengthMetres: Double
    public let model: any SurfacePressureCoefficientModel
    public init(panels: [FluidSurfacePanel], referenceLengthMetres: Double,
                model: any SurfacePressureCoefficientModel) throws {
        guard !panels.isEmpty, referenceLengthMetres.isFinite, referenceLengthMetres > 0 else { throw FluidError.invalidGeometry }
        self.panels = panels; self.referenceLengthMetres = referenceLengthMetres; self.model = model
    }
    public func evaluate(flow: FluidFlow, momentReferencePointMetres: SIMD3<Double> = .zero) throws -> SurfacePressureResult {
        try evaluate(momentReferencePointMetres: momentReferencePointMetres) { _ in flow }
    }
    /// Caller samples rotation, winds or gradients at each panel's centre.
    public func evaluate(momentReferencePointMetres: SIMD3<Double> = .zero,
                         flowAtPanel: (FluidSurfacePanel) throws -> FluidFlow) throws -> SurfacePressureResult {
        var loads = try FluidLoads(forceNewtons: .zero, momentNewtonMetres: .zero,
                                         referencePointMetres: momentReferencePointMetres)
        var coefficients = [Double](), pressures = [Double]()
        for panel in panels {
            let flow = try flowAtPanel(panel)
            let cp = try model.pressureCoefficient(panel: panel, flow: flow, referenceLengthMetres: referenceLengthMetres)
            let pressure = cp*flow.dynamicPressurePascals
            guard cp.isFinite, pressure.isFinite else { throw FluidError.nonFiniteResult }
            let local = try FluidLoads(forceNewtons: -panel.outwardNormal*(pressure*panel.areaSquareMetres),
                momentNewtonMetres: .zero, referencePointMetres: panel.centreMetres)
            loads = try loads.adding(local); coefficients.append(cp); pressures.append(pressure)
        }
        return SurfacePressureResult(pressureCoefficients: coefficients, gaugePressuresPascals: pressures, loads: loads)
    }
}
