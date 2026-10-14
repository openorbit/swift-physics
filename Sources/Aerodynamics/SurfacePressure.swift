// SPDX-License-Identifier: Apache-2.0
import Foundation

/// Flat surface quadrature element. Normal points outward; centre is the force
/// application point. All elements and flow vectors use the same body axes.
public struct AerodynamicSurfacePanel: Sendable {
    public let areaSquareMetres: Double
    public let centreMetres: SIMD3<Double>
    public let outwardNormal: SIMD3<Double>
    /// Caller-supplied exposure; this library does not trace geometry occlusion.
    public let exposedToFlow: Bool
    public init(areaSquareMetres: Double, centreMetres: SIMD3<Double>,
                outwardNormal: SIMD3<Double>, exposedToFlow: Bool = true) throws {
        let length = aeroNorm(outwardNormal)
        guard areaSquareMetres.isFinite, areaSquareMetres > 0, aeroFinite(centreMetres),
              aeroFinite(outwardNormal), length.isFinite, length > 0 else { throw AerodynamicError.invalidGeometry }
        self.areaSquareMetres = areaSquareMetres; self.centreMetres = centreMetres
        self.outwardNormal = outwardNormal/length; self.exposedToFlow = exposedToFlow
    }
}

/// Gauge pressure coefficient relative to ambient, not absolute pressure.
public protocol SurfacePressureCoefficientModel: Sendable {
    func pressureCoefficient(panel: AerodynamicSurfacePanel, flow: AerodynamicFlow,
                             referenceLengthMetres: Double) throws -> Double
}

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
              maximumKnudsenNumber <= 0.01 else { throw AerodynamicError.invalidParameter }
        self.specificHeatRatio = specificHeatRatio; self.machRange = machRange
        self.maximumKnudsenNumber = maximumKnudsenNumber
    }
    /// Normal shock followed by isentropic deceleration to stagnation.
    public func stagnationPressureCoefficient(machNumber m: Double) throws -> Double {
        guard m.isFinite, machRange.contains(m) else { throw AerodynamicError.outsideCoverage }
        let g = specificHeatRatio, m2 = m*m
        let pressureRatio = (2*g*m2-(g-1))/(g+1)
        let downstreamMachSquared = ((g-1)*m2+2)/(2*g*m2-(g-1))
        let stagnationRatio = pressureRatio*exp((g/(g-1))*log1p((g-1)*downstreamMachSquared/2))
        let result = (stagnationRatio-1)/(0.5*g*m2)
        guard result.isFinite, result > 0 else { throw AerodynamicError.nonFiniteResult }
        return result
    }
    public func pressureCoefficient(panel: AerodynamicSurfacePanel, flow: AerodynamicFlow,
                                    referenceLengthMetres: Double) throws -> Double {
        // Validate length even at rest, but rest needs no Mach/rarefaction inference.
        let kn = try flow.knudsenNumber(referenceLengthMetres: referenceLengthMetres)
        if flow.speedMetresPerSecond == 0 { return 0 }
        guard let mach = flow.machNumber, let kn else { throw AerodynamicError.missingFlowProperty }
        guard kn <= maximumKnudsenNumber else { throw AerodynamicError.outsideCoverage }
        let maximum = try stagnationPressureCoefficient(machNumber: mach)
        let incidence = max(0, min(1, (panel.outwardNormal*flow.relativeVelocityMetresPerSecond).sum()/flow.speedMetresPerSecond))
        return panel.exposedToFlow ? maximum*incidence*incidence : 0
    }
}

public struct SurfacePressureResult: Sendable {
    public let pressureCoefficients: [Double]
    public let gaugePressuresPascals: [Double]
    public let loads: AerodynamicLoads
}

/// O(N) independent-panel pressure integration. Ambient pressure is subtracted;
/// an open panel set represents aerodynamic gauge loads, not total structural pressure.
public struct SurfacePressureAerodynamics: Sendable {
    public let panels: [AerodynamicSurfacePanel]
    public let referenceLengthMetres: Double
    public let model: any SurfacePressureCoefficientModel
    public init(panels: [AerodynamicSurfacePanel], referenceLengthMetres: Double,
                model: any SurfacePressureCoefficientModel) throws {
        guard !panels.isEmpty, referenceLengthMetres.isFinite, referenceLengthMetres > 0 else { throw AerodynamicError.invalidGeometry }
        self.panels = panels; self.referenceLengthMetres = referenceLengthMetres; self.model = model
    }
    public func evaluate(flow: AerodynamicFlow, momentReferencePointMetres: SIMD3<Double> = .zero) throws -> SurfacePressureResult {
        try evaluate(momentReferencePointMetres: momentReferencePointMetres) { _ in flow }
    }
    /// Caller samples rotation, winds or gradients at each panel's centre.
    public func evaluate(momentReferencePointMetres: SIMD3<Double> = .zero,
                         flowAtPanel: (AerodynamicSurfacePanel) throws -> AerodynamicFlow) throws -> SurfacePressureResult {
        var loads = try AerodynamicLoads(forceNewtons: .zero, momentNewtonMetres: .zero,
                                         referencePointMetres: momentReferencePointMetres)
        var coefficients = [Double](), pressures = [Double]()
        for panel in panels {
            let flow = try flowAtPanel(panel)
            let cp = try model.pressureCoefficient(panel: panel, flow: flow, referenceLengthMetres: referenceLengthMetres)
            let pressure = cp*flow.dynamicPressurePascals
            guard cp.isFinite, pressure.isFinite else { throw AerodynamicError.nonFiniteResult }
            let local = try AerodynamicLoads(forceNewtons: -panel.outwardNormal*(pressure*panel.areaSquareMetres),
                momentNewtonMetres: .zero, referencePointMetres: panel.centreMetres)
            loads = try loads.adding(local); coefficients.append(cp); pressures.append(pressure)
        }
        return SurfacePressureResult(pressureCoefficients: coefficients, gaugePressuresPascals: pressures, loads: loads)
    }
}
