// SPDX-License-Identifier: Apache-2.0
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
              aeroFinite(centreOfBuoyancyMetres) else { throw AerodynamicError.invalidGeometry }
        self.displacedVolumeCubicMetres = displacedVolumeCubicMetres; self.centreOfBuoyancyMetres = centreOfBuoyancyMetres
    }
    public func evaluate(densityKilogramsPerCubicMetre density: Double, field: BuoyancyField,
                         momentReferencePointMetres: SIMD3<Double> = .zero) throws -> AerodynamicLoads {
        guard density.isFinite, density >= 0 else { throw AerodynamicError.invalidFlow }
        let force: SIMD3<Double>
        switch field {
        case .hydrostaticGravityMetresPerSecondSquared(let gravity):
            guard aeroFinite(gravity) else { throw AerodynamicError.invalidFlow }
            force = -gravity*(density*displacedVolumeCubicMetres)
        case .pressureGradientPascalsPerMetre(let gradient):
            guard aeroFinite(gradient) else { throw AerodynamicError.invalidFlow }
            force = -gradient*displacedVolumeCubicMetres
        }
        return try AerodynamicLoads(forceNewtons: force, momentNewtonMetres: .zero,
            referencePointMetres: centreOfBuoyancyMetres).referred(to: momentReferencePointMetres)
    }
}

/// Combine displaced-volume buoyancy and optional supplied-shape drag exactly once.
/// Weight of envelope, lifting gas and payload is handled by external mass/gravity.
public struct BalloonAerodynamics: Sendable {
    public let buoyancy: Buoyancy
    public let drag: DragAreaAerodynamics?
    public init(buoyancy: Buoyancy, drag: DragAreaAerodynamics? = nil) {
        self.buoyancy = buoyancy; self.drag = drag
    }
    public func evaluate(flow: AerodynamicFlow, field: BuoyancyField,
                         momentReferencePointMetres: SIMD3<Double> = .zero) throws -> AerodynamicLoads {
        let pressureLoads = try buoyancy.evaluate(densityKilogramsPerCubicMetre: flow.densityKilogramsPerCubicMetre,
            field: field, momentReferencePointMetres: momentReferencePointMetres)
        guard let drag else { return pressureLoads }
        return try pressureLoads.adding(drag.evaluate(flow: flow, momentReferencePointMetres: momentReferencePointMetres))
    }
}
