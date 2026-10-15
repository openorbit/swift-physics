// SPDX-License-Identifier: Apache-2.0
import Foundation

/// Supplied continuum/slip/transitional/free-molecular pressure and shear data.
/// Axes in order: log10(Kn), outward-normal incidence cosine, speed/sqrt(RT), Tw/T.
/// Pressure is gauge pressure/(rho RT). Shear magnitude/(rho RT) opposes tangential
/// body-relative velocity. Values must come from a consistent gas/surface dataset;
/// the table does not extrapolate either asymptotic evaluator into transitional flow.
public struct TabulatedRarefiedTraction: SurfaceTractionModel {
    public let gas: FreeMolecularGas
    public let wallTemperatureKelvin: Double
    public let pressure: AerodynamicScalarTable
    public let shear: AerodynamicScalarTable
    public init(gas: FreeMolecularGas, wallTemperatureKelvin: Double,
                pressure: AerodynamicScalarTable, shear: AerodynamicScalarTable) throws {
        guard wallTemperatureKelvin.isFinite, wallTemperatureKelvin > 0,
              pressure.axes == shear.axes, pressure.valueUnits == "1", shear.valueUnits == "1",
              pressure.axes.map(\.name) == ["log10Kn","incidenceCosine","speedRatio","wallTemperatureRatio"],
              pressure.axes.allSatisfy({ $0.units == "1" }),
              pressure.axes[1].coordinates.first! >= -1, pressure.axes[1].coordinates.last! <= 1,
              pressure.axes[2].coordinates.first! >= 0, pressure.axes[3].coordinates.first! > 0,
              shear.values.allSatisfy({ $0 >= 0 }) else { throw AerodynamicError.malformedData }
        self.gas = gas; self.wallTemperatureKelvin = wallTemperatureKelvin
        self.pressure = pressure; self.shear = shear
    }
    public func traction(panel: AerodynamicSurfacePanel, flow: AerodynamicFlow,
                         referenceLengthMetres: Double) throws -> SurfaceTraction {
        guard let temperature = flow.temperatureKelvin,
              let kn = try flow.knudsenNumber(referenceLengthMetres: referenceLengthMetres) else {
            throw AerodynamicError.missingFlowProperty
        }
        let rt = gas.specificGasConstantJoulesPerKilogramKelvin*temperature
        guard rt.isFinite, rt > 0 else { throw AerodynamicError.nonFiniteResult }
        let speed = flow.speedMetresPerSecond
        let cosine = speed == 0 ? 0 : max(-1,min(1,(flow.relativeVelocityMetresPerSecond*panel.outwardNormal).sum()/speed))
        let coordinates = [log10(kn),cosine,speed/sqrt(rt),wallTemperatureKelvin/temperature]
        let p = try pressure.value(at: coordinates), tau = try shear.value(at: coordinates)
        if !panel.exposedToFlow { return try SurfaceTraction(gaugePressurePascals: 0,shearPascals: .zero) }
        let tangent = aeroCross(panel.outwardNormal,aeroCross(flow.relativeVelocityMetresPerSecond,panel.outwardNormal))
        let tangentSpeed = aeroNorm(tangent)
        // A normal/rest query has no defined shear direction. Nonzero table shear is inconsistent.
        guard tangentSpeed > 1e-12*max(1,speed) || tau == 0 else { throw AerodynamicError.invalidCoefficients }
        let scale = flow.densityKilogramsPerCubicMetre*rt
        return try SurfaceTraction(gaugePressurePascals: p*scale,
            shearPascals: tangentSpeed > 0 ? -tangent/tangentSpeed*(tau*scale) : .zero)
    }
}
