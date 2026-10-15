// SPDX-License-Identifier: Apache-2.0
import FluidDynamics
import Foundation

/// Ideal single-species Maxwellian gas parameter, with R in J/(kg K).
/// Using a mixture-average R is an additional approximation, not a resolved mixture.
public struct FreeMolecularGas: Sendable {
    public let identifier: String
    public let specificGasConstantJoulesPerKilogramKelvin: Double
    public let source: String
    public init(identifier: String, specificGasConstantJoulesPerKilogramKelvin: Double, source: String) throws {
        guard [identifier,source].allSatisfy({ !$0.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty }),
              specificGasConstantJoulesPerKilogramKelvin.isFinite, specificGasConstantJoulesPerKilogramKelvin > 0 else {
            throw FluidError.invalidParameter
        }
        self.identifier = identifier; self.specificGasConstantJoulesPerKilogramKelvin = specificGasConstantJoulesPerKilogramKelvin
        self.source = source
    }
}
/// Maxwell reflection mixture: a specified fraction reflects specularly; the remainder
/// is fully accommodated diffuse re-emission at wall temperature, with zero net mass flux.
public struct GasSurfaceReflection: Sendable {
    public let specularFraction: Double
    public let wallTemperatureKelvin: Double
    public init(specularFraction: Double, wallTemperatureKelvin: Double) throws {
        guard specularFraction.isFinite, (0...1).contains(specularFraction),
              wallTemperatureKelvin.isFinite, wallTemperatureKelvin > 0 else { throw FluidError.invalidParameter }
        self.specularFraction = specularFraction; self.wallTemperatureKelvin = wallTemperatureKelvin
    }
}

/// Independent single-collision Maxwellian pressure and shear, including leeward
/// thermal impacts. Requires Kn≥10; no intermolecular collisions, chemistry,
/// multiple reflections, shielding calculation, ablation or continuum bridging.
public struct MaxwellianFreeMolecularTraction: SurfaceTractionModel {
    public let gas: FreeMolecularGas
    public let surface: GasSurfaceReflection
    public let minimumKnudsenNumber: Double
    public init(gas: FreeMolecularGas, surface: GasSurfaceReflection, minimumKnudsenNumber: Double = 10) throws {
        guard minimumKnudsenNumber.isFinite, minimumKnudsenNumber >= 10 else { throw FluidError.invalidParameter }
        self.gas = gas; self.surface = surface; self.minimumKnudsenNumber = minimumKnudsenNumber
    }
    public func traction(panel: FluidSurfacePanel, flow: FluidFlow,
                         referenceLengthMetres: Double) throws -> SurfaceTraction {
        guard let temperature = flow.temperatureKelvin,
              let kn = try flow.knudsenNumber(referenceLengthMetres: referenceLengthMetres) else { throw FluidError.missingFlowProperty }
        guard kn >= minimumKnudsenNumber else { throw FluidError.outsideCoverage }
        let r = gas.specificGasConstantJoulesPerKilogramKelvin
        let sigma = sqrt(r)*sqrt(temperature), wallSigma = sqrt(r)*sqrt(surface.wallTemperatureKelvin)
        guard sigma.isFinite, sigma > 0, wallSigma.isFinite else { throw FluidError.nonFiniteResult }
        let normalVelocity = (flow.relativeVelocityMetresPerSecond*panel.outwardNormal).sum()
        let a = normalVelocity/sigma
        // Explicit numerical coverage bound; not a Mach or physical gas validity limit.
        guard a.isFinite, abs(a) <= 1e6 else { throw FluidError.outsideCoverage }
        if !panel.exposedToFlow { return try SurfaceTraction(gaugePressurePascals: 0,shearPascals: .zero) }
        let phi = exp(-0.5*a*a)/sqrt(2*Double.pi)
        let massMoment: Double, pressureMoment: Double
        if a < -10 {
            // Stable Gaussian-tail expansion avoids cancellation on deep leeward faces.
            let x = -a, t = 1/(x*x)
            massMoment = phi*t*(1-3*t+15*t*t-105*t*t*t+945*t*t*t*t)
            pressureMoment = phi*(t/x)*(2-12*t+90*t*t-840*t*t*t+9450*t*t*t*t)
        } else {
            let probability = 0.5*erfc(-a/sqrt(2))
            massMoment = max(0,phi+a*probability)
            pressureMoment = max(0,(a*a+1)*probability+a*phi)
        }
        let rho = flow.densityKilogramsPerCubicMetre
        let incidentMassFlux = (rho*sigma)*massMoment
        let ambientPressure = (rho*sigma)*sigma
        let incidentPressure = ambientPressure*pressureMoment
        let specular = surface.specularFraction, diffuse = 1-specular
        let reflectedDiffusePressure = incidentMassFlux*wallSigma*sqrt(Double.pi/2)
        let pressure = (1+specular)*incidentPressure+diffuse*reflectedDiffusePressure-ambientPressure
        let tangentVelocity = fluidCross(panel.outwardNormal,
            fluidCross(flow.relativeVelocityMetresPerSecond,panel.outwardNormal))
        let shear = -tangentVelocity*(diffuse*incidentMassFlux)
        return try SurfaceTraction(gaugePressurePascals: pressure,shearPascals: shear)
    }
}
