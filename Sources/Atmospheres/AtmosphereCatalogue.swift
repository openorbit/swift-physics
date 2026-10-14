// SPDX-License-Identifier: Apache-2.0
import Foundation

/// Body centres use NAIF IDs; system barycentres deliberately have no profile.
public enum AtmosphericBody: Int32, CaseIterable, Sendable {
    case venus = 299, earth = 399, mars = 499, jupiter = 599, saturn = 699
    case uranus = 799, neptune = 899, titan = 606, triton = 801, io = 501, callisto = 504
}
public enum AtmosphericRegime: Sendable { case continuum, tenuous, exosphere }

/// Fixed nominal profiles, not forecasts or mission-grade entry models.
public struct BodyAtmosphere: AtmosphericProfile {
    public let body: AtmosphericBody
    public let regime: AtmosphericRegime
    private let profile: any AtmosphericProfile
    public var altitudeRangeMetres: ClosedRange<Double> { profile.altitudeRangeMetres }
    public var altitudeDatum: String { profile.altitudeDatum }
    public var source: String { profile.source }
    public var assumptions: String { profile.assumptions }
    private init(body: AtmosphericBody, regime: AtmosphericRegime = .continuum, profile: any AtmosphericProfile) {
        self.body = body; self.regime = regime; self.profile = profile
    }
    public func properties(at altitudeMetres: Double) throws -> AtmosphericProperties {
        try profile.properties(at: altitudeMetres)
    }
    public static func model(for body: AtmosphericBody) -> BodyAtmosphere {
        if body == .earth { return BodyAtmosphere(body: body, profile: ISAAtmosphere()) }
        let r = 8314.46261815324 // Universal gas constant divided by molar mass in g/mol.
        let fact = "https://nssdc.gsfc.nasa.gov/planetary/factsheet/"
        var pressure: Double, temperature: Double, molecularMass: Double, gravity: Double
        var gamma: Double? = 1.4, top = 100_000.0, regime = AtmosphericRegime.continuum
        var source: String
        var levels: [HydrostaticAtmosphere.TemperatureLevel]
        var datum = "Metres above nominal mean surface"
        var assumptions = "Nominal isothermal ideal gas, constant gravity and composition; no weather, winds, seasons or escape. Coverage is a chosen simulation bound, not a measured atmospheric edge. Viscosity unavailable."
        switch body {
        case .venus: (pressure,temperature,molecularMass,gravity) = (9.2e6,737,43.45,8.87); gamma = 1.3; source = fact+"venusfact.html"
        case .mars: (pressure,temperature,molecularMass,gravity) = (636,214,43.49,3.73); gamma = 1.3; source = fact+"marsfact.html"
        case .titan: (pressure,temperature,molecularMass,gravity) = (1.5e5,94,28.0,1.352); source = "https://www.nasa.gov/wp-content/uploads/2015/04/niac_sunspiral_superballbot_phaseii_finalreport_tagged.pdf"
        case .triton: (pressure,temperature,molecularMass,gravity) = (1.4,38,28.0,0.779); regime = .tenuous; source = "https://ntrs.nasa.gov/citations/19950012719"
        case .io: (pressure,temperature,molecularMass,gravity) = (1e-4,150,64.066,1.796); gamma = nil; regime = .tenuous; top = 50_000; source = "https://doi.org/10.1006/icar.1996.0186"
            assumptions += " Io SO2 pressure is order-of-magnitude only; 150 K is a chosen nominal gas temperature. Volcanism and eclipse collapse omitted; continuum transport unavailable."
        case .callisto: (pressure,temperature,molecularMass,gravity) = (7.5e-7,150,44.01,1.236); gamma = nil; regime = .exosphere; source = "https://doi.org/10.1126/science.283.5403.820"
            assumptions += " Historical CO2 component only; pressure is an equivalent number-density measure. Exponential density is a proxy, not a collisional hydrostatic atmosphere; no sound or continuum transport."
        case .jupiter: (pressure,temperature,molecularMass,gravity) = (1e5,165,2.22,24.79); source = fact+"jupiterfact.html"
        case .saturn: (pressure,temperature,molecularMass,gravity) = (1e5,134,2.07,10.44); source = fact+"saturnfact.html"
        case .uranus: (pressure,temperature,molecularMass,gravity) = (1e5,76,2.64,8.87); source = fact+"uranusfact.html"
        case .neptune: (pressure,temperature,molecularMass,gravity) = (1e5,72,2.60,11.15); source = fact+"neptunefact.html"
        case .earth: preconditionFailure("Handled above")
        }
        let gasConstant = r/molecularMass
        if [.jupiter,.saturn,.uranus,.neptune].contains(body) {
            // Dry ideal diatomic adiabat below 1 bar, limited to 100 bar.
            let lapse = -gravity*(1.4-1)/1.4/gasConstant
            let deepTemperature = temperature*pow(100,(1.4-1)/1.4)
            let bottom = (deepTemperature-temperature)/lapse
            let upperTemperature: Double
            switch body { case .jupiter: upperTemperature = 112; case .saturn: upperTemperature = 84
            case .uranus: upperTemperature = 53; default: upperTemperature = 55 }
            // Linear temperature layer calibrated to the published 0.1-bar temperature.
            let upperHeight = gasConstant/gravity*(temperature-upperTemperature)
                / log(temperature/upperTemperature)*log(10)
            levels = [(bottom,deepTemperature),(0,temperature),(upperHeight,upperTemperature),
                      (upperHeight+5*gasConstant*upperTemperature/gravity,upperTemperature)]
                .map { .init(altitudeMetres: $0.0, temperatureKelvin: $0.1) }
            datum = "Metres above the 1-bar level; negative altitudes descend into the clouds, no solid surface"
            assumptions = "Ideal gas with constant gravity/composition; dry gamma=1.4 adiabat from 1 to 100 bar, linear temperature fit from 1 to 0.1 bar, then isothermal. No cloud microphysics, winds, condensation or deep interior equation of state. Coverage is a simulation bound. Viscosity unavailable."
        } else {
            levels = [(0.0,temperature),(top,temperature)].map { .init(altitudeMetres: $0.0, temperatureKelvin: $0.1) }
        }
        return BodyAtmosphere(body: body, regime: regime, profile: try! HydrostaticAtmosphere(levels: levels,
            referencePressurePascals: pressure, specificGasConstant: gasConstant,
            gravityMetresPerSecondSquared: gravity, heatCapacityRatio: gamma,
            altitudeDatum: datum, source: source, assumptions: assumptions))
    }
    public static func model(naifBodyCode: Int32) -> BodyAtmosphere? {
        AtmosphericBody(rawValue: naifBodyCode).map { model(for: $0) }
    }
}
