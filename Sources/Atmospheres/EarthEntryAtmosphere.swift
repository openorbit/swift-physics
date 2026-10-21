// SPDX-License-Identifier: Apache-2.0
// Copyright 2026 Mattias Holm
import Foundation

/// ISA lower layers joined continuously to a sparse US Standard Atmosphere 1976
/// upper table. Log-linear pressure/density, linear kinetic temperature. Fixed
/// nominal air, not weather/solar activity, chemistry or continuum certification.
public struct EarthEntryAtmosphere: AtmosphericProfile {
    public let altitudeRangeMetres = -500.0...300000.0
    public let altitudeDatum = "Geometric metres above the simulation mean sea level"
    public let source = "https://www.ngdc.noaa.gov/stp/space-weather/online-publications/miscellaneous/us-standard-atmosphere-1976/us-standard-atmosphere_st76-1562_noaa.pdf"
    public let assumptions = "ISA to its explicit upper limit; USSA1976 Table I pp.68–70 at 90–300 km. Sparse log interpolation; match ISA at join. Upper sound speed uses gamma=1.4 and p/rho; Sutherland viscosity is a nominal transport approximation, not resolved composition."
    public init() {}
    public func properties(at altitudeMetres: Double) throws -> AtmosphericProperties {
        guard altitudeMetres.isFinite, altitudeRangeMetres.contains(altitudeMetres) else { throw AtmosphericModelError.altitudeOutsideCoverage }
        let isa = ISAAtmosphere(), join = isa.altitudeRangeMetres.upperBound
        if altitudeMetres <= join { return try isa.properties(at:altitudeMetres) }
        let base = try isa.properties(at:join)
        // Geometric km, kinetic K, pressure Pa (publication mb ×100), rho kg/m³.
        let rows: [(Double,Double,Double,Double)] = [
            (join,base.temperatureKelvin,base.pressurePascals,base.densityKilogramsPerCubicMetre),
            (90000,186.87,0.18359,3.416e-6), (95000,188.42,0.075966,1.393e-6),
            (100000,195.08,0.032011,5.604e-7), (105000,208.84,0.014477,2.325e-7),
            (110000,240,0.0071042,9.708e-8), (115000,300,0.0040075,4.061e-8),
            (120000,360,0.0025382,2.222e-8), (125000,417.23,0.0017340,1.291e-8),
            (130000,469.27,0.0012505,8.152e-9),
            (140000,559.63,0.00072028,3.831e-9), (150000,634.39,0.00045422,2.076e-9),
            (160000,696.29,0.00030395,1.233e-9), (170000,747.57,0.00021210,7.815e-10),
            (180000,790.07,0.00015271,5.194e-10), (190000,825.31,0.00011266,3.581e-10),
            (200000,854.56,0.000084736,2.541e-10), (220000,899.01,0.000050149,1.367e-10),
            (230000,915.78,0.000039276,1.029e-10), (240000,929.73,0.000031059,7.858e-11),
            (250000,941.33,0.000024767,6.073e-11), (300000,976.00,0.0000087704,1.916e-11)]
        let i = min(rows.count-2,rows.lastIndex { $0.0 <= altitudeMetres }!)
        let a = rows[i], b = rows[i+1], u = (altitudeMetres-a.0)/(b.0-a.0)
        let t = a.1+(b.1-a.1)*u, p = exp(log(a.2)+(log(b.2)-log(a.2))*u)
        let rho = exp(log(a.3)+(log(b.3)-log(a.3))*u)
        return AtmosphericProperties(pressurePascals:p,temperatureKelvin:t,densityKilogramsPerCubicMetre:rho,
            speedOfSoundMetresPerSecond:sqrt(1.4*p/rho),dynamicViscosityPascalSeconds:1.458e-6*pow(t,1.5)/(t+110.4))
    }
}
