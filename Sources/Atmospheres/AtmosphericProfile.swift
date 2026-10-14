// SPDX-License-Identifier: Apache-2.0
import Foundation

public enum AtmosphericModelError: Error, Equatable {
    case invalidParameters, altitudeOutsideCoverage, nonFiniteResult
}

/// SI properties. Optional transport quantities are absent where no continuum model is supplied.
public struct AtmosphericProperties: Sendable {
    public let pressurePascals: Double
    public let temperatureKelvin: Double
    public let densityKilogramsPerCubicMetre: Double
    public let speedOfSoundMetresPerSecond: Double?
    public let dynamicViscosityPascalSeconds: Double?
}

/// Altitude is measured in metres above the explicitly described datum.
public protocol AtmosphericProfile: Sendable {
    var altitudeRangeMetres: ClosedRange<Double> { get }
    var altitudeDatum: String { get }
    var source: String { get }
    var assumptions: String { get }
    func properties(at altitudeMetres: Double) throws -> AtmosphericProperties
}

/// Piecewise linear temperature with constant gravity and molecular composition.
/// Pressure is integrated analytically from hydrostatic balance; density uses the ideal gas law.
public struct HydrostaticAtmosphere: AtmosphericProfile {
    public struct TemperatureLevel: Sendable {
        public let altitudeMetres: Double
        public let temperatureKelvin: Double
        public init(altitudeMetres: Double, temperatureKelvin: Double) {
            self.altitudeMetres = altitudeMetres; self.temperatureKelvin = temperatureKelvin
        }
    }
    public let altitudeRangeMetres: ClosedRange<Double>
    public let altitudeDatum: String
    public let source: String
    public let assumptions: String
    public let specificGasConstant: Double
    public let gravityMetresPerSecondSquared: Double
    public let heatCapacityRatio: Double?
    private let levels: [TemperatureLevel]
    private let logPressures: [Double]

    public init(levels: [TemperatureLevel], referenceAltitudeMetres: Double = 0,
                referencePressurePascals: Double, specificGasConstant: Double,
                gravityMetresPerSecondSquared: Double, heatCapacityRatio: Double? = nil,
                altitudeDatum: String, source: String, assumptions: String) throws {
        guard levels.count >= 2,
              levels.allSatisfy({ $0.altitudeMetres.isFinite && $0.temperatureKelvin.isFinite && $0.temperatureKelvin > 0 }),
              zip(levels, levels.dropFirst()).allSatisfy({ $0.altitudeMetres < $1.altitudeMetres }),
              [referencePressurePascals, specificGasConstant, gravityMetresPerSecondSquared].allSatisfy({ $0.isFinite && $0 > 0 }),
              heatCapacityRatio.map({ $0.isFinite && $0 > 1 }) ?? true,
              !altitudeDatum.isEmpty, !source.isEmpty,
              let referenceIndex = levels.firstIndex(where: { $0.altitudeMetres == referenceAltitudeMetres }) else {
            throw AtmosphericModelError.invalidParameters
        }
        self.levels = levels; self.specificGasConstant = specificGasConstant
        self.gravityMetresPerSecondSquared = gravityMetresPerSecondSquared
        self.heatCapacityRatio = heatCapacityRatio; self.altitudeDatum = altitudeDatum
        self.source = source; self.assumptions = assumptions
        altitudeRangeMetres = levels.first!.altitudeMetres...levels.last!.altitudeMetres
        var integrals = [0.0]
        for i in 1..<levels.count {
            let lower = levels[i-1], upper = levels[i]
            integrals.append(integrals[i-1] + Self.logPressureRatio(
                height: upper.altitudeMetres-lower.altitudeMetres,
                baseTemperature: lower.temperatureKelvin, temperature: upper.temperatureKelvin,
                gravityOverR: gravityMetresPerSecondSquared/specificGasConstant))
        }
        logPressures = integrals.map { log(referencePressurePascals) + $0 - integrals[referenceIndex] }
        guard logPressures.allSatisfy(\.isFinite) else { throw AtmosphericModelError.nonFiniteResult }
    }

    private static func logPressureRatio(height: Double, baseTemperature: Double,
                                         temperature: Double, gravityOverR: Double) -> Double {
        let difference = temperature-baseTemperature
        if abs(difference) < 1e-10 * baseTemperature {
            return -gravityOverR * height/baseTemperature
        }
        return -gravityOverR * height/difference * log1p(difference/baseTemperature)
    }

    public func properties(at altitudeMetres: Double) throws -> AtmosphericProperties {
        guard altitudeMetres.isFinite, altitudeRangeMetres.contains(altitudeMetres) else {
            throw AtmosphericModelError.altitudeOutsideCoverage
        }
        let i = min(levels.lastIndex(where: { $0.altitudeMetres <= altitudeMetres })!, levels.count-2)
        let lower = levels[i], upper = levels[i+1]
        let height = altitudeMetres-lower.altitudeMetres
        let temperature = lower.temperatureKelvin + (upper.temperatureKelvin-lower.temperatureKelvin)
            * height/(upper.altitudeMetres-lower.altitudeMetres)
        let pressure = exp(logPressures[i] + Self.logPressureRatio(height: height,
            baseTemperature: lower.temperatureKelvin, temperature: temperature,
            gravityOverR: gravityMetresPerSecondSquared/specificGasConstant))
        let density = pressure/specificGasConstant/temperature
        let sound = heatCapacityRatio.map { sqrt($0*specificGasConstant*temperature) }
        guard pressure.isFinite, pressure > 0, density.isFinite, density > 0,
              sound.map(\.isFinite) ?? true else { throw AtmosphericModelError.nonFiniteResult }
        return AtmosphericProperties(pressurePascals: pressure, temperatureKelvin: temperature,
            densityKilogramsPerCubicMetre: density, speedOfSoundMetresPerSecond: sound,
            dynamicViscosityPascalSeconds: nil)
    }
}

/// ISA dry-air layers through 84.852 km geopotential altitude (about 86 km geometric).
public struct ISAAtmosphere: AtmosphericProfile {
    public let altitudeDatum = "Geometric altitude above mean sea level; internally converted to geopotential height"
    public let source = "https://ntrs.nasa.gov/citations/19770009539"
    public let assumptions = "ISA/US Standard Atmosphere lower layers; dry ideal air, no weather or thermosphere"
    public let altitudeRangeMetres = -500.0...(6_356_766 * 84_852 / (6_356_766 - 84_852))
    private static let profile = try! HydrostaticAtmosphere(
        levels: [(-1000.0,294.65),(0,288.15),(11000,216.65),(20000,216.65),
                 (32000,228.65),(47000,270.65),(51000,270.65),(71000,214.65),(84852,186.946)]
            .map { .init(altitudeMetres: $0.0, temperatureKelvin: $0.1) },
        referencePressurePascals: 101325, specificGasConstant: 287.05287,
        gravityMetresPerSecondSquared: 9.80665, heatCapacityRatio: 1.4,
        altitudeDatum: "Geopotential height", source: "US Standard Atmosphere 1976", assumptions: "Dry air")
    public init() {}
    public func properties(at altitudeMetres: Double) throws -> AtmosphericProperties {
        guard altitudeMetres.isFinite, altitudeRangeMetres.contains(altitudeMetres) else {
            throw AtmosphericModelError.altitudeOutsideCoverage
        }
        let h = 6_356_766*altitudeMetres/(6_356_766+altitudeMetres)
        let result = try Self.profile.properties(at: min(84852, h))
        let t = result.temperatureKelvin
        // Sutherland's law for dry air, Pa s.
        let viscosity = 1.458e-6 * pow(t,1.5)/(t+110.4)
        return AtmosphericProperties(pressurePascals: result.pressurePascals,
            temperatureKelvin: t, densityKilogramsPerCubicMetre: result.densityKilogramsPerCubicMetre,
            speedOfSoundMetresPerSecond: result.speedOfSoundMetresPerSecond,
            dynamicViscosityPascalSeconds: viscosity)
    }
}
