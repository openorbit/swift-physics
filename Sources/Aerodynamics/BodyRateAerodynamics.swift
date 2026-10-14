// SPDX-License-Identifier: Apache-2.0
import Foundation

/// Columns of the six-coefficient derivative with respect to pHat, qHat, rHat.
/// Signed derivatives may describe damping or excitation; no stability is inferred.
public struct BodyAerodynamicRateDerivatives: Sendable {
    public let rollRate: AerodynamicCoefficients
    public let pitchRate: AerodynamicCoefficients
    public let yawRate: AerodynamicCoefficients
    public init(rollRate: AerodynamicCoefficients, pitchRate: AerodynamicCoefficients,
                yawRate: AerodynamicCoefficients) {
        self.rollRate = rollRate; self.pitchRate = pitchRate; self.yawRate = yawRate
    }
    public func increment(nondimensionalRates rates: SIMD3<Double>) throws -> AerodynamicCoefficients {
        guard aeroFinite(rates) else { throw AerodynamicError.invalidParameter }
        return try AerodynamicCoefficients(force: rollRate.force*rates.x+pitchRate.force*rates.y+yawRate.force*rates.z,
            moment: rollRate.moment*rates.x+pitchRate.moment*rates.y+yawRate.moment*rates.z)
    }
}

/// Derivatives may depend on the static flow state, but are queried at zero rates.
public protocol BodyAerodynamicRateDerivativeModel: Sendable {
    var source: String { get }
    var version: String { get }
    func derivatives(for query: BodyAerodynamicQuery) throws -> BodyAerodynamicRateDerivatives
}

/// Prescribed derivatives within explicit Mach/Re/incidence/sideslip coverage.
/// Reference dimensions must match the dataset; values are not a capsule database.
public struct ConstantBodyRateDerivatives: BodyAerodynamicRateDerivativeModel {
    public let value: BodyAerodynamicRateDerivatives
    public let source: String
    public let version: String
    public let machRange: ClosedRange<Double>
    public let reynoldsRange: ClosedRange<Double>
    public let incidenceRangeRadians: ClosedRange<Double>
    public let sideslipRangeRadians: ClosedRange<Double>
    public let momentLengthsMetres: SIMD3<Double>
    public let reynoldsReferenceLengthMetres: Double
    public init(value: BodyAerodynamicRateDerivatives, source: String, version: String,
                machRange: ClosedRange<Double>, reynoldsRange: ClosedRange<Double>,
                incidenceRangeRadians: ClosedRange<Double>, sideslipRangeRadians: ClosedRange<Double>,
                momentLengthsMetres: SIMD3<Double>, reynoldsReferenceLengthMetres: Double) throws {
        guard [source,version].allSatisfy({ !$0.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty }),
              machRange.lowerBound.isFinite, machRange.lowerBound >= 0, machRange.upperBound.isFinite,
              reynoldsRange.lowerBound.isFinite, reynoldsRange.lowerBound > 0, reynoldsRange.upperBound.isFinite,
              [incidenceRangeRadians.lowerBound,incidenceRangeRadians.upperBound].allSatisfy({ $0.isFinite && abs($0) <= .pi }),
              [sideslipRangeRadians.lowerBound,sideslipRangeRadians.upperBound].allSatisfy({ $0.isFinite && abs($0) <= .pi/2 }),
              [momentLengthsMetres.x,momentLengthsMetres.y,momentLengthsMetres.z].allSatisfy({ $0.isFinite && $0 > 0 }),
              reynoldsReferenceLengthMetres.isFinite, reynoldsReferenceLengthMetres > 0 else { throw AerodynamicError.invalidParameter }
        self.value = value; self.source = source; self.version = version; self.machRange = machRange; self.reynoldsRange = reynoldsRange
        self.incidenceRangeRadians = incidenceRangeRadians; self.sideslipRangeRadians = sideslipRangeRadians
        self.momentLengthsMetres = momentLengthsMetres; self.reynoldsReferenceLengthMetres = reynoldsReferenceLengthMetres
    }
    public func derivatives(for query: BodyAerodynamicQuery) throws -> BodyAerodynamicRateDerivatives {
        guard query.momentLengthsMetres == momentLengthsMetres,
              query.reynoldsReferenceLengthMetres == reynoldsReferenceLengthMetres else { throw AerodynamicError.invalidGeometry }
        guard let mach = query.machNumber, let re = try query.reynoldsNumber(),
              let alpha = query.angleOfAttackRadians, let beta = query.sideslipRadians else { throw AerodynamicError.missingFlowProperty }
        guard machRange.contains(mach), reynoldsRange.contains(re), incidenceRangeRadians.contains(alpha),
              sideslipRangeRadians.contains(beta) else { throw AerodynamicError.outsideCoverage }
        return value
    }
}

/// C = Cstatic + D[pHat,qHat,rHat], with rate_i = ω_i L_i/(2V).
/// Small-rate linearization only; no alphaDot, betaDot, time-history or unsteady lift.
public struct LinearBodyRateAerodynamics: BodyAerodynamicCoefficientModel {
    public let staticModel: any BodyAerodynamicCoefficientModel
    public let derivativeModel: any BodyAerodynamicRateDerivativeModel
    public let rateReference: BodyAerodynamicRateReference
    public let minimumSpeedMetresPerSecond: Double
    public let maximumAbsoluteNondimensionalRates: SIMD3<Double>
    public var source: String { "Static: \(staticModel.source); rate derivatives: \(derivativeModel.source)" }
    public var version: String { "\(staticModel.version) / \(derivativeModel.version)" }
    public init(staticModel: any BodyAerodynamicCoefficientModel,
                derivativeModel: any BodyAerodynamicRateDerivativeModel,
                rateReference: BodyAerodynamicRateReference,
                minimumSpeedMetresPerSecond: Double,
                maximumAbsoluteNondimensionalRates: SIMD3<Double>) throws {
        guard minimumSpeedMetresPerSecond.isFinite, minimumSpeedMetresPerSecond > 0,
              aeroFinite(maximumAbsoluteNondimensionalRates),
              [maximumAbsoluteNondimensionalRates.x,maximumAbsoluteNondimensionalRates.y,
               maximumAbsoluteNondimensionalRates.z].allSatisfy({ $0 >= 0 }) else { throw AerodynamicError.invalidParameter }
        self.staticModel = staticModel; self.derivativeModel = derivativeModel; self.rateReference = rateReference
        self.minimumSpeedMetresPerSecond = minimumSpeedMetresPerSecond
        self.maximumAbsoluteNondimensionalRates = maximumAbsoluteNondimensionalRates
    }
    public func coefficients(for query: BodyAerodynamicQuery) throws -> AerodynamicCoefficients {
        guard query.angularRateReference == rateReference else { throw AerodynamicError.unsupportedQuery }
        guard query.flow.speedMetresPerSecond >= minimumSpeedMetresPerSecond,
              let rates = query.nondimensionalAngularRates else { throw AerodynamicError.outsideCoverage }
        let bounds = maximumAbsoluteNondimensionalRates
        guard abs(rates.x) <= bounds.x, abs(rates.y) <= bounds.y, abs(rates.z) <= bounds.z else { throw AerodynamicError.outsideCoverage }
        let staticQuery = try BodyAerodynamicQuery(flow: query.flow,
            reynoldsReferenceLengthMetres: query.reynoldsReferenceLengthMetres,
            momentLengthsMetres: query.momentLengthsMetres, angularRateReference: rateReference)
        let base = try staticModel.coefficients(for: staticQuery)
        let increment = try derivativeModel.derivatives(for: staticQuery).increment(nondimensionalRates: rates)
        return try AerodynamicCoefficients(force: base.force+increment.force, moment: base.moment+increment.moment)
    }
}
