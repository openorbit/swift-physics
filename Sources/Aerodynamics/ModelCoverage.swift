// SPDX-License-Identifier: Apache-2.0
import Foundation

/// Declared applicability, not a guarantee of accuracy. All conditions are conjunctive.
public struct AerodynamicFlowCoverage: Sendable {
    public let machNumber: ClosedRange<Double>
    public let reynoldsNumber: ClosedRange<Double>
    public let angleOfAttackRadians: ClosedRange<Double>
    public let sideslipRadians: ClosedRange<Double>
    public let knudsenNumber: ClosedRange<Double>
    public init(machNumber: ClosedRange<Double>, reynoldsNumber: ClosedRange<Double>,
                angleOfAttackRadians: ClosedRange<Double>, sideslipRadians: ClosedRange<Double>,
                knudsenNumber: ClosedRange<Double>) throws {
        guard [machNumber,reynoldsNumber,angleOfAttackRadians,sideslipRadians,knudsenNumber]
            .allSatisfy({ $0.lowerBound.isFinite && $0.upperBound.isFinite }),
              machNumber.lowerBound >= 0, reynoldsNumber.lowerBound > 0, knudsenNumber.lowerBound >= 0,
              angleOfAttackRadians.lowerBound >= -.pi, angleOfAttackRadians.upperBound <= .pi,
              sideslipRadians.lowerBound >= -.pi/2, sideslipRadians.upperBound <= .pi/2 else {
            throw AerodynamicError.invalidParameter
        }
        self.machNumber = machNumber; self.reynoldsNumber = reynoldsNumber
        self.angleOfAttackRadians = angleOfAttackRadians; self.sideslipRadians = sideslipRadians
        self.knudsenNumber = knudsenNumber
    }
    public func validate(flow: AerodynamicFlow, referenceLengthMetres: Double) throws {
        guard let mach = flow.machNumber, let re = try flow.reynoldsNumber(referenceLengthMetres: referenceLengthMetres),
              let kn = try flow.knudsenNumber(referenceLengthMetres: referenceLengthMetres),
              let alpha = flow.angleOfAttackRadians, let beta = flow.sideslipRadians else {
            throw AerodynamicError.missingFlowProperty
        }
        guard machNumber.contains(mach), reynoldsNumber.contains(re), knudsenNumber.contains(kn),
              angleOfAttackRadians.contains(alpha), sideslipRadians.contains(beta) else {
            throw AerodynamicError.outsideCoverage
        }
    }
}

/// Binds a coefficient provider to its reference dimensions and stated flow coverage.
/// Particularly useful for measured/CFD tables across transonic and high-incidence flow.
public struct CoveredBodyCoefficients: BodyAerodynamicCoefficientModel {
    public let model: any BodyAerodynamicCoefficientModel
    public let coverage: AerodynamicFlowCoverage
    public let reynoldsReferenceLengthMetres: Double
    public let momentLengthsMetres: SIMD3<Double>
    public var source: String { model.source }
    public var version: String { model.version }
    public init(model: any BodyAerodynamicCoefficientModel, coverage: AerodynamicFlowCoverage,
                reynoldsReferenceLengthMetres: Double, momentLengthsMetres: SIMD3<Double>) throws {
        guard reynoldsReferenceLengthMetres.isFinite, reynoldsReferenceLengthMetres > 0,
              [momentLengthsMetres.x,momentLengthsMetres.y,momentLengthsMetres.z].allSatisfy({ $0.isFinite && $0 > 0 }) else {
            throw AerodynamicError.invalidGeometry
        }
        self.model = model; self.coverage = coverage
        self.reynoldsReferenceLengthMetres = reynoldsReferenceLengthMetres; self.momentLengthsMetres = momentLengthsMetres
    }
    public func coefficients(for query: BodyAerodynamicQuery) throws -> AerodynamicCoefficients {
        guard query.reynoldsReferenceLengthMetres == reynoldsReferenceLengthMetres,
              query.momentLengthsMetres == momentLengthsMetres else { throw AerodynamicError.unsupportedQuery }
        try coverage.validate(flow: query.flow,referenceLengthMetres: reynoldsReferenceLengthMetres)
        return try model.coefficients(for: query)
    }
}

/// Explicit ordered Mach dispatch. Threshold belongs to the higher interval.
/// Gaps and other coverage checks are owned by the selected bounded provider.
/// Discontinuous databases remain discontinuous; this does not invent a transition blend.
public struct MachRegimeBodyCoefficients: BodyAerodynamicCoefficientModel {
    public let models: [any BodyAerodynamicCoefficientModel]
    public let machThresholds: [Double]
    public let source: String
    public let version: String
    public init(models: [any BodyAerodynamicCoefficientModel], machThresholds: [Double],
                source: String, version: String) throws {
        guard !models.isEmpty, models.count == machThresholds.count+1,
              machThresholds.allSatisfy({ $0.isFinite && $0 > 0 }),
              zip(machThresholds,machThresholds.dropFirst()).allSatisfy({ $0 < $1 }),
              !source.isEmpty, !version.isEmpty else { throw AerodynamicError.invalidParameter }
        self.models = models; self.machThresholds = machThresholds; self.source = source; self.version = version
    }
    public func coefficients(for query: BodyAerodynamicQuery) throws -> AerodynamicCoefficients {
        guard let mach = query.machNumber else { throw AerodynamicError.missingFlowProperty }
        return try models[machThresholds.prefix(while: { mach >= $0 }).count].coefficients(for: query)
    }
}
