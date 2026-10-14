// SPDX-License-Identifier: Apache-2.0
import Foundation

/// Whole-body inputs in X-forward/Y-right/Z-down axes. Angular rates are carried
/// for providers with calibrated damping; no damping coefficient is invented.
public struct BodyAerodynamicQuery: Sendable {
    public let flow: AerodynamicFlow
    public let reynoldsReferenceLengthMetres: Double
    public let momentLengthsMetres: SIMD3<Double>
    public let angularVelocityRadiansPerSecond: SIMD3<Double>
    public var angleOfAttackRadians: Double? { flow.angleOfAttackRadians }
    public var sideslipRadians: Double? { flow.sideslipRadians }
    public var machNumber: Double? { flow.machNumber }
    public func reynoldsNumber() throws -> Double? { try flow.reynoldsNumber(referenceLengthMetres: reynoldsReferenceLengthMetres) }
    public var nondimensionalAngularRates: SIMD3<Double>? {
        guard flow.speedMetresPerSecond > 0 else { return nil }
        return ((angularVelocityRadiansPerSecond*momentLengthsMetres)/flow.speedMetresPerSecond)*0.5
    }
    public init(flow: AerodynamicFlow, reynoldsReferenceLengthMetres: Double,
                momentLengthsMetres: SIMD3<Double>, angularVelocityRadiansPerSecond: SIMD3<Double> = .zero) throws {
        guard reynoldsReferenceLengthMetres.isFinite, reynoldsReferenceLengthMetres > 0,
              [momentLengthsMetres.x,momentLengthsMetres.y,momentLengthsMetres.z].allSatisfy({ $0.isFinite && $0 > 0 }),
              aeroFinite(angularVelocityRadiansPerSecond) else { throw AerodynamicError.invalidParameter }
        self.flow = flow; self.reynoldsReferenceLengthMetres = reynoldsReferenceLengthMetres
        self.momentLengthsMetres = momentLengthsMetres; self.angularVelocityRadiansPerSecond = angularVelocityRadiansPerSecond
        guard nondimensionalAngularRates.map(aeroFinite) ?? true else { throw AerodynamicError.nonFiniteResult }
    }
}
public protocol BodyAerodynamicCoefficientModel: Sendable {
    var source: String { get }
    var version: String { get }
    /// Six coefficients in body axes, with moments about the supplied geometry's
    /// reference point. Dataset dimensions must match that geometry.
    func coefficients(for query: BodyAerodynamicQuery) throws -> AerodynamicCoefficients
}

/// Prescribed coefficients for controlled scenarios; no capsule geometry or
/// scientific validity is implied. Coefficients are fixed in body axes, not wind axes.
public struct ConstantBodyCoefficients: BodyAerodynamicCoefficientModel {
    public let value: AerodynamicCoefficients
    public let source: String
    public let version: String
    public init(value: AerodynamicCoefficients, source: String, version: String) throws {
        guard !source.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              !version.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else { throw AerodynamicError.invalidParameter }
        self.value = value; self.source = source; self.version = version
    }
    public func coefficients(for query: BodyAerodynamicQuery) throws -> AerodynamicCoefficients { value }
}

/// Wind axes X along body-relative-to-air velocity, Y to wind-right, Z wind-down.
/// Converts drag, side force and lift into signed body force coefficients.
/// Moments remain separate body-axis coefficients supplied by the dataset.
public struct WindForceCoefficients: Sendable {
    public let drag: Double
    public let side: Double
    public let lift: Double
    public init(drag: Double, side: Double = 0, lift: Double) throws {
        guard [drag,side,lift].allSatisfy(\.isFinite), drag >= 0 else { throw AerodynamicError.invalidCoefficients }
        self.drag = drag; self.side = side; self.lift = lift
    }
    public func bodyCoefficients(angleOfAttackRadians alpha: Double, sideslipRadians beta: Double,
                                 momentBody: SIMD3<Double> = .zero) throws -> AerodynamicCoefficients {
        guard alpha.isFinite, beta.isFinite, abs(alpha) <= .pi, abs(beta) <= .pi/2 else { throw AerodynamicError.invalidParameter }
        let ca = cos(alpha), sa = sin(alpha), cb = cos(beta), sb = sin(beta)
        let x = SIMD3(ca*cb,sb,sa*cb), y = SIMD3(-ca*sb,cb,-sa*sb), z = SIMD3(-sa,0,ca)
        return try AerodynamicCoefficients(force: -drag*x+side*y-lift*z,moment: momentBody)
    }
}

/// Generic capsule/fuselage/whole-body evaluator driven by a supplied coefficient
/// model. No automatic model selection, heating, deployment or motion integration.
public struct WholeBodyAerodynamics: Sendable {
    public let geometry: AerodynamicReferenceGeometry
    public let reynoldsReferenceLengthMetres: Double
    public let model: any BodyAerodynamicCoefficientModel
    public init(geometry: AerodynamicReferenceGeometry, reynoldsReferenceLengthMetres: Double,
                model: any BodyAerodynamicCoefficientModel) throws {
        guard reynoldsReferenceLengthMetres.isFinite, reynoldsReferenceLengthMetres > 0 else { throw AerodynamicError.invalidGeometry }
        self.geometry = geometry; self.reynoldsReferenceLengthMetres = reynoldsReferenceLengthMetres; self.model = model
    }
    public func evaluate(flow: AerodynamicFlow, angularVelocityRadiansPerSecond: SIMD3<Double> = .zero,
                         momentReferencePointMetres: SIMD3<Double>? = nil) throws -> AerodynamicLoads {
        let query = try BodyAerodynamicQuery(flow: flow,reynoldsReferenceLengthMetres: reynoldsReferenceLengthMetres,
            momentLengthsMetres: geometry.momentLengthsMetres,angularVelocityRadiansPerSecond: angularVelocityRadiansPerSecond)
        let loads: AerodynamicLoads
        if flow.dynamicPressurePascals == 0 {
            loads = try AerodynamicLoads(forceNewtons: .zero,momentNewtonMetres: .zero,
                                        referencePointMetres: geometry.momentReferencePointMetres)
        } else {
            loads = try model.coefficients(for: query).loads(dynamicPressurePascals: flow.dynamicPressurePascals,geometry: geometry)
        }
        return try momentReferencePointMetres.map { try loads.referred(to: $0) } ?? loads
    }
}

/// One whole-body coefficient sample. Conditions use radians and dimensionless Mach/Reynolds.
public struct BodyCoefficientSample: Sendable, Codable {
    public let angleOfAttackRadians: Double
    public let sideslipRadians: Double
    public let machNumber: Double
    public let reynoldsNumber: Double
    public let coefficients: AerodynamicCoefficients
    public init(angleOfAttackRadians: Double, sideslipRadians: Double, machNumber: Double,
                reynoldsNumber: Double, coefficients: AerodynamicCoefficients) {
        self.angleOfAttackRadians = angleOfAttackRadians; self.sideslipRadians = sideslipRadians
        self.machNumber = machNumber; self.reynoldsNumber = reynoldsNumber; self.coefficients = coefficients
    }
}

/// Validated rectangular incidence/sideslip/Mach/Reynolds grid, linearly interpolated
/// across up to 16 corners. Static coefficients: supplied angular rates are not used.
public struct TabulatedBodyCoefficients: BodyAerodynamicCoefficientModel, Sendable, Codable {
    public let samples: [BodyCoefficientSample]
    public let boundaryPolicy: CoefficientBoundaryPolicy
    public let source: String
    public let version: String
    private let axes: [[Double]]
    private let values: [[Double]: AerodynamicCoefficients]
    public init(samples: [BodyCoefficientSample], boundaryPolicy: CoefficientBoundaryPolicy = .reject,
                source: String, version: String) throws {
        guard !samples.isEmpty, !source.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              !version.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else { throw AerodynamicError.invalidParameter }
        var records: [[Double]: AerodynamicCoefficients] = [:]
        for sample in samples {
            let key = [sample.angleOfAttackRadians,sample.sideslipRadians,sample.machNumber,sample.reynoldsNumber]
            guard key.allSatisfy(\.isFinite), abs(key[0]) <= .pi, abs(key[1]) <= .pi/2,
                  key[2] >= 0, key[3] > 0 else { throw AerodynamicError.invalidCoefficients }
            guard records[key] == nil else { throw AerodynamicError.duplicateCondition }
            records[key] = sample.coefficients
        }
        let grids = (0..<4).map { coordinate in Array(Set(records.keys.map { $0[coordinate] })).sorted() }
        var expected = 1
        for grid in grids {
            let product = expected.multipliedReportingOverflow(by: grid.count)
            guard !product.overflow, product.partialValue <= samples.count else { throw AerodynamicError.incompleteGrid }
            expected = product.partialValue
        }
        guard expected == records.count else { throw AerodynamicError.incompleteGrid }
        self.samples = samples; self.boundaryPolicy = boundaryPolicy; self.source = source; self.version = version
        axes = grids; values = records
    }
    public func coefficients(for query: BodyAerodynamicQuery) throws -> AerodynamicCoefficients {
        guard let alpha = query.angleOfAttackRadians, let beta = query.sideslipRadians else { throw AerodynamicError.unsupportedQuery }
        guard let mach = query.machNumber, let re = try query.reynoldsNumber() else { throw AerodynamicError.missingFlowProperty }
        return try coefficients(angleOfAttackRadians: alpha,sideslipRadians: beta,machNumber: mach,reynoldsNumber: re)
    }
    public func coefficients(angleOfAttackRadians: Double, sideslipRadians: Double,
                             machNumber: Double, reynoldsNumber: Double) throws -> AerodynamicCoefficients {
        let coordinates = [angleOfAttackRadians,sideslipRadians,machNumber,reynoldsNumber]
        guard coordinates.allSatisfy(\.isFinite), abs(angleOfAttackRadians) <= .pi,
              abs(sideslipRadians) <= .pi/2, machNumber >= 0, reynoldsNumber > 0 else { throw AerodynamicError.invalidParameter }
        let brackets = try zip(coordinates,axes).map { try bracket($0,grid: $1) }
        func interpolate(_ dimension: Int, key: [Double]) throws -> AerodynamicCoefficients {
            if dimension == 4 { return values[key]! }
            let (lower,upper,weight) = brackets[dimension]
            let first = try interpolate(dimension+1,key: key+[lower])
            guard lower != upper else { return first }
            let last = try interpolate(dimension+1,key: key+[upper])
            return try AerodynamicCoefficients(force: first.force*(1-weight)+last.force*weight,
                                               moment: first.moment*(1-weight)+last.moment*weight)
        }
        return try interpolate(0,key: [])
    }
    private func bracket(_ value: Double, grid: [Double]) throws -> (Double,Double,Double) {
        guard boundaryPolicy == .clamp || (grid.first!...grid.last!).contains(value) else { throw AerodynamicError.outsideCoverage }
        let v = min(grid.last!,max(grid.first!,value)), lower = grid.lastIndex(where: { $0 <= v })!
        guard lower+1 < grid.count, v != grid[lower] else { return (grid[lower],grid[lower],0) }
        return (grid[lower],grid[lower+1],(v-grid[lower])/(grid[lower+1]-grid[lower]))
    }
    private enum CodingKeys: String, CodingKey { case samples, boundaryPolicy, source, version }
    public init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        try self.init(samples: c.decode([BodyCoefficientSample].self,forKey: .samples),
            boundaryPolicy: c.decode(CoefficientBoundaryPolicy.self,forKey: .boundaryPolicy),
            source: c.decode(String.self,forKey: .source),version: c.decode(String.self,forKey: .version))
    }
    public func encode(to encoder: Encoder) throws {
        var c = encoder.container(keyedBy: CodingKeys.self)
        try c.encode(samples,forKey: .samples); try c.encode(boundaryPolicy,forKey: .boundaryPolicy)
        try c.encode(source,forKey: .source); try c.encode(version,forKey: .version)
    }
}

/// JSON I/O is separate from coefficient evaluation; decoding revalidates every sample/grid.
public enum BodyCoefficientTableLoader {
    public static func json(_ data: Data) throws -> TabulatedBodyCoefficients {
        try JSONDecoder().decode(TabulatedBodyCoefficients.self,from: data)
    }
    public static func jsonFile(_ url: URL) throws -> TabulatedBodyCoefficients { try json(Data(contentsOf: url)) }
}
