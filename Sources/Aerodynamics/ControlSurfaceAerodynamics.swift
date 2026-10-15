// SPDX-License-Identifier: Apache-2.0
import FluidDynamics
import Foundation

/// A checked incidence/Re/Mach dataset at one supplied geometric deflection.
public struct ControlSurfaceSetting: Sendable, Codable {
    public let deflectionRadians: Double
    public let airfoil: TabulatedAirfoil
    public init(deflectionRadians: Double, airfoil: TabulatedAirfoil) throws {
        guard deflectionRadians.isFinite, abs(deflectionRadians) <= .pi else { throw FluidError.invalidParameter }
        self.deflectionRadians = deflectionRadians; self.airfoil = airfoil
    }
    private enum CodingKeys: String, CodingKey { case deflectionRadians, airfoil }
    public init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        try self.init(deflectionRadians: c.decode(Double.self,forKey: .deflectionRadians),
                      airfoil: c.decode(TabulatedAirfoil.self,forKey: .airfoil))
    }
}

/// Linear interpolation across geometric control deflection, using checked polars
/// for incidence/Re/Mach. No incidence offset, thin-airfoil effectiveness factor,
/// actuator response, hinge moment or dynamic-stall behavior is inferred.
public struct TabulatedControlSurface: WingSectionCoefficientModel, Sendable, Codable {
    public let settings: [ControlSurfaceSetting]
    public let deflectionBoundaryPolicy: CoefficientBoundaryPolicy
    public let source: String
    public let version: String
    public let deflectionConvention: String
    public let momentReferenceDescription: String
    public var deflectionRangeRadians: ClosedRange<Double> {
        settings.first!.deflectionRadians...settings.last!.deflectionRadians
    }
    public init(settings: [ControlSurfaceSetting], deflectionBoundaryPolicy: CoefficientBoundaryPolicy = .reject,
                source: String, version: String, deflectionConvention: String,
                momentReferenceDescription: String) throws {
        guard !settings.isEmpty, [source,version,deflectionConvention,momentReferenceDescription]
            .allSatisfy({ !$0.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty }) else { throw FluidError.invalidParameter }
        let ordered = settings.sorted { $0.deflectionRadians < $1.deflectionRadians }
        guard zip(ordered,ordered.dropFirst()).allSatisfy({ $0.deflectionRadians < $1.deflectionRadians }) else {
            throw FluidError.duplicateCondition
        }
        // Keep the Re/Mach coordinate grid consistent across deflections. Angle
        // sampling can differ; interpolation uses the coverage of both neighbors.
        let reference = ordered[0].airfoil.polars
        let re = Set(reference.map(\.reynoldsNumber)), mach = Set(reference.map(\.machNumber))
        guard ordered.allSatisfy({ Set($0.airfoil.polars.map(\.reynoldsNumber)) == re
            && Set($0.airfoil.polars.map(\.machNumber)) == mach }) else { throw FluidError.incompleteGrid }
        self.settings = ordered; self.deflectionBoundaryPolicy = deflectionBoundaryPolicy
        self.source = source; self.version = version; self.deflectionConvention = deflectionConvention
        self.momentReferenceDescription = momentReferenceDescription
    }
    public func coefficients(for query: WingSectionQuery) throws -> AirfoilCoefficients {
        guard let alpha = query.angleOfAttackRadians, let mach = query.machNumber,
              let re = try query.reynoldsNumber() else { throw FluidError.missingFlowProperty }
        let deflection = try boundedCoefficientCoordinate(query.controlDeflectionRadians,
            range: deflectionRangeRadians,policy: deflectionBoundaryPolicy)
        let lower = settings.lastIndex { $0.deflectionRadians <= deflection }!
        func evaluate(_ index: Int) throws -> AirfoilCoefficients {
            try settings[index].airfoil.coefficients(angleOfAttackRadians: alpha,reynoldsNumber: re,machNumber: mach)
        }
        let first = try evaluate(lower)
        // Exact settings query only that dataset, never an irrelevant neighbor.
        guard settings[lower].deflectionRadians != deflection, lower+1 < settings.count else { return first }
        let next = lower+1
        return try blendCoefficients(first,evaluate(next),weight:
            (deflection-settings[lower].deflectionRadians)/(settings[next].deflectionRadians-settings[lower].deflectionRadians))
    }
    private enum CodingKeys: String, CodingKey {
        case settings, deflectionBoundaryPolicy, source, version, deflectionConvention, momentReferenceDescription
    }
    public init(from decoder: Decoder) throws {
        let c = try decoder.container(keyedBy: CodingKeys.self)
        try self.init(settings: c.decode([ControlSurfaceSetting].self,forKey: .settings),
            deflectionBoundaryPolicy: c.decode(CoefficientBoundaryPolicy.self,forKey: .deflectionBoundaryPolicy),
            source: c.decode(String.self,forKey: .source),version: c.decode(String.self,forKey: .version),
            deflectionConvention: c.decode(String.self,forKey: .deflectionConvention),
            momentReferenceDescription: c.decode(String.self,forKey: .momentReferenceDescription))
    }
}

public enum ControlSurfaceTableLoader {
    public static func load(data: Data) throws -> TabulatedControlSurface {
        do { return try JSONDecoder().decode(TabulatedControlSurface.self,from: data) }
        catch let error as FluidError { throw error }
        catch { throw FluidError.malformedData }
    }
}
