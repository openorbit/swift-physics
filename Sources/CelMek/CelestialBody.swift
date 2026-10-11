// SPDX-License-Identifier: Apache-2.0
import Foundation

public enum CelestialBodyError: Error, Equatable, Sendable {
  case invalidIdentifier, missingProvenance, invalidPhysicalValue, invalidDefinition
  case duplicateIdentifier, unknownBody
}

/// Namespaced identity, independent of display names and model availability.
/// NAIF codes preserve the distinction between system barycentres and body centres.
public enum BodyIdentifier: Hashable, Sendable {
  case naif(Int32)
  case custom(namespace: String, name: String)

  public static let solarSystemBarycentre = Self.naif(0)
  public static let mercuryBarycentre = Self.naif(1)
  public static let venusBarycentre = Self.naif(2)
  public static let earthMoonBarycentre = Self.naif(3)
  public static let marsBarycentre = Self.naif(4)
  public static let jupiterBarycentre = Self.naif(5)
  public static let saturnBarycentre = Self.naif(6)
  public static let uranusBarycentre = Self.naif(7)
  public static let neptuneBarycentre = Self.naif(8)
  public static let plutoBarycentre = Self.naif(9)
  public static let sun = Self.naif(10)
  public static let mercury = Self.naif(199)
  public static let venus = Self.naif(299)
  public static let moon = Self.naif(301)
  public static let earth = Self.naif(399)
  public static let mars = Self.naif(499)
  public static let jupiter = Self.naif(599)
  public static let saturn = Self.naif(699)
  public static let uranus = Self.naif(799)
  public static let neptune = Self.naif(899)
  public static let pluto = Self.naif(999)

  internal func validate() throws {
    if case .custom(let namespace, let name) = self {
      guard !namespace.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
            !name.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else {
        throw CelestialBodyError.invalidIdentifier
      }
    }
  }
}

/// Caller-owned source/version metadata. No implicit current-data lookup occurs.
public struct PhysicalDataProvenance: Hashable, Sendable {
  public let source: String
  public let version: String
  public let reference: String?
  public init(source: String, version: String, reference: String? = nil) throws {
    guard !source.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
          !version.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else {
      throw CelestialBodyError.missingProvenance
    }
    self.source = source; self.version = version; self.reference = reference
  }
}

/// Positive finite physical quantity in the SI units specified by its containing field.
public struct SourcedPhysicalValue: Hashable, Sendable {
  public let value: Double
  public let provenance: PhysicalDataProvenance
  public init(value: Double, provenance: PhysicalDataProvenance) throws {
    guard value.isFinite, value > 0 else { throw CelestialBodyError.invalidPhysicalValue }
    self.value = value; self.provenance = provenance
  }
}

/// Coordinate-time convention of a GM. SI unit conversion alone does not convert
/// between relativistic scale conventions; callers must select compatible dynamics.
public enum GravitationalParameterConvention: Hashable, Sendable {
  case tdbCompatible, ttCompatible, tcbCompatible
  /// Caller-defined Newtonian model, without a relativistic coordinate-time scaling.
  case newtonian
}
public struct GravitationalParameter: Hashable, Sendable {
  /// μ = GM, in m³/s². This is not a mass in kg.
  public let cubicMetresPerSecondSquared: Double
  public let convention: GravitationalParameterConvention
  public let provenance: PhysicalDataProvenance
  public init(cubicMetresPerSecondSquared: Double, convention: GravitationalParameterConvention,
              provenance: PhysicalDataProvenance) throws {
    _ = try SourcedPhysicalValue(value: cubicMetresPerSecondSquared,provenance: provenance)
    self.cubicMetresPerSecondSquared = cubicMetresPerSecondSquared
    self.convention = convention; self.provenance = provenance
  }
}

/// Geometric reference radii along caller-defined body-fixed X/Y/Z axes, in metres.
/// Metadata only: no surface, shape transform, or orientation provider is implied.
public struct BodyReferenceRadii: Hashable, Sendable {
  public let metres: SIMD3<Double>
  public let axesDescription: String
  public let provenance: PhysicalDataProvenance
  public init(metres: SIMD3<Double>, axesDescription: String, provenance: PhysicalDataProvenance) throws {
    guard [metres.x,metres.y,metres.z].allSatisfy({ $0.isFinite && $0 > 0 }),
          !axesDescription.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else {
      throw CelestialBodyError.invalidPhysicalValue
    }
    self.metres = metres; self.axesDescription = axesDescription; self.provenance = provenance
  }
}

/// Missing properties are unknown, never zero or inferred from other fields.
public struct BodyPhysicalProperties: Hashable, Sendable {
  public let gravitationalParameter: GravitationalParameter?
  /// Independently sourced mass in kilograms. No assumed μ/G conversion.
  public let massKilograms: SourcedPhysicalValue?
  public let referenceRadii: BodyReferenceRadii?
  public init(gravitationalParameter: GravitationalParameter? = nil,
              massKilograms: SourcedPhysicalValue? = nil, referenceRadii: BodyReferenceRadii? = nil) {
    self.gravitationalParameter = gravitationalParameter
    self.massKilograms = massKilograms; self.referenceRadii = referenceRadii
  }
}

public enum CelestialBodyKind: Hashable, Sendable {
  case star, planet, dwarfPlanet, naturalSatellite, minorBody, barycentre, artificial, custom
}

/// Immutable body metadata, separate from ephemeris, orientation, shape and gravity models.
public struct CelestialBody: Hashable, Sendable {
  public let identifier: BodyIdentifier
  public let name: String
  public let kind: CelestialBodyKind
  public let properties: BodyPhysicalProperties
  public let provenance: PhysicalDataProvenance
  public init(identifier: BodyIdentifier, name: String, kind: CelestialBodyKind,
              properties: BodyPhysicalProperties = BodyPhysicalProperties(),
              provenance: PhysicalDataProvenance) throws {
    try identifier.validate()
    guard !name.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
          kind != .barycentre || properties.referenceRadii == nil else {
      throw CelestialBodyError.invalidDefinition
    }
    self.identifier = identifier; self.name = name; self.kind = kind
    self.properties = properties; self.provenance = provenance
  }
}

/// Immutable injectable catalogue. Duplicate IDs fail; display names need not be unique.
public struct BodyCatalogue: Sendable {
  public let bodies: [CelestialBody]
  private let byIdentifier: [BodyIdentifier: CelestialBody]
  public init(bodies: [CelestialBody]) throws {
    var lookup: [BodyIdentifier: CelestialBody] = [:]
    for body in bodies {
      guard lookup.updateValue(body,forKey: body.identifier) == nil else { throw CelestialBodyError.duplicateIdentifier }
    }
    self.bodies = bodies; byIdentifier = lookup
  }
  public func body(_ identifier: BodyIdentifier) throws -> CelestialBody {
    guard let body = byIdentifier[identifier] else { throw CelestialBodyError.unknownBody }
    return body
  }
}
