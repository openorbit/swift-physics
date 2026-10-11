import Foundation
import Testing
@testable import CelMek

@Suite struct CelestialBodyTests {
  private func provenance() throws -> PhysicalDataProvenance {
    try PhysicalDataProvenance(source: "game configuration",version: "1")
  }
  @Test func catalogueMatchesIndependentKernelSubset() throws {
    struct Reference: Decodable { let naif: Int32; let km3PerSecond2: Double }
    let url = try #require(Bundle.module.url(forResource: "jpl-gm-2022",withExtension: "json",subdirectory: "Fixtures"))
    let rows = try JSONDecoder().decode([Reference].self,from: Data(contentsOf: url))
    let catalogue = BodyCatalogue.jplGM2022
    #expect(rows.count == 20)
    #expect(catalogue.bodies.count == 21)
    for row in rows {
      let body = try catalogue.body(.naif(row.naif))
      let gm = try #require(body.properties.gravitationalParameter)
      #expect(gm.cubicMetresPerSecondSquared == row.km3PerSecond2*1e9)
      #expect(gm.convention == .tdbCompatible)
      #expect(gm.provenance.version.contains("2022-12-14"))
      #expect(body.properties.massKilograms == nil)
      #expect(body.properties.referenceRadii == nil)
    }
    #expect(try catalogue.body(.solarSystemBarycentre).properties.gravitationalParameter == nil)
    let earth = try #require(catalogue.body(.earth).properties.gravitationalParameter)
    let moon = try #require(catalogue.body(.moon).properties.gravitationalParameter)
    let emb = try #require(catalogue.body(.earthMoonBarycentre).properties.gravitationalParameter)
    #expect(abs(earth.cubicMetresPerSecondSquared+moon.cubicMetresPerSecondSquared-emb.cubicMetresPerSecondSquared) < 0.125)
    let jupiter = try catalogue.body(.jupiter), system = try catalogue.body(.jupiterBarycentre)
    #expect(jupiter.identifier != system.identifier)
    #expect(jupiter.kind == .planet && system.kind == .barycentre)
    #expect(jupiter.properties.gravitationalParameter!.cubicMetresPerSecondSquared < system.properties.gravitationalParameter!.cubicMetresPerSecondSquared)
    #expect(jupiter.properties.gravitationalParameter!.provenance.source.contains("satellite"))
    #expect(system.properties.gravitationalParameter!.provenance.source.contains("DE440"))
  }
  @Test func injectableCustomBodyAndNamespaceIdentity() throws {
    let source = try provenance()
    let mu = try GravitationalParameter(cubicMetresPerSecondSquared: 1e12,convention: .newtonian,provenance: source)
    let radii = try BodyReferenceRadii(metres: SIMD3(1000,900,800),axesDescription: "game body-fixed XYZ",provenance: source)
    let mass = try SourcedPhysicalValue(value: 2e20,provenance: source)
    let properties = BodyPhysicalProperties(gravitationalParameter: mu,massKilograms: mass,referenceRadii: radii)
    let first = try CelestialBody(identifier: .custom(namespace: "world-a",name: "terra"),name: "Terra",kind: .custom,
      properties: properties,provenance: source)
    let second = try CelestialBody(identifier: .custom(namespace: "world-b",name: "terra"),name: "Terra",kind: .custom,provenance: source)
    let catalogue = try BodyCatalogue(bodies: [first,second])
    #expect(try catalogue.body(first.identifier).properties.referenceRadii?.metres == SIMD3(1000,900,800))
    #expect(try catalogue.body(second.identifier).properties.massKilograms == nil)
    #expect(first.identifier != .earth)
    #expect(throws: CelestialBodyError.unknownBody) { try catalogue.body(.earth) }
    #expect(throws: CelestialBodyError.duplicateIdentifier) { try BodyCatalogue(bodies: [first,first]) }
  }
  @Test func invalidPropertiesAndMetadataFailExplicitly() throws {
    let source = try provenance()
    for value in [0.0,-1,.nan,.infinity] {
      #expect(throws: CelestialBodyError.invalidPhysicalValue) {
        try GravitationalParameter(cubicMetresPerSecondSquared: value,convention: .newtonian,provenance: source)
      }
      #expect(throws: CelestialBodyError.invalidPhysicalValue) { try SourcedPhysicalValue(value: value,provenance: source) }
    }
    #expect(throws: CelestialBodyError.missingProvenance) { try PhysicalDataProvenance(source: " ",version: "1") }
    #expect(throws: CelestialBodyError.invalidIdentifier) {
      try CelestialBody(identifier: .custom(namespace: "",name: "body"),name: "Body",kind: .custom,provenance: source)
    }
    #expect(throws: CelestialBodyError.invalidPhysicalValue) {
      try BodyReferenceRadii(metres: SIMD3(1,0,1),axesDescription: "XYZ",provenance: source)
    }
    let radii = try BodyReferenceRadii(metres: SIMD3(1,1,1),axesDescription: "XYZ",provenance: source)
    #expect(throws: CelestialBodyError.invalidDefinition) {
      try CelestialBody(identifier: .earthMoonBarycentre,name: "EMB",kind: .barycentre,
        properties: BodyPhysicalProperties(referenceRadii: radii),provenance: source)
    }
    #expect(throws: CelestialBodyError.invalidDefinition) {
      try CelestialBody(identifier: .earth,name: "",kind: .planet,provenance: source)
    }
  }
}
