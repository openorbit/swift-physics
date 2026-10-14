// SPDX-License-Identifier: Apache-2.0
import Testing
import Math
import MagneticFields
@testable import CelMek
@Suite struct BodyMagneticFieldTests {
    private struct CustomField: MagneticFieldModel {
        let metadata: MagneticFieldMetadata
        func field(at query: MagneticFieldQuery) throws -> MagneticFluxDensity {
            try MagneticFluxDensity(vectorTesla: SIMD3(query.positionMetres.x,query.decimalYear,0))
        }
    }
    private func metadata() throws -> MagneticFieldMetadata {
        try MagneticFieldMetadata(name: "Test",source: "Synthetic",version: "1",axesDescription: "XYZ",
            assumptions: "Test",radiusRangeMetres: 0...100,decimalYearRange: 2025...2030)
    }
    @Test func harmonicModelsPlugIntoBodyProvider() throws {
        let position = SIMD3(7_000_000.0,0,0)
        let models: [any MagneticFieldModel] = [try WMMMagneticField.wmm2025(), try IGRFMagneticField.igrf14()]
        for model in models {
            let provider = try ModelBodyMagneticField(bodyIdentifier: .earth,model: model,decimalYear: { _ in 2026.5 })
            let expected = try model.checkedField(at: MagneticFieldQuery(positionMetres: position,decimalYear: 2026.5))
            let actual = try provider.magneticField(at: .j2000,displacementMetres: position)
            #expect(actual.vectorTesla == expected.vectorTesla)
            #expect(actual.epoch == .j2000)
        }
    }
    @Test func customProviderCalendarAndBodyValidation() throws {
        let custom = try ModelBodyMagneticField(bodyIdentifier: .earth,model: CustomField(metadata: metadata()),
                                                decimalYear: { _ in 2026.5 })
        let field = try custom.magneticField(at: .j2000,displacementMetres: SIMD3(5,0,0))
        #expect(field.epoch == .j2000)
        #expect(field.vectorTesla == SIMD3(5,2026.5,0))
        let invalidDate = try ModelBodyMagneticField(bodyIdentifier: .earth,model: CustomField(metadata: metadata()),
                                                     decimalYear: { _ in 2031 })
        #expect(throws: MagneticFieldError.outsideCoverage) {
            try invalidDate.magneticField(at: .j2000,displacementMetres: .zero)
        }
        #expect(BodyIdentifier.venus.magneticModel == nil)
        #expect(throws: BodyModelError.bodyMismatch) {
            try ModelBodyMagneticField(bodyIdentifier: .mars,model: BodyMagneticField.model(for: .earth),decimalYear: { _ in 2026 })
        }
        let earth = try BodyCatalogue.jplGM2022.body(.earth)
        #expect(try CelestialBodyModels(body: earth,magneticField: custom).magneticField != nil)
        #expect(throws: BodyModelError.bodyMismatch) {
            try CelestialBodyModels(body: BodyCatalogue.jplGM2022.body(.mars),magneticField: custom)
        }
    }
    @Test func rotatesPositionAndVectorAndChecksIdentity() throws {
        let frame = Frame(name: "Magnetic reference",pseudoInertial: true)
        let provenance = try PhysicalDataProvenance(source: "Test",version: "1")
        let orientation = try IdentifiedBodyOrientationProvider(bodyIdentifier: .earth,
            provider: FixedAxisOrientationProvider(referenceFrame: frame,axis: SIMD3(0,0,1),rate: 0,
                initialRotation: Quaternion(angle: .pi/2,axis: SIMD3(0,0,1))),provenance: provenance)
        let model = try ModelBodyMagneticField(bodyIdentifier: .earth,model: CustomField(metadata: metadata()),
                                               decimalYear: { _ in 2026 })
        let oriented = try OrientedBodyMagneticField(provider: model,orientation: orientation)
        let field = try oriented.magneticField(at: .j2000,displacementMetres: SIMD3(0,10,0),in: frame)
        // Reference position +Y becomes body -X. Body field (-10,2026,0) becomes (2026,10,0).
        #expect(abs(field.fieldTeslas.x-2026) < 1e-10)
        #expect(abs(field.fieldTeslas.y-10) < 1e-10)
        #expect(field.frame === frame)
        #expect(field.origin == .earth)
        let spacecraft = SpacecraftState(cartesian: try CartesianState(epoch: .j2000,frame: frame,
            position: SIMD3(0,10,0),velocity: .zero),mass: try SpacecraftMassProperties(dryMassKilograms: 1))
        let provider: any MagneticFieldProvider = oriented
        #expect(try provider.field(at: spacecraft).fieldTeslas == field.fieldTeslas)
        #expect(throws: BodyFixedFrameError.frameMismatch) {
            try oriented.magneticField(at: .j2000,displacementMetres: .zero,in: Frame(name: "Other",pseudoInertial: true))
        }
        let other = try IdentifiedBodyOrientationProvider(bodyIdentifier: .mars,provider: orientation,provenance: provenance)
        #expect(throws: BodyModelError.bodyMismatch) { try OrientedBodyMagneticField(provider: model,orientation: other) }
    }
}
