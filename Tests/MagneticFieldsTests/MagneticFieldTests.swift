// SPDX-License-Identifier: Apache-2.0
import Testing
@testable import MagneticFields
@Suite struct MagneticFieldTests {
    private func metadata(radius: ClosedRange<Double> = 0...100,
                          years: ClosedRange<Double>? = nil) throws -> MagneticFieldMetadata {
        try MagneticFieldMetadata(name: "Test", source: "Synthetic", version: "1", axesDescription: "XYZ",
                                  assumptions: "Test", radiusRangeMetres: radius,decimalYearRange: years)
    }
    private func query(_ p: SIMD3<Double>, year: Double = 2026) throws -> MagneticFieldQuery {
        try MagneticFieldQuery(positionMetres: p,decimalYear: year)
    }
    @Test func analyticDipoleAndInverseCube() throws {
        let model = try DipoleMagneticField(referenceRadiusMetres: 10,momentAxisFieldTesla: SIMD3(0,0,2),metadata: metadata())
        #expect(try model.field(at: query(SIMD3(10,0,0))).vectorTesla == SIMD3(0,0,-2))
        #expect(try model.field(at: query(SIMD3(0,0,10))).vectorTesla == SIMD3(0,0,4))
        #expect(try model.field(at: query(SIMD3(0,0,-10))).vectorTesla == SIMD3(0,0,4))
        #expect(try model.field(at: query(SIMD3(20,0,0))).vectorTesla == SIMD3(0,0,-0.25))
        let tilted = try DipoleMagneticField(referenceRadiusMetres: 10,momentAxisFieldTesla: SIMD3(2,0,0),
                                             offsetMetres: SIMD3(1,2,3),metadata: metadata())
        #expect(try tilted.field(at: query(SIMD3(11,2,3))).vectorTesla == SIMD3(4,0,0))
        #expect(throws: MagneticFieldError.singularField) { try tilted.field(at: query(SIMD3(1,2,3))) }
    }
    @Test func coverageAndValidation() throws {
        let model = try ConstantMagneticField(vectorTesla: SIMD3(1,2,3),metadata: metadata(radius: 10...20,years: 2025...2030))
        #expect(try model.checkedField(at: query(SIMD3(10,0,0),year: 2025)).vectorTesla == SIMD3(1,2,3))
        #expect(throws: MagneticFieldError.outsideCoverage) { try model.checkedField(at: query(SIMD3(21,0,0))) }
        #expect(throws: MagneticFieldError.outsideCoverage) { try model.checkedField(at: query(SIMD3(10,0,0),year: 2031)) }
        #expect(throws: MagneticFieldError.invalidPosition) { try query(SIMD3(.nan,0,0)) }
        #expect(throws: MagneticFieldError.invalidDate) { try query(.zero,year: .infinity) }
        #expect(throws: MagneticFieldError.invalidParameters) {
            try DipoleMagneticField(referenceRadiusMetres: -1,momentAxisFieldTesla: .zero,metadata: metadata())
        }
        #expect(throws: MagneticFieldError.nonFiniteResult) { try MagneticFluxDensity(vectorTesla: SIMD3(.infinity,0,0)) }
    }
    @Test func catalogueStrengthsAndMissingModels() throws {
        for body in MagneticBody.allCases {
            let model = BodyMagneticField.model(for: body), radius = model.dipole.referenceRadiusMetres
            let equator = try model.field(at: query(SIMD3(radius,0,0)))
            let pole = try model.field(at: query(SIMD3(0,0,radius)))
            #expect(abs(pole.magnitudeTesla/equator.magnitudeTesla-2) < 1e-14)
            #expect(model.metadata.decimalYearRange == nil)
        }
        let earth = BodyMagneticField.model(for: .earth)
        let north = try earth.field(at: query(SIMD3(0,0,earth.dipole.referenceRadiusMetres)))
        #expect(abs(north.vectorTesla.z + 61.2e-6) < 1e-15)
        #expect(BodyMagneticField.model(naifBodyCode: 299) == nil)
        #expect(BodyMagneticField.model(naifBodyCode: 3) == nil)
    }
}
