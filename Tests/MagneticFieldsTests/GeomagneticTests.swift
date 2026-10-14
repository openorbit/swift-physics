// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
@testable import MagneticFields
@Suite struct GeomagneticTests {
    @Test func officialWMM2025Vectors() throws {
        let model = try WMMMagneticField.wmm2025()
        let file = Bundle.module.url(forResource: "WMM2025_TEST_VALUES",withExtension: "txt",subdirectory: "Fixtures")!
        let text = try String(contentsOf: file,encoding: .utf8)
        var count = 0
        for line in text.split(whereSeparator: \.isNewline) where !line.hasPrefix("#") {
            let v = line.split(whereSeparator: \.isWhitespace).map { Double($0)! }
            let query = try GeodeticMagneticQuery(latitudeRadians: v[2] * .pi/180,longitudeRadians: v[3] * .pi/180,
                ellipsoidHeightMetres: v[1]*1000,decimalYear: v[0])
            let result = try model.elements(at: query), ned = result.northEastDownTesla*1e9
            #expect(abs(ned.x-v[4]) < 0.11)
            #expect(abs(ned.y-v[5]) < 0.11)
            #expect(abs(ned.z-v[6]) < 0.11)
            #expect(abs(result.totalIntensityTesla*1e9-v[8]) < 0.11)
            count += 1
        }
        #expect(count > 10)
    }
    @Test func officialIGRF14ReferenceVectors() throws {
        let model = try IGRFMagneticField.igrf14()
        let file = Bundle.module.url(forResource: "IGRF14_REFERENCE",withExtension: "txt",subdirectory: "Fixtures")!
        let text = try String(contentsOf: file,encoding: .utf8)
        var count = 0
        for line in text.split(whereSeparator: \.isNewline) where !line.hasPrefix("#") {
            let v = line.split(whereSeparator: \.isWhitespace).map { Double($0)! }
            let result = try model.elements(at: GeodeticMagneticQuery(latitudeRadians: v[2] * .pi/180,
                longitudeRadians: v[3] * .pi/180,ellipsoidHeightMetres: v[1]*1000,decimalYear: v[0]))
            let ned = result.northEastDownTesla*1e9
            #expect(abs(ned.x-v[4]) < 0.2)
            #expect(abs(ned.y-v[5]) < 0.2)
            #expect(abs(ned.z-v[6]) < 0.2)
            #expect(abs(result.totalIntensityTesla*1e9-v[7]) < 0.2)
            count += 1
        }
        #expect(count == 48)
    }
    @Test func igrfInterpolationAndSecularVariation() throws {
        let text = """
        g/h n m 2000.0 2005.0 2005-10
        g 1 0 -30000 -29000 10
        g 1 1 0 0 0
        h 1 1 0 0 0
        """
        let model = try IGRFMagneticField(coefficientText: text)
        for (year,g) in [(2000.0,-30000.0),(2002.5,-29500),(2005,-29000),(2010,-28950)] {
            let field = try model.field(at: MagneticFieldQuery(positionMetres: SIMD3(6371200,0,0),decimalYear: year))
            #expect(abs(field.vectorTesla.z + g*1e-9) < 1e-15)
        }
        let actual = try IGRFMagneticField.igrf14()
        #expect(actual.metadata.decimalYearRange == 1900...2030)
        #expect(try actual.elements(at: GeodeticMagneticQuery(latitudeRadians: 0,longitudeRadians: 0,
            ellipsoidHeightMetres: 0,decimalYear: 2025)).totalIntensityTesla > 1e-5)
        #expect(throws: MagneticFieldError.outsideCoverage) {
            try actual.field(at: MagneticFieldQuery(positionMetres: SIMD3(6371200,0,0),decimalYear: 2031))
        }
    }
    @Test func coefficientFilesAndNewEditions() throws {
        let directory = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: directory,withIntermediateDirectories: true)
        defer { try? FileManager.default.removeItem(at: directory) }
        let file = directory.appendingPathComponent("WMM.COF")
        try "2030 NEXT date\n1 0 -30000 0 10 0\n1 1 0 0 0 0\n9999999999".write(to: file,atomically: true,encoding: .utf8)
        let model = try WMMMagneticField(coefficientFile: file)
        #expect(model.metadata.decimalYearRange == 2030...2035)
        let field = try model.field(at: MagneticFieldQuery(positionMetres: SIMD3(6371200,0,0),decimalYear: 2035))
        #expect(abs(field.vectorTesla.z-29950e-9) < 1e-15)
        let igrfFile = directory.appendingPathComponent("igrf.txt")
        try "g/h n m 2090.0 2095.0 2095-00\ng 1 0 -30000 -29000 10\ng 1 1 0 0 0\nh 1 1 0 0 0".write(to: igrfFile,atomically: true,encoding: .utf8)
        #expect(try IGRFMagneticField(coefficientFile: igrfFile).metadata.decimalYearRange == 2090...2100)
    }
    @Test func coefficientValidationAndPoles() throws {
        let valid = "2025 TEST date\n1 0 -30000 0 0 0\n1 1 1 2 0 0\n9999999999"
        #expect(throws: GeomagneticCoefficientError.incompleteCoefficients) {
            try WMMMagneticField(coefficientText: "2025 TEST date\n1 0 -30000 0 0 0\n9999999999")
        }
        #expect(throws: GeomagneticCoefficientError.duplicateCoefficient) {
            try WMMMagneticField(coefficientText: valid.replacingOccurrences(of: "9999999999",with: "1 0 1 0 0 0\n9999999999"))
        }
        #expect(throws: GeomagneticCoefficientError.invalidHeader) { try IGRFMagneticField(coefficientText: "g/h n m rubbish") }
        let model = try WMMMagneticField(coefficientText: valid)
        for sign in [-1.0,1.0] {
            let exact = try model.field(at: MagneticFieldQuery(positionMetres: SIMD3(0,0,sign*6371200),decimalYear: 2025))
            let near = try model.field(at: MagneticFieldQuery(positionMetres: SIMD3(0.001,0,sign*6371200),decimalYear: 2025))
            #expect(abs(exact.vectorTesla.x-near.vectorTesla.x) < 1e-12)
            #expect(abs(exact.vectorTesla.y-near.vectorTesla.y) < 1e-12)
            #expect(abs(exact.vectorTesla.z-near.vectorTesla.z) < 1e-12)
        }
    }
}
