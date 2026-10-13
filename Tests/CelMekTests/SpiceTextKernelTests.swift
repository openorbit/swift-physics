// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
@testable import CelMek

@Suite struct SpiceTextKernelTests {
    func snapshot(_ text: String) throws -> SpiceKernelSnapshot {
        let helper = SpiceKernelSetTests(), dir = try helper.directory()
        defer { try? FileManager.default.removeItem(at: dir) }
        let url = try helper.write(text, "kernel.txt", in: dir)
        let set = try SpiceKernelSet(); try set.load(url)
        return set.snapshot()
    }
    let lsk = """
    KPL/LSK
    \\begindata
    DELTET/DELTA_T_A = 32.184
    DELTET/K = 1.657D-3
    DELTET/EB = 1.671D-2
    DELTET/M = (6.239996 1.99096871D-7)
    DELTET/DELTA_AT = (10 @1972-JAN-1 11 @1972-JUL-1)
    """
    @Test func datesAreUnscaledAndValidated() throws {
        #expect(try spiceDateSeconds("@2000-JAN-1-12:00:00") == 0)
        #expect(try spiceDateSeconds("@2000-JAN-1") == -43200)
        #expect(try spiceDateSeconds("@2000-JAN-2-12:00:00.125") == 86400.125)
        #expect(try spiceDateSeconds("@1972-JAN-1") == -883656000)
        for date in ["@2001-FEB-29", "@2000-JAN-1-24:00:00", "@2000-JAN-1-23:59:60", "@00-JAN-1", "@2000-13-1"] {
            #expect(throws: (any Error).self) { try spiceDateSeconds(date) }
        }
        let pool = try snapshot("KPL/PCK\n\\begindata\nDATES=(1 @2000-JAN-1 @2000-JAN-1-12:00:00)")
        #expect(pool["DATES"]?.value == .numbers([1,-43200,0]))
    }
    @Test func leapSecondsUseExistingUTCConversionAndExplicitExpiry() throws {
        let pool = try snapshot(lsk), data = try pool.leapSecondData()
        let table = try data.table(validUntilModifiedJulianDay: 41683, version: "synthetic")
        #expect(data.entries.map(\.effectiveModifiedJulianDay) == [41317,41499])
        #expect(data.sources.count == 1)
        #expect(try table.dayLength(on: 41498) == 86401)
        let leap = try Epoch(utc: UTCDate(modifiedJulianDay: 41498, secondOfDay: Duration(seconds: 86400)), using: table)
        let midnight = try Epoch(utc: UTCDate(modifiedJulianDay: 41499, secondOfDay: .zero), using: table)
        #expect(try midnight.utc(using: table).secondOfDay == .zero)
        #expect(try leap.utc(using: table).secondOfDay.seconds == 86400)
        #expect(throws: UTCError.outsideCoverage) { try table.taiMinusUTC(on: 41683) }
        for malformed in [lsk.replacingOccurrences(of: "11 @1972", with: "13 @1972"),
                          lsk.replacingOccurrences(of: "32.184", with: "32.0"),
                          lsk.replacingOccurrences(of: "@1972-JUL-1", with: "@1972-JUL-1-01:00:00")] {
            #expect(throws: (any Error).self) { try snapshot(malformed).leapSecondData() }
        }
    }
    @Test func pckRetainsOrientationAndConvertsRadii() throws {
        let pool = try snapshot("KPL/PCK\n\\begindata\nBODY399_RADII=(6378 6378 6357)\nBODY399_POLE_RA=(0 -.641 0)\nBODY399_POLE_DEC=(90 -.557 0)\nBODY399_PM=(190.16 360.9856235 0)\nBODY399_NUT_PREC_RA=(1 2)")
        #expect(try pool.bodyRadiiMetres(bodyID: 399) == [6378000,6378000,6357000])
        let parameters = try pool.bodyOrientationParameters(bodyID: 399)
        #expect(parameters.primeMeridian == [190.16,360.9856235,0])
        #expect(parameters.additionalParameters["NUT_PREC_RA"] == [1,2])
        #expect(parameters.sources.count == 1)
        #expect(throws: SpiceTextKernelError.self) { try pool.bodyRadiiMetres(bodyID: 301) }
        #expect(throws: SpiceTextKernelError.self) { try snapshot("KPL/PCK\n\\begindata\nBODY1_RADII=(1 0 2)").bodyRadiiMetres(bodyID: 1) }
    }
    @Test func fixedMatrixDirectionAndUnsupportedDefinitions() throws {
        let fk = "KPL/FK\n\\begindata\nFRAME_123_NAME='TEST'\nFRAME_123_CLASS=4\nFRAME_123_CLASS_ID=123\nFRAME_123_CENTER=-10\nTKFRAME_123_RELATIVE='J2000'\nTKFRAME_123_SPEC='MATRIX'\nTKFRAME_123_MATRIX=(0 1 0 -1 0 0 0 0 1)"
        let frame = try snapshot(fk).fixedFrameDefinition(frameID: 123)
        // First column maps +X to relative +Y; no transpose or handedness flip.
        #expect(Array(frame.matrix[0..<3]) == [0,1,0])
        #expect(frame.centreID == -10)
        #expect(frame.relativeFrameName == "J2000")
        #expect(throws: SpiceTextKernelError.unsupportedFrameClass(5)) { try snapshot(fk.replacingOccurrences(of: "CLASS=4", with: "CLASS=5")).fixedFrameDefinition(frameID: 123) }
        #expect(throws: SpiceTextKernelError.unsupportedFrameRepresentation("UNKNOWN")) { try snapshot(fk.replacingOccurrences(of: "'MATRIX'", with: "'UNKNOWN'")).fixedFrameDefinition(frameID: 123) }
        for bad in ["(1 0 0 0 1 0 0 0 -1)", "(2 0 0 0 1 0 0 0 1)", "(1 0)"] {
            #expect(throws: SpiceTextKernelError.self) { try snapshot(fk.replacingOccurrences(of: "(0 1 0 -1 0 0 0 0 1)", with: bad)).fixedFrameDefinition(frameID: 123) }
        }
    }
}
