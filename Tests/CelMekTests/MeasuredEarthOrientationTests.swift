import Foundation
import Testing
@testable import CelMek

@Suite struct MeasuredEarthOrientationTests {
  private struct Reference: Decodable {
    let mjd: Double, old: Bool, dx: Double, dy: Double, matrix: [Double]
  }
  private struct Provider: EarthOrientationProvider {
    let dx: Double, dy: Double
    let reference: CelestialPoleReference
    let provenance = "independent reference test"
    func sample(at epoch: Epoch) throws -> EarthOrientationSample {
      try EarthOrientationSample(epoch: epoch, ut1MinusTT: .zero, polarMotionX: 0, polarMotionY: 0,
        celestialPoleX: dx, celestialPoleY: dy, poleReference: reference)
    }
  }
  @Test func observedPolesMatchIndependentERFAAndBothFrameChains() throws {
    let url = try #require(Bundle.module.url(forResource: "observed-pole-erfa", withExtension: "json", subdirectory: "Fixtures"))
    let rows = try JSONDecoder().decode([Reference].self, from: Data(contentsOf: url))
    #expect(rows.count == 36)
    for row in rows {
      let epoch = try Epoch(julianDate: row.mjd+2400000.5, scale: .tt)
      let data = Provider(dx: row.dx, dy: row.dy, reference: row.old ? .iau2000A : .iau2006A)
      func matrix(_ edge: EarthFrameTransformProvider.Edge) throws -> [Double] {
        try EarthFrameTransformProvider(edge: edge, earthOrientation: data).matrix(at: epoch)
      }
      let cirf = try matrix(.celestialIntermediate), tod = try matrix(.trueOfDate)
      let era = try matrix(.earthRotation), gast = try matrix(.apparentSidereal)
      for i in 0..<9 { #expect(abs(cirf[i]-row.matrix[i]) < 2e-14) }
      // Every terrestrial coordinate must agree through the CIO and equinox routes.
      for i in 0..<3 { for j in 0..<3 {
        let a = (0..<3).reduce(0.0) { $0+era[i*3+$1]*cirf[$1*3+j] }
        let b = (0..<3).reduce(0.0) { $0+gast[i*3+$1]*tod[$1*3+j] }
        #expect(abs(a-b) < 3e-14)
      }}
      for i in 6..<9 { #expect(abs(tod[i]-cirf[i]) < 3e-16) }
    }
  }
  private func record(_ day: Int64, dut1: Double, x: Double = 0) throws -> EarthOrientationRecord {
    try EarthOrientationRecord(utc: UTCDate(modifiedJulianDay: day, secondOfDay: .zero),
      dut1: CelMek.Duration(seconds: dut1), polarMotionX: x, polarMotionY: 0)
  }
  @Test func leapSecondInterpolationAndUT1Conversion() throws {
    let leap = LeapSecondTable.iersBulletinC72
    let records = try [record(57753,dut1: -0.4), record(57754,dut1: 0.6), record(57755,dut1: 0.6)]
    for method in [TabulatedEarthOrientationProvider.Interpolation.linear, .naturalCubic] {
      let table = try TabulatedEarthOrientationProvider(records: records, leapSeconds: leap,
        interpolation: method, maximumGap: CelMek.Duration(seconds: 86401), provenance: "synthetic leap case")
      let utc = try UTCDate(modifiedJulianDay: 57753, secondOfDay: CelMek.Duration(seconds: 86400.5))
      let epoch = try Epoch(utc: utc, using: leap)
      let sample = try table.sample(at: epoch)
      #expect(abs(sample.ut1MinusTT.seconds + 68.584) < 2e-14)
      #expect(abs(try sample.dut1(using: leap).seconds + 0.4) < 2e-14)
      let next = try table.sample(at: records[1].sample(using: leap).epoch)
      #expect(abs(try next.dut1(using: leap).seconds - 0.6) < 2e-14)
      let jd = try epoch.ut1JulianDate(using: table)
      #expect(jd.wholeDay == 2457754)
      #expect(abs(jd.fractionalDay - (0.5+0.1/86400)) < 2e-16)
      #expect(throws: EarthOrientationError.outsideCoverage) {
        try table.sample(at: table.coverage.lowerBound.advanced(by: CelMek.Duration(seconds: -1)))
      }
    }
  }
  @Test func tableValidationAndCubicInterpolation() throws {
    let leap = LeapSecondTable.iersBulletinC72
    let rows = try [record(60000,dut1: 0,x: 0), record(60001,dut1: 0,x: 1e-6), record(60002,dut1: 0,x: 0)]
    let table = try TabulatedEarthOrientationProvider(records: rows, leapSeconds: leap,
      interpolation: .naturalCubic, maximumGap: CelMek.Duration(seconds: 86400), provenance: "cubic")
    let epoch = try table.coverage.lowerBound.advanced(by: CelMek.Duration(seconds: 43200))
    #expect(abs(try table.sample(at: epoch).polarMotionX - 0.6875e-6) < 1e-20)
    #expect(throws: EarthOrientationError.invalidTable) {
      try TabulatedEarthOrientationProvider(records: [rows[0],rows[0]], leapSeconds: leap,
        interpolation: .linear, maximumGap: CelMek.Duration(seconds: 86401), provenance: "duplicates")
    }
    #expect(throws: EarthOrientationError.excessiveGap) {
      try TabulatedEarthOrientationProvider(records: rows, leapSeconds: leap,
        interpolation: .linear, maximumGap: CelMek.Duration(seconds: 1), provenance: "gap")
    }
    #expect(throws: EarthOrientationError.invalidParameters) {
      try EarthOrientationSample(epoch: .j2000, ut1MinusTT: .zero, polarMotionX: 0, polarMotionY: 0, celestialPoleX: .nan)
    }
  }
  private func line(flag: String = "I", dx: String = "0.123", year: String = "23") -> String {
    var bytes = Array(repeating: UInt8(32), count: 185)
    func put(_ column: Int, _ text: String) {
      for (i,b) in text.utf8.enumerated() { bytes[column-1+i] = b }
    }
    put(1,year); put(3," 2"); put(5,"25"); put(8,"60000.00")
    put(17,flag); put(19," 0.100000"); put(38," 0.200000")
    put(58,flag); put(59,"-0.1234567"); put(96,flag)
    put(98,dx); put(117,"-0.456")
    return String(decoding: bytes, as: UTF8.self)
  }
  @Test func finalsUnitsPredictionSelectionAndMalformedRows() throws {
    let rows = try Finals2000A.records(from: line()+"\n"+line(flag: "P"), predictions: .observedOnly)
    #expect(rows.count == 1)
    let r = try #require(rows.first)
    #expect(r.utc.modifiedJulianDay == 60000)
    #expect(r.poleReference == .iau2000A)
    #expect(abs(r.polarMotionX - 0.1 * .pi/(180*3600)) < 1e-20)
    #expect(abs(r.celestialPoleX - 0.123 * .pi/(180*3600*1000)) < 1e-23)
    #expect(abs(r.dut1.seconds + 0.1234567) < 1e-16)
    #expect(try Finals2000A.records(from: line(flag: "P"), predictions: .includePredictions).count == 1)
    #expect(try Finals2000A.records(from: line(), predictions: .observedOnly, modifiedJulianDays: 60001...60002).isEmpty)
    #expect(throws: EarthOrientationError.invalidFile(line: 1)) {
      try Finals2000A.records(from: line(dx: "NaN"), predictions: .observedOnly)
    }
    #expect(throws: EarthOrientationError.invalidFile(line: 1)) {
      try Finals2000A.records(from: line(year: "22"), predictions: .observedOnly)
    }
  }
}
