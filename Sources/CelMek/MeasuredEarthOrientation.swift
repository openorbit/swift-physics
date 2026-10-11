// SPDX-License-Identifier: Apache-2.0
import Foundation

/// Measured EOP at a UTC label. All angles are radians; DUT1 is UT1−UTC in seconds.
public struct EarthOrientationRecord {
  public let utc: UTCDate
  public let dut1: Duration
  public let polarMotionX: Double
  public let polarMotionY: Double
  public let celestialPoleX: Double
  public let celestialPoleY: Double
  public let poleReference: CelestialPoleReference
  public init(utc: UTCDate, dut1: Duration, polarMotionX: Double, polarMotionY: Double,
              celestialPoleX: Double = 0, celestialPoleY: Double = 0,
              poleReference: CelestialPoleReference = .iau2006A) {
    self.utc = utc; self.dut1 = dut1
    self.polarMotionX = polarMotionX; self.polarMotionY = polarMotionY
    self.celestialPoleX = celestialPoleX; self.celestialPoleY = celestialPoleY
    self.poleReference = poleReference
  }
  public func sample(using leapSeconds: LeapSecondTable) throws -> EarthOrientationSample {
    let epoch = try leapSeconds.epoch(from: utc)
    let offset = try leapSeconds.taiMinusUTC(on: utc.modifiedJulianDay)
    return try EarthOrientationSample(epoch: epoch,
      ut1MinusTT: dut1.subtracting(Duration(wholeSeconds: Int64(offset)+32, fractionalSecond: 0.184)),
      polarMotionX: polarMotionX, polarMotionY: polarMotionY,
      celestialPoleX: celestialPoleX, celestialPoleY: celestialPoleY, poleReference: poleReference)
  }
}

extension EarthOrientationSample {
  /// Reconstruct UT1−UTC with the caller's versioned leap-second table.
  public func dut1(using leapSeconds: LeapSecondTable) throws -> Duration {
    let utc = try epoch.utc(using: leapSeconds)
    let offset = try leapSeconds.taiMinusUTC(on: utc.modifiedJulianDay)
    return try ut1MinusTT.adding(Duration(wholeSeconds: Int64(offset)+32, fractionalSecond: 0.184))
  }
}

/// Two-part UT1 Julian date. UT1 is a measured scale rather than a TimeScale enum case.
public struct UT1JulianDateComponents: Sendable {
  public let wholeDay: Int64
  public let fractionalDay: Double
}
extension Epoch {
  public func ut1JulianDate(using earthOrientation: any EarthOrientationProvider) throws -> UT1JulianDateComponents {
    let sample = try earthOrientation.sample(at: self)
    guard sample.epoch == self else { throw FrameTransformError.epochMismatch }
    let shifted = try advanced(by: sample.ut1MinusTT).julianDate(in: .tt)
    return UT1JulianDateComponents(wholeDay: shifted.wholeDay, fractionalDay: shifted.fractionalDay)
  }
}

/// Explicit interpolation on atomic seconds. UT1−UTC jumps are removed before interpolation.
/// Natural cubic is C2 within coverage; linear interpolation has derivative discontinuities at knots.
/// Neither policy extrapolates. Natural endpoint curvature is an interpolation assumption, not LOD data.
public struct TabulatedEarthOrientationProvider: EarthOrientationProvider {
  public enum Interpolation: Sendable { case linear, naturalCubic }
  public let provenance: String
  public let coverage: ClosedRange<Epoch>
  public let interpolation: Interpolation
  private let epochs: [Epoch]
  private let times: [Double]
  private let values: [[Double]]
  private let secondDerivatives: [[Double]]
  private let reference: CelestialPoleReference

  public init(records: [EarthOrientationRecord], leapSeconds: LeapSecondTable,
              interpolation: Interpolation, maximumGap: Duration, provenance: String) throws {
    guard records.count >= 2, maximumGap.seconds > 0,
          !provenance.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else {
      throw EarthOrientationError.invalidTable
    }
    let samples = try records.map { try $0.sample(using: leapSeconds) }
    let epochs = samples.map(\.epoch)
    let times = try epochs.map { try $0.duration(since: epochs[0]).seconds }
    for i in 1..<times.count {
      guard times[i] > times[i-1], samples[i].poleReference == samples[0].poleReference else {
        throw EarthOrientationError.invalidTable
      }
      guard times[i]-times[i-1] <= maximumGap.seconds else { throw EarthOrientationError.excessiveGap }
    }
    let values = (0..<5).map { field in samples.map { s in
      [s.ut1MinusTT.seconds,s.polarMotionX,s.polarMotionY,s.celestialPoleX,s.celestialPoleY][field]
    }}
    self.epochs = epochs; self.times = times; self.values = values
    self.interpolation = interpolation; self.reference = samples[0].poleReference
    self.coverage = epochs[0]...epochs.last!
    self.provenance = "\(provenance); leap seconds: \(leapSeconds.version); interpolation: \(interpolation)"
    self.secondDerivatives = values.map { v in
      guard interpolation == .naturalCubic, times.count > 2 else { return Array(repeating: 0, count: times.count) }
      var m = Array(repeating: 0.0, count: times.count)
      var diagonal = m, rhs = m
      for i in 1..<(times.count-1) {
        let left = times[i]-times[i-1], right = times[i+1]-times[i]
        let factor = i == 1 ? 0 : left/diagonal[i-1]
        diagonal[i] = 2*(left+right)-factor*(times[i]-times[i-1])
        rhs[i] = 6*((v[i+1]-v[i])/right-(v[i]-v[i-1])/left)-factor*rhs[i-1]
      }
      for i in stride(from: times.count-2, through: 1, by: -1) {
        m[i] = (rhs[i]-(times[i+1]-times[i])*m[i+1])/diagonal[i]
      }
      return m
    }
  }
  public func sample(at epoch: Epoch) throws -> EarthOrientationSample {
    guard coverage.contains(epoch) else { throw EarthOrientationError.outsideCoverage }
    let t = try epoch.duration(since: epochs[0]).seconds
    var low = 0, high = times.count-1
    while high-low > 1 {
      let mid = (low+high)/2
      if times[mid] <= t { low = mid } else { high = mid }
    }
    let h = times[high]-times[low], a = (times[high]-t)/h, b = (t-times[low])/h
    let v = (0..<5).map { f in
      a*values[f][low]+b*values[f][high]
        + ((a*a*a-a)*secondDerivatives[f][low]+(b*b*b-b)*secondDerivatives[f][high])*h*h/6
    }
    return try EarthOrientationSample(epoch: epoch, ut1MinusTT: Duration(seconds: v[0]),
      polarMotionX: v[1], polarMotionY: v[2], celestialPoleX: v[3], celestialPoleY: v[4], poleReference: reference)
  }
}

/// USNO/IERS finals2000A Bulletin A columns. Offsets refer to IAU 2000A; file angles
/// are arcseconds for polar motion and milliarcseconds for dX/dY. Bulletin B and errors
/// are not consumed. Blank forecast-tail rows are skipped; incomplete selected rows fail.
public enum Finals2000A {
  public enum PredictionPolicy { case observedOnly, includePredictions }
  public static func records(from text: String, predictions: PredictionPolicy,
                             modifiedJulianDays: ClosedRange<Int64>? = nil) throws -> [EarthOrientationRecord] {
    let arcsecond = Double.pi/(180*3600)
    var result: [EarthOrientationRecord] = []
    for (index, raw) in text.split(separator: "\n", omittingEmptySubsequences: false).enumerated() {
      let row = Array(raw.utf8)
      if row.allSatisfy({ $0 == 32 || $0 == 13 || $0 == 9 }) { continue }
      func field(_ start: Int, _ end: Int) -> String {
        guard row.count >= end else { return "" }
        return String(decoding: row[(start-1)..<end], as: UTF8.self).trimmingCharacters(in: .whitespacesAndNewlines)
      }
      func number(_ start: Int, _ end: Int) throws -> Double {
        guard let v = Double(field(start,end)), v.isFinite else { throw EarthOrientationError.invalidFile(line: index+1) }
        return v
      }
      let mjd = try number(8,15)
      guard mjd.rounded(.towardZero) == mjd, mjd >= Double(Int64.min), mjd < Double(Int64.max) else {
        throw EarthOrientationError.invalidFile(line: index+1)
      }
      let day = Int64(mjd)
      if let range = modifiedJulianDays, !range.contains(day) { continue }
      let flags = [field(17,17),field(58,58),field(96,96)]
      if flags.allSatisfy(\.isEmpty), [field(19,27),field(38,46),field(59,68),field(98,106),field(117,125)].allSatisfy(\.isEmpty) { continue }
      if predictions == .observedOnly, flags.contains("P") { continue }
      guard flags.allSatisfy({ $0 == "I" || $0 == "P" }) else { throw EarthOrientationError.invalidFile(line: index+1) }
      // Date columns are redundant with MJD; reject corrupted or wrong-format files.
      guard let shortYear = Int(field(1,2)), let month = Int(field(3,4)), let date = Int(field(5,6)) else {
        throw EarthOrientationError.invalidFile(line: index+1)
      }
      guard (0...99).contains(shortYear) else { throw EarthOrientationError.invalidFile(line: index+1) }
      let year = shortYear + (day <= 51543 ? 1900 : 2000)
      let calendar = try GregorianDay(year: year, month: month, day: date)
      guard calendar.modifiedJulianDay == day else { throw EarthOrientationError.invalidFile(line: index+1) }
      result.append(try EarthOrientationRecord(utc: UTCDate(modifiedJulianDay: day, secondOfDay: .zero),
        dut1: Duration(seconds: number(59,68)), polarMotionX: number(19,27)*arcsecond,
        polarMotionY: number(38,46)*arcsecond, celestialPoleX: number(98,106)*arcsecond/1000,
        celestialPoleY: number(117,125)*arcsecond/1000, poleReference: .iau2000A))
    }
    return result
  }
}
