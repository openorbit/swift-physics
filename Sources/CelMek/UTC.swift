// SPDX-License-Identifier: Apache-2.0
import Foundation

public enum UTCError: Error, Equatable, Sendable {
  case invalidTable, outsideCoverage, invalidSecondOfDay, malformedData
}

/// UTC civil day (integer Modified Julian Date at midnight) and labelled seconds
/// since midnight. The leap second is 86400..<86401; it is not the next midnight.
/// Validation of day length belongs to the selected leap-second table.
public struct UTCDate: Hashable, Sendable {
  public let modifiedJulianDay: Int64
  public let secondOfDay: Duration
  public init(modifiedJulianDay: Int64, secondOfDay: Duration) throws {
    guard secondOfDay >= .zero, secondOfDay.wholeSeconds < 86401 else { throw UTCError.invalidSecondOfDay }
    self.modifiedJulianDay = modifiedJulianDay
    self.secondOfDay = secondOfDay
  }
}

/// Immutable, versioned post-1971 UTC history. Coverage ends at an exclusive
/// UTC midnight. Unknown dates never extrapolate the last known offset.
public struct LeapSecondTable: Sendable {
  public struct Entry: Hashable, Sendable {
    public let effectiveModifiedJulianDay: Int64
    public let taiMinusUTC: Int
    public init(effectiveModifiedJulianDay: Int64, taiMinusUTC: Int) {
      self.effectiveModifiedJulianDay = effectiveModifiedJulianDay
      self.taiMinusUTC = taiMinusUTC
    }
  }
  public let version: String
  public let source: String
  public let entries: [Entry]
  public let validUntilModifiedJulianDay: Int64

  public init(entries: [Entry], validUntilModifiedJulianDay: Int64, version: String, source: String) throws {
    guard let first = entries.first, let last = entries.last,
          first.effectiveModifiedJulianDay >= 41317,
          last.effectiveModifiedJulianDay < validUntilModifiedJulianDay,
          !version.isEmpty, !source.isEmpty else { throw UTCError.invalidTable }
    for (index, entry) in entries.enumerated() {
      guard (-86400...86400).contains(entry.taiMinusUTC) else { throw UTCError.invalidTable }
      if index > 0 {
        let previous = entries[index-1]
        guard entry.effectiveModifiedJulianDay > previous.effectiveModifiedJulianDay,
              abs(entry.taiMinusUTC - previous.taiMinusUTC) == 1 else { throw UTCError.invalidTable }
      }
    }
    self.entries = entries; self.validUntilModifiedJulianDay = validUntilModifiedJulianDay
    self.version = version; self.source = source
    // Reject tables outside Epoch's supported range, including integer overflow.
    _ = try midnight(first.effectiveModifiedJulianDay, offset: first.taiMinusUTC)
    _ = try midnight(validUntilModifiedJulianDay, offset: last.taiMinusUTC)
  }

  /// Load IERS Leap_Second.dat rows (MJD, day, month, year, TAI−UTC).
  /// Caller supplies version, provenance and exclusive coverage from metadata.
  public init(iersData: String, validUntilModifiedJulianDay: Int64, version: String, source: String) throws {
    var entries: [Entry] = []
    for line in iersData.split(whereSeparator: \.isNewline) {
      let row = line.split(separator: "#", omittingEmptySubsequences: false)[0]
      let fields = row.split(whereSeparator: \.isWhitespace)
      if fields.isEmpty { continue }
      guard fields.count == 5, let day = Double(fields[0]), day.isFinite,
            day.rounded() == day, day >= 41317, day < Double(Int64.max),
            let offset = Int(fields[4]), fields[1...3].allSatisfy({ Int($0) != nil }) else {
        throw UTCError.malformedData
      }
      entries.append(Entry(effectiveModifiedJulianDay: Int64(day), taiMinusUTC: offset))
    }
    try self.init(entries: entries, validUntilModifiedJulianDay: validUntilModifiedJulianDay, version: version, source: source)
  }

  public func taiMinusUTC(on modifiedJulianDay: Int64) throws -> Int {
    guard modifiedJulianDay >= entries[0].effectiveModifiedJulianDay,
          modifiedJulianDay < validUntilModifiedJulianDay else { throw UTCError.outsideCoverage }
    return entries.last(where: { $0.effectiveModifiedJulianDay <= modifiedJulianDay })!.taiMinusUTC
  }
  public func dayLength(on modifiedJulianDay: Int64) throws -> Int {
    let offset = try taiMinusUTC(on: modifiedJulianDay)
    let next = entries.first(where: { $0.effectiveModifiedJulianDay == modifiedJulianDay+1 })
    return 86400 + (next.map { $0.taiMinusUTC-offset } ?? 0)
  }
  public func epoch(from utc: UTCDate) throws -> Epoch {
    let offset = try taiMinusUTC(on: utc.modifiedJulianDay)
    guard utc.secondOfDay.wholeSeconds < Int64(try dayLength(on: utc.modifiedJulianDay)) else {
      throw UTCError.invalidSecondOfDay
    }
    return try midnight(utc.modifiedJulianDay, offset: offset).advanced(by: utc.secondOfDay)
  }
  public func utc(from epoch: Epoch) throws -> UTCDate {
    var low = entries[0].effectiveModifiedJulianDay
    var high = validUntilModifiedJulianDay
    guard epoch >= (try midnight(low, offset: entries[0].taiMinusUTC)),
          epoch < (try midnight(high, offset: entries.last!.taiMinusUTC)) else { throw UTCError.outsideCoverage }
    // Search actual atomic midnights; positive and negative leap days are unambiguous.
    while high-low > 1 {
      let mid = low+(high-low)/2
      if try epoch < midnight(mid, offset: taiMinusUTC(on: mid)) { high = mid } else { low = mid }
    }
    let start = try midnight(low, offset: taiMinusUTC(on: low))
    var elapsed = try epoch.duration(since: start)
    let length = Int64(try dayLength(on: low))
    // Subtraction can round a fraction immediately below midnight up to 1.
    // The midnight search already established that this instant is in this day.
    if elapsed.wholeSeconds >= length {
      elapsed = try Duration(wholeSeconds: length-1, fractionalSecond: Double(1).nextDown)
    }
    return try UTCDate(modifiedJulianDay: low, secondOfDay: elapsed)
  }
  private func midnight(_ day: Int64, offset: Int) throws -> Epoch {
    let (jd, overflow) = day.addingReportingOverflow(2400000)
    guard !overflow else { throw TimeArithmeticError.outOfRange }
    let nominal = try Epoch(julianDate: JulianDateComponents(wholeDay: jd, fractionalDay: 0.5, scale: .tt))
    return try nominal.advanced(by: Duration(wholeSeconds: Int64(offset)+32, fractionalSecond: 0.184))
  }

  /// Frozen IERS Bulletin C 72 snapshot, July 2026; strict expiry at 2027-06-28.
  /// No automatic network access. Replace/inject a table to extend coverage.
  public static let iersBulletinC72: LeapSecondTable = try! LeapSecondTable(
    entries: zip([41317,41499,41683,42048,42413,42778,43144,43509,43874,44239,
                  44786,45151,45516,46247,47161,47892,48257,48804,49169,49534,
                  50083,50630,51179,53736,54832,56109,57204,57754], 10...37)
      .map { Entry(effectiveModifiedJulianDay: Int64($0.0), taiMinusUTC: $0.1) },
    validUntilModifiedJulianDay: 61584, version: "IERS Bulletin C 72 (2026-07)",
    source: "https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat")
}

extension Epoch {
  public init(utc: UTCDate, using table: LeapSecondTable) throws { self = try table.epoch(from: utc) }
  public func utc(using table: LeapSecondTable) throws -> UTCDate { try table.utc(from: self) }
}
