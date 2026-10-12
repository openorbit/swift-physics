//
// SPDX-License-Identifier: Apache-2.0
//
// Copyright 2022 Mattias Holm
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

import Foundation

public struct EphemerisVector {
  public let position: SIMD3<Double>
  public let velocity: SIMD3<Double>
}

public struct NutationAngles {
  public let longitude: Double
  public let obliquity: Double
  public let longitudeRate: Double
  public let obliquityRate: Double
}

public enum DE441Error: Error, CustomStringConvertible {
  case invalidFile(String)
  case dateOutOfRange(Double)
  case unsupportedBody(Int)

  public var description: String {
    switch self {
    case .invalidFile(let reason): "Invalid DE441 file: \(reason)"
    case .dateOutOfRange(let date): "Julian date \(date) is outside this DE441 file"
    case .unsupportedBody(let body): "DE441 body \(body) is not supported"
    }
  }
}

/// Random-access reader for JPL's little-endian DE441 binary export.
///
/// Body numbers follow `testpo.441`: 1...11 are Mercury through Sun, 12 is
/// the solar-system barycentre, and 13 is the Earth-Moon barycentre.
public final class DE441Ephemeris {
  public static let astronomicalUnit = 149_597_870.7
  public static let earthMoonMassRatio = 81.3005682214972154

  private struct CoefficientLayout {
    let start: Int
    let count: Int
    let subdivisions: Int
  }

  private static let recordSize = 8_144
  private static let dataRecordOffset = 2
  private static let recordDuration = 32.0
  private static let layouts = [
    CoefficientLayout(start: 3, count: 14, subdivisions: 4),
    CoefficientLayout(start: 171, count: 10, subdivisions: 2),
    CoefficientLayout(start: 231, count: 13, subdivisions: 2),
    CoefficientLayout(start: 309, count: 11, subdivisions: 1),
    CoefficientLayout(start: 342, count: 8, subdivisions: 1),
    CoefficientLayout(start: 366, count: 7, subdivisions: 1),
    CoefficientLayout(start: 387, count: 6, subdivisions: 1),
    CoefficientLayout(start: 405, count: 6, subdivisions: 1),
    CoefficientLayout(start: 423, count: 6, subdivisions: 1),
    CoefficientLayout(start: 441, count: 13, subdivisions: 8),
    CoefficientLayout(start: 753, count: 11, subdivisions: 2),
  ]
  private static let nutationLayout = CoefficientLayout(start: 819, count: 10, subdivisions: 4)

  public let metadata: DE441FileMetadata
  private let handle: FileHandle
  private let dataRecordCount: Int
  private var cachedRecordIndex: Int?
  private var cachedRecord: [Double] = []
  public private(set) var startDate = 0.0
  public private(set) var endDate = 0.0

  public init(url: URL) throws {
    let openedHandle = try FileHandle(forReadingFrom: url)
    handle = openedHandle
    var completed = false
    defer { if !completed { try? openedHandle.close() } }
    let size = try handle.seekToEnd()
    guard size % UInt64(Self.recordSize) == 0 else {
      throw DE441Error.invalidFile("size is not a multiple of \(Self.recordSize) bytes")
    }
    guard size <= UInt64(Int.max) else { throw DE441Error.invalidFile("file size exceeds platform range") }
    let records = Int(size) / Self.recordSize
    guard records > Self.dataRecordOffset else {
      throw DE441Error.invalidFile("no coefficient records")
    }
    dataRecordCount = records - Self.dataRecordOffset
    try handle.seek(toOffset: 0)
    guard let header = try handle.read(upToCount: Self.recordSize), header.count == Self.recordSize,
          let constants = try handle.read(upToCount: Self.recordSize), constants.count == Self.recordSize else {
      throw DE441Error.invalidFile("short binary header")
    }
    metadata = try DE441FileMetadata(header: header, constants: constants,
                                    fileSizeBytes: size, dataRecordCount: dataRecordCount)
    let first = try readRecord(0)
    let last = try readRecord(dataRecordCount - 1)
    startDate = first[0]
    endDate = last[1]
    guard startDate == metadata.startJulianDateTDB, endDate == metadata.endJulianDateTDB,
          startDate.isFinite, endDate.isFinite,
          endDate - startDate == Double(dataRecordCount) * Self.recordDuration,
          first[1] - first[0] == Self.recordDuration,
          last[1] - last[0] == Self.recordDuration else {
      throw DE441Error.invalidFile("unexpected coefficient record span")
    }
    completed = true
  }

  deinit {
    try? handle.close()
  }

  /// Position in AU and velocity in AU/day, evaluated at a TDB Julian date.
  public func vector(target: Int, center: Int, at julianDate: Double) throws -> EphemerisVector {
    guard julianDate.isFinite, julianDate >= startDate, julianDate <= endDate else {
      throw DE441Error.dateOutOfRange(julianDate)
    }
    let targetVector = try barycentricVector(body: target, at: julianDate)
    let centerVector = try barycentricVector(body: center, at: julianDate)
    return EphemerisVector(
      position: targetVector.position - centerVector.position,
      velocity: targetVector.velocity - centerVector.velocity
    )
  }

  /// Nutation in longitude and obliquity, in radians and radians/day.
  func nutation(at julianDate: Double) throws -> NutationAngles {
    let values = try storedValues(layout: Self.nutationLayout, coordinates: 2, at: julianDate)
    return NutationAngles(
      longitude: values[0].value,
      obliquity: values[1].value,
      longitudeRate: values[0].derivative,
      obliquityRate: values[1].derivative
    )
  }

  private func barycentricVector(body: Int, at julianDate: Double) throws -> EphemerisVector {
    switch body {
    case 1, 2, 4...9, 11:
      return try storedVector(index: body - 1, at: julianDate)
    case 3:
      let emb = try storedVector(index: 2, at: julianDate)
      let moon = try storedVector(index: 9, at: julianDate)
      return emb - moon / (1 + metadata.earthMoonMassRatio)
    case 10:
      let earth = try barycentricVector(body: 3, at: julianDate)
      return earth + (try storedVector(index: 9, at: julianDate))
    case 12:
      return EphemerisVector(position: .zero, velocity: .zero)
    case 13:
      return try storedVector(index: 2, at: julianDate)
    default:
      throw DE441Error.unsupportedBody(body)
    }
  }

  private func storedVector(index: Int, at julianDate: Double) throws -> EphemerisVector {
    let values = try storedValues(layout: Self.layouts[index], coordinates: 3, at: julianDate)
    return EphemerisVector(
      position: SIMD3(values[0].value, values[1].value, values[2].value) / metadata.astronomicalUnitKilometres,
      velocity: SIMD3(values[0].derivative, values[1].derivative, values[2].derivative)
        / metadata.astronomicalUnitKilometres
    )
  }

  private func storedValues(
    layout: CoefficientLayout,
    coordinates: Int,
    at julianDate: Double
  ) throws -> [(value: Double, derivative: Double)] {
    guard julianDate.isFinite, julianDate >= startDate, julianDate <= endDate else {
      throw DE441Error.dateOutOfRange(julianDate)
    }
    var recordIndex = Int(floor((julianDate - startDate) / Self.recordDuration))
    if recordIndex == dataRecordCount { recordIndex -= 1 }
    let record = try readRecord(recordIndex)
    guard julianDate >= record[0], julianDate <= record[1] else {
      throw DE441Error.invalidFile("non-contiguous coefficient records")
    }

    let subdivisionDuration = Self.recordDuration / Double(layout.subdivisions)
    var subdivision = Int(floor((julianDate - record[0]) / subdivisionDuration))
    if subdivision == layout.subdivisions { subdivision -= 1 }
    let subdivisionStart = record[0] + Double(subdivision) * subdivisionDuration
    let time = 2 * (julianDate - subdivisionStart) / subdivisionDuration - 1
    let scale = 2 / subdivisionDuration
    let base = layout.start - 1 + subdivision * layout.count * coordinates

    var result: [(value: Double, derivative: Double)] = []
    for coordinate in 0..<coordinates {
      let start = base + coordinate * layout.count
      let values = evaluate(
        coefficients: record[start ..< start + layout.count],
        at: time
      )
      guard values.value.isFinite, values.derivative.isFinite else {
        throw DE441Error.invalidFile("nonfinite interpolated value")
      }
      result.append((values.value, values.derivative * scale))
    }
    return result
  }

  private func evaluate(coefficients: ArraySlice<Double>, at x: Double) -> (value: Double, derivative: Double) {
    var tPrevious = 1.0
    var tCurrent = x
    var dPrevious = 0.0
    var dCurrent = 1.0
    var value = coefficients[coefficients.startIndex]
    var derivative = 0.0
    for order in 1..<coefficients.count {
      let coefficient = coefficients[coefficients.startIndex + order]
      value += coefficient * tCurrent
      derivative += coefficient * dCurrent
      let nextT = 2 * x * tCurrent - tPrevious
      let nextD = 2 * tCurrent + 2 * x * dCurrent - dPrevious
      tPrevious = tCurrent
      tCurrent = nextT
      dPrevious = dCurrent
      dCurrent = nextD
    }
    return (value, derivative)
  }

  private func readRecord(_ index: Int) throws -> [Double] {
    guard index >= 0, index < dataRecordCount else { throw DE441Error.invalidFile("coefficient record index out of range") }
    if cachedRecordIndex == index { return cachedRecord }
    let offset = UInt64(index + Self.dataRecordOffset) * UInt64(Self.recordSize)
    try handle.seek(toOffset: offset)
    guard let data = try handle.read(upToCount: Self.recordSize), data.count == Self.recordSize else {
      throw DE441Error.invalidFile("short read at coefficient record \(index)")
    }
    let record: [Double] = data.withUnsafeBytes { bytes in
      (0..<(Self.recordSize / 8)).map { item in
        let bits = bytes.loadUnaligned(fromByteOffset: item * 8, as: UInt64.self)
        return Double(bitPattern: UInt64(littleEndian: bits))
      }
    }
    guard record.allSatisfy(\.isFinite) else { throw DE441Error.invalidFile("nonfinite coefficient record") }
    let expectedStart = metadata.startJulianDateTDB + Double(index) * Self.recordDuration
    guard record[0] == expectedStart, record[1] == expectedStart + Self.recordDuration else {
      throw DE441Error.invalidFile("invalid span at coefficient record \(index)")
    }
    cachedRecordIndex = index
    cachedRecord = record
    return record
  }
}

private func + (lhs: EphemerisVector, rhs: EphemerisVector) -> EphemerisVector {
  EphemerisVector(position: lhs.position + rhs.position, velocity: lhs.velocity + rhs.velocity)
}

private func - (lhs: EphemerisVector, rhs: EphemerisVector) -> EphemerisVector {
  EphemerisVector(position: lhs.position - rhs.position, velocity: lhs.velocity - rhs.velocity)
}

private func / (lhs: EphemerisVector, rhs: Double) -> EphemerisVector {
  EphemerisVector(position: lhs.position / rhs, velocity: lhs.velocity / rhs)
}
