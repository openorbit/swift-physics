// SPDX-License-Identifier: Apache-2.0
import Foundation
import XCTest
@testable import CelMek

final class DE441EphemerisTests: XCTestCase {
  private let epoch = 2_451_545.0

  // Synthetic records are deliberately independent of astronomical datasets.
  // Every component uses p(x) = a + 2x + 3(2x² - 1), p'(x) = 2 + 12x.
  private func fixture(spans: [(Double, Double)]? = nil) throws -> URL {
    let spans = spans ?? [(epoch, epoch + 32), (epoch + 32, epoch + 64)]
    let layouts = [(3,14,4,3), (171,10,2,3), (231,13,2,3),
                   (309,11,1,3), (342,8,1,3), (366,7,1,3),
                   (387,6,1,3), (405,6,1,3), (423,6,1,3),
                   (441,13,8,3), (753,11,2,3), (819,10,4,2)]
    var data = Data(repeating: 0, count: 2 * 8_144)
    for (recordIndex, span) in spans.enumerated() {
      var record = [Double](repeating: 0, count: 1_018)
      record[0] = span.0
      record[1] = span.1
      for (body, layout) in layouts.enumerated() {
        let (start, count, subdivisions, coordinates) = layout
        for subdivision in 0..<subdivisions {
          for axis in 0..<coordinates {
            let base = start - 1 + (subdivision * coordinates + axis) * count
            record[base] = Double(1000 * recordIndex + 100 * body + 10 * subdivision + axis)
            record[base + 1] = 2
            record[base + 2] = 3
          }
        }
      }
      for value in record {
        var bits = value.bitPattern.littleEndian
        withUnsafeBytes(of: &bits) { data.append(contentsOf: $0) }
      }
    }
    let url = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
    try data.write(to: url)
    addTeardownBlock { try FileManager.default.removeItem(at: url) }
    return url
  }

  func testPolynomialValuesAndVelocities() throws {
    let reader = try DE441Ephemeris(url: fixture())
    XCTAssertEqual(reader.startDate, epoch)
    XCTAssertEqual(reader.endDate, epoch + 64)
    // Mercury: first subdivision has duration 8 days; t=6 gives x=0.5.
    let state = try reader.vector(target: 1, center: 12, at: epoch + 6)
    for axis in 0..<3 {
      XCTAssertEqual(state.position[axis], (Double(axis) - 0.5) / DE441Ephemeris.astronomicalUnit, accuracy: 1e-15)
      XCTAssertEqual(state.velocity[axis], 2 / DE441Ephemeris.astronomicalUnit, accuracy: 1e-15)
    }
  }

  func testSubdivisionRecordAndCoverageBoundaries() throws {
    let reader = try DE441Ephemeris(url: fixture())
    // Exact internal boundaries select the next subdivision or record.
    for (offset, expected) in [(0.0, 1.0), (8, 11), (32, 1001), (64, 1035), (0, 1)] {
      let state = try reader.vector(target: 1, center: 12, at: epoch + offset)
      XCTAssertEqual(state.position.x, expected / DE441Ephemeris.astronomicalUnit, accuracy: 1e-15)
    }
  }

  func testEarthMoonReconstructionAndRelativeVelocity() throws {
    let reader = try DE441Ephemeris(url: fixture())
    let date = epoch + 2
    let earth = try reader.vector(target: 3, center: 12, at: date)
    let moon = try reader.vector(target: 10, center: 12, at: date)
    let emb = try reader.vector(target: 13, center: 12, at: date)
    let relative = try reader.vector(target: 10, center: 3, at: date)
    let ratio = DE441Ephemeris.earthMoonMassRatio
    // Moon first subdivision midpoint: x=0, derivative scale=1/2.
    for axis in 0..<3 {
      XCTAssertEqual(relative.position[axis], Double(897 + axis) / DE441Ephemeris.astronomicalUnit, accuracy: 1e-15)
      XCTAssertEqual(relative.velocity[axis], 1 / DE441Ephemeris.astronomicalUnit, accuracy: 1e-15)
      XCTAssertEqual((earth.position[axis] * ratio + moon.position[axis]) / (ratio + 1), emb.position[axis], accuracy: 1e-15)
      XCTAssertEqual((earth.velocity[axis] * ratio + moon.velocity[axis]) / (ratio + 1), emb.velocity[axis], accuracy: 1e-15)
    }
    let reverse = try reader.vector(target: 3, center: 10, at: date)
    XCTAssertEqual(reverse.position, -relative.position)
    XCTAssertEqual(reverse.velocity, -relative.velocity)
  }

  func testNutationUsesRadiansAndDailyDerivatives() throws {
    let reader = try DE441Ephemeris(url: fixture())
    let angles = try reader.nutation(at: epoch + 6)
    XCTAssertEqual(angles.longitude, 1099.5)
    XCTAssertEqual(angles.obliquity, 1100.5)
    XCTAssertEqual(angles.longitudeRate, 2)
    XCTAssertEqual(angles.obliquityRate, 2)
  }

  func testInvalidQueriesThrowIncludingBarycentreOnly() throws {
    let reader = try DE441Ephemeris(url: fixture())
    for date in [epoch - 1, epoch + 65, .nan, .infinity, -.infinity] {
      XCTAssertThrowsError(try reader.vector(target: 12, center: 12, at: date)) {
        guard case DE441Error.dateOutOfRange = $0 else { return XCTFail("Unexpected error: \($0)") }
      }
      XCTAssertThrowsError(try reader.nutation(at: date))
    }
    for body in [0, 14, -1] {
      XCTAssertThrowsError(try reader.vector(target: body, center: 12, at: epoch))
      XCTAssertThrowsError(try reader.vector(target: 12, center: body, at: epoch))
    }
  }

  func testMalformedFilesThrow() throws {
    for spans in [[], [(epoch, epoch + 31)], [(epoch, epoch + 32), (epoch + 33, epoch + 65)], [(.nan, .nan)]] as [[(Double, Double)]] {
      XCTAssertThrowsError(try DE441Ephemeris(url: fixture(spans: spans)))
    }
    let url = try fixture()
    let handle = try FileHandle(forWritingTo: url)
    try handle.truncate(atOffset: 123)
    try handle.close()
    XCTAssertThrowsError(try DE441Ephemeris(url: url))
  }
}
