// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
@testable import CelMek

enum LocalDE441Assets {
    static let directory: URL = {
        if let override = ProcessInfo.processInfo.environment["CELMEK_DE441_DATA_DIRECTORY"] {
            return URL(fileURLWithPath: override, isDirectory: true)
        }
        return URL(fileURLWithPath: #filePath).deletingLastPathComponent().deletingLastPathComponent()
            .deletingLastPathComponent().appendingPathComponent("Data/DE441", isDirectory: true)
    }()
    static var available: Bool {
        FileManager.default.fileExists(atPath: directory.appendingPathComponent("linux_m13000p17000.441").path)
        && FileManager.default.fileExists(atPath: directory.appendingPathComponent("testpo.441").path)
    }
}

@Suite(.enabled(if: LocalDE441Assets.available, "Optional DE441 binary and testpo.441 are not installed"))
struct DE441LocalDataTests {
    @Test func commonAdapterAtInteriorAndRecordBoundaryEpochs() throws {
        let reader = try DE441Ephemeris(url: LocalDE441Assets.directory.appendingPathComponent("linux_m13000p17000.441"))
        let model = try FairheadTDBModel()
        let start = try Epoch.j2000.advanced(by: Duration(seconds: -86400*40))
        let end = try Epoch.j2000.advanced(by: Duration(seconds: 86400*40))
        let provider = try DE441StateProvider(reader: reader, target: .earth, coverage: start...end, tdbModel: model)
        let index = floor((2451545-reader.startDate)/32)
        let boundary = reader.startDate+index*32
        for jd in [boundary, boundary+4, boundary+32, 2451545] {
            let date = try JulianDateComponents(value: jd, scale: .tdb)
            let epoch = try Epoch(julianDate: date, using: model)
            let actual = try provider.checkedState(at: epoch)
            let native = try reader.vector(target: .earth, center: .solarSystemBarycentre,
                                           at: epoch.julianDate(in: .tdb, using: model))
            let rate = try epoch.tdbCorrection(using: model).rate
            let scale = reader.metadata.astronomicalUnitKilometres*1000
            for axis in 0..<3 {
                #expect(abs(actual.position[axis]-native.position[axis]*scale) < 1e-5)
                #expect(abs(actual.velocity[axis]-native.velocity[axis]*(scale/86400*rate)) < 1e-9)
            }
            let before = try provider.state(at: epoch.advanced(by: Duration(seconds: -0.5)))
            let after = try provider.state(at: epoch.advanced(by: Duration(seconds: 0.5)))
            for axis in 0..<3 { #expect(abs(after.position[axis]-before.position[axis]-actual.velocity[axis]) < 0.002) }
        }
    }

    @Test func binaryHeaderAndJPLTestpoReferences() throws {
        let reader = try DE441Ephemeris(url: LocalDE441Assets.directory.appendingPathComponent("linux_m13000p17000.441"))
        #expect(reader.metadata.constants.count == 645)
        #expect(reader.metadata.constants["DENUM"] == 441)
        #expect(reader.metadata.recordSizeBytes == 8144)
        let text = try String(contentsOf: LocalDE441Assets.directory.appendingPathComponent("testpo.441"), encoding: .ascii)
        guard let eot = text.range(of: "EOT") else { throw DE441Error.invalidFile("testpo has no EOT marker") }
        let fullValidation = ProcessInfo.processInfo.environment["CELMEK_DE441_FULL_VALIDATION"] == "1"
        var checked = 0, outside = 0, supported = 0
        var targets: Set<Int> = []
        for line in text[eot.upperBound...].split(whereSeparator: \.isNewline) {
            let fields = line.dropFirst(15).split(whereSeparator: \.isWhitespace)
            guard fields.count == 5, let date = Double(fields[0]),
                  let target = Int(fields[1]), let center = Int(fields[2]),
                  let component = Int(fields[3]), let expected = Double(fields[4]) else {
                throw DE441Error.invalidFile("malformed testpo row")
            }
            if date < reader.startDate || date > reader.endDate { outside += 1; continue }
            // Native reader currently exposes vectors and nutation, not libration.
            if target == 15 { continue }
            supported += 1
            // Normal local validation samples throughout the whole file. Full
            // validation is opt-in because it reads nearly the entire 2.6 GB export.
            if !fullValidation && supported % 83 != 1 { continue }
            targets.insert(target)
            let actual: Double
            if target == 14 {
                let value = try reader.nutation(at: date)
                switch component {
                case 1: actual = value.longitude
                case 2: actual = value.obliquity
                case 3: actual = value.longitudeRate
                case 4: actual = value.obliquityRate
                default: throw DE441Error.invalidFile("invalid nutation component")
                }
            } else {
                guard (1...6).contains(component) else { throw DE441Error.invalidFile("invalid vector component") }
                let value = try reader.vector(target: target, center: center, at: date)
                actual = component <= 3 ? value.position[component-1] : value.velocity[component-4]
            }
            #expect(abs(actual-expected) <= 1e-12, "JPL reference at JD \(date), target \(target), component \(component)")
            checked += 1
        }
        #expect(checked > (fullValidation ? 100_000 : 3_000))
        #expect(targets.contains(14))
        #expect((1...13).allSatisfy { targets.contains($0) })
        #expect(outside > 0) // ASCII master header/testpo coverage exceeds this binary export.
    }
}
