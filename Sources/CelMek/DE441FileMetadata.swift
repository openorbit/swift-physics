// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

/// JPL testpo identifiers. Mars through Pluto denote planetary-system
/// barycentres, not resolved planet centres. Values are not NAIF identifiers.
public enum DE441Body: Int, CaseIterable, Sendable {
    case mercury = 1, venus, earth, marsBarycentre, jupiterBarycentre
    case saturnBarycentre, uranusBarycentre, neptuneBarycentre, plutoBarycentre
    case moon, sun, solarSystemBarycentre, earthMoonBarycentre
}

/// Validated binary-header metadata for the supported native DE441 profile.
/// Coverage is TDB Julian date; constants retain their JPL-defined native units.
public struct DE441FileMetadata {
    public let ephemerisNumber: Int
    public let titles: [String]
    public let startJulianDateTDB: Double
    public let endJulianDateTDB: Double
    public let recordDurationDays: Double
    public let coefficientRecordCount: Int
    public let fileSizeBytes: UInt64
    public let recordSizeBytes = 8144
    public let astronomicalUnitKilometres: Double
    public let earthMoonMassRatio: Double
    public let constants: [String: Double]
    public var frame: Frame { CelestialFrames.icrf }
    public let modelTimeScale = TimeScale.tdb
    public var availableBodies: [DE441Body] { DE441Body.allCases }
    public let hasNutationCoefficients = true
    public let hasLibrationCoefficients = true

    init(header: Data, constants values: Data, fileSizeBytes: UInt64, dataRecordCount: Int) throws {
        func int(_ offset: Int) -> Int {
            header.withUnsafeBytes { Int(Int32(littleEndian: $0.loadUnaligned(fromByteOffset: offset, as: Int32.self))) }
        }
        func double(_ data: Data, _ offset: Int) -> Double {
            data.withUnsafeBytes { Double(bitPattern: UInt64(littleEndian: $0.loadUnaligned(fromByteOffset: offset, as: UInt64.self))) }
        }
        func text(_ offset: Int, _ length: Int) throws -> String {
            guard let value = String(data: header.subdata(in: offset..<offset+length), encoding: .ascii) else {
                throw DE441Error.invalidFile("non-ASCII header text")
            }
            return value.trimmingCharacters(in: .whitespacesAndNewlines)
        }
        guard header.count == 8144, values.count == 8144 else { throw DE441Error.invalidFile("short header") }
        guard int(2840) == 441 else { throw DE441Error.invalidFile("unsupported ephemeris or byte order (expected little-endian DE441)") }
        let count = int(2676)
        guard count > 0, count <= 1018 else { throw DE441Error.invalidFile("invalid constant count") }
        // This initial implementation deliberately recognizes one documented profile.
        let expected = [3,14,4,171,10,2,231,13,2,309,11,1,342,8,1,366,7,1,
                        387,6,1,405,6,1,423,6,1,441,13,8,753,11,2,819,10,4,899,10,4]
        let actual = (0..<36).map { int(2696+4*$0) } + (0..<3).map { int(2844+4*$0) }
        guard actual == expected else { throw DE441Error.invalidFile("unsupported DE441 coefficient layout") }
        let extraOffset = 2856 + max(0,count-400)*6
        guard extraOffset+24 <= header.count,
              (0..<6).map({ int(extraOffset+4*$0) }) == [1019,0,0,1019,0,0] else {
            throw DE441Error.invalidFile("unsupported extended coefficient layout")
        }
        let start = double(header,2652), end = double(header,2660), step = double(header,2668)
        let au = double(header,2680), ratio = double(header,2688)
        guard start.isFinite, end.isFinite, start < end, step == 32,
              end-start == Double(dataRecordCount)*step,
              au.isFinite, au > 0, ratio.isFinite, ratio > 0 else {
            throw DE441Error.invalidFile("invalid binary coverage or physical constants")
        }
        var constants: [String:Double] = [:]
        for index in 0..<count {
            let name = try text(index < 400 ? 252+index*6 : 2856+(index-400)*6,6)
            let value = double(values,index*8)
            guard !name.isEmpty, !name.contains("\0"), value.isFinite, constants[name] == nil else {
                throw DE441Error.invalidFile("invalid or duplicate constant")
            }
            constants[name] = value
        }
        guard constants["DENUM"] == 441, constants["AU"] == au, constants["EMRAT"] == ratio else {
            throw DE441Error.invalidFile("inconsistent header constants")
        }
        self.ephemerisNumber = 441; self.titles = try (0..<3).map { try text($0*84,84) }
        self.startJulianDateTDB = start; self.endJulianDateTDB = end; self.recordDurationDays = step
        self.coefficientRecordCount = dataRecordCount; self.fileSizeBytes = fileSizeBytes
        self.astronomicalUnitKilometres = au; self.earthMoonMassRatio = ratio; self.constants = constants
    }
}

extension DE441Ephemeris {
    /// Low-level native AU/AU-day query; a common SI ephemeris adapter is separate.
    public func vector(target: DE441Body, center: DE441Body, at julianDateTDB: Double) throws -> EphemerisVector {
        try vector(target: target.rawValue, center: center.rawValue, at: julianDateTDB)
    }
}
