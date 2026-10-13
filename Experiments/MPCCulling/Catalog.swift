// SPDX-License-Identifier: Apache-2.0
import Foundation
import CelMek
import Math

let au = 149_597_870_700.0
// Explicit experimental Newtonian convention, not a precision ephemeris.
let solarMu = 1.32712440018e20
let catalogFrame = Frame(name: "MPC heliocentric ecliptic J2000 experiment", pseudoInertial: true)
let catalogGM = try! GravitationalParameter(cubicMetresPerSecondSquared: solarMu, convention: .newtonian,
    provenance: PhysicalDataProvenance(source: "Experimental nominal solar GM", version: "1"))
enum ExperimentError: Error { case malformed(String), invalidQuery, outsideIndexWindow }

func calendarEpoch(_ year: Int, _ month: Int, _ day: Double) throws -> Epoch {
    guard (1...9999).contains(year), (1...12).contains(month), day.isFinite else {
        throw ExperimentError.malformed("calendar date")
    }
    let leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0)
    let days = [31, leap ? 29 : 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    guard day >= 1, day < Double(days[month - 1] + 1) else {
        throw ExperimentError.malformed("calendar date")
    }
    let y = month <= 2 ? year - 1 : year
    let m = month <= 2 ? month + 12 : month
    let a = floor(Double(y) / 100)
    let jd = floor(365.25 * Double(y + 4716)) + floor(30.6001 * Double(m + 1)) + day + 2 - a + floor(a / 4) - 1524.5
    return try Epoch(julianDate: jd, scale: .tt)
}
func packedEpoch(_ value: String) throws -> Epoch {
    let chars = Array(value.utf8)
    func digit(_ c: UInt8) throws -> Int {
        if (48...57).contains(c) { return Int(c - 48) }
        if (65...90).contains(c) { return Int(c - 65) + 10 }
        throw ExperimentError.malformed("packed epoch \(value)")
    }
    guard chars.count == 5 else { throw ExperimentError.malformed("packed epoch \(value)") }
    return try calendarEpoch(digit(chars[0]) * 100 + digit(chars[1]) * 10 + digit(chars[2]),
                             digit(chars[3]), Double(digit(chars[4])))
}
struct OrbitRecord {
    let id: String
    let elementEpoch: Epoch? // Comets can omit this; perihelion time remains available.
    let orbit: TwoBodyPropagator
    let perihelion: Double
    let aphelion: Double
    let maximumSpeed: Double
    init(id: String, epoch: Epoch, elementEpoch: Epoch?, q: Double, e: Double,
         inclination: Double, node: Double, peri: Double, anomaly: Double) throws {
        let elements = try OsculatingOrbitalElements(epoch: epoch, frame: catalogFrame, centralBody: .sun,
            gravitationalParameter: catalogGM, semiLatusRectumMetres: q * (1 + e), eccentricity: e,
            inclination: inclination, longitudeOfAscendingNode: node, argumentOfPeriapsis: peri, trueAnomaly: anomaly)
        self.id = id; self.elementEpoch = elementEpoch
        orbit = try TwoBodyPropagator(initialState: elements.cartesianState(), centralBody: .sun, gravitationalParameter: catalogGM)
        perihelion = q; aphelion = e < 1 ? q * (1 + e) / (1 - e) : .infinity
        maximumSpeed = sqrt(solarMu / q * (1 + e))
    }
}

// Bounded-buffer line reader: never creates a String for the entire MPCORB file.
func readLines(_ url: URL, consume: (String, Int) throws -> Void) throws {
    let file = try FileHandle(forReadingFrom: url)
    defer { try? file.close() }
    var pending = Data(), number = 0
    while let chunk = try file.read(upToCount: 1 << 20), !chunk.isEmpty {
        pending.append(chunk)
        var start = pending.startIndex
        while let end = pending[start...].firstIndex(of: 10) {
            number += 1
            try consume(String(decoding: pending[start..<end], as: UTF8.self), number)
            start = end + 1
        }
        pending = Data(pending[start...])
    }
    if !pending.isEmpty { try consume(String(decoding: pending, as: UTF8.self), number + 1) }
}
func loadAsteroids(_ url: URL, limit: Int?) throws -> [OrbitRecord] {
    var records = [OrbitRecord](), inData = false
    let radians = Double.pi / 180
    try readLines(url) { line, number in
        if !inData { if line.hasPrefix("----------") { inData = true }; return }
        if line.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty { return }
        if let limit, records.count >= limit { return }
        let bytes = Array(line.utf8)
        guard bytes.count >= 103 else { throw ExperimentError.malformed("MPCORB line \(number): short record") }
        func field(_ lo: Int, _ hi: Int) -> String {
            String(decoding: bytes[(lo-1)..<hi], as: UTF8.self).trimmingCharacters(in: .whitespaces)
        }
        func numeric(_ lo: Int, _ hi: Int) throws -> Double {
            guard let x = Double(field(lo, hi)), x.isFinite else { throw ExperimentError.malformed("MPCORB line \(number), columns \(lo)-\(hi)") }
            return x
        }
        let epoch = try packedEpoch(field(21, 25))
        let e = try numeric(71, 79), a = try numeric(93, 103) * au
        guard a > 0, e >= 0, e < 1 else { throw ExperimentError.malformed("unsupported asteroid conic at line \(number)") }
        let anomaly = try ConicAnomaly.from(meanAnomaly: numeric(27, 35) * radians, eccentricity: e).trueAnomaly(eccentricity: e)
        records.append(try OrbitRecord(id: "mpc:asteroid:" + field(1, 7), epoch: epoch, elementEpoch: epoch,
            q: a * (1 - e), e: e, inclination: numeric(60, 68) * radians,
            node: numeric(49, 57) * radians, peri: numeric(38, 46) * radians, anomaly: anomaly))
    }
    guard inData else { throw ExperimentError.malformed("missing MPCORB header separator") }
    return records
}
struct CometInput: Decodable {
    let Orbit_type: String
    let Designation_and_name: String
    let Comet_num: Int?
    let Provisional_packed_desig: String?
    let Year_of_perihelion: Int, Month_of_perihelion: Int
    let Day_of_perihelion: Double
    let Perihelion_dist: Double, e: Double, Peri: Double, Node: Double, i: Double
    let Epoch_year: Int?, Epoch_month: Int?, Epoch_day: Int?
}
func loadComets(_ url: URL) throws -> [OrbitRecord] {
    try JSONDecoder().decode([CometInput].self, from: Data(contentsOf: url)).map { c in
        guard let designation = c.Comet_num.flatMap({ _ in c.Designation_and_name.split(separator: "/").first.map(String.init) }) ?? c.Provisional_packed_desig else {
            throw ExperimentError.malformed("comet without designation")
        }
        let epoch = try calendarEpoch(c.Year_of_perihelion, c.Month_of_perihelion, c.Day_of_perihelion)
        let elementEpoch: Epoch?
        if let y = c.Epoch_year, let m = c.Epoch_month, let d = c.Epoch_day {
            elementEpoch = try calendarEpoch(y, m, Double(d))
        } else {
            guard c.Epoch_year == nil, c.Epoch_month == nil, c.Epoch_day == nil else {
                throw ExperimentError.malformed("partial comet element epoch")
            }
            elementEpoch = nil
        }
        return try OrbitRecord(id: "mpc:comet:" + c.Orbit_type + ":" + designation,
            epoch: epoch, elementEpoch: elementEpoch, q: c.Perihelion_dist * au, e: c.e,
            inclination: c.i * .pi / 180, node: c.Node * .pi / 180, peri: c.Peri * .pi / 180, anomaly: 0)
    }
}
