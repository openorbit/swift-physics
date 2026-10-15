//
// SPDX-License-Identifier: Apache-2.0
//
// Copyright 2026 Mattias Holm
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
import Math

public enum SpiceTextKernelError: Error, Equatable {
    case missingVariable(String)
    case invalidVariable(String)
    case unsupportedFrameClass(Int32)
    case unsupportedFrameRepresentation(String)
}

// Bounded Gregorian date subset for ordinary LSK and SCLK metadata.
// Pool dates carry no implied time scale.
internal func spiceDateSeconds(_ token: String) throws -> Double {
    guard token.utf8.count <= 80 else { throw SpiceTextKernelError.invalidVariable(token) }
    let pattern = #"^@(?:([0-9]{4})-([A-Za-z]{3})-([0-9]{1,2})|([0-9]{1,2})-([A-Za-z]{3})-([0-9]{4}))(?:(?:-|//?)([0-9]{1,2})(?::([0-9]{2}))?(?::([0-9]{2}(?:\.[0-9]+)?))?)?$"#
    let regex = try NSRegularExpression(pattern: pattern)
    let ns = token as NSString
    guard let match = regex.firstMatch(in: token, range: NSRange(location: 0, length: ns.length)) else {
        throw SpiceTextKernelError.invalidVariable(token)
    }
    func field(_ i: Int) -> String? { let r = match.range(at: i); return r.location == NSNotFound ? nil : ns.substring(with: r) }
    let months = ["JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"]
    guard let month = months.firstIndex(of: (field(2) ?? field(5)!).uppercased()),
          let year = Int(field(1) ?? field(6)!), year >= 1,
          let day = Int(field(3) ?? field(4)!) else { throw SpiceTextKernelError.invalidVariable(token) }
    let date = try GregorianDay(year: year, month: month + 1, day: day)
    let hour = Int(field(7) ?? "0")!, minute = Int(field(8) ?? "0")!, second = Double(field(9) ?? "0")!
    guard hour < 24, minute < 60, second < 60 else { throw SpiceTextKernelError.invalidVariable(token) }
    return Double(date.modifiedJulianDay - 51544) * 86400 - 43200 + Double(hour * 3600 + minute * 60) + second
}

/// Validated LSK parameters. The periodic coefficients are retained, not silently
/// installed as a replacement for a caller's TDB model.
public struct SpiceLeapSecondData: Sendable {
    public let ttMinusTAI: Double
    public let periodicAmplitude: Double
    public let eccentricity: Double
    public let meanAnomaly: [Double]
    public let entries: [LeapSecondTable.Entry]
    public let sources: [SpiceKernelSource]
    /// LSKs do not certify a future validity date; callers must supply one.
    public func table(validUntilModifiedJulianDay: Int64, version: String) throws -> LeapSecondTable {
        try LeapSecondTable(entries: entries, validUntilModifiedJulianDay: validUntilModifiedJulianDay,
                            version: version, source: sources.map(\.path).joined(separator: "; "))
    }
}

/// Raw orientation parameters retain their documented degree/day/century units.
/// Evaluation and frame construction are separate from validated pool decoding.
public struct SpiceBodyOrientationParameters: Sendable {
    public let poleRightAscension: [Double]
    public let poleDeclination: [Double]
    public let primeMeridian: [Double]
    public let additionalParameters: [String: [Double]]
    public let sources: [SpiceKernelSource]
}

public struct SpiceFixedFrameDefinition: Sendable {
    public let frameID: Int32
    public let name: String
    public let centreID: Int32
    public let classID: Int32
    public let relativeFrameName: String
    /// Column-major proper rotation, mapping this frame's vectors to relative axes.
    public let matrix: [Double]
    public let sources: [SpiceKernelSource]
}

extension SpiceKernelSnapshot {
    internal func numeric(_ name: String, count: Int? = nil) throws -> [Double] {
        guard let variable = self[name] else { throw SpiceTextKernelError.missingVariable(name) }
        guard case .numbers(let values) = variable.value, !values.isEmpty,
              values.allSatisfy(\.isFinite), count == nil || values.count == count else {
            throw SpiceTextKernelError.invalidVariable(name)
        }
        return values
    }
    private func string(_ name: String) throws -> String {
        guard let variable = self[name] else { throw SpiceTextKernelError.missingVariable(name) }
        guard case .strings(let values) = variable.value, values.count == 1, !values[0].isEmpty else {
            throw SpiceTextKernelError.invalidVariable(name)
        }
        return values[0]
    }
    private func integer(_ name: String) throws -> Int32 {
        let value = try numeric(name, count: 1)[0]
        guard let integer = Int32(exactly: value) else { throw SpiceTextKernelError.invalidVariable(name) }
        return integer
    }
    internal func sources(for names: [String]) -> [SpiceKernelSource] {
        var result: [SpiceKernelSource] = []
        for name in names { for source in self[name]?.sources ?? [] where !result.contains(source) { result.append(source) } }
        return result
    }
    public func leapSecondData() throws -> SpiceLeapSecondData {
        let names = ["DELTET/DELTA_T_A", "DELTET/K", "DELTET/EB", "DELTET/M", "DELTET/DELTA_AT"]
        let offset = try numeric(names[0], count: 1)[0]
        // Existing Epoch has the standard fixed TT-TAI relation.
        guard offset == 32.184 else { throw SpiceTextKernelError.invalidVariable(names[0]) }
        let k = try numeric(names[1], count: 1)[0], eb = try numeric(names[2], count: 1)[0]
        let m = try numeric(names[3], count: 2), pairs = try numeric(names[4])
        guard pairs.count % 2 == 0 else { throw SpiceTextKernelError.invalidVariable(names[4]) }
        var entries: [LeapSecondTable.Entry] = []
        for i in stride(from: 0, to: pairs.count, by: 2) {
            let day = (pairs[i+1] + 43200) / 86400 + 51544
            guard let mjd = Int64(exactly: day), let tai = Int(exactly: pairs[i]) else {
                throw SpiceTextKernelError.invalidVariable(names[4])
            }
            entries.append(.init(effectiveModifiedJulianDay: mjd, taiMinusUTC: tai))
        }
        guard let last = entries.last, last.effectiveModifiedJulianDay < Int64.max else {
            throw SpiceTextKernelError.invalidVariable(names[4])
        }
        let provenance = sources(for: names)
        _ = try LeapSecondTable(entries: entries, validUntilModifiedJulianDay: last.effectiveModifiedJulianDay + 1,
                               version: "validation", source: provenance.map(\.path).joined(separator: "; "))
        return SpiceLeapSecondData(ttMinusTAI: offset, periodicAmplitude: k, eccentricity: eb,
                                   meanAnomaly: m, entries: entries, sources: provenance)
    }
    /// Caller confirms conventional PCK kilometres; returns metres.
    public func bodyRadiiMetres(bodyID: Int32) throws -> [Double] {
        let name = "BODY\(bodyID)_RADII"
        let radii = try numeric(name, count: 3).map { $0 * 1000 }
        guard radii.allSatisfy({ $0.isFinite && $0 > 0 }) else { throw SpiceTextKernelError.invalidVariable(name) }
        return radii
    }
    public func bodyOrientationParameters(bodyID: Int32) throws -> SpiceBodyOrientationParameters {
        let prefix = "BODY\(bodyID)_"
        let names = [prefix + "POLE_RA", prefix + "POLE_DEC", prefix + "PM"]
        let ra = try numeric(names[0], count: 3), dec = try numeric(names[1], count: 3), pm = try numeric(names[2], count: 3)
        var additional: [String: [Double]] = [:]
        for suffix in ["NUT_PREC_RA", "NUT_PREC_DEC", "NUT_PREC_PM", "LONG_AXIS", "CONSTANTS_JED_EPOCH", "CONSTANTS_REF_FRAME"] {
            if self[prefix + suffix] != nil { additional[suffix] = try numeric(prefix + suffix) }
        }
        return SpiceBodyOrientationParameters(poleRightAscension: ra, poleDeclination: dec, primeMeridian: pm,
            additionalParameters: additional, sources: sources(for: names + additional.keys.sorted().map { prefix + $0 }))
    }
    /// Explicit ID lookup; built-in frames and name aliases are not synthesized.
    public func fixedFrameDefinition(frameID: Int32) throws -> SpiceFixedFrameDefinition {
        guard frameID != 0 else { throw SpiceTextKernelError.invalidVariable("frameID") }
        let prefix = "FRAME_\(frameID)_"
        let name = try string(prefix + "NAME"), frameClass = try integer(prefix + "CLASS")
        guard frameClass == 4 else { throw SpiceTextKernelError.unsupportedFrameClass(frameClass) }
        let classID = try integer(prefix + "CLASS_ID"), centre = try integer(prefix + "CENTER")
        guard classID == frameID else { throw SpiceTextKernelError.invalidVariable(prefix + "CLASS_ID") }
        let tk = "TKFRAME_\(frameID)_"
        let relative = try string(tk + "RELATIVE"), spec = try string(tk + "SPEC")
        var usedNames = [tk + "RELATIVE", tk + "SPEC"]
        let matrix: [Double]
        func columns(_ q: Quaternion) -> [Double] {
            let m = q.matrix
            return (0..<9).map { m[$0 % 3, $0 / 3] }
        }
        switch spec {
        case "MATRIX":
            matrix = try numeric(tk + "MATRIX", count: 9); usedNames.append(tk + "MATRIX")
        case "ANGLES":
            let angles = try numeric(tk + "ANGLES", count: 3), axes = try numeric(tk + "AXES", count: 3)
            let units = try string(tk + "UNITS")
            let scales: [String: Double] = ["RADIANS": 1, "DEGREES": .pi/180, "ARCSECONDS": .pi/648000,
                "ARCMINUTES": .pi/10800, "HOURANGLE": .pi/12, "MINUTEANGLE": .pi/720, "SECONDANGLE": .pi/43200]
            guard let scale = scales[units] else { throw SpiceTextKernelError.invalidVariable(tk + "UNITS") }
            var q = Quaternion.identity
            for i in 0..<3 {
                guard let axis = Int(exactly: axes[i]), (1...3).contains(axis) else {
                    throw SpiceTextKernelError.invalidVariable(tk + "AXES")
                }
                let angle = angles[i] * scale
                guard angle.isFinite else { throw SpiceTextKernelError.invalidVariable(tk + "ANGLES") }
                var direction = SIMD3<Double>.zero; direction[axis-1] = 1
                // Public SPICE Euler matrices rotate coordinates with negative Hamilton angles.
                q = q * Quaternion(angle: -angle.truncatingRemainder(dividingBy: 2 * .pi), axis: direction)
            }
            matrix = columns(q); usedNames += [tk + "ANGLES", tk + "AXES", tk + "UNITS"]
        case "QUATERNION":
            let q = try numeric(tk + "Q", count: 4)
            let norm = q.reduce(0.0) { $0 + $1 * $1 }
            guard norm.isFinite, abs(norm-1) <= 1e-10 else { throw SpiceTextKernelError.invalidVariable(tk + "Q") }
            matrix = columns(Quaternion(x: q[1], y: q[2], z: q[3], w: q[0])); usedNames.append(tk + "Q")
        default: throw SpiceTextKernelError.unsupportedFrameRepresentation(spec)
        }
        for a in 0..<3 { for b in 0..<3 {
            let dot = (0..<3).reduce(0.0) { $0 + matrix[a*3+$1] * matrix[b*3+$1] }
            guard abs(dot - (a == b ? 1 : 0)) <= 1e-10 else { throw SpiceTextKernelError.invalidVariable(tk + "MATRIX") }
        } }
        let det = matrix[0]*(matrix[4]*matrix[8]-matrix[7]*matrix[5]) - matrix[3]*(matrix[1]*matrix[8]-matrix[7]*matrix[2]) + matrix[6]*(matrix[1]*matrix[5]-matrix[4]*matrix[2])
        guard abs(det - 1) <= 1e-10 else { throw SpiceTextKernelError.invalidVariable(tk + "MATRIX") }
        return SpiceFixedFrameDefinition(frameID: frameID, name: name, centreID: centre, classID: classID,
            relativeFrameName: relative, matrix: matrix, sources: sources(for: [prefix+"NAME",prefix+"CLASS",prefix+"CLASS_ID",prefix+"CENTER"] + usedNames))
    }
}
