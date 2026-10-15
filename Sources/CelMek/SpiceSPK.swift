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

public enum SpiceSPKError: Error, Equatable {
    case invalidEpoch
    case invalidSegment
    case unsupportedType(Int32)
    case outOfCoverage
    case noSegment(targetID: Int32, secondsTDB: Double)
    case malformedSegment(path: String, array: Int, reason: String)
    case numericalRangeExceeded
}

/// Primitive geometric state in native SPK axes, relative to the descriptor's centre.
/// Velocity is per TDB argument second; this is not yet an atomic-Epoch/SI adapter.
public struct SpiceSPKState: Sendable {
    public let secondsTDB: Double
    public let positionKilometres: SIMD3<Double>
    public let velocityKilometresPerTDBSecond: SIMD3<Double>
    public let descriptor: SpiceSPKDescriptor
    public let sourcePath: String
    public let segmentName: String
    public let arrayIndex: Int
    public let recordIndex: Int
}

/// Validated type 2/3/21 directory and bounded on-demand record evaluation.
/// This value retains its reader; the underlying binary asset must remain unchanged.
public struct SpiceSPKSegment: Sendable {
    public let file: SpiceDAFFile
    public let arrayIndex: Int
    public let descriptor: SpiceSPKDescriptor
    public let initialSecondsTDB: Double
    /// Zero for variable-interval type 21 records.
    public let intervalSeconds: Double
    public let recordSizeWords: Int
    public let recordCount: Int
    /// Chebyshev degree for types 2/3; -1 for type 21 difference lines.
    public let polynomialDegree: Int
    private let type21: SpiceSPKType21Directory?

    public init(file: SpiceDAFFile, array index: Int) throws {
        guard file.metadata.kind == "SPK", file.metadata.arrays.indices.contains(index),
              let descriptor = file.metadata.arrays[index].spk else { throw SpiceSPKError.invalidSegment }
        guard [2, 3, 21].contains(descriptor.dataType) else { throw SpiceSPKError.unsupportedType(descriptor.dataType) }
        if descriptor.dataType == 21 {
            let directory = try file.spk21Directory(array: index)
            self.file = file; self.arrayIndex = index; self.descriptor = descriptor
            self.initialSecondsTDB = descriptor.startSecondsTDB; self.intervalSeconds = 0
            self.recordSizeWords = directory.size; self.recordCount = directory.count
            self.polynomialDegree = -1; self.type21 = directory
            return
        }
        self.type21 = nil
        let array = file.metadata.arrays[index], words = array.finalAddress - array.initialAddress + 1
        func malformed(_ reason: String) -> SpiceSPKError { .malformedSegment(path: file.url.path, array: index, reason: reason) }
        guard words >= 4 else { throw malformed("missing directory") }
        let directory = try file.readWords(array: index, offset: words - 4, count: 4)
        func integer(_ value: Double) throws -> Int {
            // DAF addresses already bound the possible count/size below Int32.max.
            guard value.isFinite, value >= 1, value <= Double(words - 4), value.rounded(.towardZero) == value else {
                throw malformed("invalid directory integer")
            }
            return Int(value)
        }
        let initial = directory[0], interval = directory[1]
        let size = try integer(directory[2]), count = try integer(directory[3])
        let components = descriptor.dataType == 2 ? 3 : 6
        guard size >= 2 + components, (size - 2) % components == 0,
              (words - 4) % size == 0, (words - 4) / size == count else { throw malformed("inconsistent record layout") }
        guard size <= file.limits.maximumReadWords else { throw SpiceDAFError.resourceLimit("SPK record words") }
        let end = initial + interval * Double(count)
        guard initial.isFinite, interval.isFinite, interval > 0, (interval/2) > 0,
              end.isFinite, end > initial, initial + interval > initial, end - interval < end else {
            throw malformed("invalid or unresolvable directory time grid")
        }
        guard descriptor.startSecondsTDB >= initial - spkTimeSlack(initial, interval),
              descriptor.endSecondsTDB <= end + spkTimeSlack(end, interval) else { throw malformed("descriptor exceeds record coverage") }
        self.file = file; self.arrayIndex = index; self.descriptor = descriptor
        self.initialSecondsTDB = initial; self.intervalSeconds = interval
        self.recordSizeWords = size; self.recordCount = count; self.polynomialDegree = (size - 2)/components - 1
    }

    /// Types 2/3 select the later record at an internal boundary; type 21 selects the earlier.
    /// The final endpoint uses the last covering record.
    /// Only the selected record is read/validated. No extrapolation beyond descriptor coverage.
    public func state(atSecondsTDB seconds: Double) throws -> SpiceSPKState {
        guard seconds.isFinite else { throw SpiceSPKError.invalidEpoch }
        guard seconds >= descriptor.startSecondsTDB, seconds <= descriptor.endSecondsTDB else { throw SpiceSPKError.outOfCoverage }
        func malformed(_ reason: String) -> SpiceSPKError { .malformedSegment(path: file.url.path, array: arrayIndex, reason: reason) }
        if let directory = type21 {
            let (record, position, velocity) = try directory.evaluate(file: file, array: arrayIndex, seconds: seconds)
            return SpiceSPKState(secondsTDB: seconds, positionKilometres: position, velocityKilometresPerTDBSecond: velocity,
                descriptor: descriptor, sourcePath: file.url.path, segmentName: file.metadata.arrays[arrayIndex].name,
                arrayIndex: arrayIndex, recordIndex: record)
        }
        let ratio = (seconds - initialSecondsTDB)/intervalSeconds
        guard ratio.isFinite else { throw SpiceSPKError.numericalRangeExceeded }
        var record = Int(max(0, min(Double(recordCount - 1), floor(ratio))))
        // Correct a quotient rounded across one boundary using the directory's actual grid.
        if record + 1 < recordCount, seconds >= initialSecondsTDB + Double(record + 1)*intervalSeconds { record += 1 }
        if record > 0, seconds < initialSecondsTDB + Double(record)*intervalSeconds { record -= 1 }
        let values = try file.readWords(array: arrayIndex, offset: record * recordSizeWords, count: recordSizeWords)
        guard values.allSatisfy(\.isFinite), values[1] > 0 else { throw malformed("non-finite coefficients or invalid record radius") }
        let start = initialSecondsTDB + Double(record)*intervalSeconds
        let expectedMid = start + intervalSeconds/2
        guard abs(values[0] - expectedMid) <= spkTimeSlack(expectedMid, intervalSeconds),
              abs(values[1] - intervalSeconds/2) <= spkTimeSlack(intervalSeconds/2, intervalSeconds) else {
            throw malformed("record interval disagrees with directory")
        }
        let normalized = (seconds - values[0])/values[1]
        let slack = spkTimeSlack(values[0], intervalSeconds)/values[1] + 16*Double.ulpOfOne
        guard normalized.isFinite, normalized >= -1 - slack, normalized <= 1 + slack else { throw malformed("epoch outside selected coefficient interval") }
        let x = max(-1, min(1, normalized)), coefficientCount = polynomialDegree + 1
        var position = SIMD3<Double>.zero, velocity = SIMD3<Double>.zero
        for axis in 0..<3 {
            let polynomial = try spkPolynomial(values, offset: 2 + axis*coefficientCount, count: coefficientCount, x: x, differentiate: descriptor.dataType == 2)
            position[axis] = polynomial.value
            if descriptor.dataType == 2 { velocity[axis] = polynomial.derivative/values[1] }
            else {
                velocity[axis] = try spkPolynomial(values, offset: 2 + (axis + 3)*coefficientCount, count: coefficientCount, x: x, differentiate: false).value
            }
        }
        guard (0..<3).allSatisfy({ position[$0].isFinite && velocity[$0].isFinite }) else { throw SpiceSPKError.numericalRangeExceeded }
        return SpiceSPKState(secondsTDB: seconds, positionKilometres: position, velocityKilometresPerTDBSecond: velocity,
            descriptor: descriptor, sourcePath: file.url.path, segmentName: file.metadata.arrays[arrayIndex].name,
            arrayIndex: arrayIndex, recordIndex: record)
    }
}

// Allow only storage-rounding discrepancies in record/directory consistency checks.
private func spkTimeSlack(_ time: Double, _ interval: Double) -> Double { 8 * max(time.ulp, interval.ulp) }

/// Clenshaw summation of c0*T0 + c1*T1 + ... and its differentiated recurrence.
/// In particular c0 is used in full; some mathematical series use a different convention.
private func spkPolynomial(_ coefficients: [Double], offset: Int, count: Int, x: Double, differentiate: Bool) throws -> (value: Double, derivative: Double) {
    var b1 = 0.0, b2 = 0.0, d1 = 0.0, d2 = 0.0
    if count > 1 {
        for k in stride(from: count - 1, through: 1, by: -1) {
            let b = coefficients[offset + k] + 2*x*b1 - b2
            let d = differentiate ? 2*b1 + 2*x*d1 - d2 : 0
            guard b.isFinite, d.isFinite else { throw SpiceSPKError.numericalRangeExceeded }
            b2 = b1; b1 = b; d2 = d1; d1 = d
        }
    }
    let value = coefficients[offset] + x*b1 - b2, derivative = differentiate ? b1 + x*d1 - d2 : 0
    guard value.isFinite, derivative.isFinite else { throw SpiceSPKError.numericalRangeExceeded }
    return (value, derivative)
}

extension SpiceDAFFile {
    /// Direct primitive state from a caller-selected SPK array.
    public func spkState(array index: Int, atSecondsTDB seconds: Double) throws -> SpiceSPKState {
        guard seconds.isFinite else { throw SpiceSPKError.invalidEpoch }
        return try SpiceSPKSegment(file: self, array: index).state(atSecondsTDB: seconds)
    }
}
extension SpiceKernelSnapshot {
    /// Selects the highest-priority covering target segment and evaluates its native primitive state.
    /// A corrupt or unsupported selected segment fails; lower-priority data is not substituted.
    /// Does not yet chain centres, transform axes, or apply aberration corrections.
    public func spkState(targetID: Int32, atSecondsTDB seconds: Double) throws -> SpiceSPKState {
        guard seconds.isFinite else { throw SpiceSPKError.invalidEpoch }
        var checks = Int.max
        if case .selected(let segment) = try lookupSPK(targetID: targetID, seconds: seconds, remainingChecks: &checks, checkLimit: Int.max) {
            return try segment.file.spkState(array: segment.arrayIndex, atSecondsTDB: seconds)
        }
        throw SpiceSPKError.noSegment(targetID: targetID, secondsTDB: seconds)
    }
}
