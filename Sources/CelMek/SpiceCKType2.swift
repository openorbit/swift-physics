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
@preconcurrency import Math

/// Continuous constant-rate CK intervals. Seconds/tick belongs to each pointing record.
public struct SpiceCKType2Segment {
    public let descriptor: SpiceCKDescriptor
    private let file: SpiceDAFFile
    private let array: Int
    private let starts: [Double]
    private let stops: [Double]
    public init(file: SpiceDAFFile, array: Int, maximumRecords: Int = 100_000) throws {
        guard file.metadata.arrays.indices.contains(array), let d = file.metadata.arrays[array].ck else { throw SpiceCKError.malformedSegment }
        guard d.dataType == 2 else { throw SpiceCKError.unsupportedType(d.dataType) }
        guard d.hasAngularVelocity else { throw SpiceCKError.malformedSegment }
        guard maximumRecords > 0, maximumRecords <= Int.max/11 else { throw SpiceCKError.resourceLimit }
        let a = file.metadata.arrays[array], length = a.finalAddress-a.initialAddress+1
        // The format has no count footer. Invert its strictly increasing length formula.
        var lo = 1, hi = min(maximumRecords,length/10)
        guard hi >= lo else { throw SpiceCKError.malformedSegment }
        while lo < hi {
            let mid = lo+(hi-lo)/2
            if 10*mid+(mid-1)/100 < length { lo = mid+1 } else { hi = mid }
        }
        let n = lo
        guard 10*n+(n-1)/100 == length else {
            if length > 10*maximumRecords+(maximumRecords-1)/100 { throw SpiceCKError.resourceLimit }
            throw SpiceCKError.malformedSegment
        }
        func read(_ offset: Int) throws -> [Double] {
            var result: [Double] = [], i = 0
            while i < n {
                let count = min(n-i,file.limits.maximumReadWords)
                result += try file.readWords(array: array,offset: offset+i,count: count); i += count
            }
            return result
        }
        let starts = try read(8*n), stops = try read(9*n)
        for i in 0..<n {
            guard starts[i].isFinite, stops[i].isFinite, starts[i] >= d.startEncodedTicks,
                  stops[i] <= d.endEncodedTicks, starts[i] < stops[i],
                  i == 0 || stops[i-1] <= starts[i] else { throw SpiceCKError.malformedSegment }
        }
        for i in 0..<((n-1)/100) {
            let j = (i+1)*100
            let expected = stops[j-1]+(starts[j]-stops[j-1])/2
            let value = try file.readWords(array: array,offset: 10*n+i,count: 1)[0]
            guard value.isFinite, abs(value-expected) <= 8*max(value.ulp,expected.ulp) else { throw SpiceCKError.malformedSegment }
        }
        self.file = file; self.array = array; self.descriptor = d; self.starts = starts; self.stops = stops
    }
    /// Inclusive intervals; at a shared endpoint the interval starting there wins.
    /// Gaps remain gaps even when the caller uses tolerance for other CK types.
    public func pointing(atEncodedTicks ticks: Double) throws -> SpiceCKPointing? {
        guard ticks.isFinite else { throw SpiceCKError.invalidQuery }
        var lo = 0, hi = starts.count
        while lo < hi { let mid = lo+(hi-lo)/2; if starts[mid] <= ticks { lo = mid+1 } else { hi = mid } }
        guard lo > 0 else { return nil }
        let i = lo-1
        guard ticks <= stops[i] else { return nil }
        let v = try file.readWords(array: array,offset: 8*i,count: 8)
        guard v.allSatisfy(\.isFinite), v[7] > 0, abs(v.prefix(4).reduce(0) { $0+$1*$1 }-1) <= 1e-10 else { throw SpiceCKError.malformedSegment }
        let initial = try attitudeQuaternion(Quaternion(x: v[1],y: v[2],z: v[3],w: v[0]))
        let omega = SIMD3(v[4],v[5],v[6]), magnitude = hypot(hypot(v[4],v[5]),v[6])
        let angle = magnitude*((ticks-starts[i])*v[7])
        guard angle.isFinite else { throw SpiceCKError.numericalRangeExceeded }
        // Physical reference-axis spin: inverse coordinate rotation multiplies C on the right.
        let delta = magnitude == 0 ? Quaternion.identity : Quaternion(angle: -angle.truncatingRemainder(dividingBy: 2 * .pi),axis: omega)
        return SpiceCKPointing(requestedEncodedTicks: ticks,encodedTicks: ticks,referenceToInstrument: try attitudeQuaternion(initial*delta),
            angularVelocityReferenceRadiansPerSecond: omega, secondsPerEncodedTick: v[7],
            descriptor: descriptor,sourcePath: file.url.path,recordIndex: i)
    }
}
