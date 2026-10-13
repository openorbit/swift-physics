// SPDX-License-Identifier: Apache-2.0
import Foundation
@preconcurrency import Math

/// CK type-3 fixed-axis interpolation, limited to declared interpolation intervals.
/// Stored angular velocities are interpolated independently of pointing.
public struct SpiceCKType3Segment {
    public let descriptor: SpiceCKDescriptor
    private let file: SpiceDAFFile
    private let array: Int
    private let width: Int
    private let times: [Double]
    private let intervalIndices: [Int]
    public init(file: SpiceDAFFile, array: Int, maximumRecords: Int = 100_000) throws {
        guard file.metadata.arrays.indices.contains(array), let d = file.metadata.arrays[array].ck else { throw SpiceCKError.malformedSegment }
        guard d.dataType == 3 else { throw SpiceCKError.unsupportedType(d.dataType) }
        guard maximumRecords > 0 else { throw SpiceCKError.resourceLimit }
        let a = file.metadata.arrays[array], length = a.finalAddress-a.initialAddress+1
        guard length >= 2 else { throw SpiceCKError.malformedSegment }
        let footer = try file.readWords(array: array,offset: length-2,count: 2)
        guard footer.allSatisfy({ $0.isFinite && $0 >= 1 && $0 < Double(Int.max) && $0.rounded() == $0 }) else { throw SpiceCKError.malformedSegment }
        guard footer[1] <= Double(maximumRecords) else { throw SpiceCKError.resourceLimit }
        let n = Int(footer[1]), m = Int(footer[0]), width = d.hasAngularVelocity ? 7 : 4
        guard m <= n, n <= (length-2)/(width+1) else { throw SpiceCKError.malformedSegment }
        let timeDirectory = (n-1)/100, startDirectory = (m-1)/100
        guard length == n*(width+1)+timeDirectory+m+startDirectory+2 else { throw SpiceCKError.malformedSegment }
        func read(_ offset: Int, _ count: Int) throws -> [Double] {
            var values: [Double] = [], i = 0
            while i < count {
                let size = min(count-i,file.limits.maximumReadWords)
                values += try file.readWords(array: array,offset: offset+i,count: size); i += size
            }
            return values
        }
        let times = try read(n*width,n), starts = try read(n*(width+1)+timeDirectory,m)
        guard times.allSatisfy({ $0.isFinite && $0 >= d.startEncodedTicks && $0 <= d.endEncodedTicks }),
              zip(times,times.dropFirst()).allSatisfy({ $0 < $1 }),
              starts.allSatisfy(\.isFinite), starts.first == times.first,
              zip(starts,starts.dropFirst()).allSatisfy({ $0 < $1 }) else { throw SpiceCKError.malformedSegment }
        var indices: [Int] = [], cursor = 0
        for start in starts {
            while cursor < n, times[cursor] < start { cursor += 1 }
            guard cursor < n, times[cursor] == start else { throw SpiceCKError.malformedSegment }
            indices.append(cursor)
        }
        for (values,offset,count) in [(times,n*(width+1),timeDirectory),(starts,n*(width+1)+timeDirectory+m,startDirectory)] {
            for i in 0..<count {
                let value = try file.readWords(array: array,offset: offset+i,count: 1)[0]
                guard value == values[(i+1)*100-1] else { throw SpiceCKError.malformedSegment }
            }
        }
        self.file = file; self.array = array; self.width = width; self.times = times
        self.intervalIndices = indices; self.descriptor = d
    }
    /// Nearest gap endpoints preserve their own timestamps. Equal-distance ties use the later endpoint.
    /// Tolerance applies within descriptor bounds, matching the library's bounded search policy.
    public func pointing(atEncodedTicks ticks: Double, toleranceTicks: Double = 0) throws -> SpiceCKPointing? {
        guard ticks.isFinite, toleranceTicks.isFinite, toleranceTicks >= 0 else { throw SpiceCKError.invalidQuery }
        guard ticks >= descriptor.startEncodedTicks, ticks <= descriptor.endEncodedTicks else { return nil }
        func upperBound(_ values: [Double], _ value: Double) -> Int {
            var lo = 0, hi = values.count
            while lo < hi { let mid = lo+(hi-lo)/2; if values[mid] <= value { lo = mid+1 } else { hi = mid } }
            return lo
        }
        // Locate the interval by its start record; no interpolation across a gap.
        var lo = 0, hi = intervalIndices.count
        while lo < hi { let mid = lo+(hi-lo)/2; if times[intervalIndices[mid]] <= ticks { lo = mid+1 } else { hi = mid } }
        let interval = lo-1
        var selected = ticks, first: Int, second: Int
        if interval >= 0 {
            let start = intervalIndices[interval]
            let end = interval+1 < intervalIndices.count ? intervalIndices[interval+1]-1 : times.count-1
            if ticks <= times[end] {
                first = max(start,min(end,upperBound(times,ticks)-1))
                second = min(first+1,end)
            } else {
                first = end; second = end
                if interval+1 < intervalIndices.count {
                    let next = intervalIndices[interval+1]
                    if times[next]-ticks <= ticks-times[end] { first = next; second = next }
                }
                selected = times[first]
            }
        } else { first = 0; second = 0; selected = times[0] }
        guard abs(selected-ticks) <= toleranceTicks else { return nil }
        func record(_ i: Int) throws -> (Quaternion,SIMD3<Double>?) {
            let v = try file.readWords(array: array,offset: i*width,count: width)
            guard v.allSatisfy(\.isFinite), abs(v.prefix(4).reduce(0) { $0+$1*$1 }-1) <= 1e-10 else { throw SpiceCKError.malformedSegment }
            return (try attitudeQuaternion(Quaternion(x: v[1],y: v[2],z: v[3],w: v[0])),width == 7 ? SIMD3(v[4],v[5],v[6]) : nil)
        }
        let a = try record(first)
        var q = a.0, rate = a.1
        if second != first, selected != times[first] {
            let b = try record(second), weight = (selected-times[first])/(times[second]-times[first])
            // Physical reference-axis relative rotation, shortest arc. Canonicalize exact half turns.
            var relative = try attitudeQuaternion(b.0.conj*a.0)
            let components = [relative.w,relative.x,relative.y,relative.z]
            if components.first(where: { $0 != 0 }).map({ $0 < 0 }) == true {
                relative = Quaternion(x: -relative.x,y: -relative.y,z: -relative.z,w: -relative.w)
            }
            let magnitude = hypot(hypot(relative.x,relative.y),relative.z)
            if magnitude > 0 {
                let angle = 2*atan2(magnitude,relative.w)*weight
                q = try attitudeQuaternion(a.0*Quaternion(angle: -angle,axis: relative.vector))
            }
            if let av = a.1, let bv = b.1 {
                rate = av*(1-weight)+bv*weight
                guard attitudeFinite(rate!) else { throw SpiceCKError.numericalRangeExceeded }
            }
        }
        return SpiceCKPointing(requestedEncodedTicks: ticks,encodedTicks: selected,referenceToInstrument: q,
            angularVelocityReferenceRadiansPerSecond: rate,secondsPerEncodedTick: nil,
            descriptor: descriptor,sourcePath: file.url.path,recordIndex: first)
    }
}
