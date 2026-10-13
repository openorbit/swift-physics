// SPDX-License-Identifier: Apache-2.0
import Foundation

public enum SpiceSCLKError: Error, Equatable {
    case invalidSpacecraftID, invalidFormat, invalidPartition, invalidFields
    case unsupportedClockType(Double), unsupportedTimeSystem(Double)
    case invalidCorrelation, outsideCoverage, correlationGap, ambiguousCorrelation
    case numericalRangeExceeded, tdbModelRequired
}

public struct SpiceSCLKReading: Equatable, Sendable {
    /// Explicit one-based partition; reset readings are never resolved by guessing.
    public let partition: Int
    public let fields: [Int64]
    public init(partition: Int, fields: [Int64]) { self.partition = partition; self.fields = fields }
}
public struct SpiceSCLKPartition: Sendable {
    public let startRawTicks: Double
    public let endRawTicks: Double
    public let startEncodedTicks: Double
    public var endEncodedTicks: Double { startEncodedTicks + (endRawTicks - startRawTicks) }
}
public struct SpiceSCLKCorrelation: Sendable {
    public let encodedTicks: Double
    /// Seconds past J2000 noon in the declared parallel time system.
    public let parallelSeconds: Double
    public let parallelSecondsPerTick: Double
}
public struct SpiceSCLKSample: Sendable {
    public let encodedTicks: Double
    public let parallelSeconds: Double
    public let parallelSecondsPerTick: Double
    public let correlationIndex: Int
}

/// Immutable type-1 numeric/string clock and correlation data from a snapshot.
/// Counts use Double encoded ticks as specified by SPICE; precision is finite.
public struct SpiceSCLKClock: Sendable {
    public let spacecraftID: Int32
    public let parallelTimeScale: TimeScale
    public let moduli: [Int64]
    public let offsets: [Int64]
    public let outputDelimiter: String
    public let partitions: [SpiceSCLKPartition]
    public let correlations: [SpiceSCLKCorrelation]
    public let encodedTickCoverage: ClosedRange<Double>
    public let provenance: PhysicalDataProvenance
    public let cacheScope: EphemerisCacheScope

    public init(snapshot: SpiceKernelSnapshot, spacecraftID: Int32,
                encodedTickCoverage: ClosedRange<Double>, provenance: PhysicalDataProvenance) throws {
        guard spacecraftID < 0 else { throw SpiceSCLKError.invalidSpacecraftID }
        let suffix = String(-Int64(spacecraftID)), p = "SCLK01_"
        var names: [String] = []
        func values(_ name: String, count: Int? = nil) throws -> [Double] {
            names.append(name); return try snapshot.numeric(name,count: count)
        }
        let type = try values("SCLK_DATA_TYPE_" + suffix,count: 1)[0]
        guard type == 1 else { throw SpiceSCLKError.unsupportedClockType(type) }
        let countValue = try values(p + "N_FIELDS_" + suffix,count: 1)[0]
        guard let count = Int(exactly: countValue), (1...16).contains(count) else { throw SpiceSCLKError.invalidFormat }
        let m = try values(p + "MODULI_" + suffix,count: count)
        let o = try values(p + "OFFSETS_" + suffix,count: count)
        guard m.allSatisfy({ Int64(exactly: $0) != nil && $0 > 0 }), o.allSatisfy({ Int64(exactly: $0) != nil && $0 >= 0 }) else { throw SpiceSCLKError.invalidFormat }
        let moduli = m.map { Int64($0) }, offsets = o.map { Int64($0) }
        var total = 1.0, mostSignificantTicks = 1.0
        for (i, modulus) in m.enumerated() {
            total *= modulus
            if i > 0 { mostSignificantTicks *= modulus }
            guard total.isFinite, total <= 9_007_199_254_740_992 else { throw SpiceSCLKError.invalidFormat }
        }
        let delimiterCode = try values(p + "OUTPUT_DELIM_" + suffix,count: 1)[0]
        guard let delimiter = Int(exactly: delimiterCode), (1...5).contains(delimiter) else { throw SpiceSCLKError.invalidFormat }
        let timeName = p + "TIME_SYSTEM_" + suffix
        let timeCode = try snapshot[timeName].map { _ in try values(timeName,count: 1)[0] } ?? 1
        guard timeCode == 1 || timeCode == 2 else { throw SpiceSCLKError.unsupportedTimeSystem(timeCode) }
        let starts = try values("SCLK_PARTITION_START_" + suffix), ends = try values("SCLK_PARTITION_END_" + suffix)
        guard starts.count == ends.count else { throw SpiceSCLKError.invalidPartition }
        var partitions: [SpiceSCLKPartition] = [], encodedEnd = 0.0
        for (start,end) in zip(starts,ends) {
            guard start >= 0, start.rounded() == start, end.rounded() == end, end > start, end < total else { throw SpiceSCLKError.invalidPartition }
            let partition = SpiceSCLKPartition(startRawTicks: start,endRawTicks: end,startEncodedTicks: encodedEnd)
            let next = partition.endEncodedTicks
            guard next > encodedEnd, next <= 9_007_199_254_740_992 else { throw SpiceSCLKError.invalidPartition }
            encodedEnd = next; partitions.append(partition)
        }
        let raw = try values(p + "COEFFICIENTS_" + suffix)
        guard raw.count % 3 == 0 else { throw SpiceSCLKError.invalidCorrelation }
        var correlations: [SpiceSCLKCorrelation] = []
        for i in stride(from: 0,to: raw.count,by: 3) {
            let tick = raw[i], time = raw[i+1], rate = raw[i+2]/mostSignificantTicks
            guard tick >= 0, tick <= encodedEnd, rate.isFinite, rate > 0,
                  abs(time) < 1e15, correlations.last.map({ tick > $0.encodedTicks }) ?? true else { throw SpiceSCLKError.invalidCorrelation }
            correlations.append(SpiceSCLKCorrelation(encodedTicks: tick,parallelSeconds: time,parallelSecondsPerTick: rate))
        }
        guard encodedTickCoverage.lowerBound.isFinite, encodedTickCoverage.upperBound.isFinite,
              encodedTickCoverage.lowerBound >= correlations[0].encodedTicks,
              encodedTickCoverage.upperBound <= encodedEnd else { throw SpiceSCLKError.outsideCoverage }
        self.spacecraftID = spacecraftID; self.parallelTimeScale = timeCode == 1 ? .tdb : .tt
        self.moduli = moduli; self.offsets = offsets; self.outputDelimiter = [".",":","-",","," "][delimiter-1]
        self.partitions = partitions; self.correlations = correlations; self.encodedTickCoverage = encodedTickCoverage
        self.provenance = try snapshot.physicalProvenance(names: names,supplied: provenance); cacheScope = snapshot.cacheScope
        _ = try sample(atEncodedTicks: encodedTickCoverage.lowerBound)
        _ = try sample(atEncodedTicks: encodedTickCoverage.upperBound)
    }
    public func encodedTicks(for reading: SpiceSCLKReading) throws -> Double {
        guard (1...partitions.count).contains(reading.partition) else { throw SpiceSCLKError.invalidPartition }
        guard reading.fields.count == moduli.count else { throw SpiceSCLKError.invalidFields }
        var raw = 0.0
        for i in moduli.indices {
            let (value, overflow) = reading.fields[i].subtractingReportingOverflow(offsets[i])
            guard !overflow, value >= 0, value < moduli[i] else { throw SpiceSCLKError.invalidFields }
            raw = raw * Double(moduli[i]) + Double(value)
        }
        let partition = partitions[reading.partition-1]
        guard raw >= partition.startRawTicks, raw <= partition.endRawTicks else { throw SpiceSCLKError.invalidPartition }
        let ticks = partition.startEncodedTicks + (raw-partition.startRawTicks)
        guard encodedTickCoverage.contains(ticks) else { throw SpiceSCLKError.outsideCoverage }
        return ticks
    }
    /// Continuous input is preserved by correlation APIs. Civil readings require
    /// integral ticks, and a caller selects the partition at shared reset boundaries.
    public func reading(atEncodedTicks ticks: Double, partition number: Int) throws -> SpiceSCLKReading {
        guard ticks.isFinite, encodedTickCoverage.contains(ticks), ticks.rounded() == ticks else { throw SpiceSCLKError.outsideCoverage }
        guard (1...partitions.count).contains(number) else { throw SpiceSCLKError.invalidPartition }
        let partition = partitions[number-1]
        guard ticks >= partition.startEncodedTicks, ticks <= partition.endEncodedTicks else { throw SpiceSCLKError.invalidPartition }
        var raw = Int64(partition.startRawTicks + (ticks-partition.startEncodedTicks)), fields = Array(repeating: Int64(0),count: moduli.count)
        for i in moduli.indices.reversed() {
            let value = raw % moduli[i]; raw /= moduli[i]
            let (field, overflow) = value.addingReportingOverflow(offsets[i])
            guard !overflow else { throw SpiceSCLKError.invalidFields }; fields[i] = field
        }
        guard raw == 0 else { throw SpiceSCLKError.invalidFields }
        return SpiceSCLKReading(partition: number,fields: fields)
    }
    /// Strict canonical subset: explicit partition and all integral fields.
    public func reading(from string: String) throws -> SpiceSCLKReading {
        guard string.utf8.count <= 512 else { throw SpiceSCLKError.invalidFields }
        let parts = string.split(separator: "/",omittingEmptySubsequences: false)
        guard parts.count == 2, let partition = Int(parts[0].trimmingCharacters(in: .whitespaces)) else { throw SpiceSCLKError.invalidPartition }
        let fieldString = parts[1].trimmingCharacters(in: .whitespaces)
        guard fieldString.range(of: #"^[0-9]+(?:\s*(?:[.:,-]\s*|\s+)[0-9]+)*$"#, options: .regularExpression) != nil else {
            throw SpiceSCLKError.invalidFields
        }
        let tokens = fieldString.split(whereSeparator: { " .:,-".contains($0) })
        guard tokens.count == moduli.count else { throw SpiceSCLKError.invalidFields }
        let fields = try tokens.map { token -> Int64 in
            guard token.allSatisfy(\.isNumber), let value = Int64(token) else { throw SpiceSCLKError.invalidFields }; return value
        }
        let reading = SpiceSCLKReading(partition: partition,fields: fields)
        _ = try encodedTicks(for: reading); return reading
    }
    public func string(for reading: SpiceSCLKReading) throws -> String {
        _ = try encodedTicks(for: reading)
        return String(reading.partition) + "/" + reading.fields.map(String.init).joined(separator: outputDelimiter)
    }
    public func sample(atEncodedTicks ticks: Double) throws -> SpiceSCLKSample {
        guard ticks.isFinite, encodedTickCoverage.contains(ticks) else { throw SpiceSCLKError.outsideCoverage }
        var low = 0, high = correlations.count
        while low < high { let mid = low+(high-low)/2; if correlations[mid].encodedTicks <= ticks { low = mid+1 } else { high = mid } }
        let index = low-1, correlation = correlations[index]
        let time = correlation.parallelSeconds + correlation.parallelSecondsPerTick*(ticks-correlation.encodedTicks)
        guard time.isFinite, abs(time) < 1e15 else { throw SpiceSCLKError.numericalRangeExceeded }
        return SpiceSCLKSample(encodedTicks: ticks,parallelSeconds: time,parallelSecondsPerTick: correlation.parallelSecondsPerTick,correlationIndex: index)
    }
    /// Inverse requires a unique covered branch. Discontinuous correlations may
    /// leave gaps or duplicate times; neither is silently resolved.
    public func sample(atParallelSeconds time: Double) throws -> SpiceSCLKSample {
        guard time.isFinite else { throw SpiceSCLKError.outsideCoverage }
        var result: SpiceSCLKSample?
        for i in correlations.indices {
            let c = correlations[i], lower = max(c.encodedTicks,encodedTickCoverage.lowerBound)
            let upper = min(i+1 < correlations.count ? correlations[i+1].encodedTicks : encodedTickCoverage.upperBound,encodedTickCoverage.upperBound)
            guard lower <= upper else { continue }
            let ticks = c.encodedTicks + (time-c.parallelSeconds)/c.parallelSecondsPerTick
            let exclusiveUpper = i+1 < correlations.count && upper == correlations[i+1].encodedTicks
            guard ticks.isFinite, ticks >= lower, ticks <= upper, !(exclusiveUpper && ticks == upper) else { continue }
            let sample = try sample(atEncodedTicks: ticks)
            if result != nil { throw SpiceSCLKError.ambiguousCorrelation }; result = sample
        }
        guard let result else { throw SpiceSCLKError.correlationGap }; return result
    }
    public func epoch(atEncodedTicks ticks: Double, using tdbModel: (any TDBTimeModel)? = nil) throws -> Epoch {
        let time = try sample(atEncodedTicks: ticks).parallelSeconds
        let nominal = Epoch(offsetFromJ2000: try Duration(seconds: time))
        if parallelTimeScale == .tt { return nominal }
        guard let tdbModel else { throw SpiceSCLKError.tdbModelRequired }
        return try Epoch.invertTDB(try nominal.julianDate(in: .tt).withScale(.tdb),using: tdbModel)
    }
    public func sample(at epoch: Epoch, using tdbModel: (any TDBTimeModel)? = nil) throws -> SpiceSCLKSample {
        let time: Double
        if parallelTimeScale == .tt { time = epoch.offsetFromJ2000.seconds }
        else {
            guard let tdbModel else { throw SpiceSCLKError.tdbModelRequired }
            let correction = try epoch.tdbCorrection(using: tdbModel)
            time = try epoch.offsetFromJ2000.adding(correction.tdbMinusTT).seconds
        }
        return try sample(atParallelSeconds: time)
    }
}

private extension JulianDateComponents {
    func withScale(_ scale: TimeScale) throws -> JulianDateComponents {
        try JulianDateComponents(wholeDay: wholeDay,fractionalDay: fractionalDay,scale: scale)
    }
}
