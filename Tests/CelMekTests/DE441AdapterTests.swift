// SPDX-License-Identifier: Apache-2.0
import Foundation
import Synchronization
import Math
import Testing
@testable import CelMek

@Suite struct DE441AdapterTests {
    private let jd = 2451545.0
    private func date(_ offset: Double) throws -> JulianDateComponents {
        try JulianDateComponents(wholeDay: 2451545, fractionalDay: offset, scale: .tdb)
    }
    private struct LinearTDB: TDBTimeModel {
        let coverage: ClosedRange<Epoch>
        let provenance = "Test affine TDB model"
        func correction(at epoch: Epoch) throws -> TDBCorrection {
            let t = try epoch.duration(since: .j2000).seconds
            return try TDBCorrection(epoch: epoch, tdbMinusTT: Duration(seconds: 0.1*t), rate: 1.1, acceleration: 0)
        }
    }

    @Test func splitDatesRetainFinePrecisionAndEnforceScale() throws {
        let url = try DE441EphemerisTests().fixture()
        defer { try? FileManager.default.removeItem(at: url) }
        let reader = try DE441Ephemeris(url: url)
        let first = try date(6), second = try date(6 + 1e-10)
        #expect(first.value == second.value) // A combined JD loses the distinction.
        let a = try reader.vector(target: .mercury, center: .solarSystemBarycentre, at: first)
        let b = try reader.vector(target: .mercury, center: .solarSystemBarycentre, at: second)
        #expect(b.position.x > a.position.x)
        let deltaKm = (b.position.x-a.position.x)*reader.metadata.astronomicalUnitKilometres
        #expect(abs(deltaKm-2e-10) < 1e-15)
        #expect(throws: DE441Error.self) {
            try reader.vector(target: .sun, center: .sun, at: JulianDateComponents(value: jd, scale: .tt))
        }
        let end = try JulianDateComponents(wholeDay: 2451609, fractionalDay: 1e-10, scale: .tdb)
        #expect(throws: DE441Error.self) { try reader.vector(target: .sun, center: .sun, at: end) }
    }

    @Test func commonSIAdapterUsesCoordinateTimeRate() throws {
        let url = try DE441EphemerisTests().fixture()
        defer { try? FileManager.default.removeItem(at: url) }
        let reader = try DE441Ephemeris(url: url)
        let upper = try Epoch.j2000.advanced(by: Duration(seconds: 1e6))
        let model = LinearTDB(coverage: .j2000...upper)
        let provider = try DE441StateProvider(reader: reader, target: .mercury,
            coverage: .j2000...upper, tdbModel: model, accuracy: .low)
        let epoch = try Epoch.j2000.advanced(by: Duration(seconds: 6*86400/1.1))
        let state = try provider.checkedState(at: epoch)
        #expect(abs(state.position.x+500) < 1e-7)
        #expect(abs(state.velocity.x-2000/86400*1.1) < 1e-12)
        #expect(state.acceleration == nil)
        #expect(provider.metadata.body == .mercury)
        #expect(provider.metadata.origin == .solarSystemBarycentre)
        #expect(state.frame === CelestialFrames.icrf)
        #expect(provider.accuracy == .low && provider.effectiveAccuracy == .full)
        let collection = try EphemerisCollection(providers: [provider])
        #expect(try collection.state(of: .mercury, at: epoch) == state)
        #expect(throws: EphemerisError.outsideCoverage) {
            try provider.state(at: upper.advanced(by: Duration(seconds: 1)))
        }
        #expect(throws: EphemerisError.invalidMetadata) {
            try DE441StateProvider(reader: reader, target: .moon, center: .earth, coverage: .j2000...upper, tdbModel: model)
        }
        let frame = FrameTree(name: "Earth-centred ICRF axes").root
        let lunar = try DE441StateProvider(reader: reader, target: .moon, center: .earth, frame: frame,
                                           coverage: .j2000...upper, tdbModel: model)
        #expect(lunar.metadata.origin == .earth)
        #expect(try lunar.state(at: .j2000).frame === frame)
    }

    @Test func reconstructionUsesFileAUAndMassRatio() throws {
        let url = try DE441EphemerisTests().fixture()
        defer { try? FileManager.default.removeItem(at: url) }
        var data = try Data(contentsOf: url)
        for (offset,value) in [(2680,1000.0),(8152,1000.0),(2688,40.0),(8160,40.0)] {
            var bits = value.bitPattern.littleEndian
            withUnsafeBytes(of: &bits) { data.replaceSubrange(offset..<offset+8, with: $0) }
        }
        try data.write(to: url)
        let reader = try DE441Ephemeris(url: url)
        let earth = try reader.vector(target: .earth, center: .solarSystemBarycentre, at: jd+2)
        let relative = try reader.vector(target: .moon, center: .earth, at: jd+2)
        #expect(abs(earth.position.x-(198.875-897.0/41)/1000) < 1e-15)
        #expect(abs(relative.position.x-0.897) < 1e-15)
        #expect(abs(relative.velocity.x-0.001) < 1e-15)
    }

    @Test func nonfiniteInteriorCoefficientsAndArithmeticOverflowFail() throws {
        let url = try DE441EphemerisTests().fixture(spans: [(jd,jd+32),(jd+32,jd+64),(jd+64,jd+96)])
        defer { try? FileManager.default.removeItem(at: url) }
        let original = try Data(contentsOf: url)
        for coefficients in [[Double.nan,2,3],[1e308,1e308,1e308]] {
            var data = original
            for (index,value) in coefficients.enumerated() {
                var bits = value.bitPattern.littleEndian
                let offset = 3*8144+(2+index)*8
                withUnsafeBytes(of: &bits) { data.replaceSubrange(offset..<offset+8, with: $0) }
            }
            try data.write(to: url)
            let reader = try DE441Ephemeris(url: url)
            #expect(throws: DE441Error.self) {
                try reader.vector(target: .mercury, center: .solarSystemBarycentre, at: jd+38)
            }
        }
    }

    @Test func boundedLRUAndDisabledCacheHaveIdenticalResults() throws {
        let url = try DE441EphemerisTests().fixture(spans: [(jd,jd+32),(jd+32,jd+64),(jd+64,jd+96)])
        defer { try? FileManager.default.removeItem(at: url) }
        let cached = try DE441Ephemeris(url: url, maximumCachedRecords: 2)
        let uncached = try DE441Ephemeris(url: url, maximumCachedRecords: 0)
        cached.clearRecordCache()
        let before = cached.cacheStatistics.diskReads
        for offset in [1.0,33,1,65,1,33] {
            let a = try cached.vector(target: .mercury, center: .sun, at: jd+offset)
            let b = try uncached.vector(target: .mercury, center: .sun, at: jd+offset)
            #expect(a.position == b.position && a.velocity == b.velocity)
            #expect(cached.cacheStatistics.cachedRecords <= 2)
        }
        #expect(cached.cacheStatistics.diskReads-before == 4)
        #expect(uncached.cacheStatistics.cachedRecords == 0)
        #expect(throws: DE441Error.self) { try DE441Ephemeris(url: url, maximumCachedRecords: -1) }
    }

    // Sendability of the reader relies on its internal file/cache lock; the
    // failure collection is protected independently by Mutex.
    private final class ConcurrentQueries: @unchecked Sendable {
        let reader: DE441Ephemeris
        private let failedIndices = Mutex<[Int]>([])

        init(reader: DE441Ephemeris) { self.reader = reader }

        func recordFailure(_ index: Int) {
            failedIndices.withLock { $0.append(index) }
        }

        var failures: [Int] {
            failedIndices.withLock { $0 }
        }
    }

    @Test func sharedReaderConcurrentQueriesAreRepeatable() throws {
        let url = try DE441EphemerisTests().fixture()
        defer { try? FileManager.default.removeItem(at: url) }
        let reader = try DE441Ephemeris(url: url, maximumCachedRecords: 1)
        let queries = ConcurrentQueries(reader: reader)
        let julianDate = jd
        DispatchQueue.concurrentPerform(iterations: 100) { index in
            do {
                let offset = index % 2 == 0 ? 6.0 : 38.0
                let value = try queries.reader.vector(target: .mercury, center: .solarSystemBarycentre, at: julianDate+offset)
                let expected = (index % 2 == 0 ? -0.5 : 999.5)/queries.reader.metadata.astronomicalUnitKilometres
                if abs(value.position.x-expected) > 1e-15 { queries.recordFailure(index) }
            } catch { queries.recordFailure(index) }
        }
        #expect(queries.failures.isEmpty)
    }
}
