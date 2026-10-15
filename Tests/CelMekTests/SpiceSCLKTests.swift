// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
@testable import CelMek

@Suite struct SpiceSCLKTests {
    let helpers = SpiceTextKernelTests()
    var provenance: PhysicalDataProvenance { try! PhysicalDataProvenance(source: "Synthetic documented SCLK",version: "clock-1") }
    let text = """
    KPL/SCLK
    \\begindata
    SCLK_KERNEL_ID=@04-SEP-1990//4:23:00
    SCLK_DATA_TYPE_77=1
    SCLK01_N_FIELDS_77=2
    SCLK01_MODULI_77=(100 10)
    SCLK01_OFFSETS_77=(0 1)
    SCLK01_OUTPUT_DELIM_77=2
    SCLK_PARTITION_START_77=(100 50)
    SCLK_PARTITION_END_77=(200 150)
    SCLK01_COEFFICIENTS_77=(0 0 20 100 200 10)
    """
    func clock(_ input: String? = nil, coverage: ClosedRange<Double> = 0...200) throws -> SpiceSCLKClock {
        try SpiceSCLKClock(snapshot: helpers.snapshot(input ?? text),spacecraftID: -77,encodedTickCoverage: coverage,provenance: provenance)
    }
    func epoch(_ seconds: Double) throws -> Epoch { Epoch(offsetFromJ2000: try CelMek.Duration(seconds: seconds)) }
    @Test func partitionResetsOffsetsAndCanonicalStrings() throws {
        let model = try clock()
        let first = SpiceSCLKReading(partition: 1,fields: [15,1]), second = SpiceSCLKReading(partition: 2,fields: [15,1])
        #expect(try model.encodedTicks(for: first) == 50)
        #expect(try model.encodedTicks(for: second) == 200)
        #expect(try model.string(for: first) == "1/15:1")
        #expect(try model.reading(from: " 1 / 15 : 1 ") == first)
        #expect(try model.reading(atEncodedTicks: 50,partition: 1) == first)
        #expect(try model.encodedTicks(for: SpiceSCLKReading(partition: 1,fields: [20,1])) == 100)
        #expect(try model.encodedTicks(for: SpiceSCLKReading(partition: 2,fields: [5,1])) == 100)
        #expect(try model.reading(atEncodedTicks: 100,partition: 1).fields == [20,1])
        #expect(try model.reading(atEncodedTicks: 100,partition: 2).fields == [5,1])
        #expect(model.cacheScope != EphemerisCacheScope())
        #expect(model.provenance.source.contains("kernel.txt") && model.provenance.version == "clock-1")
        for string in ["15:1","1/15::1","1/15:0","1/15:11","3/15:1","1/99:1","1/15:1:2","1/-15:1"] {
            #expect(throws: (any Error).self) { try model.reading(from: string) }
        }
        #expect(throws: SpiceSCLKError.invalidPartition) { try model.encodedTicks(for: SpiceSCLKReading(partition: Int.min,fields: [15,1])) }
        #expect(throws: SpiceSCLKError.invalidPartition) { try model.reading(atEncodedTicks: 50,partition: Int.min) }
        #expect(throws: SpiceSCLKError.self) { try model.reading(atEncodedTicks: 50.5,partition: 1) }
    }
    @Test func correlationRatesUseMostSignificantFieldAndPreserveFractions() throws {
        let model = try clock()
        for (ticks,time,rate,index) in [(0.0,0.0,2.0,0),(50.25,100.5,2.0,0),(100.0,200.0,1.0,1),(150.125,250.125,1.0,1),(200.0,300.0,1.0,1)] {
            let sample = try model.sample(atEncodedTicks: ticks)
            #expect(sample.parallelSeconds == time)
            #expect(sample.parallelSecondsPerTick == rate && sample.correlationIndex == index)
            #expect(try model.sample(atParallelSeconds: time).encodedTicks == ticks)
        }
        #expect(try model.sample(atEncodedTicks: 100.0.nextDown).correlationIndex == 0)
        #expect(throws: SpiceSCLKError.outsideCoverage) { try model.sample(atEncodedTicks: 200.0.nextUp) }
        #expect(throws: SpiceSCLKError.outsideCoverage) { try model.sample(atEncodedTicks: .nan) }
        let limited = try clock(coverage: 20...180)
        #expect(throws: SpiceSCLKError.outsideCoverage) { try limited.sample(atEncodedTicks: 19) }
        #expect(throws: SpiceSCLKError.correlationGap) { try limited.sample(atParallelSeconds: 0) }
    }
    @Test func discontinuityGapsAndAmbiguityAreNotHidden() throws {
        let gap = try clock(text.replacingOccurrences(of: "100 200 10", with: "100 250 10"))
        #expect(try gap.sample(atEncodedTicks: 100).parallelSeconds == 250)
        #expect(throws: SpiceSCLKError.correlationGap) { try gap.sample(atParallelSeconds: 225) }
        let overlap = try clock(text.replacingOccurrences(of: "100 200 10", with: "100 150 10"))
        #expect(throws: SpiceSCLKError.ambiguousCorrelation) { try overlap.sample(atParallelSeconds: 175) }
        #expect(try overlap.sample(atParallelSeconds: 100).encodedTicks == 50)
    }
    @Test func tdbAndTTCorrelationUseExplicitEpochModels() throws {
        let model = try clock(), coverage = try epoch(-100)...epoch(1000)
        let time = SpiceSPKAdapterTests.Clock(coverage: coverage,offset: 0.25,drift: 0.001)
        let e = try model.epoch(atEncodedTicks: 50,using: time)
        #expect(abs(e.offsetFromJ2000.seconds - 99.75/1.001) < 1e-10)
        #expect(abs(try model.sample(at: e,using: time).encodedTicks - 50) < 1e-10)
        #expect(throws: SpiceSCLKError.tdbModelRequired) { try model.epoch(atEncodedTicks: 50) }
        #expect(throws: SpiceSCLKError.tdbModelRequired) { try model.sample(at: e) }
        let tt = try clock(text + "\nSCLK01_TIME_SYSTEM_77=2")
        #expect(tt.parallelTimeScale == .tt)
        #expect(try tt.epoch(atEncodedTicks: 50).offsetFromJ2000.seconds == 100)
        #expect(try tt.sample(at: epoch(100)).encodedTicks == 50)
        let stale = SpiceSPKAdapterTests.Clock(coverage: coverage,stale: true)
        #expect(throws: TDBModelError.epochMismatch) { try model.sample(at: e,using: stale) }
    }
    @Test func malformedConstantsAndBoundedRepresentation() throws {
        for (old,new) in [("TYPE_77=1","TYPE_77=2"),("FIELDS_77=2","FIELDS_77=17"),
                          ("(100 10)","(100 0)"),("(100 10)","(1e12 1e12)"),
                          ("(0 1)","(0 -1)"),("DELIM_77=2","DELIM_77=6"),
                          ("(200 150)","(99 150)"),("100 200 10","0 200 10"),
                          ("100 200 10","100 200 0"),("100 200 10","100 200") ] {
            #expect(throws: (any Error).self) { try clock(text.replacingOccurrences(of: old,with: new)) }
        }
        #expect(throws: SpiceSCLKError.unsupportedTimeSystem(3)) { try clock(text+"\nSCLK01_TIME_SYSTEM_77=3") }
        #expect(throws: SpiceSCLKError.outsideCoverage) { try clock(coverage: 0...201) }
        #expect(throws: SpiceSCLKError.invalidSpacecraftID) {
            try SpiceSCLKClock(snapshot: helpers.snapshot(text),spacecraftID: 77,encodedTickCoverage: 0...200,provenance: provenance)
        }
    }
    @Test func dateMetadataSupportsDocumentedSCLKSpelling() throws {
        #expect(try spiceDateSeconds("@04-SEP-1990//4:23:00") == spiceDateSeconds("@1990-SEP-4-04:23:00"))
        #expect(try spiceDateSeconds("@04-SEP-1990/4:23") == spiceDateSeconds("@1990-SEP-4-04:23:00"))
        #expect(try spiceDateSeconds("@2000-JAN-1-12") == 0)
        #expect(throws: (any Error).self) { try spiceDateSeconds("@04-SEP-1990//4:99:00") }
    }
    @Test func kernelPrecedenceUnloadAndRetainedClock() throws {
        let helper = SpiceKernelSetTests(), dir = try helper.directory(); defer { try? FileManager.default.removeItem(at: dir) }
        let first = try helper.write(text,"first.tsc",in: dir)
        let second = try helper.write("KPL/SCLK\n\\begindata\nSCLK01_COEFFICIENTS_77=(0 10 20 100 210 10)","second.tsc",in: dir)
        let manifest = try helper.write("KPL/MK\n\\begindata\nKERNELS_TO_LOAD=('first.tsc')","clock.tm",in: dir)
        let set = try SpiceKernelSet(); try set.load(manifest)
        let snapshot = set.snapshot(), old = try SpiceSCLKClock(snapshot: snapshot,spacecraftID: -77,encodedTickCoverage: 0...200,provenance: provenance)
        #expect(snapshot.loads[0].sources.contains { $0.path == resolvedSpiceTestPath(first) && $0.kind == "SCLK" })
        let id = try set.load(second)
        let current = try SpiceSCLKClock(snapshot: set.snapshot(),spacecraftID: -77,encodedTickCoverage: 0...200,provenance: provenance)
        #expect(try current.sample(atEncodedTicks: 0).parallelSeconds == 10)
        #expect(current.provenance.source.contains("first.tsc") && current.provenance.source.contains("second.tsc"))
        try set.unload(id)
        #expect(try old.sample(atEncodedTicks: 0).parallelSeconds == 0)
        #expect(old.cacheScope != current.cacheScope)
        #expect(set.snapshot()["SCLK01_COEFFICIENTS_77"]?.value == snapshot["SCLK01_COEFFICIENTS_77"]?.value)
    }
}
