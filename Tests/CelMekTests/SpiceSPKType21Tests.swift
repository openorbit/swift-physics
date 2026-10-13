// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
#if canImport(CryptoKit)
import CryptoKit
#endif
@testable import CelMek

private struct SPK21Fixture {
    var daf: DAFFixture
    let size = 71
    let records: Int
    var epochOffset: Int { 4096 + size * records * 8 }
    var tailOffset: Int { epochOffset + (records + records/100) * 8 }
    init(big: Bool = false, records: Int = 2) {
        self.records = records; daf = DAFFixture(big: big)
        let words = size * records + records + records/100 + 2, last = 512 + words
        daf.bytes = daf.bytes.prefix(4096) + Data(repeating: 0, count: ((last+127)/128)*1024 - 4096)
        daf.summary(record: 3, previous: 0, next: 0, count: 1); daf.int(3,80); daf.int(Int32(last+1),84)
        daf.array(record: 3, slot: 0, target: 900021, type: 21, start: 513, end: Int32(last))
        daf.double(0, 2048+24); daf.double(Double(records*10),2048+32)
        for r in 0..<records {
            let base = 4096 + r * size * 8
            daf.double(Double(r*10),base)
            // Reference p/v are interleaved; constant acceleration has order 1.
            for (i,value) in [10.0+Double(r)*100,1,20,2,30,3].enumerated() { daf.double(value,base+(16+i)*8) }
            for axis in 0..<3 { daf.double(Double(axis+1),base+(22+axis*15)*8) }
            daf.double(3,base+67*8)
            for axis in 0..<3 { daf.double(1,base+(68+axis)*8) }
            daf.double(Double((r+1)*10),epochOffset+r*8)
        }
        for i in 0..<records/100 { daf.double(Double((i+1)*1000),epochOffset+(records+i)*8) }
        daf.double(15,tailOffset); daf.double(Double(records),tailOffset+8)
    }
}
private struct SPK21References: Decodable {
    struct Asset: Decodable { let name,sha256: String; let dimension,count: Int }
    struct State: Decodable { let file: String; let t: Double; let record: Int; let state: [Double] }
    let toolkit,wrapper,basis: String
    let files: [Asset]
    let states: [State]
}
@Suite struct SpiceSPKType21Tests {
    @Test func binaryOracleHighOrdersAndHundredRecordDirectories() throws {
        let root = try #require(Bundle.module.url(forResource: "reference", withExtension: "json", subdirectory: "Fixtures/spk21"))
        let refs = try JSONDecoder().decode(SPK21References.self, from: Data(contentsOf: root))
        #expect(refs.toolkit == "CSPICE_N0067" && refs.wrapper == "8.3.0")
        for asset in refs.files {
            let url = root.deletingLastPathComponent().appendingPathComponent(asset.name)
            #if canImport(CryptoKit)
            let hash = SHA256.hash(data: try Data(contentsOf: url)).map { String(format: "%02x", $0) }.joined()
            #expect(hash == asset.sha256)
            #endif
            let file = try SpiceDAFFile(url: url), segment = try SpiceSPKSegment(file: file, array: 0)
            #expect(segment.recordSizeWords == 4*asset.dimension+11 && segment.recordCount == asset.count)
            let kernels = try SpiceKernelSet(); try kernels.load(url)
            for row in refs.states where row.file == asset.name {
                let state = try segment.state(atSecondsTDB: row.t)
                #expect(state.recordIndex == row.record)
                let resolved = try kernels.snapshot().spkState(targetID: segment.descriptor.targetID, centreID: 0, atSecondsTDB: row.t)
                for axis in 0..<3 {
                    #expect(abs(state.positionKilometres[axis]-row.state[axis]) <= 6e-8)
                    #expect(abs(state.velocityKilometresPerTDBSecond[axis]-row.state[axis+3]) <= 2e-13)
                    #expect(resolved.positionKilometres[axis] == state.positionKilometres[axis])
                }
            }
        }
    }
    @Test func realHorizonsCeresAgainstBinaryOracle() throws {
        // Real references use the same numeric state shape, but omit the per-row file field.
        struct RealReference: Decodable {
            struct State: Decodable { let t: Double; let record: Int; let state: [Double] }
            let file,sha256: String
            let body,center,frame: Int32
            let states: [State]
        }
        let referenceURL = try #require(Bundle.module.url(forResource: "horizons-reference", withExtension: "json", subdirectory: "Fixtures/spk21"))
        let refs = try JSONDecoder().decode(RealReference.self,from: Data(contentsOf: referenceURL))
        let url = referenceURL.deletingLastPathComponent().appendingPathComponent(refs.file)
        #if canImport(CryptoKit)
        #expect(SHA256.hash(data: try Data(contentsOf: url)).map { String(format: "%02x",$0) }.joined() == refs.sha256)
        #endif
        let kernels = try SpiceKernelSet(); try kernels.load(url)
        let helper = SpiceSPKAdapterTests()
        let first = try #require(refs.states.first), last = try #require(refs.states.last)
        let coverage = try helper.epoch(first.t)...helper.epoch(last.t)
        let frame = Frame(name: "Ceres oracle heliocentric J2000", pseudoInertial: true)
        let provider = try SpiceSPKStateProvider(snapshot: kernels.snapshot(), body: .naif(refs.body),
            origin: .naif(refs.center), frame: frame, frameBindings: helper.bindings(frame, id: refs.frame),
            coverage: coverage, tdbModel: SpiceSPKAdapterTests.Clock(coverage: coverage), dataProvenance: helper.provenance)
        for row in refs.states {
            let state = try kernels.snapshot().spkState(targetID: refs.body,centreID: refs.center,frameID: refs.frame,atSecondsTDB: row.t)
            let si = try provider.checkedState(at: helper.epoch(row.t))
            #expect(si.acceleration == nil && si.frame === frame)
            #expect(state.contributions[0].state.recordIndex == row.record)
            for axis in 0..<3 {
                #expect(abs(state.positionKilometres[axis]-row.state[axis]) < 1e-7)
                #expect(abs(state.velocityKilometresPerTDBSecond[axis]-row.state[axis+3]) < 1e-12)
                #expect(abs(si.position[axis]-1000*row.state[axis]) < 1e-4)
                #expect(abs(si.velocity[axis]-1000*row.state[axis+3]) < 1e-9)
            }
        }
    }
    @Test func analyticConstantAccelerationBothByteOrdersAndEndpoints() throws {
        for big in [false,true] {
            let fixture = SPK21Fixture(big: big)
            try SpiceDAFTests().withFile(fixture.daf.bytes) { url in
                let s = try SpiceSPKSegment(file: SpiceDAFFile(url: url), array: 0)
                for t in [0.0,5,10,10.0.nextUp,15,20] {
                    let record = t <= 10 ? 0 : 1, dt = t - Double(record*10)
                    let state = try s.state(atSecondsTDB: t)
                    #expect(state.recordIndex == record)
                    for axis in 0..<3 {
                        let speed = Double(axis+1), initial = Double((axis+1)*10) + (axis == 0 ? Double(record*100) : 0)
                        #expect(abs(state.positionKilometres[axis] - (initial + speed*dt + speed*dt*dt/2)) < 2e-12)
                        #expect(abs(state.velocityKilometresPerTDBSecond[axis] - speed*(1+dt)) < 2e-13)
                    }
                }
                #expect(throws: SpiceSPKError.outOfCoverage) { try s.state(atSecondsTDB: -Double.leastNonzeroMagnitude) }
                #expect(throws: SpiceSPKError.outOfCoverage) { try s.state(atSecondsTDB: 20.0.nextUp) }
                #expect(throws: SpiceSPKError.invalidEpoch) { try s.state(atSecondsTDB: .nan) }
            }
        }
    }
    @Test func malformedLayoutEpochsOrdersStepsAndResourceLimits() throws {
        var cases = [SPK21Fixture]()
        for (offset,value) in [(0,14.0),(0,26),(0,15.5),(8,0),(8,1.5)] {
            var f = SPK21Fixture(); f.daf.double(value,f.tailOffset+offset); cases.append(f)
        }
        var f = SPK21Fixture(); f.daf.double(10,f.epochOffset+8); cases.append(f)
        f = SPK21Fixture(); f.daf.double(.nan,f.epochOffset); cases.append(f)
        f = SPK21Fixture(); f.daf.double(19,f.epochOffset+8); cases.append(f)
        f = SPK21Fixture(); f.daf.double(.infinity,4096); cases.append(f)
        f = SPK21Fixture(); f.daf.double(1.5,4096+68*8); cases.append(f)
        f = SPK21Fixture(); f.daf.double(16,4096+68*8); cases.append(f)
        f = SPK21Fixture(); f.daf.double(2,4096+68*8); cases.append(f) // required step is zero
        f = SPK21Fixture(); f.daf.double(1,4096+67*8); cases.append(f)
        f = SPK21Fixture(records: 100); f.daf.double(999,f.epochOffset+100*8); cases.append(f)
        for fixture in cases {
            try SpiceDAFTests().withFile(fixture.daf.bytes) { url in
                #expect(throws: SpiceSPKError.self) { try SpiceDAFFile(url: url).spkState(array: 0,atSecondsTDB: 5) }
            }
        }
        try SpiceDAFTests().withFile(SPK21Fixture().daf.bytes) { url in
            let file = try SpiceDAFFile(url: url, limits: SpiceDAFLimits(maximumReadWords: 70))
            #expect(throws: SpiceDAFError.self) { try SpiceSPKSegment(file: file,array: 0) }
        }
    }
}
