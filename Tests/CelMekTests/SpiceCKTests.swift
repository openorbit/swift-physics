// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

struct CKFixture {
    var daf: DAFFixture
    init(big: Bool = false, velocity: Bool = true, type: Int32 = 1, times: [Double] = [0,50,100]) {
        daf = DAFFixture(big: big)
        let width = velocity ? 7 : 4
        let size = times.count*(width+1)+(times.count-1)/100+1
        daf.bytes = Data(repeating: 0,count: max(8,(512+size+127)/128)*1024)
        // Recreate the documented file record and integrity marker before replacing its summary.
        let header = DAFFixture(big: big).bytes
        daf.bytes.replaceSubrange(0..<4096,with: header.prefix(4096))
        daf.text("DAF/CK",0,8); daf.int(3,80); daf.int(Int32(513+size),84)
        daf.summary(record: 3,previous: 0,next: 0,count: 1)
        daf.double(0,2072); daf.double(100,2080)
        for (i,v) in [Int32(-77001),1,type,velocity ? 1 : 0,513,Int32(512+size)].enumerated() { daf.int(v,2088+i*4) }
        for i in times.indices {
            let angle = Double(i)*Double.pi/4
            var values = [cos(angle/2),0,0,sin(angle/2)]
            if velocity { values += [1,0,0] }
            for (j,value) in values.enumerated() { daf.double(value,4096+(i*width+j)*8) }
            daf.double(times[i],4096+(times.count*width+i)*8)
        }
        for i in 1...max(1,(times.count-1)/100) where i <= (times.count-1)/100 {
            daf.double((times[i*100-1]+times[i*100])/2,4096+(times.count*(width+1)+i-1)*8)
        }
        daf.double(Double(times.count),4096+(size-1)*8)
    }
}
@Suite struct SpiceCKTests {
    let files = SpiceDAFTests()
    @Test func discreteSamplesBothEndianAndMissingRates() throws {
        for big in [false,true] { for velocity in [false,true] {
            let f = CKFixture(big: big,velocity: velocity)
            try files.withFile(f.daf.bytes) { url in
                let file = try SpiceDAFFile(url: url), s = try SpiceCKSegment(file: file,array: 0)
                #expect(file.metadata.arrays[0].ck?.hasAngularVelocity == velocity)
                #expect(try s.pointing(atEncodedTicks: 25) == nil)
                let p = try #require(try s.pointing(atEncodedTicks: 25,toleranceTicks: 25))
                #expect(p.encodedTicks == 50 && p.requestedEncodedTicks == 25 && p.recordIndex == 1)
                #expect((p.angularVelocityReferenceRadiansPerSecond != nil) == velocity)
                let v = p.referenceToInstrument.matrix*SIMD3<Double>(1,0,0)
                #expect(abs(v.x-sqrt(0.5)) < 1e-12 && abs(v.y-sqrt(0.5)) < 1e-12)
                #expect(try s.pointing(atEncodedTicks: 100)?.recordIndex == 2)
                #expect(throws: SpiceCKError.invalidQuery) { try s.pointing(atEncodedTicks: .nan) }
            }
        } }
    }
    @Test func directoriesMalformedDataAndUnsupportedTypes() throws {
        var f = CKFixture(times: (0..<101).map(Double.init))
        try files.withFile(f.daf.bytes) { url in
            let segment = try SpiceCKSegment(file: SpiceDAFFile(url: url),array: 0)
            #expect(try segment.pointing(atEncodedTicks: 100)?.recordIndex == 100)
        }
        f.daf.double(99,4096+101*8*8)
        try files.withFile(f.daf.bytes) { url in
            #expect(throws: SpiceCKError.malformedSegment) { try SpiceCKSegment(file: SpiceDAFFile(url: url),array: 0) }
        }
        for type in [Int32(4),5] {
            try files.withFile(CKFixture(type: type).daf.bytes) { url in
                #expect(throws: SpiceCKError.unsupportedType(type)) { try SpiceCKSegment(file: SpiceDAFFile(url: url),array: 0) }
            }
        }
        var bad = CKFixture(); bad.daf.double(2,4096)
        try files.withFile(bad.daf.bytes) { url in
            #expect(throws: SpiceCKError.malformedSegment) { try SpiceCKSegment(file: SpiceDAFFile(url: url),array: 0).pointing(atEncodedTicks: 0) }
        }
    }
    @Test func priorityGapsAndSelectedFailures() throws {
        try files.withFile(CKFixture().daf.bytes) { oldURL in
            try files.withFile(CKFixture(times: [0,100]).daf.bytes) { newURL in
                let set = try SpiceKernelSet(); try set.load(oldURL); let newest = try set.load(newURL)
                let retained = set.snapshot()
                #expect(try retained.ckPointing(instrumentID: -77001,atEncodedTicks: 50).sourcePath == oldURL.path)
                #expect(try retained.ckPointing(instrumentID: -77001,atEncodedTicks: 100).sourcePath == newURL.path)
                #expect(throws: SpiceCKError.resourceLimit) {
                    try retained.ckPointing(instrumentID: -77001,atEncodedTicks: 50,maximumDescriptorChecks: 1)
                }
                try set.unload(newest)
                #expect(try retained.ckPointing(instrumentID: -77001,atEncodedTicks: 100).sourcePath == newURL.path)
                #expect(try set.snapshot().ckPointing(instrumentID: -77001,atEncodedTicks: 100).sourcePath == oldURL.path)
                try files.withFile(CKFixture(type: 4).daf.bytes) { unsupported in
                    try set.load(unsupported)
                    #expect(throws: SpiceCKError.unsupportedType(4)) { try set.snapshot().ckPointing(instrumentID: -77001,atEncodedTicks: 50) }
                }
            }
        }
    }
    @Test func clockAdapterPreservesUnavailableAccelerationAndSnapshot() throws {
        for velocity in [false,true] {
            try files.withFile(CKFixture(velocity: velocity).daf.bytes) { url in
                let set = try SpiceKernelSet(), id = try set.load(url), snapshot = set.snapshot()
                let helper = SpiceSCLKTests()
                let clock = try helper.clock(helper.text + "\nSCLK01_TIME_SYSTEM_77=2\n")
                let ref = CelestialFrames.icrf, body = Frame(parent: ref,name: "CK body",pseudoInertial: false)
                let provider = try SpiceCKAttitudeProvider(snapshot: snapshot,instrumentID: -77001,clock: clock,
                    referenceFrame: ref,bodyFrame: body,frameBindings: SpiceFrameBindings([.icrfJ2000(provenance: helper.provenance)]),angularRateTimeScale: .atomic)
                try set.unload(id)
                if velocity {
                    let state = try provider.attitude(at: helper.epoch(100))
                    #expect(state.angularAccelerationBodyRadiansPerSecondSquared == nil)
                    #expect(abs(state.angularVelocityBodyRadiansPerSecond.x-sqrt(0.5)) < 1e-12)
                    #expect(abs(state.angularVelocityBodyRadiansPerSecond.y-sqrt(0.5)) < 1e-12)
                    #expect(throws: AttitudeError.angularAccelerationUnavailable) { try state.orientation() }
                } else {
                    #expect(throws: SpiceCKError.angularVelocityUnavailable) { try provider.attitude(at: helper.epoch(100)) }
                }
                #expect(throws: SpiceCKError.coverageGap) { try provider.attitude(at: helper.epoch(50)) }
                #expect(throws: SpiceCKError.missingInstrument) { try snapshot.ckPointing(instrumentID: -1,atEncodedTicks: 0) }
            }
        }
    }
}
