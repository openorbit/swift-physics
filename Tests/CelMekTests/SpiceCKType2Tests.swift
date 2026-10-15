// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

struct CK2Fixture {
    var daf: DAFFixture
    init(big: Bool = false, starts: [Double] = [0,50], stops: [Double] = [50,100], rate: Double = 3) {
        daf = DAFFixture(big: big)
        let n = starts.count, length = 10*n+(n-1)/100
        daf.bytes.append(Data(repeating: 0,count: max(0,((512+length+127)/128)*1024-daf.bytes.count)))
        daf.text("DAF/CK",0,8); daf.int(3,80); daf.int(Int32(513+length),84)
        daf.summary(record: 3,previous: 0,next: 0,count: 1)
        daf.double(starts[0],2072); daf.double(stops[n-1],2080)
        for (i,v) in [Int32(-77001),1,2,1,513,Int32(512+length)].enumerated() { daf.int(v,2088+i*4) }
        for i in 0..<n {
            for (j,v) in [1.0,0,0,0,0,0,0.01,rate].enumerated() { daf.double(v,4096+(8*i+j)*8) }
            daf.double(starts[i],4096+(8*n+i)*8); daf.double(stops[i],4096+(9*n+i)*8)
        }
        for i in 0..<((n-1)/100) {
            let j = (i+1)*100
            daf.double((stops[j-1]+starts[j])/2,4096+(10*n+i)*8)
        }
    }
}
@Suite struct SpiceCKType2Tests {
    let files = SpiceDAFTests()
    @Test func constantSpinDirectionSharedEndpointsAndGaps() throws {
        for big in [false,true] {
            try files.withFile(CK2Fixture(big: big).daf.bytes) { url in
                let segment = try SpiceCKType2Segment(file: SpiceDAFFile(url: url),array: 0)
                let sample = try #require(try segment.pointing(atEncodedTicks: 25))
                let x = sample.referenceToInstrument.matrix*SIMD3<Double>(1,0,0)
                #expect(abs(x.x-cos(0.75)) < 1e-12 && abs(x.y+sin(0.75)) < 1e-12)
                #expect(sample.secondsPerEncodedTick == 3 && sample.encodedTicks == 25)
                #expect(try segment.pointing(atEncodedTicks: 50)?.recordIndex == 1)
                #expect(try segment.pointing(atEncodedTicks: 100)?.recordIndex == 1)
                #expect(try segment.pointing(atEncodedTicks: 100.0.nextUp) == nil)
            }
        }
        try files.withFile(CK2Fixture(starts: [0,60],stops: [40,100]).daf.bytes) { url in
            let set = try SpiceKernelSet(); try set.load(url)
            #expect(throws: SpiceCKError.coverageGap) { try set.snapshot().ckPointing(instrumentID: -77001,atEncodedTicks: 50,toleranceTicks: 20) }
        }
    }
    @Test func directoryBoundsAndMalformedIntervals() throws {
        let starts = (0..<101).map { Double($0*2) }, stops = starts.map { $0+1 }
        var f = CK2Fixture(starts: starts,stops: stops)
        try files.withFile(f.daf.bytes) { url in
            let segment = try SpiceCKType2Segment(file: SpiceDAFFile(url: url),array: 0)
            #expect(try segment.pointing(atEncodedTicks: 200.5)?.recordIndex == 100)
            #expect(throws: SpiceCKError.resourceLimit) { try SpiceCKType2Segment(file: SpiceDAFFile(url: url),array: 0,maximumRecords: 100) }
        }
        f.daf.double(1,4096+1010*8)
        try files.withFile(f.daf.bytes) { url in
            #expect(throws: SpiceCKError.malformedSegment) { try SpiceCKType2Segment(file: SpiceDAFFile(url: url),array: 0) }
        }
        for fixture in [CK2Fixture(starts: [0,30],stops: [40,100]), CK2Fixture(starts: [0,50],stops: [0,100]), CK2Fixture(rate: 0)] {
            try files.withFile(fixture.daf.bytes) { url in
                #expect(throws: SpiceCKError.malformedSegment) { try SpiceCKType2Segment(file: SpiceDAFFile(url: url),array: 0).pointing(atEncodedTicks: 0) }
            }
        }
    }
    @Test func atomicDerivativesFollowClockAndRecordRate() throws {
        for tdb in [false,true] {
            try files.withFile(CK2Fixture(starts: [0],stops: [100]).daf.bytes) { url in
                let set = try SpiceKernelSet(); try set.load(url)
                let helper = SpiceSCLKTests()
                let clock = try helper.clock(tdb ? helper.text : helper.text+"\nSCLK01_TIME_SYSTEM_77=2\n")
                let model = SpiceSPKAdapterTests.Clock(coverage: try helper.epoch(0)...helper.epoch(150),drift: 0.01,quadratic: 0.0001)
                let ref = CelestialFrames.icrf, body = Frame(parent: ref,name: "CK2 body",pseudoInertial: false)
                let provider = try SpiceCKAttitudeProvider(snapshot: set.snapshot(),instrumentID: -77001,clock: clock,
                    referenceFrame: ref,bodyFrame: body,frameBindings: SpiceFrameBindings([.icrfJ2000(provenance: helper.provenance)]),
                    angularRateTimeScale: .atomic,tdbModel: model)
                let t = 40.0, h = 0.001, state = try provider.attitude(at: helper.epoch(t))
                let before = try provider.attitude(at: helper.epoch(t-h)), after = try provider.attitude(at: helper.epoch(t+h))
                let expectedRate = 0.015*(tdb ? 1.01+0.0002*t : 1)
                #expect(abs(state.angularVelocityBodyRadiansPerSecond.z-expectedRate) < 1e-12)
                let alpha = try #require(state.angularAccelerationBodyRadiansPerSecondSquared)
                #expect(abs(alpha.z-(tdb ? 0.000003 : 0)) < 1e-12)
                let derivative = (after.angularVelocityBodyRadiansPerSecond-before.angularVelocityBodyRadiansPerSecond)/(2*h)
                #expect(abs(derivative.z-alpha.z) < 1e-10)
                let x = SIMD3<Double>(1,0,0)
                let dx = (after.bodyToReference.matrix*x-before.bodyToReference.matrix*x)/(2*h)
                let expected = state.bodyToReference.matrix*SIMD3<Double>(0,expectedRate,0)
                #expect(abs(dx.x-expected.x) < 1e-9 && abs(dx.y-expected.y) < 1e-9)
                let orientation = try state.orientation()
                #expect(abs(orientation.angularVelocity.z+expectedRate) < 1e-12)
            }
        }
    }
}
