// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

struct CK3Fixture {
    var daf: DAFFixture
    init(big: Bool = false, velocity: Bool = true, times: [Double] = [0,20,60,80,100], starts: [Double] = [0,60,100], angles: [Double]? = nil) {
        daf = DAFFixture(big: big)
        let n = times.count, m = starts.count, width = velocity ? 7 : 4
        let length = n*(width+1)+(n-1)/100+m+(m-1)/100+2
        daf.bytes.append(Data(repeating: 0,count: max(0,((512+length+127)/128)*1024-daf.bytes.count)))
        daf.text("DAF/CK",0,8); daf.int(3,80); daf.int(Int32(513+length),84)
        daf.summary(record: 3,previous: 0,next: 0,count: 1)
        daf.double(times.first!,2072); daf.double(times.last!,2080)
        for (i,v) in [Int32(-77001),1,3,velocity ? 1 : 0,513,Int32(512+length)].enumerated() { daf.int(v,2088+i*4) }
        for i in 0..<n {
            let a = angles?[i] ?? Double(i)*0.4
            var values = [cos(a/2),0,0,sin(a/2)]
            if velocity { values += [Double(i),0,9] }
            for (j,v) in values.enumerated() { daf.double(v,4096+(width*i+j)*8) }
            daf.double(times[i],4096+(n*width+i)*8)
        }
        for i in 0..<((n-1)/100) { daf.double(times[(i+1)*100-1],4096+(n*(width+1)+i)*8) }
        let offset = n*(width+1)+(n-1)/100
        for (i,v) in starts.enumerated() { daf.double(v,4096+(offset+i)*8) }
        for i in 0..<((m-1)/100) { daf.double(starts[(i+1)*100-1],4096+(offset+m+i)*8) }
        daf.double(Double(m),4096+(length-2)*8); daf.double(Double(n),4096+(length-1)*8)
    }
}
@Suite struct SpiceCKType3Tests {
    let files = SpiceDAFTests()
    @Test func independentPointingAndRatesBothEndianAndGapEndpoints() throws {
        for big in [false,true] { for velocity in [false,true] {
            try files.withFile(CK3Fixture(big: big,velocity: velocity).daf.bytes) { url in
                let s = try SpiceCKType3Segment(file: SpiceDAFFile(url: url),array: 0)
                let p = try #require(try s.pointing(atEncodedTicks: 10))
                let x = p.referenceToInstrument.matrix*SIMD3<Double>(1,0,0)
                #expect(abs(x.x-cos(0.2)) < 1e-12 && abs(x.y-sin(0.2)) < 1e-12)
                #expect(p.encodedTicks == 10 && p.secondsPerEncodedTick == nil)
                if velocity { #expect(p.angularVelocityReferenceRadiansPerSecond == SIMD3(0.5,0,9)) }
                else { #expect(p.angularVelocityReferenceRadiansPerSecond == nil) }
                #expect(try s.pointing(atEncodedTicks: 40) == nil)
                let endpoint = try #require(try s.pointing(atEncodedTicks: 40,toleranceTicks: 20))
                #expect(endpoint.encodedTicks == 60 && endpoint.requestedEncodedTicks == 40)
                #expect(try s.pointing(atEncodedTicks: 90) == nil)
                #expect(try s.pointing(atEncodedTicks: 100)?.recordIndex == 4)
                #expect(throws: SpiceCKError.invalidQuery) { try s.pointing(atEncodedTicks: 10,toleranceTicks: -1) }
            }
        } }
    }
    @Test func shortestArcAntipodalAndTinyRotations() throws {
        for (angles,expected) in [([0.0,2*Double.pi-0.4],-0.2),([0,2*Double.pi],0),([0,1e-12],5e-13)] {
            try files.withFile(CK3Fixture(times: [0,100],starts: [0],angles: angles).daf.bytes) { url in
                let s = try SpiceCKType3Segment(file: SpiceDAFFile(url: url),array: 0)
                let p = try #require(try s.pointing(atEncodedTicks: 50))
                let x = p.referenceToInstrument.matrix*SIMD3<Double>(1,0,0)
                #expect(abs(x.x-cos(expected)) < 1e-14 && abs(x.y-sin(expected)) < 1e-14)
            }
        }
    }
    @Test func interpolationWithNoncommutingInitialAttitude() throws {
        var f = CK3Fixture(times: [0,100],starts: [0])
        let initial = Quaternion(angle: 0.6,axis: SIMD3<Double>(1,0,0))
        let end = initial*Quaternion(angle: 0.8,axis: SIMD3<Double>(0,0,1))
        for (i,q) in [initial,end].enumerated() {
            for (j,v) in [q.w,q.x,q.y,q.z].enumerated() { f.daf.double(v,4096+(7*i+j)*8) }
        }
        try files.withFile(f.daf.bytes) { url in
            let s = try SpiceCKType3Segment(file: SpiceDAFFile(url: url),array: 0)
            let p = try #require(try s.pointing(atEncodedTicks: 25))
            let expected = (initial*Quaternion(angle: 0.2,axis: SIMD3<Double>(0,0,1))).matrix*SIMD3<Double>(1,0,0)
            let actual = p.referenceToInstrument.matrix*SIMD3<Double>(1,0,0)
            #expect((0..<3).allSatisfy { abs(actual[$0]-expected[$0]) < 1e-12 })
        }
    }
    @Test func directoriesAndMalformedIntervalMembership() throws {
        let times = (0..<101).map(Double.init)
        var f = CK3Fixture(times: times,starts: times)
        try files.withFile(f.daf.bytes) { url in
            let s = try SpiceCKType3Segment(file: SpiceDAFFile(url: url),array: 0)
            #expect(try s.pointing(atEncodedTicks: 100)?.recordIndex == 100)
            #expect(try s.pointing(atEncodedTicks: 50.5) == nil)
            #expect(throws: SpiceCKError.resourceLimit) { try SpiceCKType3Segment(file: SpiceDAFFile(url: url),array: 0,maximumRecords: 100) }
        }
        f.daf.double(98,4096+808*8)
        try files.withFile(f.daf.bytes) { url in
            #expect(throws: SpiceCKError.malformedSegment) { try SpiceCKType3Segment(file: SpiceDAFFile(url: url),array: 0) }
        }
        for f in [CK3Fixture(starts: [0,61]),CK3Fixture(starts: [20,60]),CK3Fixture(times: [0,20,20,80,100])] {
            try files.withFile(f.daf.bytes) { url in
                #expect(throws: SpiceCKError.malformedSegment) { try SpiceCKType3Segment(file: SpiceDAFFile(url: url),array: 0) }
            }
        }
    }
    @Test func adapterRatesAreExplicitAndAccelerationUnavailable() throws {
        for velocity in [false,true] {
            try files.withFile(CK3Fixture(velocity: velocity).daf.bytes) { url in
                let set = try SpiceKernelSet(); try set.load(url)
                let helper = SpiceSCLKTests(), clock = try helper.clock(helper.text+"\nSCLK01_TIME_SYSTEM_77=2\n")
                let ref = CelestialFrames.icrf, body = Frame(parent: ref,name: "CK3 body",pseudoInertial: false)
                let p = try SpiceCKAttitudeProvider(snapshot: set.snapshot(),instrumentID: -77001,clock: clock,
                    referenceFrame: ref,bodyFrame: body,frameBindings: SpiceFrameBindings([.icrfJ2000(provenance: helper.provenance)]),angularRateTimeScale: .atomic)
                if velocity {
                    let state = try p.attitude(at: helper.epoch(20))
                    #expect(state.angularVelocityBodyRadiansPerSecond.z == 9)
                    #expect(state.angularAccelerationBodyRadiansPerSecondSquared == nil)
                    #expect(throws: AttitudeError.angularAccelerationUnavailable) { try state.orientation() }
                } else {
                    #expect(throws: SpiceCKError.angularVelocityUnavailable) { try p.attitude(at: helper.epoch(20)) }
                }
                #expect(throws: SpiceCKError.coverageGap) { try p.attitude(at: helper.epoch(80)) }
            }
        }
    }
}
