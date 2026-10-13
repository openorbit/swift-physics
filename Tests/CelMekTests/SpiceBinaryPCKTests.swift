// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

struct PCKFixture {
    var spk: SPKFixture
    var bytes: Data { spk.daf.bytes }
    init(big: Bool = false, type: Int32 = 2, classID: Int32 = 301, frameID: Int32 = 1,
         coefficients: [[Double]]? = nil, initial: Double = -20, interval: Double = 10, records: Int = 2) {
        spk = SPKFixture(big: big, initial: initial, interval: interval, coefficients: coefficients, records: records)
        spk.daf.text("DAF/PCK",0,8); spk.daf.int(5,12)
        let end = Int32(512 + spk.size*records + 4)
        for (i,v) in [classID,frameID,type,Int32(513),end,Int32(123456)].enumerated() { spk.daf.int(v,2048+40+i*4) }
    }
}
@Suite struct SpiceBinaryPCKTests {
    let files = SpiceDAFTests()
    var provenance: PhysicalDataProvenance { try! PhysicalDataProvenance(source: "Independent public-format PCK",version: "fixture-1") }
    func epoch(_ t: Double) throws -> Epoch { Epoch(offsetFromJ2000: try CelMek.Duration(seconds: t)) }
    func near(_ a: SIMD3<Double>, _ b: SIMD3<Double>, tolerance: Double = 1e-11) -> Bool {
        (0..<3).allSatisfy { abs(a[$0]-b[$0]) < tolerance }
    }
    @Test func independentAngleDerivativesBothEndianAndBoundaries() throws {
        for big in [false,true] {
            let f = PCKFixture(big: big)
            try files.withFile(f.bytes) { url in
                let file = try SpiceDAFFile(url: url), segment = try SpiceBinaryPCKSegment(file: file,array: 0)
                #expect(file.metadata.integerComponentCount == 5)
                #expect(file.metadata.arrays[0].spk == nil)
                #expect(segment.descriptor.frameClassID == 301 && segment.descriptor.frameID == 1)
                for t in [-20.0,-18.75,-15,-12.5,-10,-7.5,-5,-1.25,0] {
                    let record = t < -10 ? 0 : 1, mid = record == 0 ? -15.0 : -5.0, x = (t-mid)/5
                    let sample = try segment.angles(atSecondsTDB: t)
                    let reference = SpiceSPKTests().reference(x,record: record,type: 2,radius: 5)
                    let second = SIMD3<Double>((-8+12*x+24*x*x)/25,0,(24-96*x*x)/25)
                    #expect(near(sample.eulerAnglesRadians,reference.0))
                    #expect(near(sample.ratesRadiansPerTDBSecond,reference.1))
                    #expect(near(sample.accelerationsRadiansPerTDBSecondSquared,second))
                    #expect(sample.recordIndex == record && sample.sourcePath == url.path)
                }
                #expect(try segment.angles(atSecondsTDB: (-10.0).nextDown).recordIndex == 0)
                #expect(try segment.angles(atSecondsTDB: (-10.0).nextUp).recordIndex == 1)
                #expect(throws: SpiceBinaryPCKError.outOfCoverage) { try segment.angles(atSecondsTDB: 0.0.nextUp) }
                #expect(throws: SpiceBinaryPCKError.invalidEpoch) { try segment.angles(atSecondsTDB: .nan) }
            }
        }
    }
    @Test func degreeZeroAndHighDegreeAcrossStorageRecords() throws {
        for degree in [0,60] {
            var coefficients = Array(repeating: 0.0,count: degree+1); coefficients[degree] = 1
            let f = PCKFixture(coefficients: [coefficients,coefficients,coefficients])
            try files.withFile(f.bytes) { url in
                let segment = try SpiceBinaryPCKSegment(file: SpiceDAFFile(url: url),array: 0)
                let x = 0.4, theta = acos(x), n = Double(degree)
                let sample = try segment.angles(atSecondsTDB: -15+5*x)
                #expect(abs(sample.eulerAnglesRadians.x - cos(n*theta)) < 1e-12)
                let first = degree == 0 ? 0 : n*sin(n*theta)/sqrt(1-x*x)/5
                let second = degree == 0 ? 0 : (x*(first*5)-n*n*cos(n*theta))/(1-x*x)/25
                #expect(abs(sample.ratesRadiansPerTDBSecond.x-first) < 1e-10)
                #expect(abs(sample.accelerationsRadiansPerTDBSecondSquared.x-second) < 1e-8)
            }
        }
    }
    @Test func descriptorsDirectoryAndRecordFailuresAreExplicit() throws {
        for type: Int32 in [3,20,99] {
            try files.withFile(PCKFixture(type: type).bytes) { url in
                #expect(throws: SpiceBinaryPCKError.unsupportedType(type)) { try SpiceDAFFile(url: url).pckAngles(array: 0,atSecondsTDB: -15) }
            }
        }
        for field in [8,12,2048+40,2048+44,2048+48] {
            var f = PCKFixture(); f.spk.daf.int(0,field)
            try files.withFile(f.bytes) { url in #expect(throws: SpiceDAFError.self) { try SpiceDAFFile(url: url) } }
        }
        let f = PCKFixture()
        for (offset,value) in [(f.spk.directoryOffset+8,0.0),(f.spk.directoryOffset+16,4.5),(f.spk.directoryOffset+24,3.0),
                               (4096,123.0),(4104,-1.0),(4112,Double.nan)] {
            var bad = f; bad.spk.daf.double(value,offset)
            try files.withFile(bad.bytes) { url in
                #expect(throws: (any Error).self) { try SpiceDAFFile(url: url).pckAngles(array: 0,atSecondsTDB: -15) }
            }
        }
        try files.withFile(f.bytes) { url in
            #expect(throws: SpiceDAFError.self) { try SpiceBinaryPCKSegment(file: SpiceDAFFile(url: url,limits: SpiceDAFLimits(maximumReadWords: 10)),array: 0) }
        }
    }
    @Test func priorityGapsUnloadAndMixedSPKPool() throws {
        let directory = try SpiceKernelSetTests().directory(); defer { try? FileManager.default.removeItem(at: directory) }
        let a = directory.appendingPathComponent("a.bpc"), b = directory.appendingPathComponent("b.bpc"), spk = directory.appendingPathComponent("state.bsp")
        try PCKFixture().bytes.write(to: a); try PCKFixture(type: 20).bytes.write(to: b); try SPKFixture().daf.bytes.write(to: spk)
        let set = try SpiceKernelSet(); try set.load(a); try set.load(spk)
        let old = set.snapshot(); #expect(try old.pckAngles(frameClassID: 301,atSecondsTDB: -15).sourcePath == a.path)
        #expect(try old.spkState(targetID: 399,atSecondsTDB: -15).sourcePath == spk.path)
        let id = try set.load(b)
        #expect(throws: SpiceBinaryPCKError.unsupportedType(20)) { try set.snapshot().pckAngles(frameClassID: 301,atSecondsTDB: -15) }
        try set.unload(id)
        #expect(try set.snapshot().pckAngles(frameClassID: 301,atSecondsTDB: -15).sourcePath == a.path)
        #expect(old.cacheScope != set.snapshot().cacheScope)
        #expect(throws: SpiceBinaryPCKError.coverageGap(301)) { try old.pckAngles(frameClassID: 301,atSecondsTDB: 1) }
        #expect(throws: SpiceBinaryPCKError.missingFrameClass(399)) { try old.pckAngles(frameClassID: 399,atSecondsTDB: -15) }
        #expect(throws: SpiceBinaryPCKError.descriptorLimit) { try old.pckAngles(frameClassID: 301,atSecondsTDB: -15,maximumDescriptorChecks: 0) }
    }
    @Test func providerEulerConventionAndNonUnitAtomicClock() throws {
        let f = PCKFixture(coefficients: [[0,0,0],[0,0,0],[0.5,0.2,0.01]])
        try files.withFile(f.bytes) { url in
            let set = try SpiceKernelSet(), id = try set.load(url), coverage = try epoch(-18)...epoch(-12)
            let clock = SpiceSPKAdapterTests.Clock(coverage: coverage,offset: 0.1,drift: 0.01,quadratic: 0.0001)
            let provider = try SpiceBinaryPCKOrientationProvider(snapshot: set.snapshot(),frameClassID: 301,referenceFrame: CelestialFrames.icrf,
                frameBindings: SpiceFrameBindings([.icrfJ2000(provenance: provenance)]),coverage: coverage,tdbModel: clock,provenance: provenance)
            let t = -15.0, sample = try provider.evaluation(at: epoch(t))
            let et = t+0.1+0.01*t+0.0001*t*t, x = (et+15)/5
            let w = 0.49+0.2*x+0.02*x*x, first = (0.2+0.04*x)/5, second = 0.04/25
            #expect(near(sample.orientation.rotation.matrix * SIMD3<Double>(1,0,0),SIMD3(cos(w),-sin(w),0)))
            let rate = 1.01+0.0002*t
            #expect(near(sample.orientation.angularVelocity,SIMD3(0,0,-first*rate)))
            #expect(near(sample.orientation.angularAcceleration,SIMD3(0,0,-second*rate*rate-first*0.0002)))
            #expect(sample.nativeAngles.descriptor.frameClassID == 301)
            try set.unload(id)
            #expect(try provider.evaluation(at: epoch(t)).nativeAngles.eulerAnglesRadians == sample.nativeAngles.eulerAnglesRadians)
            #expect(throws: SpiceBinaryPCKError.outOfCoverage) { try provider.orientation(at: epoch(-11)) }
        }
    }
    @Test func coupledEulerRatesAndFrameComposition() throws {
        let f = PCKFixture(frameID: 2)
        try files.withFile(f.bytes) { url in
            let set = try SpiceKernelSet(); try set.load(url)
            var tree: FrameTree? = FrameTree()
            let root = tree!.root
            let q = Quaternion(angle: 0.3,axis: SIMD3(0,0,1))
            let native = try tree!.makeFrame(parent: root,name: "alternate inertial axes",pseudoInertial: true,
                provider: FixedTransformProvider(rotation: q,translation: SIMD3(123,456,789)))
            let bindings = try SpiceFrameBindings([SpiceFrameBinding(frameID: 2,axesFrame: native,provenance: provenance)])
            let coverage = try epoch(-18)...epoch(-12), clock = SpiceSPKAdapterTests.Clock(coverage: coverage)
            let provider = try SpiceBinaryPCKOrientationProvider(snapshot: set.snapshot(),frameClassID: 301,referenceFrame: root,
                frameBindings: bindings,frameTree: tree,coverage: coverage,tdbModel: clock,provenance: provenance)
            let t = -15.0, h = 0.00001, input = SIMD3<Double>(0.3,0.5,0.7)
            let sample = try provider.evaluation(at: epoch(t)), angles = sample.nativeAngles.eulerAnglesRadians
            func z(_ a: Double, _ v: SIMD3<Double>) -> SIMD3<Double> { SIMD3(cos(a)*v.x-sin(a)*v.y,sin(a)*v.x+cos(a)*v.y,v.z) }
            func x(_ a: Double, _ v: SIMD3<Double>) -> SIMD3<Double> { SIMD3(v.x,cos(a)*v.y-sin(a)*v.z,sin(a)*v.y+cos(a)*v.z) }
            let expected = z(-angles.z,x(-angles.y,z(-angles.x,z(0.3,input))))
            let r = sample.orientation.rotation.matrix * input
            #expect(near(r,expected))
            let plus = try provider.orientation(at: epoch(t+h)), minus = try provider.orientation(at: epoch(t-h))
            let numerical = (plus.rotation.matrix*input-minus.rotation.matrix*input)/(2*h)
            let omega = sample.orientation.angularVelocity
            #expect(near(numerical,SIMD3(omega.y*r.z-omega.z*r.y,omega.z*r.x-omega.x*r.z,omega.x*r.y-omega.y*r.x),tolerance: 1e-8))
            #expect(near((plus.angularVelocity-minus.angularVelocity)/(2*h),sample.orientation.angularAcceleration,tolerance: 1e-8))
            tree = nil
            #expect(throws: SpiceBinaryPCKError.frameTreeUnavailable) { try provider.orientation(at: epoch(t)) }
        }
    }
    @Test func laterArrayPriorityAndBoundedDescriptorSearch() throws {
        var f = PCKFixture()
        let recordWords = f.spk.size*2+4, secondStart = 513+recordWords, secondEnd = secondStart+recordWords-1
        let payload = f.spk.daf.bytes.subdata(in: 4096..<(4096+recordWords*8))
        f.spk.daf.bytes.replaceSubrange((secondStart-1)*8..<secondEnd*8,with: payload)
        f.spk.daf.int(Int32(secondEnd+1),84)
        f.spk.daf.summary(record: 3,previous: 0,next: 0,count: 2)
        let descriptor = 2048+24+40
        f.spk.daf.double(-20,descriptor); f.spk.daf.double(0,descriptor+8)
        for (i,v) in [Int32(301),Int32(1),Int32(2),Int32(secondStart),Int32(secondEnd)].enumerated() {
            f.spk.daf.int(v,descriptor+16+i*4)
        }
        f.spk.daf.double(777,(secondStart-1)*8+16)
        try files.withFile(f.bytes) { url in
            let set = try SpiceKernelSet(); try set.load(url)
            let snapshot = set.snapshot()
            #expect(try snapshot.pckAngles(frameClassID: 301,atSecondsTDB: -15,maximumDescriptorChecks: 1).arrayIndex == 1)
            #expect(try snapshot.pckAngles(frameClassID: 301,atSecondsTDB: -15).eulerAnglesRadians.x == 778.25)
            #expect(throws: SpiceBinaryPCKError.descriptorLimit) { try snapshot.pckAngles(frameClassID: 399,atSecondsTDB: -15,maximumDescriptorChecks: 1) }
        }
    }

}
