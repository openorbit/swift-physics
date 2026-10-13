// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
@testable import CelMek

struct SPKFixture {
    var daf: DAFFixture
    let size: Int
    let directoryOffset: Int
    init(type: Int32 = 2, big: Bool = false, initial: Double = -20, interval: Double = 10,
         coefficients: [[Double]]? = nil, records: Int = 2, target: Int32 = 399) {
        daf = DAFFixture(big: big)
        let position = [[3.0,2,-1,0.5,0.25],[-7,0,0,0,0],[1,0,2,0,-1]]
        let velocity = [[11.0,3,0,0,0],[4,0,-1,0,0],[0,0,0,1,0]]
        let all = coefficients ?? (type == 3 ? position + velocity : position)
        size = 2 + all.count*all[0].count
        directoryOffset = 4096 + size*records*8
        let last = 512 + size*records + 4
        if daf.bytes.count < last*8 { daf.bytes.append(Data(repeating: 0,count: ((last+127)/128)*1024-daf.bytes.count)) }
        daf.summary(record: 3,previous: 0,next: 0,count: 1); daf.int(3,80); daf.int(Int32(last+1),84)
        daf.array(record: 3,slot: 0,target: target,type: type,start: 513,end: Int32(last))
        daf.double(initial,2048+24); daf.double(initial+interval*Double(records),2048+32)
        for r in 0..<records {
            let offset = 4096+r*size*8
            daf.double(initial+interval*(Double(r)+0.5),offset); daf.double(interval/2,offset+8)
            for (axis,cs) in all.enumerated() {
                for (k,c) in cs.enumerated() {
                    daf.double(c + (axis == 0 && k == 0 ? Double(r)*100 : 0),offset+16+(axis*cs.count+k)*8)
                }
            }
        }
        for (i,v) in [initial,interval,Double(size),Double(records)].enumerated() { daf.double(v,directoryOffset+i*8) }
    }
}
@Suite struct SpiceSPKTests {
    let files = SpiceDAFTests()
    func close(_ a: SIMD3<Double>, _ b: SIMD3<Double>, tolerance: Double = 1e-12) -> Bool {
        (0..<3).allSatisfy { abs(a[$0]-b[$0]) <= tolerance }
    }
    // Independent expansion in ordinary powers, rather than the evaluation recurrence.
    func reference(_ x: Double, record: Int, type: Int32, radius: Double) -> (SIMD3<Double>,SIMD3<Double>) {
        let p = SIMD3(4.25+0.5*x-4*x*x+2*x*x*x+2*pow(x,4)+Double(record)*100,
                      -7,-2+12*x*x-8*pow(x,4))
        let v = type == 2 ? SIMD3((0.5-8*x+6*x*x+8*x*x*x)/radius,0,(24*x-32*x*x*x)/radius)
            : SIMD3(11+3*x,5-2*x*x,4*x*x*x-3*x)
        return (p,v)
    }
    @Test func type2And3IndependentPolynomialStatesInBothByteOrders() throws {
        for type: Int32 in [2,3] { for big in [false,true] {
            let fixture = SPKFixture(type: type,big: big)
            try files.withFile(fixture.daf.bytes) { url in
                let file = try SpiceDAFFile(url: url), segment = try SpiceSPKSegment(file: file,array: 0)
                #expect(segment.polynomialDegree == 4)
                #expect(segment.recordCount == 2)
                for t in [-20.0,-18.75,-15,-12.5,-10,-7.5,-5,-1.25,0] {
                    let record = t < -10 ? 0 : 1, mid = record == 0 ? -15.0 : -5.0
                    let expected = reference((t-mid)/5,record: record,type: type,radius: 5)
                    let state = try segment.state(atSecondsTDB: t)
                    #expect(close(state.positionKilometres,expected.0))
                    #expect(close(state.velocityKilometresPerTDBSecond,expected.1))
                    #expect(state.recordIndex == record)
                    #expect(state.secondsTDB == t)
                    #expect(state.descriptor.targetID == 399)
                    #expect(state.descriptor.centreID == 0)
                    #expect(state.descriptor.frameID == 1)
                    #expect(state.sourcePath == url.path)
                }
            }
        } }
    }
    @Test func boundarySelectionAndDescriptorSubinterval() throws {
        var fixture = SPKFixture(initial: 10,interval: 4)
        fixture.daf.double(11,2048+24); fixture.daf.double(17,2048+32)
        try files.withFile(fixture.daf.bytes) { url in
            let segment = try SpiceSPKSegment(file: SpiceDAFFile(url: url),array: 0)
            #expect(try segment.state(atSecondsTDB: 14.0.nextDown).recordIndex == 0)
            #expect(try segment.state(atSecondsTDB: 14).recordIndex == 1)
            #expect(try segment.state(atSecondsTDB: 14.0.nextUp).recordIndex == 1)
            #expect(throws: SpiceSPKError.outOfCoverage) { try segment.state(atSecondsTDB: 11.0.nextDown) }
            #expect(throws: SpiceSPKError.outOfCoverage) { try segment.state(atSecondsTDB: 17.0.nextUp) }
            #expect(throws: SpiceSPKError.invalidEpoch) { try segment.state(atSecondsTDB: .nan) }
        }
        for initial in [0.1,1e9,-1e9] {
            let interval = 0.3, f = SPKFixture(initial: initial,interval: interval)
            try files.withFile(f.daf.bytes) { url in
                let s = try SpiceSPKSegment(file: SpiceDAFFile(url: url),array: 0)
                let boundary = initial+interval, end = initial+interval*2
                #expect(try s.state(atSecondsTDB: boundary.nextDown).recordIndex == 0)
                #expect(try s.state(atSecondsTDB: boundary).recordIndex == 1)
                #expect(try s.state(atSecondsTDB: end).recordIndex == 1)
            }
        }
    }
    @Test func constantsHigherDegreeAndCrossRecordReads() throws {
        for type: Int32 in [2,3] {
            let coefficients = type == 2 ? [[2.0],[-3],[5]] : [[2.0],[-3],[5],[7],[8],[-9]]
            let fixture = SPKFixture(type: type,coefficients: coefficients,records: 1)
            try files.withFile(fixture.daf.bytes) { url in
                let s = try SpiceDAFFile(url: url).spkState(array: 0,atSecondsTDB: -15)
                #expect(s.positionKilometres == SIMD3(2,-3,5))
                #expect(s.velocityKilometresPerTDBSecond == (type == 2 ? .zero : SIMD3(7,8,-9)))
            }
        }
        var cs = Array(repeating: Array(repeating: 0.0,count: 21),count: 3)
        cs[0][20] = 1
        let f = SPKFixture(coefficients: cs)
        try files.withFile(f.daf.bytes) { url in
            let segment = try SpiceSPKSegment(file: SpiceDAFFile(url: url),array: 0)
            for x in [-1.0,-0.75,0,0.5,1] {
                let t = -5 + 5*x, state = try segment.state(atSecondsTDB: t)
                #expect(abs(state.positionKilometres.x - (100+cos(20*acos(x)))) < 1e-12)
                let d = abs(x) == 1 ? (x > 0 ? 400.0 : -400.0) : 20*sin(20*acos(x))/sqrt(1-x*x)
                #expect(abs(state.velocityKilometresPerTDBSecond.x-d/5) < 1e-11)
            }
        }
    }
    @Test func unsupportedOrCorruptHigherPrioritySegmentDoesNotFallBack() throws {
        try files.withFile(SPKFixture().daf.bytes) { url in
            let set = try SpiceKernelSet(), first = try set.load(url)
            #expect(try set.snapshot().spkState(targetID: 399,atSecondsTDB: -15).positionKilometres.x == 4.25)
            var unsupported = SPKFixture()
            unsupported.daf.int(99,2048+24+28)
            try files.withFile(unsupported.daf.bytes) { badURL in
                let second = try set.load(badURL)
                #expect(throws: SpiceSPKError.unsupportedType(99)) { try set.snapshot().spkState(targetID: 399,atSecondsTDB: -15) }
                try set.unload(second)
                #expect(try set.snapshot().spkState(targetID: 399,atSecondsTDB: -15).sourcePath == url.path)
            }
            var corrupt = SPKFixture(); corrupt.daf.double(.nan,corrupt.directoryOffset)
            try files.withFile(corrupt.daf.bytes) { badURL in
                let second = try set.load(badURL)
                #expect(throws: SpiceSPKError.self) { try set.snapshot().spkState(targetID: 399,atSecondsTDB: -15) }
                try set.unload(second)
            }
            #expect(throws: SpiceSPKError.noSegment(targetID: 3,secondsTDB: -15)) { try set.snapshot().spkState(targetID: 3,atSecondsTDB: -15) }
            #expect(throws: SpiceSPKError.noSegment(targetID: 399,secondsTDB: 1)) { try set.snapshot().spkState(targetID: 399,atSecondsTDB: 1) }
            #expect(throws: SpiceSPKError.invalidEpoch) { try set.snapshot().spkState(targetID: 399,atSecondsTDB: .infinity) }
            try set.unload(first)
        }
    }
    @Test func malformedDirectoryAndRecords() throws {
        var variants: [SPKFixture] = []
        for (offset,value) in [(0,Double.nan),(8,0.0),(8,-1.0),(16,17.5),(16,18.0),(24,0.0),(24,3.0),(24,Double.greatestFiniteMagnitude)] {
            var f = SPKFixture(); f.daf.double(value,f.directoryOffset+offset); variants.append(f)
        }
        var f = SPKFixture(); f.daf.double(-21,2048+24); variants.append(f)
        f = SPKFixture(); f.daf.double(1,2048+32); variants.append(f)
        f = SPKFixture(); f.daf.double(0,4096+8); variants.append(f)
        f = SPKFixture(); f.daf.double(-14,4096); variants.append(f)
        f = SPKFixture(); f.daf.double(4,4096+8); variants.append(f)
        f = SPKFixture(); f.daf.double(.infinity,4096+16); variants.append(f)
        f = SPKFixture(); f.daf.double(Double.greatestFiniteMagnitude,4096+16+4*8); variants.append(f)
        for fixture in variants {
            try files.withFile(fixture.daf.bytes) { url in
                #expect(throws: SpiceSPKError.self) { try SpiceDAFFile(url: url).spkState(array: 0,atSecondsTDB: -14) }
            }
        }
        try files.withFile(SPKFixture().daf.bytes) { url in
            let bounded = try SpiceDAFFile(url: url,limits: SpiceDAFLimits(maximumReadWords: 4))
            #expect(throws: SpiceDAFError.self) { try SpiceSPKSegment(file: bounded,array: 0) }
            #expect(throws: SpiceSPKError.invalidSegment) { try SpiceSPKSegment(file: bounded,array: -1) }
        }
    }
    @Test func laterArrayWinsOnlyInsideItsCoverage() throws {
        var f = SPKFixture()
        let words = f.size*2+4, start = 513+words, end = start+words-1
        let offset = (start-1)*8
        let payload = f.daf.bytes.subdata(in: 4096..<4096+words*8)
        f.daf.bytes.replaceSubrange(offset..<offset+words*8,with: payload)
        f.daf.summary(record: 3,previous: 0,next: 0,count: 2)
        f.daf.array(record: 3,slot: 1,target: 399,start: Int32(start),end: Int32(end))
        f.daf.double(-18,2048+24+40); f.daf.double(-12,2048+32+40)
        f.daf.double(50,offset+16); f.daf.int(Int32(end+1),84)
        try files.withFile(f.daf.bytes) { url in
            let set = try SpiceKernelSet(); _ = try set.load(url)
            let view = set.snapshot()
            #expect(try view.spkState(targetID: 399,atSecondsTDB: -15).arrayIndex == 1)
            #expect(try view.spkState(targetID: 399,atSecondsTDB: -15).positionKilometres.x == 51.25)
            #expect(try view.spkState(targetID: 399,atSecondsTDB: -19).arrayIndex == 0)
            #expect(try view.spkState(targetID: 399,atSecondsTDB: -18).arrayIndex == 1)
            #expect(try view.spkState(targetID: 399,atSecondsTDB: -12).arrayIndex == 1)
            #expect(try view.spkState(targetID: 399,atSecondsTDB: -11).arrayIndex == 0)
        }
    }

}
