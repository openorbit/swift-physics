// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
@testable import CelMek

/// Generated directly from the published DAF header/summary layout; no toolkit writer.
struct DAFFixture {
    var bytes = Data(repeating: 0, count: 8*1024)
    let big: Bool
    mutating func text(_ value: String, _ offset: Int, _ width: Int) {
        bytes.replaceSubrange(offset..<offset+width, with: Array(value.utf8) + Array(repeating: 32, count: width-value.utf8.count))
    }
    mutating func int(_ value: Int32, _ offset: Int) {
        let bits = UInt32(bitPattern: value)
        let encoded = (0..<4).map { UInt8(truncatingIfNeeded: bits >> (8*(big ? 3-$0 : $0))) }
        bytes.replaceSubrange(offset..<offset+4, with: encoded)
    }
    mutating func double(_ value: Double, _ offset: Int) {
        let bits = value.bitPattern
        bytes.replaceSubrange(offset..<offset+8, with: (0..<8).map { UInt8(truncatingIfNeeded: bits >> (8*(big ? 7-$0 : $0))) })
    }
    mutating func summary(record: Int, previous: Int, next: Int, count: Int) {
        let offset = (record-1)*1024
        double(Double(next),offset); double(Double(previous),offset+8); double(Double(count),offset+16)
    }
    mutating func array(record: Int, slot: Int, target: Int32, type: Int32 = 2, start: Int32, end: Int32) {
        let offset = (record-1)*1024+24+slot*40
        double(-100,offset); double(200,offset+8)
        for (i,v) in [target,Int32(0),Int32(1),type,start,end].enumerated() { int(v,offset+16+4*i) }
        text("Synthetic \(target)", record*1024+slot*40,40)
    }
    init(big: Bool = false) {
        self.big = big
        text("DAF/SPK",0,8); int(2,8); int(6,12); text("Independent synthetic DAF",16,60)
        int(3,76); int(6,80); int(1025,84); text(big ? "BIG-IEEE" : "LTL-IEEE",88,8)
        let ftp = Array("FTPSTR:".utf8)+[13,58,10,58,13,10,58,13,0,58,129,58,16,206,58]+Array("ENDFTP".utf8)
        bytes.replaceSubrange(699..<727,with: ftp)
        bytes.replaceSubrange(1024..<1032,with: Array("Comments".utf8))
        summary(record: 3,previous: 0,next: 6,count: 2)
        array(record: 3,slot: 0,target: -123,start: 513,end: 515)
        array(record: 3,slot: 1,target: 399,type: 99,start: 516,end: 518)
        summary(record: 6,previous: 3,next: 0,count: 1)
        array(record: 6,slot: 0,target: 301,type: 3,start: 897,end: 1024)
        for (i,v) in [1.25,-2.5,3.75,4,5,6].enumerated() { double(v,4096+i*8) }
        double(123,8192-8)
    }
}
@Suite struct SpiceDAFTests {
    func withFile(_ bytes: Data, _ action: (URL) throws -> Void) throws {
        let url = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString + ".bsp")
        try bytes.write(to: url); defer { try? FileManager.default.removeItem(at: url) }
        try action(url)
    }
    @Test func bothByteOrdersAndBoundedArrayReads() throws {
        for big in [false,true] {
            try withFile(DAFFixture(big: big).bytes) { url in
                let file = try SpiceDAFFile(url: url), meta = file.metadata
                #expect(meta.byteOrder == (big ? .bigEndian : .littleEndian))
                #expect(meta.internalName == "Independent synthetic DAF")
                #expect(meta.arrays.count == 3)
                #expect(meta.arrays[0].spk?.targetID == -123)
                #expect(meta.arrays[1].spk?.dataType == 99)
                #expect(meta.arrays[2].spk?.startSecondsTDB == -100)
                #expect(meta.arrays[2].spk?.endSecondsTDB == 200)
                #expect(meta.reservedRecords.count == 1024)
                #expect(String(decoding: meta.reservedRecords.prefix(8),as: UTF8.self) == "Comments")
                #expect(try file.readWords(array: 0,count: 3) == [1.25,-2.5,3.75])
                #expect(try file.readWords(array: 1,offset: 1,count: 2) == [5,6])
                #expect(try file.readWords(array: 2,offset: 127,count: 1) == [123])
                #expect(throws: SpiceDAFError.invalidRange) { try file.readWords(array: 0,offset: 2,count: 2) }
                #expect(throws: SpiceDAFError.invalidRange) { try file.readWords(array: -1,count: 1) }
                #expect(throws: SpiceDAFError.self) { try file.readWords(array: 0,count: Int.max) }
            }
        }
    }
    @Test func malformedHeadersLinksDescriptorsAndRanges() throws {
        var variants: [Data] = []
        func changed(_ edit: (inout DAFFixture) -> Void) {
            var f = DAFFixture(); edit(&f); variants.append(f.bytes)
        }
        changed { $0.int(-1,8) }
        changed { $0.int(5,12) }
        changed { $0.int(0,76) }
        changed { $0.int(3,80) }
        changed { $0.int(100,84) }
        changed { $0.bytes[699] = 0 }
        changed { $0.double(.nan,2048) }
        changed { $0.double(3.5,2048) }
        changed { $0.double(3,5120) }
        changed { $0.double(4,5120+8) }
        changed { $0.double(26,2048+16) }
        changed { $0.double(300,2048+24) }
        changed { $0.int(512,2048+24+36) }
        changed { $0.int(0,2048+24+32) }
        changed { $0.int(300,2048+24+32) }
        changed { $0.int(514,2048+24+40+32) }
        changed { $0.bytes[3072] = 255 }
        variants.append(Data(DAFFixture().bytes.dropLast()))
        for bytes in variants { try withFile(bytes) { url in
            #expect(throws: SpiceDAFError.self) { try SpiceDAFFile(url: url) }
        } }
    }
    @Test func unsupportedFormatsAndLimits() throws {
        for format in ["VAX-GFLT", "        "] {
            var f = DAFFixture(); f.text(format,88,8)
            try withFile(f.bytes) { url in #expect(throws: SpiceDAFError.self) { try SpiceDAFFile(url: url) } }
        }
        try withFile(DAFFixture().bytes) { url in
            for limits in [try SpiceDAFLimits(maximumSummaryRecords: 1),try SpiceDAFLimits(maximumArrays: 2),try SpiceDAFLimits(maximumReservedBytes: 10)] {
                #expect(throws: SpiceDAFError.self) { try SpiceDAFFile(url: url,limits: limits) }
            }
            let limited = try SpiceDAFFile(url: url,limits: SpiceDAFLimits(maximumReadWords: 2))
            #expect(throws: SpiceDAFError.self) { try limited.readWords(array: 0,count: 3) }
        }
    }
    @Test func kernelSetBinaryOwnershipAndManifestRollback() throws {
        try withFile(DAFFixture().bytes) { url in
            let set = try SpiceKernelSet(limits: SpiceKernelLimits(maximumFileBytes: 1))
            let id = try set.load(url), snapshot = set.snapshot()
            #expect(snapshot.binaryKernels.count == 1)
            #expect(snapshot.loads[0].sources[0].kind == "SPK")
            try set.unload(id)
            #expect(set.snapshot().binaryKernels.isEmpty)
            #expect(try snapshot.binaryKernels[0].file.readWords(array: 0,count: 1) == [1.25])
            let manifest = url.deletingPathExtension().appendingPathExtension("tm")
            defer { try? FileManager.default.removeItem(at: manifest) }
            try Data("KPL/MK\n\\begindata\nKERNELS_TO_LOAD=('\(url.lastPathComponent)' 'missing.bsp')".utf8).write(to: manifest)
            let ordinary = try SpiceKernelSet(), before = ordinary.snapshot()
            #expect(throws: SpiceKernelError.self) { try ordinary.load(manifest) }
            #expect(ordinary.snapshot().cacheScope == before.cacheScope)
            #expect(ordinary.snapshot().binaryKernels.isEmpty)
        }
    }
    @Test func readsAcrossRecordBoundariesAndDetectsTruncation() throws {
        var fixture = DAFFixture()
        fixture.bytes.append(Data(repeating: 0,count: 1024))
        fixture.int(1153,84)
        fixture.int(1152,5120+24+36)
        fixture.double(456,8192)
        try withFile(fixture.bytes) { url in
            let reader = try SpiceDAFFile(url: url)
            #expect(try reader.readWords(array: 2,offset: 127,count: 2) == [123,456])
            let writer = try FileHandle(forWritingTo: url)
            try writer.truncate(atOffset: 4096); try writer.close()
            #expect(throws: SpiceDAFError.self) { try reader.readWords(array: 2,count: 1) }
        }
    }
    @Test func concurrentFileReadsUseIndependentPositions() async throws {
        let url = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try DAFFixture(big: true).bytes.write(to: url)
        defer { try? FileManager.default.removeItem(at: url) }
        let reader = try SpiceDAFFile(url: url)
        try await withThrowingTaskGroup(of: Bool.self) { group in
            for i in 0..<100 {
                group.addTask {
                    let index = i % 2
                    return try reader.readWords(array: index,count: 3) == (index == 0 ? [1.25,-2.5,3.75] : [4,5,6])
                }
            }
            for try await correct in group { #expect(correct) }
        }
    }

    @Test func invalidCKDimensionsAndEmptyRecords() throws {
        for big in [false,true] {
            var fixture = DAFFixture(big: big)
            fixture.summary(record: 3,previous: 0,next: 6,count: 1)
            fixture.summary(record: 6,previous: 3,next: 0,count: 0)
            try withFile(fixture.bytes) { url in
                let file = try SpiceDAFFile(url: url)
                #expect(file.metadata.arrays.count == 1)
                #expect(try file.readWords(array: 0,count: 3) == [1.25,-2.5,3.75])
            }
            fixture.text("DAF/CK",0,8); fixture.int(0,8); fixture.int(3,12)
            try withFile(fixture.bytes) { url in
                #expect(throws: SpiceDAFError.malformed("CK requires ND=2 NI=6")) { try SpiceDAFFile(url: url) }
            }
        }
    }

}
