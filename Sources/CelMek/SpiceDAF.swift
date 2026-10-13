// SPDX-License-Identifier: Apache-2.0
import Foundation

public enum SpiceDAFError: Error, Equatable {
    case invalidLimits
    case unsupportedFormat(String)
    case malformed(String)
    case resourceLimit(String)
    case invalidRange
    case io(String)
}
public enum SpiceDAFByteOrder: String, Sendable { case littleEndian = "LTL-IEEE", bigEndian = "BIG-IEEE" }
public struct SpiceDAFLimits: Sendable {
    public let maximumSummaryRecords: Int
    public let maximumArrays: Int
    public let maximumReadWords: Int
    public let maximumReservedBytes: Int
    public init(maximumSummaryRecords: Int = 10_000, maximumArrays: Int = 250_000,
                maximumReadWords: Int = 131_072, maximumReservedBytes: Int = 8_388_608) throws {
        guard maximumSummaryRecords > 0, maximumArrays > 0, maximumReadWords > 0,
              maximumReadWords <= Int.max / 8, maximumReservedBytes >= 0 else { throw SpiceDAFError.invalidLimits }
        self.maximumSummaryRecords = maximumSummaryRecords; self.maximumArrays = maximumArrays
        self.maximumReadWords = maximumReadWords; self.maximumReservedBytes = maximumReservedBytes
    }
}
/// Inclusive ephemeris-second bounds in the kernel's TDB argument, not atomic Epoch seconds.
public struct SpiceSPKDescriptor: Equatable, Sendable {
    public let startSecondsTDB: Double
    public let endSecondsTDB: Double
    public let targetID: Int32
    public let centreID: Int32
    public let frameID: Int32
    public let dataType: Int32
}
public struct SpicePCKDescriptor: Equatable, Sendable {
    public let startSecondsTDB: Double
    public let endSecondsTDB: Double
    public let frameClassID: Int32
    public let frameID: Int32
    public let dataType: Int32
}
/// Bounds are encoded spacecraft-clock ticks, not ephemeris seconds.
public struct SpiceCKDescriptor: Equatable, Sendable {
    public let startEncodedTicks: Double
    public let endEncodedTicks: Double
    public let instrumentID: Int32
    public let frameID: Int32
    public let dataType: Int32
    public let hasAngularVelocity: Bool
}
public struct SpiceDAFArray: Equatable, Sendable {
    public let name: String
    public let doubleComponents: [Double]
    public let integerComponents: [Int32]
    public let initialAddress: Int
    public let finalAddress: Int
    /// Present for DAF/SPK only. Unsupported types remain visible to inspection.
    public let spk: SpiceSPKDescriptor?
    public let pck: SpicePCKDescriptor?
    public let ck: SpiceCKDescriptor?
}
public struct SpiceDAFMetadata: Sendable {
    public let kind: String
    public let internalName: String
    public let byteOrder: SpiceDAFByteOrder
    public let fileSizeBytes: UInt64
    public let doubleComponentCount: Int
    public let integerComponentCount: Int
    public let firstFreeAddress: Int
    public let arrays: [SpiceDAFArray]
    /// Preserved verbatim. Byte-level comment delimiters are not inferred from toolkit code.
    public let reservedRecords: Data
}

/// Read-only modern IEEE DAF reader. File seek/read operations are serialized.
/// The asset must remain unchanged for the lifetime of this reader and all retained snapshots.
public final class SpiceDAFFile: @unchecked Sendable {
    public let url: URL
    public let limits: SpiceDAFLimits
    public let metadata: SpiceDAFMetadata
    private let handle: FileHandle
    private let lock = NSLock()
    public init(url: URL, limits: SpiceDAFLimits? = nil) throws {
        let options = try limits ?? SpiceDAFLimits()
        let opened: FileHandle
        do { opened = try FileHandle(forReadingFrom: url) } catch { throw SpiceDAFError.io(String(describing: error)) }
        do {
            let size = try opened.seekToEnd()
            let reader = DAFReader(handle: opened, size: size)
            let parsed = try reader.inspect(limits: options)
            self.url = url; self.limits = options; self.metadata = parsed; self.handle = opened
        } catch { try? opened.close(); throw error }
    }
    deinit { try? handle.close() }
    /// Zero-based array index and element offset; addresses in metadata are one-based DAF words.
    public func readWords(array index: Int, offset: Int = 0, count: Int) throws -> [Double] {
        guard metadata.arrays.indices.contains(index), offset >= 0, count > 0 else { throw SpiceDAFError.invalidRange }
        guard count <= limits.maximumReadWords else { throw SpiceDAFError.resourceLimit("read words") }
        let array = metadata.arrays[index], length = array.finalAddress - array.initialAddress + 1
        guard offset <= length, count <= length - offset else { throw SpiceDAFError.invalidRange }
        lock.lock(); defer { lock.unlock() }
        let reader = DAFReader(handle: handle, size: metadata.fileSizeBytes)
        let data = try reader.read(offset: UInt64(array.initialAddress - 1 + offset) * 8, count: count * 8)
        return (0..<count).map { dafDouble(data, $0 * 8, metadata.byteOrder) }
    }
}
private func dafInt(_ data: Data, _ offset: Int, _ order: SpiceDAFByteOrder) -> Int32 {
    let raw = data.withUnsafeBytes { $0.loadUnaligned(fromByteOffset: offset, as: UInt32.self) }
    return Int32(bitPattern: order == .littleEndian ? UInt32(littleEndian: raw) : UInt32(bigEndian: raw))
}
private func dafDouble(_ data: Data, _ offset: Int, _ order: SpiceDAFByteOrder) -> Double {
    let raw = data.withUnsafeBytes { $0.loadUnaligned(fromByteOffset: offset, as: UInt64.self) }
    return Double(bitPattern: order == .littleEndian ? UInt64(littleEndian: raw) : UInt64(bigEndian: raw))
}
private func dafText(_ data: Data, _ offset: Int, _ count: Int) throws -> String {
    let bytes = data.subdata(in: offset..<offset + count)
    guard bytes.allSatisfy({ (32...126).contains($0) }) else { throw SpiceDAFError.malformed("non-printable header/name text") }
    return String(decoding: bytes, as: UTF8.self).trimmingCharacters(in: .whitespaces)
}
private struct DAFReader {
    let handle: FileHandle
    let size: UInt64
    func read(offset: UInt64, count: Int) throws -> Data {
        guard count >= 0, offset <= size, UInt64(count) <= size - offset else { throw SpiceDAFError.malformed("read outside file") }
        do {
            try handle.seek(toOffset: offset)
            var result = Data()
            while result.count < count {
                let chunk = try handle.read(upToCount: count - result.count) ?? Data()
                guard !chunk.isEmpty else { throw SpiceDAFError.malformed("truncated file") }
                result.append(chunk)
            }
            return result
        } catch let error as SpiceDAFError { throw error }
        catch { throw SpiceDAFError.io(String(describing: error)) }
    }
    func inspect(limits: SpiceDAFLimits) throws -> SpiceDAFMetadata {
        guard size >= 1024, size % 1024 == 0 else { throw SpiceDAFError.malformed("incomplete DAF records") }
        let header = try read(offset: 0, count: 1024)
        let identifier = try dafText(header, 0, 8)
        guard ["DAF/SPK", "DAF/CK", "DAF/PCK"].contains(identifier) else { throw SpiceDAFError.unsupportedFormat(identifier) }
        let format = try dafText(header, 88, 8)
        guard let order = SpiceDAFByteOrder(rawValue: format) else { throw SpiceDAFError.unsupportedFormat("DAF binary representation: \(format)") }
        // Exact integrity bytes are published in the DAF file-record octal example.
        let ftp: [UInt8] = Array("FTPSTR:".utf8) + [13,58,10,58,13,10,58,13,0,58,129,58,16,206,58] + Array("ENDFTP".utf8)
        guard Array(header[699..<727]) == ftp else { throw SpiceDAFError.malformed("missing or damaged FTP integrity string") }
        let nd = Int(dafInt(header, 8, order)), ni = Int(dafInt(header, 12, order))
        guard (0...124).contains(nd), (2...250).contains(ni), nd + (ni + 1)/2 <= 125 else { throw SpiceDAFError.malformed("invalid summary dimensions") }
        if identifier == "DAF/SPK", nd != 2 || ni != 6 { throw SpiceDAFError.malformed("SPK requires ND=2 NI=6") }
        if identifier == "DAF/CK", nd != 2 || ni != 6 { throw SpiceDAFError.malformed("CK requires ND=2 NI=6") }
        if identifier == "DAF/PCK", nd != 2 || ni != 5 { throw SpiceDAFError.malformed("PCK requires ND=2 NI=5") }
        let first = Int(dafInt(header, 76, order)), last = Int(dafInt(header, 80, order)), free = Int(dafInt(header, 84, order))
        let records = size / 1024
        guard first >= 2, last >= 2, UInt64(first + 1) <= records, UInt64(last + 1) <= records,
              free >= 1, UInt64(free) <= size/8 + 1 else { throw SpiceDAFError.malformed("invalid file pointers") }
        let reservedSize = (first - 2) * 1024
        guard reservedSize <= limits.maximumReservedBytes else { throw SpiceDAFError.resourceLimit("reserved records") }
        let reserved = try read(offset: 1024, count: reservedSize)
        let summaryWords = nd + (ni + 1)/2, nameBytes = summaryWords * 8
        var visited: Set<Int> = [], structural: Set<Int> = [], arrays: [SpiceDAFArray] = []
        var current = first, previous = 0
        func control(_ value: Double, maximum: UInt64) throws -> Int {
            guard value.isFinite, value >= 0, value.rounded(.towardZero) == value,
                  value <= Double(maximum), value < Double(Int.max) else { throw SpiceDAFError.malformed("invalid summary control value") }
            return Int(value)
        }
        while current != 0 {
            guard visited.count < limits.maximumSummaryRecords else { throw SpiceDAFError.resourceLimit("summary records") }
            guard current >= first, UInt64(current + 1) <= records,
                  !structural.contains(current), !structural.contains(current + 1) else { throw SpiceDAFError.malformed("cyclic or overlapping summary/name records") }
            visited.insert(current); structural.insert(current); structural.insert(current + 1)
            let summary = try read(offset: UInt64(current - 1)*1024, count: 1024)
            let names = try read(offset: UInt64(current)*1024, count: 1024)
            let next = try control(dafDouble(summary, 0, order), maximum: records)
            let prev = try control(dafDouble(summary, 8, order), maximum: records)
            let count = try control(dafDouble(summary, 16, order), maximum: UInt64(125/summaryWords))
            guard prev == previous else { throw SpiceDAFError.malformed("broken backward link") }
            guard arrays.count <= limits.maximumArrays - count else { throw SpiceDAFError.resourceLimit("array count") }
            for index in 0..<count {
                let offset = 24 + index * nameBytes
                let doubles = (0..<nd).map { dafDouble(summary, offset + $0*8, order) }
                let integers = (0..<ni).map { dafInt(summary, offset + nd*8 + $0*4, order) }
                let start = Int(integers[ni-2]), end = Int(integers[ni-1])
                guard start > (first - 1)*128, end >= start, end < free,
                      UInt64(end) <= size/8 else { throw SpiceDAFError.malformed("invalid array addresses") }
                guard doubles.allSatisfy(\.isFinite) else { throw SpiceDAFError.malformed("non-finite summary") }
                var spk: SpiceSPKDescriptor?
                var pck: SpicePCKDescriptor?
                var ck: SpiceCKDescriptor?
                if identifier == "DAF/SPK" {
                    guard doubles[0] <= doubles[1], integers[2] != 0, integers[3] > 0 else { throw SpiceDAFError.malformed("invalid SPK descriptor") }
                    spk = SpiceSPKDescriptor(startSecondsTDB: doubles[0], endSecondsTDB: doubles[1], targetID: integers[0], centreID: integers[1], frameID: integers[2], dataType: integers[3])
                }
                if identifier == "DAF/PCK" {
                    guard doubles[0] <= doubles[1], integers[0] != 0, integers[1] != 0, integers[2] > 0 else {
                        throw SpiceDAFError.malformed("invalid PCK descriptor")
                    }
                    pck = SpicePCKDescriptor(startSecondsTDB: doubles[0], endSecondsTDB: doubles[1],
                        frameClassID: integers[0], frameID: integers[1], dataType: integers[2])
                }
                if identifier == "DAF/CK" {
                    guard doubles[0] >= 0, doubles[0] <= doubles[1], integers[0] != 0,
                          integers[1] != 0, integers[2] > 0, [0,1].contains(integers[3]) else {
                        throw SpiceDAFError.malformed("invalid CK descriptor")
                    }
                    ck = SpiceCKDescriptor(startEncodedTicks: doubles[0], endEncodedTicks: doubles[1],
                        instrumentID: integers[0], frameID: integers[1], dataType: integers[2], hasAngularVelocity: integers[3] == 1)
                }
                arrays.append(SpiceDAFArray(name: try dafText(names, index*nameBytes, nameBytes), doubleComponents: doubles,
                    integerComponents: integers, initialAddress: start, finalAddress: end, spk: spk, pck: pck, ck: ck))
            }
            if next == 0, current != last { throw SpiceDAFError.malformed("tail pointer mismatch") }
            previous = current; current = next
        }
        // Detect element ranges aliasing structural records or other arrays, regardless of list order.
        let structuralSorted = structural.sorted()
        let ranges = arrays.sorted { $0.initialAddress < $1.initialAddress }
        var priorEnd = 0, recordIndex = 0
        for array in ranges {
            guard array.initialAddress > priorEnd else { throw SpiceDAFError.malformed("overlapping arrays") }
            priorEnd = array.finalAddress
            let low = (array.initialAddress - 1)/128 + 1, high = (array.finalAddress - 1)/128 + 1
            while recordIndex < structuralSorted.count, structuralSorted[recordIndex] < low { recordIndex += 1 }
            guard recordIndex == structuralSorted.count || structuralSorted[recordIndex] > high else { throw SpiceDAFError.malformed("array aliases summary/name records") }
        }
        guard structural.allSatisfy({ $0*128 < free }) else { throw SpiceDAFError.malformed("structural records beyond free address") }
        return SpiceDAFMetadata(kind: String(identifier.dropFirst(4)), internalName: try dafText(header, 16, 60), byteOrder: order,
            fileSizeBytes: size, doubleComponentCount: nd, integerComponentCount: ni, firstFreeAddress: free,
            arrays: arrays, reservedRecords: reserved)
    }
}
