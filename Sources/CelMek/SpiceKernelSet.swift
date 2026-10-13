// SPDX-License-Identifier: Apache-2.0
import Foundation

public enum SpiceKernelError: Error, Equatable {
    case invalidLimits
    case resourceLimit(String)
    case unsupportedFormat(String)
    case syntax(source: String, line: Int, reason: String)
    case dependency(source: String, path: String, reason: String)
    case recursiveManifest([String])
    case unknownLoad
}

/// Resource limits are library policy, not NAIF toolkit capacity constants.
public struct SpiceKernelLimits: Sendable {
    public let maximumFileBytes: Int
    public let maximumFiles: Int
    public let maximumDepth: Int
    public let maximumAssignments: Int
    public let maximumValues: Int
    public init(maximumFileBytes: Int = 8_388_608, maximumFiles: Int = 256,
                maximumDepth: Int = 16, maximumAssignments: Int = 100_000,
                maximumValues: Int = 1_000_000) throws {
        guard maximumFileBytes > 0, maximumFiles > 0, maximumDepth > 0,
              maximumAssignments > 0, maximumValues > 0 else { throw SpiceKernelError.invalidLimits }
        self.maximumFileBytes = maximumFileBytes; self.maximumFiles = maximumFiles
        self.maximumDepth = maximumDepth; self.maximumAssignments = maximumAssignments
        self.maximumValues = maximumValues
    }
}

public enum SpiceKernelValue: Equatable, Sendable {
    case numbers([Double])
    case strings([String])
}
public struct SpiceKernelLoadID: Hashable, Sendable {
    fileprivate let rawValue: UUID
}
public struct SpiceKernelSource: Equatable, Sendable {
    public let path: String
    public let kind: String
    /// Immediate manifest, if this file was loaded as a dependency.
    public let parentManifest: String?
}
public struct SpiceKernelVariable: Equatable, Sendable {
    public let value: SpiceKernelValue
    /// Sources contributing to the effective value, in assignment order.
    public let sources: [SpiceKernelSource]
}
public struct SpiceKernelLoad: Sendable {
    public let id: SpiceKernelLoadID
    public let rootPath: String
    public let sources: [SpiceKernelSource]
}

/// Immutable evaluation view. Retained views remain valid after set mutation.
public struct SpiceKernelSnapshot: Sendable {
    /// Changes after every successful mutation; suitable as an external cache namespace.
    public let cacheScope: EphemerisCacheScope
    public let loads: [SpiceKernelLoad]
    public let variables: [String: SpiceKernelVariable]
    public subscript(name: String) -> SpiceKernelVariable? { variables[name] }
}

private struct SpiceAssignment {
    let name: String
    let append: Bool
    let value: SpiceKernelValue
    let line: Int
}
private struct SpiceLoadedGroup {
    let metadata: SpiceKernelLoad
    let assignments: [(SpiceAssignment, SpiceKernelSource)]
}

/// Caller-owned kernel pool. Mutations are serialized; each load commits atomically.
/// Only the documented numeric/string text subset is currently supported.
public final class SpiceKernelSet: @unchecked Sendable {
    public let limits: SpiceKernelLimits
    private let lock = NSLock()
    private var groups: [SpiceLoadedGroup] = []
    private var scope = EphemerisCacheScope()
    public init(limits: SpiceKernelLimits? = nil) throws {
        self.limits = try limits ?? SpiceKernelLimits()
    }
    public func snapshot() -> SpiceKernelSnapshot {
        lock.lock(); defer { lock.unlock() }
        return SpiceKernelSnapshot(cacheScope: scope, loads: groups.map(\.metadata), variables: pool(groups))
    }
    /// Paths are resolved relative to the containing manifest, independent of process cwd.
    /// A repeated load is a new ownership group and raises its data priority.
    @discardableResult public func load(_ url: URL) throws -> SpiceKernelLoadID {
        lock.lock(); defer { lock.unlock() }
        var sources: [SpiceKernelSource] = []
        var assignments: [(SpiceAssignment, SpiceKernelSource)] = []
        var assignmentCount = 0, valueCount = 0
        func visit(_ input: URL, parent: String?, stack: [String]) throws {
            let file = input.standardizedFileURL.resolvingSymlinksInPath()
            let path = file.path
            guard !stack.contains(path) else { throw SpiceKernelError.recursiveManifest(stack + [path]) }
            guard stack.count < limits.maximumDepth else { throw SpiceKernelError.resourceLimit("dependency depth") }
            guard sources.count < limits.maximumFiles else { throw SpiceKernelError.resourceLimit("file count") }
            let data: Data
            do {
                let handle = try FileHandle(forReadingFrom: file)
                defer { try? handle.close() }
                data = try handle.read(upToCount: limits.maximumFileBytes + (limits.maximumFileBytes < Int.max ? 1 : 0)) ?? Data()
            } catch {
                throw SpiceKernelError.dependency(source: parent ?? path, path: path, reason: "cannot read file: \(error)")
            }
            guard data.count <= limits.maximumFileBytes else { throw SpiceKernelError.resourceLimit("file bytes: \(path)") }
            var parser = try SpiceTextParser(data: data, source: path, limits: limits)
            let parsed = try parser.parse()
            let source = SpiceKernelSource(path: path, kind: parser.kind, parentManifest: parent)
            sources.append(source)
            assignmentCount += parsed.count
            valueCount += parsed.reduce(0) { $0 + $1.value.count }
            guard assignmentCount <= limits.maximumAssignments, valueCount <= limits.maximumValues else {
                throw SpiceKernelError.resourceLimit("assignment/value count")
            }
            if source.kind == "MK" {
                let local = pool([SpiceLoadedGroup(metadata: SpiceKernelLoad(id: SpiceKernelLoadID(rawValue: UUID()), rootPath: path, sources: [source]), assignments: parsed.map { ($0, source) })])
                func strings(_ name: String) throws -> [String] {
                    guard let item = local[name] else { return [] }
                    guard case .strings(let result) = item.value else {
                        throw SpiceKernelError.dependency(source: path, path: name, reason: "manifest variable must contain strings")
                    }
                    return try continued(result, source: path)
                }
                let symbols = try strings("PATH_SYMBOLS"), replacements = try strings("PATH_VALUES")
                guard symbols.count == replacements.count, Set(symbols).count == symbols.count,
                      symbols.allSatisfy({ !$0.isEmpty && $0.utf8.allSatisfy { ($0 >= 65 && $0 <= 90) || ($0 >= 97 && $0 <= 122) || ($0 >= 48 && $0 <= 57) || $0 == 95 } }) else {
                    throw SpiceKernelError.dependency(source: path, path: "PATH_SYMBOLS", reason: "invalid or mismatched substitutions")
                }
                for raw in try strings("KERNELS_TO_LOAD") {
                    var expanded = "", i = raw.startIndex
                    while i < raw.endIndex {
                        if raw[i] != "$" { expanded.append(raw[i]); i = raw.index(after: i); continue }
                        let start = raw.index(after: i)
                        var end = start
                        while end < raw.endIndex, raw[end].isASCII, raw[end].isLetter || raw[end].isNumber || raw[end] == "_" { end = raw.index(after: end) }
                        let name = String(raw[start..<end])
                        guard let index = symbols.firstIndex(of: name) else {
                            throw SpiceKernelError.dependency(source: path, path: raw, reason: "undefined path symbol \(name)")
                        }
                        expanded += replacements[index]; i = end
                    }
                    guard !expanded.isEmpty, expanded.utf8.count <= 255 else {
                        throw SpiceKernelError.dependency(source: path, path: expanded, reason: "empty or oversized dependency path")
                    }
                    let child = expanded.hasPrefix("/") ? URL(fileURLWithPath: expanded) : file.deletingLastPathComponent().appendingPathComponent(expanded)
                    try visit(child, parent: path, stack: stack + [path])
                }
                // Manifest assignments configure loading, not the persistent kernel pool.
            } else { assignments.append(contentsOf: parsed.map { ($0, source) }) }
        }
        try visit(url, parent: nil, stack: [])
        let id = SpiceKernelLoadID(rawValue: UUID())
        let candidate = SpiceLoadedGroup(metadata: SpiceKernelLoad(id: id, rootPath: url.standardizedFileURL.path, sources: sources), assignments: assignments)
        let proposed = groups + [candidate]
        guard proposed.reduce(0, { $0 + $1.metadata.sources.count }) <= limits.maximumFiles,
              proposed.reduce(0, { $0 + $1.assignments.count }) <= limits.maximumAssignments,
              proposed.reduce(0, { $0 + $1.assignments.reduce(0, { $0 + $1.0.value.count }) }) <= limits.maximumValues else {
            throw SpiceKernelError.resourceLimit("loaded kernel set capacity")
        }
        groups = proposed; scope = EphemerisCacheScope()
        return id
    }
    public func unload(_ id: SpiceKernelLoadID) throws {
        lock.lock(); defer { lock.unlock() }
        guard let index = groups.firstIndex(where: { $0.metadata.id == id }) else { throw SpiceKernelError.unknownLoad }
        groups.remove(at: index); scope = EphemerisCacheScope()
    }
    public func removeAll() {
        lock.lock(); defer { lock.unlock() }
        groups.removeAll(); scope = EphemerisCacheScope()
    }
}
private extension SpiceKernelValue {
    var count: Int { switch self { case .numbers(let a): return a.count; case .strings(let a): return a.count } }
}
private func pool(_ groups: [SpiceLoadedGroup]) -> [String: SpiceKernelVariable] {
    var result: [String: SpiceKernelVariable] = [:]
    for group in groups { for (assignment, source) in group.assignments {
        var value = assignment.value, sources = [source]
        if assignment.append, let previous = result[assignment.name] {
            switch (previous.value, value) {
            case (.numbers(let a), .numbers(let b)): value = .numbers(a + b); sources = previous.sources + sources
            case (.strings(let a), .strings(let b)): value = .strings(a + b); sources = previous.sources + sources
            default: break // A different type replaces the previous pool value.
            }
        }
        result[assignment.name] = SpiceKernelVariable(value: value, sources: sources)
    } }
    return result
}
private func continued(_ elements: [String], source: String) throws -> [String] {
    var result: [String] = [], pending = "", awaitingContinuation = false
    for element in elements {
        let trimmed = element.replacingOccurrences(of: #" +$"#, with: "", options: .regularExpression)
        if trimmed.hasSuffix("+") { pending += trimmed.dropLast(); awaitingContinuation = true }
        else { result.append(pending + element); pending = ""; awaitingContinuation = false }
    }
    guard !awaitingContinuation else { throw SpiceKernelError.dependency(source: source, path: pending, reason: "unterminated path continuation") }
    return result
}

private struct SpiceTextParser {
    let source: String
    let limits: SpiceKernelLimits
    let kind: String
    var bytes: [UInt8] = []
    var lines: [Int] = []
    var index = 0
    init(data: Data, source: String, limits: SpiceKernelLimits) throws {
        self.source = source; self.limits = limits
        guard data.allSatisfy({ (32...126).contains($0) || [9,10,13].contains($0) }), let text = String(data: data, encoding: .ascii) else {
            throw SpiceKernelError.unsupportedFormat("non-ASCII text or binary kernel: \(source)")
        }
        let rows = text.replacingOccurrences(of: "\r\n", with: "\n").replacingOccurrences(of: "\r", with: "\n").components(separatedBy: "\n")
        let header = rows[0].trimmingCharacters(in: .whitespaces)
        guard ["KPL/PCK", "KPL/FK", "KPL/LSK", "KPL/MK"].contains(header) else { throw SpiceKernelError.unsupportedFormat(header) }
        kind = String(header.dropFirst(4))
        var active = false
        for (line, row) in rows.enumerated().dropFirst() {
            let control = row.trimmingCharacters(in: .whitespaces)
            if control == "\\begindata" { active = true; continue }
            if control == "\\begintext" { active = false; continue }
            if active {
                bytes += row.utf8; bytes.append(10)
                lines += Array(repeating: line + 1, count: row.utf8.count + 1)
            }
        }
    }
    func error(_ reason: String) -> SpiceKernelError {
        .syntax(source: source, line: lines.isEmpty ? 1 : lines[min(index, lines.count - 1)], reason: reason)
    }
    mutating func whitespace() { while index < bytes.count, [9,10,32].contains(bytes[index]) { index += 1 } }
    mutating func parse() throws -> [SpiceAssignment] {
        var result: [SpiceAssignment] = [], count = 0
        while true {
            whitespace(); if index == bytes.count { return result }
            let start = index, line = lines[index]
            while index < bytes.count, ![9,10,32,44,40,41,61].contains(bytes[index]), !(bytes[index] == 43 && index + 1 < bytes.count && bytes[index + 1] == 61) { index += 1 }
            let name = String(decoding: bytes[start..<index], as: UTF8.self)
            guard !name.isEmpty, name.utf8.count <= 32 else { throw error("invalid variable name") }
            whitespace()
            var append = false
            if index < bytes.count, bytes[index] == 43 { append = true; index += 1 }
            guard index < bytes.count, bytes[index] == 61 else { throw error("expected assignment operator") }
            index += 1; whitespace()
            let vector = index < bytes.count && bytes[index] == 40
            if vector { index += 1 }
            var numbers: [Double] = [], strings: [String] = []
            while true {
                whitespace()
                if vector, index < bytes.count, bytes[index] == 41 { index += 1; break }
                guard index < bytes.count else { throw error("missing value or closing parenthesis") }
                if bytes[index] == 39 {
                    index += 1; var string: [UInt8] = []; var closed = false
                    while index < bytes.count {
                        let byte = bytes[index]; index += 1
                        if byte == 39 {
                            if index < bytes.count, bytes[index] == 39 { string.append(39); index += 1 }
                            else { closed = true; break }
                        } else { guard byte != 10 else { throw error("string crosses line") }; string.append(byte) }
                    }
                    guard closed, string.count <= 80 else { throw error("unterminated or oversized string") }
                    strings.append(String(decoding: string, as: UTF8.self))
                } else {
                    let start = index
                    while index < bytes.count, ![9,10,32,44,40,41].contains(bytes[index]) { index += 1 }
                    let token = String(decoding: bytes[start..<index], as: UTF8.self)
                    guard !token.hasPrefix("@") else { throw error("@ date literals are not supported yet") }
                    let normalized = token.replacingOccurrences(of: "D", with: "E").replacingOccurrences(of: "d", with: "e")
                    guard normalized.range(of: #"^[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$"#, options: .regularExpression) != nil,
                          let number = Double(normalized), number.isFinite else { throw error("invalid finite number") }
                    numbers.append(number)
                }
                count += 1
                guard count <= limits.maximumValues else { throw SpiceKernelError.resourceLimit("parsed values") }
                guard numbers.isEmpty || strings.isEmpty else { throw error("mixed string/numeric vector") }
                if !vector { break }
                let before = index; whitespace()
                if index < bytes.count, bytes[index] == 44 { index += 1 }
                else if index < bytes.count, bytes[index] != 41, index == before { throw error("missing value separator") }
            }
            guard !numbers.isEmpty || !strings.isEmpty else { throw error("empty vector") }
            result.append(SpiceAssignment(name: name, append: append, value: strings.isEmpty ? .numbers(numbers) : .strings(strings), line: line))
            guard result.count <= limits.maximumAssignments else { throw SpiceKernelError.resourceLimit("parsed assignments") }
        }
    }
}
