// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
@testable import CelMek

@Suite struct SpiceKernelSetTests {
    func directory() throws -> URL {
        let url = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: url, withIntermediateDirectories: true)
        return url
    }
    func write(_ text: String, _ name: String, in directory: URL) throws -> URL {
        let url = directory.appendingPathComponent(name)
        try Data(text.utf8).write(to: url)
        return url
    }
    @Test func precedenceUnloadSnapshotsAndIsolation() throws {
        let dir = try directory(); defer { try? FileManager.default.removeItem(at: dir) }
        let a = try write("KPL/PCK\ncomments ignored\n\\begindata\nX = (1, 2D1)\nX += 3\nMessage = 'You can''t'\n\\begintext\nX=100", "a.tpc", in: dir)
        let b = try write("KPL/FK\n\\begindata\nX = 9\nX += 10", "b.tf", in: dir)
        let set = try SpiceKernelSet(), other = try SpiceKernelSet()
        let aID = try set.load(a), old = set.snapshot()
        #expect(old["X"]?.value == .numbers([1,20,3]))
        #expect(old["Message"]?.value == .strings(["You can't"]))
        #expect(old["MESSAGE"] == nil)
        let bID = try set.load(b), current = set.snapshot()
        #expect(current["X"]?.value == .numbers([9,10]))
        #expect(current["X"]?.sources.allSatisfy { $0.path == b.path } == true)
        #expect(current.cacheScope != old.cacheScope)
        try set.unload(bID)
        #expect(set.snapshot()["X"]?.value == old["X"]?.value)
        #expect(current["X"]?.value == .numbers([9,10]))
        #expect(other.snapshot().variables.isEmpty)
        try set.unload(aID)
        #expect(set.snapshot().variables.isEmpty)
        #expect(throws: SpiceKernelError.unknownLoad) { try set.unload(aID) }
    }
    @Test func manifestPathsDependenciesReloadAndRollback() throws {
        let dir = try directory(); defer { try? FileManager.default.removeItem(at: dir) }
        let child = try write("KPL/PCK\n\\begindata\nGM = 123", "child.tpc", in: dir)
        let manifest = try write("KPL/MK\n\\begindata\nPATH_SYMBOLS = ('ROOT')\nPATH_VALUES = ('.')\nKERNELS_TO_LOAD = ('$ROOT/chi+' 'ld.tpc')", "set.tm", in: dir)
        let set = try SpiceKernelSet(), id = try set.load(manifest), old = set.snapshot()
        #expect(old["GM"]?.value == .numbers([123]))
        #expect(old["PATH_SYMBOLS"] == nil)
        #expect(old.loads[0].sources.count == 2)
        #expect(old["GM"]?.sources[0].parentManifest == manifest.path)
        _ = try write("KPL/PCK\n\\begindata\nGM = 456", "child.tpc", in: dir)
        let second = try set.load(manifest)
        #expect(set.snapshot()["GM"]?.value == .numbers([456]))
        try set.unload(second)
        #expect(set.snapshot()["GM"]?.value == .numbers([123]))
        let bad = try write("KPL/MK\n\\begindata\nKERNELS_TO_LOAD = ('child.tpc' 'absent.tpc')", "bad.tm", in: dir)
        let before = set.snapshot()
        #expect(throws: SpiceKernelError.self) { try set.load(bad) }
        #expect(set.snapshot().cacheScope == before.cacheScope)
        #expect(set.snapshot().loads.count == 1)
        try set.unload(id)
        #expect(set.snapshot().variables.isEmpty)
        #expect(FileManager.default.fileExists(atPath: child.path))
    }
    @Test func malformedUnsupportedAndBoundedWork() throws {
        let dir = try directory(); defer { try? FileManager.default.removeItem(at: dir) }
        let badValues = ["A=(1 'two')", "A=(1", "A='unterminated", "A=nan", "A=inf", "A=@2000-FEB-30", "A=()", "A=(1(2))", "A=(1'bad')"]
        let set = try SpiceKernelSet()
        for (i, value) in badValues.enumerated() {
            let file = try write("KPL/PCK\n\\begindata\n" + value, "bad\(i).tpc", in: dir)
            #expect(throws: SpiceKernelError.self) { try set.load(file) }
        }
        let binary = try write("DAF/SPK ", "x.bsp", in: dir)
        #expect(throws: SpiceDAFError.self) { try set.load(binary) }
        let cycle = try write("KPL/MK\n\\begindata\nKERNELS_TO_LOAD='cycle.tm'", "cycle.tm", in: dir)
        #expect(throws: SpiceKernelError.self) { try set.load(cycle) }
        let symbols = try write("KPL/MK\n\\begindata\nKERNELS_TO_LOAD='$UNKNOWN/a'", "symbol.tm", in: dir)
        #expect(throws: SpiceKernelError.self) { try set.load(symbols) }
        let file = try write("KPL/PCK\n\\begindata\nA=(1 2 3)", "limit.tpc", in: dir)
        let bounded = try SpiceKernelSet(limits: SpiceKernelLimits(maximumValues: 2))
        #expect(throws: SpiceKernelError.self) { try bounded.load(file) }
        let byteBounded = try SpiceKernelSet(limits: SpiceKernelLimits(maximumFileBytes: 10))
        #expect(throws: SpiceKernelError.self) { try byteBounded.load(file) }
        #expect(set.snapshot().variables.isEmpty)
    }
    @Test func setCapacityTypeChangesAndClear() throws {
        let dir = try directory(); defer { try? FileManager.default.removeItem(at: dir) }
        let file = try write("KPL/PCK\r\n\\begindata\r\nA=.5d+2\r\nA += 'changed'\r\nA += 'again'", "types.tpc", in: dir)
        let set = try SpiceKernelSet(limits: SpiceKernelLimits(maximumFiles: 1))
        _ = try set.load(file)
        #expect(set.snapshot()["A"]?.value == .strings(["changed", "again"]))
        let old = set.snapshot()
        #expect(throws: SpiceKernelError.self) { try set.load(file) }
        #expect(set.snapshot().cacheScope == old.cacheScope)
        set.removeAll()
        #expect(set.snapshot().cacheScope != old.cacheScope)
        #expect(old["A"]?.value == .strings(["changed", "again"]))
        #expect(set.snapshot().loads.isEmpty)
        let dangling = try write("KPL/MK\n\\begindata\nKERNELS_TO_LOAD=('+')", "dangling.tm", in: dir)
        #expect(throws: SpiceKernelError.self) { try set.load(dangling) }
    }

}
