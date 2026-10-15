// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
@testable import CelMek

@Suite struct SpiceSPKResolutionTests {
    struct Link {
        let target: Int32
        let centre: Int32
        var p = SIMD3<Double>.zero
        var v = SIMD3<Double>.zero
        var frame: Int32 = 1
        var type: Int32 = 3
        var initial = -20.0
        var corrupt = false
    }
    func withKernels(_ links: [Link], _ action: (SpiceKernelSet,[URL]) throws -> Void) throws {
        let directory = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: directory,withIntermediateDirectories: true)
        defer { try? FileManager.default.removeItem(at: directory) }
        let set = try SpiceKernelSet()
        var urls: [URL] = []
        for (index,link) in links.enumerated() {
            var fixture = SPKFixture(type: 3,initial: link.initial,coefficients: [[link.p.x],[link.p.y],[link.p.z],[link.v.x],[link.v.y],[link.v.z]],records: 1,target: link.target)
            fixture.daf.int(link.centre,2048+24+20); fixture.daf.int(link.frame,2048+24+24)
            fixture.daf.int(link.type,2048+24+28)
            if link.corrupt { fixture.daf.double(.nan,fixture.directoryOffset) }
            let url = directory.appendingPathComponent("\(index).bsp")
            try fixture.daf.bytes.write(to: url); _ = try set.load(url); urls.append(url)
        }
        try action(set,urls)
    }
    func near(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> Bool { (0..<3).allSatisfy { abs(a[$0]-b[$0]) < 1e-12 } }
    var solarLinks: [Link] {
        [Link(target: 3,centre: 0,p: SIMD3(1000,2000,3000),v: SIMD3(10,20,30)),
         Link(target: 399,centre: 3,p: SIMD3(-4,5,6),v: SIMD3(0.1,0.2,0.3)),
         Link(target: 301,centre: 3,p: SIMD3(380000,-10,11),v: SIMD3(1,2,3)),
         Link(target: -123,centre: 301,p: SIMD3(1,2,3),v: SIMD3(4,5,6))]
    }
    @Test func commonCentreForwardReverseAndBodyBarycentreIdentity() throws {
        try withKernels(solarLinks) { set,_ in
            let view = set.snapshot()
            let moon = try view.spkState(targetID: 301,centreID: 399,atSecondsTDB: -15)
            #expect(moon.positionKilometres == SIMD3(380004,-15,5))
            #expect(near(moon.velocityKilometresPerTDBSecond,SIMD3(0.9,1.8,2.7)))
            #expect(moon.commonCentreID == 3)
            #expect(moon.contributions.map { $0.state.descriptor.targetID } == [301,399])
            #expect(moon.contributions.map(\.side) == [.target,.centre])
            #expect(moon.cacheScope == view.cacheScope)
            let reverse = try view.spkState(targetID: 399,centreID: 301,atSecondsTDB: -15)
            #expect(reverse.positionKilometres == -moon.positionKilometres)
            #expect(near(reverse.velocityKilometresPerTDBSecond,-moon.velocityKilometresPerTDBSecond))
            let craft = try view.spkState(targetID: -123,centreID: 399,atSecondsTDB: -15)
            #expect(craft.positionKilometres == SIMD3(380005,-13,8))
            #expect(near(craft.velocityKilometresPerTDBSecond,SIMD3(4.9,6.8,8.7)))
            #expect(craft.contributions.map { $0.state.descriptor.targetID } == [-123,301,399])
            let earth = try view.spkState(targetID: 399,centreID: 0,atSecondsTDB: -15)
            #expect(earth.positionKilometres == SIMD3(996,2005,3006))
            #expect(near(earth.velocityKilometresPerTDBSecond,SIMD3(10.1,20.2,30.3)))
            let emb = try view.spkState(targetID: 3,centreID: 0,atSecondsTDB: -15)
            #expect(earth.positionKilometres != emb.positionKilometres)
            let origin = try view.spkState(targetID: 0,centreID: 399,atSecondsTDB: -15)
            #expect(origin.positionKilometres == -earth.positionKilometres)
        }
    }
    @Test func localQueriesDoNotRequireAbsoluteOrUnusedAncestorData() throws {
        var links = Array(solarLinks.dropFirst())
        links.append(Link(target: 3,centre: 0,type: 99))
        try withKernels(links) { set,_ in
            let view = set.snapshot()
            let moon = try view.spkState(targetID: 301,centreID: 399,atSecondsTDB: -15)
            #expect(moon.positionKilometres == SIMD3(380004,-15,5))
            #expect(moon.contributions.count == 2)
            #expect(throws: SpiceSPKError.unsupportedType(99)) { try view.spkState(targetID: 399,centreID: 0,atSecondsTDB: -15) }
            let verified1 = try view.spkState(targetID: -123,centreID: 301,atSecondsTDB: -15).contributions.count == 1
            #expect(verified1)
        }
        try withKernels(Array(solarLinks.dropFirst())) { set,_ in
            let verified2 = try set.snapshot().spkState(targetID: 301,centreID: 399,atSecondsTDB: -15).commonCentreID == 3
            #expect(verified2)
            #expect(throws: SpiceSPKResolutionError.missingLink(bodyID: 3,secondsTDB: -15)) {
                try set.snapshot().spkState(targetID: 399,centreID: 0,atSecondsTDB: -15)
            }
        }
        // Reverse a direct link even when the desired target has no own segment.
        try withKernels([Link(target: 301,centre: 399,p: SIMD3(1,2,3))]) { set,_ in
            let verified3 = try set.snapshot().spkState(targetID: 399,centreID: 301,atSecondsTDB: -15).positionKilometres == SIMD3(-1,-2,-3)
            #expect(verified3)
        }
    }
    @Test func precedenceChangesCentresAndSnapshotOwnership() throws {
        try withKernels(solarLinks + [Link(target: 399,centre: 0,p: SIMD3(900,1900,2900))]) { set,urls in
            let old = set.snapshot(), query = try old.spkState(targetID: 301,centreID: 399,atSecondsTDB: -15)
            #expect(query.commonCentreID == 0)
            #expect(query.positionKilometres == SIMD3(380100,90,111))
            #expect(query.contributions.map(\.kernelIndex) == [2,0,4])
            #expect(query.contributions.last?.state.sourcePath == resolvedSpiceTestPath(urls[4]))
            let last = old.loads.last!.id
            try set.unload(last)
            let verified4 = try set.snapshot().spkState(targetID: 301,centreID: 399,atSecondsTDB: -15).commonCentreID == 3
            #expect(verified4)
            let verified5 = try old.spkState(targetID: 301,centreID: 399,atSecondsTDB: -15).positionKilometres == query.positionKilometres
            #expect(verified5)
            #expect(set.snapshot().cacheScope != old.cacheScope)
        }
        try withKernels(solarLinks + [Link(target: 399,centre: 0,type: 99)]) { set,_ in
            #expect(throws: SpiceSPKError.unsupportedType(99)) { try set.snapshot().spkState(targetID: 301,centreID: 399,atSecondsTDB: -15) }
        }
        try withKernels(solarLinks + [Link(target: 399,centre: 0,corrupt: true)]) { set,_ in
            #expect(throws: SpiceSPKError.self) { try set.snapshot().spkState(targetID: 301,centreID: 399,atSecondsTDB: -15) }
        }
        try withKernels(solarLinks + [Link(target: 399,centre: 0,type: 99,initial: -40)]) { set,_ in
            let verified6 = try set.snapshot().spkState(targetID: 301,centreID: 399,atSecondsTDB: -15).commonCentreID == 3
            #expect(verified6)
        }
    }
    @Test func missingLinksGapsAndCyclesHaveDistinctFailures() throws {
        try withKernels([Link(target: 399,centre: 9)]) { set,_ in
            #expect(throws: SpiceSPKResolutionError.missingLink(bodyID: 9,secondsTDB: -15)) { try set.snapshot().spkState(targetID: 399,centreID: 0,atSecondsTDB: -15) }
        }
        try withKernels([Link(target: 399,centre: 3),Link(target: 3,centre: 0,initial: -40)]) { set,_ in
            #expect(throws: SpiceSPKResolutionError.coverageGap(bodyID: 3,secondsTDB: -15)) { try set.snapshot().spkState(targetID: 399,centreID: 0,atSecondsTDB: -15) }
        }
        try withKernels([Link(target: 399,centre: 3),Link(target: 3,centre: 399)]) { set,_ in
            #expect(throws: SpiceSPKResolutionError.cycle([399,3,399])) { try set.snapshot().spkState(targetID: 399,centreID: 0,atSecondsTDB: -15) }
            // A direct query does not require or evaluate the cyclic ancestor above its common centre.
            let verified7 = try set.snapshot().spkState(targetID: 399,centreID: 3,atSecondsTDB: -15).contributions.count == 1
            #expect(verified7)
        }
        try withKernels([Link(target: 399,centre: 399)]) { set,_ in
            #expect(throws: SpiceSPKResolutionError.cycle([399,399])) { try set.snapshot().spkState(targetID: 399,centreID: 0,atSecondsTDB: -15) }
        }
    }
    @Test func frameContractsWorkLimitsAndIdentityQueries() throws {
        try withKernels([Link(target: 301,centre: 399,p: SIMD3(1,2,3),frame: 17)]) { set,_ in
            let view = set.snapshot()
            #expect(throws: SpiceSPKResolutionError.frameMismatch(bodyID: 301,expected: 1,actual: 17)) { try view.spkState(targetID: 301,centreID: 399,atSecondsTDB: -15) }
            let verified8 = try view.spkState(targetID: 301,centreID: 399,frameID: 17,atSecondsTDB: -15).positionKilometres == SIMD3(1,2,3)
            #expect(verified8)
        }
        try withKernels(solarLinks) { set,_ in
            let view = set.snapshot()
            #expect(throws: SpiceSPKResolutionError.chainLimit(1)) { try view.spkState(targetID: -123,centreID: 0,atSecondsTDB: -15,options: SpiceSPKResolutionOptions(maximumChainLinks: 1)) }
            #expect(throws: SpiceSPKResolutionError.descriptorLimit(1)) { try view.spkState(targetID: 399,centreID: 0,atSecondsTDB: -15,options: SpiceSPKResolutionOptions(maximumDescriptorChecks: 1)) }
            #expect(throws: SpiceSPKResolutionError.invalidFrame) { try view.spkState(targetID: 301,centreID: 301,frameID: 0,atSecondsTDB: -15) }
            #expect(throws: SpiceSPKError.invalidEpoch) { try view.spkState(targetID: 301,centreID: 301,atSecondsTDB: .nan) }
        }
        let empty = try SpiceKernelSet().snapshot()
        let same = try empty.spkState(targetID: -123,centreID: -123,atSecondsTDB: 10)
        #expect(same.positionKilometres == .zero && same.velocityKilometresPerTDBSecond == .zero)
        #expect(same.contributions.isEmpty && same.exploredChainLinks == 0 && same.inspectedDescriptors == 0)
        #expect(throws: SpiceSPKResolutionError.invalidOptions) { try SpiceSPKResolutionOptions(maximumChainLinks: 0) }
        #expect(throws: SpiceSPKResolutionError.invalidOptions) { try SpiceSPKResolutionOptions(maximumDescriptorChecks: -1) }
    }
    @Test func compensatedCompositionAndOverflowRejection() throws {
        let large = 1e16
        try withKernels([Link(target: -1,centre: 3,p: SIMD3(large,0,0),v: SIMD3(large,0,0)),
                         Link(target: 3,centre: 0,p: SIMD3(1,0,0),v: SIMD3(1,0,0)),
                         Link(target: -2,centre: 0,p: SIMD3(large,0,0),v: SIMD3(large,0,0))]) { set,_ in
            let s = try set.snapshot().spkState(targetID: -1,centreID: -2,atSecondsTDB: -15)
            #expect(s.positionKilometres.x == 1 && s.velocityKilometresPerTDBSecond.x == 1)
        }
        let huge = Double.greatestFiniteMagnitude*0.75
        try withKernels([Link(target: -1,centre: 3,p: SIMD3(huge,0,0)),Link(target: 3,centre: 0,p: SIMD3(huge,0,0))]) { set,_ in
            #expect(throws: SpiceSPKError.numericalRangeExceeded) { try set.snapshot().spkState(targetID: -1,centreID: 0,atSecondsTDB: -15) }
        }
    }
}
