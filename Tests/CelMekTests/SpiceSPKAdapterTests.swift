// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct SpiceSPKAdapterTests {
    let fixtures = SpiceSPKResolutionTests()
    var provenance: PhysicalDataProvenance { try! PhysicalDataProvenance(source: "Synthetic public-format SPK",version: "test-1") }
    func epoch(_ seconds: Double) throws -> Epoch { Epoch(offsetFromJ2000: try CelMek.Duration(seconds: seconds)) }
    func range(_ start: Double = -18, _ end: Double = -12) throws -> ClosedRange<Epoch> { try epoch(start)...epoch(end) }
    func bindings(_ frame: Frame = CelestialFrames.icrf, id: Int32 = 1) throws -> SpiceFrameBindings {
        try SpiceFrameBindings([SpiceFrameBinding(frameID: id,axesFrame: frame,provenance: provenance)])
    }
    struct Clock: TDBTimeModel {
        let coverage: ClosedRange<Epoch>
        var offset = 0.0
        var drift = 0.0
        var quadratic = 0.0
        var stale = false
        let provenance = "Test clock, specified TT offset and derivatives"
        func correction(at epoch: Epoch) throws -> TDBCorrection {
            let t = epoch.offsetFromJ2000.seconds
            return try TDBCorrection(epoch: stale ? .j2000 : epoch,
                tdbMinusTT: CelMek.Duration(seconds: offset+drift*t+quadratic*t*t),
                rate: 1+drift+2*quadratic*t,acceleration: 2*quadratic)
        }
    }
    @Test func type2SIAtomicRateAndEpochPrecision() throws {
        let f = SPKFixture()
        try SpiceDAFTests().withFile(f.daf.bytes) { url in
            let set = try SpiceKernelSet(); _ = try set.load(url)
            let coverage = try range(), clock = Clock(coverage: try range(-100,100),offset: 0.125,drift: 0.001,quadratic: 0.00001)
            let provider = try SpiceSPKStateProvider(snapshot: set.snapshot(),body: .earth,
                frameBindings: bindings(),coverage: coverage,tdbModel: clock,dataProvenance: provenance,accuracy: .low)
            let t = try epoch(-15), sample = try provider.evaluation(at: t)
            let et = -15+0.125+0.001*(-15)+0.00001*225
            let reference = SpiceSPKTests().reference((et+15)/5,record: 0,type: 2,radius: 5)
            #expect(fixtures.near(sample.state.position,reference.0*1000))
            #expect(fixtures.near(sample.state.velocity,reference.1*1000*sample.tdbCorrection.rate))
            #expect(sample.state.epoch == t && sample.state.frame === CelestialFrames.icrf)
            #expect(sample.state.acceleration == nil)
            #expect(sample.chain.secondsTDB == et)
            #expect(provider.accuracy == .low && provider.effectiveAccuracy == .full)
            #expect(provider.metadata.modelTimeScale == .tdb)
            #expect(provider.metadata.accelerationMethod == .unavailable)
            #expect(provider.metadata.origin == .solarSystemBarycentre && provider.metadata.body == .earth)
            #expect(provider.metadata.provenance.source.contains(clock.provenance))
            let h = 0.001
            let plus = try provider.checkedState(at: epoch(-15+h)), minus = try provider.checkedState(at: epoch(-15-h))
            let derivative = (plus.position-minus.position)/(2*h)
            #expect((0..<3).allSatisfy { abs(derivative[$0]-sample.state.velocity[$0]) < 1e-4 })
            let collection = try EphemerisCollection(providers: [provider])
            let collected = try collection.state(of: .earth,at: t)
            #expect(collected.position == sample.state.position)
        }
        let large = 1e9, far = SPKFixture(initial: large,interval: 1)
        try SpiceDAFTests().withFile(far.daf.bytes) { url in
            let set = try SpiceKernelSet(); _ = try set.load(url)
            let coverage = try range(large,large+1), clock = Clock(coverage: coverage,offset: 0.02)
            let provider = try SpiceSPKStateProvider(snapshot: set.snapshot(),body: .earth,
                frameBindings: bindings(),coverage: coverage,tdbModel: clock,dataProvenance: provenance)
            let at = Epoch(offsetFromJ2000: try CelMek.Duration(wholeSeconds: 1_000_000_000,fractionalSecond: 0.123456))
            let result = try provider.evaluation(at: at)
            #expect(result.chain.secondsTDB == large+0.143456)
            #expect(result.state.epoch == at)
        }
    }
    @Test func type3PreservesSuppliedVelocitiesAndSnapshotAfterUnload() throws {
        try fixtures.withKernels([.init(target: -123,centre: 0,p: SIMD3(1,2,3),v: SIMD3(4,5,6))]) { set,_ in
            let view = set.snapshot(), coverage = try range(), clock = Clock(coverage: try range(-100,100),drift: 0.1)
            let provider = try SpiceSPKStateProvider(snapshot: view,body: .naif(-123),
                frameBindings: bindings(),coverage: coverage,tdbModel: clock,dataProvenance: provenance)
            try set.unload(view.loads[0].id)
            let s = try provider.checkedState(at: epoch(-15))
            #expect(s.position == SIMD3(1000,2000,3000))
            #expect(fixtures.near(s.velocity,SIMD3(4400,5500,6600)))
            #expect(s.acceleration == nil)
            #expect(set.snapshot().binaryKernels.isEmpty)
            #expect(provider.metadata.provenance.version == provenance.version)
        }
    }
    @Test func mixedNativeAxesUseFrameSystemBeforeComposition() throws {
        let root = Frame(name: "SPK test inertial",pseudoInertial: true), tree = try FrameTree(root: root)
        let native = try tree.makeFrame(parent: root,name: "Native axes",pseudoInertial: true,
            provider: FixedTransformProvider(rotation: Quaternion(angle: -.pi/2,axis: SIMD3(0,0,1))))
        let axes = try SpiceFrameBindings([SpiceFrameBinding(frameID: 1,axesFrame: root,provenance: provenance),
                                          SpiceFrameBinding(frameID: 17,axesFrame: native,provenance: provenance)])
        let links = [SpiceSPKResolutionTests.Link(target: 301,centre: 0,p: SIMD3(5,6,7),v: SIMD3(0.5,0.6,0.7)),
                     .init(target: 399,centre: 0,p: SIMD3(1,2,3),v: SIMD3(0.1,0.2,0.3),frame: 17)]
        try fixtures.withKernels(links) { set,_ in
            let coverage = try range(), clock = Clock(coverage: coverage)
            // root is explicitly asserted to represent Earth-centred output axes in this synthetic test.
            let provider = try SpiceSPKStateProvider(snapshot: set.snapshot(),body: .moon,origin: .earth,frame: root,
                frameBindings: axes,frameTree: tree,coverage: coverage,tdbModel: clock,dataProvenance: provenance)
            let sample = try provider.evaluation(at: epoch(-15))
            #expect(fixtures.near(sample.state.position,SIMD3(7000,5000,4000)))
            #expect(fixtures.near(sample.state.velocity,SIMD3(700,500,400)))
            #expect(sample.chain.contributions.map { $0.state.descriptor.frameID } == [1,17])
            #expect(sample.state.frame === root)
        }
    }
    @Test func rotatingAxesIncludeRateAndOriginTranslationIsNotDoubleCounted() throws {
        let root = Frame(name: "Rotating test origin",pseudoInertial: true), tree = try FrameTree(root: root)
        let rotating = try tree.makeFrame(parent: root,name: "Rotating native",
            provider: UniformRotationProvider(referenceEpoch: .j2000,axis: SIMD3(0,0,1),rate: 0.2))
        try fixtures.withKernels([.init(target: 399,centre: 0,p: SIMD3(1,0,0),frame: 17)]) { set,_ in
            let coverage = try range(), clock = Clock(coverage: coverage)
            let provider = try SpiceSPKStateProvider(snapshot: set.snapshot(),body: .earth,frame: root,
                frameBindings: bindings(rotating,id: 17),frameTree: tree,coverage: coverage,tdbModel: clock,dataProvenance: provenance)
            let state = try provider.checkedState(at: epoch(-15)), angle = -3.0
            #expect(fixtures.near(state.position,SIMD3(1000*cos(angle),1000*sin(angle),0)))
            #expect(fixtures.near(state.velocity,SIMD3(-200*sin(angle),200*cos(angle),0)))
            let h = 0.001, plus = try provider.state(at: epoch(-15+h)), minus = try provider.state(at: epoch(-15-h))
            #expect((0..<3).allSatisfy { abs((plus.position[$0]-minus.position[$0])/(2*h)-state.velocity[$0]) < 1e-5 })
        }
        try fixtures.withKernels([.init(target: 399,centre: 0,p: SIMD3(10,20,30),v: SIMD3(1,2,3)),
                                 .init(target: 301,centre: 0,p: SIMD3(12,23,34),v: SIMD3(1.5,2.3,3.1))]) { set,_ in
            let coverage = try range(), clock = Clock(coverage: coverage), axes = try bindings()
            let earth = try SpiceSPKStateProvider(snapshot: set.snapshot(),body: .earth,
                frameBindings: axes,coverage: coverage,tdbModel: clock,dataProvenance: provenance)
            let celestialTree = try CelestialFrames.makeTree()
            let gcrf = try celestialTree.makeGCRF(earth: earth,accelerationPolicy: .assumeZero)
            let moon = try SpiceSPKStateProvider(snapshot: set.snapshot(),body: .moon,origin: .earth,frame: gcrf,
                frameBindings: axes,frameTree: celestialTree,coverage: coverage,tdbModel: clock,dataProvenance: provenance)
            let relative = try moon.checkedState(at: epoch(-15))
            #expect(relative.position == SIMD3(2000,3000,4000))
            #expect(fixtures.near(relative.velocity,SIMD3(500,300,100)))
            #expect(relative.frame === gcrf && moon.metadata.origin == .earth)
        }
    }
    @Test func bindingsCoverageTimeContractsAndWeakTreeLifetime() throws {
        let coverage = try range(), clock = Clock(coverage: coverage)
        #expect(throws: SpiceSPKAdapterError.duplicateFrameBinding(1)) { try SpiceFrameBindings([SpiceFrameBinding.icrfJ2000(provenance: provenance),SpiceFrameBinding.icrfJ2000(provenance: provenance)]) }
        #expect(throws: SpiceSPKAdapterError.invalidFrameBinding) { try SpiceFrameBinding(frameID: 0,axesFrame: CelestialFrames.icrf,provenance: provenance) }
        try fixtures.withKernels([.init(target: 399,centre: 0,p: SIMD3(1,2,3))]) { set,_ in
            let view = set.snapshot()
            #expect(throws: SpiceSPKAdapterError.originFrameRequired) { try SpiceSPKStateProvider(snapshot: view,body: .moon,origin: .earth,frameBindings: bindings(),coverage: coverage,tdbModel: clock,dataProvenance: provenance) }
            #expect(throws: SpiceSPKAdapterError.unsupportedBodyIdentifier) { try SpiceSPKStateProvider(snapshot: view,body: .custom(namespace: "test",name: "craft"),frameBindings: bindings(),coverage: coverage,tdbModel: clock,dataProvenance: provenance) }
            #expect(throws: SpiceSPKAdapterError.unknownNativeFrame(1)) { try SpiceSPKStateProvider(snapshot: view,body: .earth,frameBindings: SpiceFrameBindings([]),coverage: coverage,tdbModel: clock,dataProvenance: provenance) }
            #expect(throws: TDBModelError.epochMismatch) { try SpiceSPKStateProvider(snapshot: view,body: .earth,frameBindings: bindings(),coverage: coverage,tdbModel: Clock(coverage: coverage,stale: true),dataProvenance: provenance) }
            #expect(throws: EphemerisError.outsideCoverage) { try SpiceSPKStateProvider(snapshot: view,body: .earth,frameBindings: bindings(),coverage: range(-19,-11),tdbModel: clock,dataProvenance: provenance) }
            let provider = try SpiceSPKStateProvider(snapshot: view,body: .earth,frameBindings: bindings(),coverage: coverage,tdbModel: clock,dataProvenance: provenance)
            #expect(throws: EphemerisError.outsideCoverage) { try provider.state(at: epoch(-19)) }
            var tree: FrameTree? = try CelestialFrames.makeTree()
            let destination = CelestialFrames.meanEquatorJ2000
            #expect(throws: SpiceSPKAdapterError.frameTreeRequired) { try SpiceSPKStateProvider(snapshot: view,body: .earth,frame: destination,frameBindings: bindings(),coverage: coverage,tdbModel: clock,dataProvenance: provenance) }
            let rotated = try SpiceSPKStateProvider(snapshot: view,body: .earth,frame: destination,
                frameBindings: bindings(),frameTree: tree!,coverage: coverage,tdbModel: clock,dataProvenance: provenance)
            tree = nil
            #expect(throws: SpiceSPKAdapterError.frameTreeUnavailable) { try rotated.state(at: epoch(-15)) }
        }
        try fixtures.withKernels([.init(target: 399,centre: 0,initial: -20),.init(target: 399,centre: 0,initial: 0)]) { set,_ in
            let broad = try range(-19,9)
            let provider = try SpiceSPKStateProvider(snapshot: set.snapshot(),body: .earth,
                frameBindings: bindings(),coverage: broad,tdbModel: Clock(coverage: broad),dataProvenance: provenance)
            #expect(throws: SpiceSPKResolutionError.coverageGap(bodyID: 399,secondsTDB: -5)) { try provider.state(at: epoch(-5)) }
        }
    }
}
