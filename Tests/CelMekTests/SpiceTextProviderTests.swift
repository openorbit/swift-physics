// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct SpiceTextProviderTests {
    let helper = SpiceTextKernelTests()
    var provenance: PhysicalDataProvenance { try! PhysicalDataProvenance(source: "Synthetic public-format data", version: "fixture-1") }
    func epoch(_ t: Double) throws -> Epoch { Epoch(offsetFromJ2000: try CelMek.Duration(seconds: t)) }
    var base: String { "KPL/PCK\n\\begindata\nBODY301_POLE_RA=(270 0 0)\nBODY301_POLE_DEC=(90 0 0)\nBODY301_PM=(0 10 0)\n" }
    func provider(_ text: String, clock: (any TDBTimeModel)? = nil, bindings: SpiceFrameBindings? = nil) throws -> SpicePCKOrientationProvider {
        let coverage = try clock?.coverage ?? (epoch(-1e10)...epoch(1e10))
        return try SpicePCKOrientationProvider(snapshot: helper.snapshot(text), bodyID: 301, systemID: 3,
            frameBindings: bindings ?? SpiceFrameBindings([.icrfJ2000(provenance: provenance)]), coverage: coverage,
            tdbModel: clock ?? SpiceSPKAdapterTests.Clock(coverage: coverage), provenance: provenance)
    }
    func near(_ a: SIMD3<Double>, _ b: SIMD3<Double>, _ tolerance: Double = 1e-10) -> Bool {
        (0..<3).allSatisfy { abs(a[$0]-b[$0]) < tolerance }
    }
    @Test func sourcedSIConstantsAndConventions() throws {
        let pool = try helper.snapshot("KPL/PCK\n\\begindata\nBODY301_RADII=(1 2 3)\nBODY301_GM=4.5")
        let radii = try pool.bodyReferenceRadii(bodyID: 301, axesDescription: "Declared ellipsoid axes", provenance: provenance)
        #expect(radii.metres == SIMD3(1000,2000,3000))
        let gm = try pool.bodyGravitationalParameter(bodyID: 301, convention: .tdbCompatible, provenance: provenance)
        #expect(gm.cubicMetresPerSecondSquared == 4.5e9 && gm.convention == .tdbCompatible)
        #expect(gm.provenance.version == "fixture-1" && gm.provenance.source.contains("kernel.txt"))
        #expect(radii.provenance.source.contains("kernel.txt"))
        #expect(throws: (any Error).self) { try helper.snapshot("KPL/PCK\n\\begindata\nBODY301_GM=-1").bodyGravitationalParameter(bodyID: 301, convention: .newtonian, provenance: provenance) }
    }
    @Test func spinDirectionAndAtomicDerivatives() throws {
        let coverage = try epoch(-10000)...epoch(10000)
        let clock = SpiceSPKAdapterTests.Clock(coverage: coverage, offset: 20, drift: 0.01, quadratic: 1e-8)
        let model = try provider(base, clock: clock), t = 1000.0, sample = try model.orientation(at: epoch(t))
        let et = t+20+0.01*t+1e-8*t*t, angle = -10 * et/86400 * .pi/180
        let v = sample.rotation.matrix * SIMD3<Double>(1,0,0)
        #expect(near(v, SIMD3(cos(angle),sin(angle),0)))
        #expect(abs(sample.angularVelocity.z - (-10 * .pi/180/86400)*(1.01+2e-8*t)) < 1e-16)
        #expect(abs(sample.angularAcceleration.z - (-10 * .pi/180/86400)*2e-8) < 1e-20)
        #expect(sample.referenceFrame === CelestialFrames.icrf)
        #expect(model.provenance.source.contains(clock.provenance))
    }
    @Test func periodicTermsHigherDegreeAndCoupledRates() throws {
        let text = base.replacingOccurrences(of: "(270 0 0)", with: "(30 .2 .01)").replacingOccurrences(of: "(90 0 0)", with: "(60 -.3 .02)") + "BODY3_MAX_PHASE_DEGREE=3\nBODY3_NUT_PREC_ANGLES=(10 30 2 .5 80 -20 1 .1)\nBODY301_NUT_PREC_RA=(2 3)\nBODY301_NUT_PREC_DEC=(4 5)\nBODY301_NUT_PREC_PM=(6 7)"
        let model = try provider(text), t = 2e8, c = t/(36525*86400)
        let p1 = (10+30*c+2*c*c+0.5*c*c*c)*Double.pi/180
        let p2 = (80-20*c+c*c+0.1*c*c*c)*Double.pi/180
        let ra = (30+0.2*c+0.01*c*c+2*sin(p1)+3*sin(p2))*Double.pi/180
        let dec = (60-0.3*c+0.02*c*c+4*cos(p1)+5*cos(p2))*Double.pi/180
        let w = (10*t/86400+6*sin(p1)+7*sin(p2))*Double.pi/180
        // Independent explicit basis rotations from public PCK Euler definitions.
        func z(_ a: Double, _ v: SIMD3<Double>) -> SIMD3<Double> { SIMD3(cos(a)*v.x-sin(a)*v.y,sin(a)*v.x+cos(a)*v.y,v.z) }
        func x(_ a: Double, _ v: SIMD3<Double>) -> SIMD3<Double> { SIMD3(v.x,cos(a)*v.y-sin(a)*v.z,sin(a)*v.y+cos(a)*v.z) }
        let sample = try model.orientation(at: epoch(t)), input = SIMD3<Double>(0.3,0.5,0.7)
        let r = sample.rotation.matrix * input
        #expect(near(r, z(-w,x(dec - .pi/2,z(-ra - .pi/2,input))), 1e-11))
        let h = 10.0, plus = try model.orientation(at: epoch(t+h)), minus = try model.orientation(at: epoch(t-h))
        let numerical = (plus.rotation.matrix * input - minus.rotation.matrix * input)/(2*h)
        let omega = sample.angularVelocity
        let expected = SIMD3(omega.y*r.z-omega.z*r.y,omega.z*r.x-omega.x*r.z,omega.x*r.y-omega.y*r.x)
        #expect(near(numerical,expected,1e-11))
        #expect(near((plus.angularVelocity-minus.angularVelocity)/(2*h),sample.angularAcceleration,1e-16))
    }
    @Test func alternateEpochAxesAndMalformedDependencies() throws {
        let tree = FrameTree(name: "explicit alternate inertial axes")
        let bindings = try SpiceFrameBindings([SpiceFrameBinding(frameID: 2, axesFrame: tree.root, provenance: provenance)])
        let text = base + "BODY3_CONSTS_REF_FRAME=2\nBODY3_CONSTS_JED_EPOCH=2451546"
        let model = try provider(text, bindings: bindings), sample = try model.orientation(at: epoch(86400))
        #expect(sample.referenceFrame === tree.root && model.referenceFrameID == 2)
        #expect(near(sample.rotation.matrix * SIMD3<Double>(1,0,0), SIMD3(1,0,0)))
        #expect(throws: SpiceTextProviderError.conflictingReferenceParameters) { try provider(text + "\nBODY3_CONSTANTS_REF_FRAME=1", bindings: bindings) }
        for suffix in ["BODY301_NUT_PREC_PM=(1)", "BODY3_MAX_PHASE_DEGREE=4", "BODY3_NUT_PREC_ANGLES=(1 2 3)", "BODY3_CONSTS_REF_FRAME=1.5"] {
            #expect(throws: (any Error).self) { try provider(base+suffix) }
        }
        #expect(throws: SpiceTextProviderError.outsideCoverage) { try model.orientation(at: epoch(1e10+1)) }
    }
    struct Centre: StateProvider {
        let frame: Frame
        func state(at epoch: Epoch) throws -> CartesianState {
            try CartesianState(epoch: epoch, frame: frame, position: SIMD3(100+2*epoch.offsetFromJ2000.seconds,0,0), velocity: SIMD3(2,0,0), acceleration: .zero)
        }
    }
    @Test func fixedFramesInvertLocalMatrixAndTranslateCentre() throws {
        let prefix = "KPL/FK\n\\begindata\nFRAME_123_NAME='TEST'\nFRAME_123_CLASS=4\nFRAME_123_CLASS_ID=123\nFRAME_123_CENTER=301\nTKFRAME_123_RELATIVE='J2000'\nTKFRAME_123_SPEC='MATRIX'\nTKFRAME_123_MATRIX="
        for matrix in ["(1 0 0 0 1 0 0 0 1)", "(1 0 0 0 -1 0 0 0 -1)", "(-1 0 0 0 1 0 0 0 -1)", "(-1 0 0 0 -1 0 0 0 1)", "(0 1 0 -1 0 0 0 0 1)"] {
            let definition = try helper.snapshot(prefix+matrix).fixedFrameDefinition(frameID: 123)
            let tree = FrameTree(), centre = Centre(frame: tree.root)
            let frame = try definition.makeFrame(in: tree, relativeAxes: tree.root, centre: centre, provenance: provenance)
            let state = try CartesianState(epoch: .j2000, frame: tree.root, position: SIMD3(101,2,3), velocity: SIMD3(3,2,3), acceleration: .zero)
            let actual = try tree.transform(from: tree.root,to: frame,at: .j2000).applying(to: state,destination: frame)
            let m = definition.matrix, input = SIMD3<Double>(1,2,3)
            func dotColumn(_ column: Int) -> Double {
                var result = 0.0
                for row in 0..<3 { result += m[column*3+row]*input[row] }
                return result
            }
            let expected = SIMD3<Double>(dotColumn(0),dotColumn(1),dotColumn(2))
            #expect(near(actual.position,expected))
            #expect(near(actual.velocity,expected))
            #expect(actual.acceleration == .zero)
        }
    }
    @Test func eulerUnitsOrderAndScalarFirstQuaternion() throws {
        let prefix = "KPL/FK\n\\begindata\nFRAME_123_NAME='TEST'\nFRAME_123_CLASS=4\nFRAME_123_CLASS_ID=123\nFRAME_123_CENTER=301\nTKFRAME_123_RELATIVE='J2000'\n"
        let units: [(String, Double)] = [("DEGREES",90),("RADIANS",.pi/2),("ARCSECONDS",324000),
            ("ARCMINUTES",5400),("HOURANGLE",6),("MINUTEANGLE",360),("SECONDANGLE",21600)]
        for (unit, angle) in units {
            let d = try helper.snapshot(prefix + "TKFRAME_123_SPEC='ANGLES'\nTKFRAME_123_ANGLES=(\(angle) 0 0)\nTKFRAME_123_AXES=(3 2 1)\nTKFRAME_123_UNITS='\(unit)'").fixedFrameDefinition(frameID: 123)
            let model = try SpiceFixedOrientationProvider(definition: d, relativeAxes: CelestialFrames.icrf, provenance: provenance)
            #expect(near(try model.orientation(at: .j2000).rotation.matrix * SIMD3<Double>(1,0,0), SIMD3(0,1,0)))
        }
        let a = try helper.snapshot(prefix + "TKFRAME_123_SPEC='ANGLES'\nTKFRAME_123_ANGLES=(90 90 90)\nTKFRAME_123_AXES=(3 2 1)\nTKFRAME_123_UNITS='DEGREES'").fixedFrameDefinition(frameID: 123)
        // M=[90]z[90]y[90]x: local +X maps to relative +Z.
        #expect(near(SIMD3(a.matrix[0],a.matrix[1],a.matrix[2]), SIMD3(0,0,1)))
        let q = try helper.snapshot(prefix + "TKFRAME_123_SPEC='QUATERNION'\nTKFRAME_123_Q=(0.7071067811865476 0 0 0.7071067811865476)").fixedFrameDefinition(frameID: 123)
        #expect(near(SIMD3(q.matrix[0],q.matrix[1],q.matrix[2]), SIMD3(0,1,0)))
        for tail in ["TKFRAME_123_SPEC='QUATERNION'\nTKFRAME_123_Q=(2 0 0 0)",
                     "TKFRAME_123_SPEC='ANGLES'\nTKFRAME_123_ANGLES=(1 2 3)\nTKFRAME_123_AXES=(1 4 3)\nTKFRAME_123_UNITS='DEGREES'"] {
            #expect(throws: SpiceTextKernelError.self) { try helper.snapshot(prefix+tail).fixedFrameDefinition(frameID: 123) }
        }
    }
    @Test func retainedProviderSurvivesUnloadAndBodyFixedComposition() throws {
        let fixture = SpiceKernelSetTests(), dir = try fixture.directory()
        defer { try? FileManager.default.removeItem(at: dir) }
        let file = try fixture.write(base, "orientation.tpc", in: dir)
        let set = try SpiceKernelSet(), id = try set.load(file), snapshot = set.snapshot()
        let tree = FrameTree(), coverage = try epoch(-100)...epoch(100)
        let bindings = try SpiceFrameBindings([SpiceFrameBinding(frameID: 1, axesFrame: tree.root, provenance: provenance)])
        let model = try SpicePCKOrientationProvider(snapshot: snapshot, bodyID: 301, systemID: 3, frameBindings: bindings,
            coverage: coverage, tdbModel: SpiceSPKAdapterTests.Clock(coverage: coverage), provenance: provenance)
        let before = try model.orientation(at: epoch(10))
        try set.unload(id)
        #expect(model.cacheScope == snapshot.cacheScope && model.cacheScope != set.snapshot().cacheScope)
        #expect(try model.orientation(at: epoch(10)).rotation.matrix * SIMD3<Double>(1,0,0) == before.rotation.matrix * SIMD3<Double>(1,0,0))
        let frame = try tree.makeBodyFixedFrame(name: "loaded lunar axes", centre: Centre(frame: tree.root), orientation: model)
        let state = try Centre(frame: tree.root).state(at: epoch(10))
        let transformed = try tree.transform(from: tree.root, to: frame, at: epoch(10)).applying(to: state, destination: frame)
        #expect(near(transformed.position,.zero))
        #expect(near(transformed.velocity,.zero))
        #expect(near(transformed.acceleration!,.zero))
    }

}
