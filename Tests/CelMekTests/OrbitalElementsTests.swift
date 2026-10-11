import Foundation
import Testing
@testable import CelMek

@Suite struct OrbitalElementsTests {
    private func gm(_ value: Double = 3.986004418e14) throws -> GravitationalParameter {
        try GravitationalParameter(cubicMetresPerSecondSquared: value, convention: .newtonian,
            provenance: PhysicalDataProvenance(source: "Test model", version: "1", reference: "analytic fixture"))
    }
    private func norm(_ x: SIMD3<Double>) -> Double { hypot(hypot(x.x, x.y), x.z) }

    @Test func analyticallyKnownPolarState() throws {
        // p=6 Mm, e=0.5, i=Omega=omega=nu=pi/2. Perifocal axes
        // are +Z and -Y. mu=p gives a unit velocity scale.
        let tree = FrameTree()
        let state = try CartesianState(epoch: .j2000, frame: tree.root,
            position: SIMD3(0, -6e6, 0), velocity: SIMD3(0, -0.5, -1))
        let elements = try OsculatingOrbitalElements(state: state, centralBody: .earth,
                                                      gravitationalParameter: gm(6e6))
        #expect(abs(elements.semiLatusRectumMetres - 6e6) < 1e-8)
        #expect(abs(elements.eccentricity - 0.5) < 1e-14)
        for angle in [elements.inclination, elements.longitudeOfAscendingNode,
                      elements.argumentOfPeriapsis, elements.trueAnomaly] {
            #expect(abs(angle - .pi/2) < 1e-14)
        }
        let rebuilt = try elements.cartesianState()
        #expect(norm(rebuilt.position - state.position) < 1e-8)
        #expect(norm(rebuilt.velocity - state.velocity) < 1e-14)
        #expect(rebuilt.acceleration == nil)
        #expect(rebuilt.epoch == state.epoch)
        #expect(rebuilt.frame === state.frame)
    }

    @Test func conicsAndSingularitiesRoundTrip() throws {
        let tree = FrameTree(), parameter = try gm()
        for e in [0.0, 1e-8, 0.3, 0.999999999, 1, 1.000000001, 2.0] {
            for i in [0.0, 1e-8, 0.7, .pi/2, .pi-1e-8, .pi] {
                for nu in [-0.8, 0.0, 0.8] {
                    let original = try OsculatingOrbitalElements(epoch: .j2000, frame: tree.root,
                        centralBody: .earth, gravitationalParameter: parameter,
                        semiLatusRectumMetres: 7e6, eccentricity: e, inclination: i,
                        longitudeOfAscendingNode: 1.2, argumentOfPeriapsis: 0.6, trueAnomaly: nu)
                    let state = try original.cartesianState()
                    let extracted = try OsculatingOrbitalElements(state: state, centralBody: .earth,
                                                                  gravitationalParameter: parameter)
                    let rebuilt = try extracted.cartesianState()
                    #expect(norm(rebuilt.position - state.position) < 1e-6)
                    #expect(norm(rebuilt.velocity - state.velocity) < 1e-8)
                    #expect(abs(extracted.eccentricity - e) < 2e-14)
                    #expect(extracted.centralBody == .earth)
                    #expect(extracted.gravitationalParameter == parameter)
                    if e == 0 { #expect(extracted.argumentOfPeriapsis == 0) }
                    if i == 0 || i == .pi { #expect(extracted.longitudeOfAscendingNode == 0) }
                }
            }
        }
    }

    @Test func circularEquatorialCartesianFixtures() throws {
        let tree = FrameTree(), mu = try gm(1)
        for sign in [-1.0, 1.0] {
            let state = try CartesianState(epoch: .j2000, frame: tree.root,
                position: SIMD3(0, 1, 0), velocity: SIMD3(-sign, 0, 0))
            let elements = try OsculatingOrbitalElements(state: state, centralBody: .earth,
                                                          gravitationalParameter: mu)
            #expect(elements.angleConvention == .circularEquatorial)
            #expect(elements.inclination == (sign > 0 ? 0 : .pi))
            #expect(abs(elements.trueAnomaly - sign * .pi/2) < 1e-14)
            #expect(elements.semiMajorAxisMetres == 1)
            #expect(norm(try elements.cartesianState().velocity - state.velocity) < 1e-14)
        }
        let parabola = try OsculatingOrbitalElements(epoch: .j2000, frame: tree.root, centralBody: .earth,
            gravitationalParameter: mu, semiLatusRectumMetres: 2, eccentricity: 1, inclination: 0,
            longitudeOfAscendingNode: 0, argumentOfPeriapsis: 0, trueAnomaly: 0)
        #expect(parabola.conic == .parabolic)
        #expect(parabola.semiMajorAxisMetres == nil)
        #expect(parabola.periapsisDistanceMetres == 1)
    }

    @Test func explicitSingularityThresholds() throws {
        let tree = FrameTree(), mu = try gm()
        let original = try OsculatingOrbitalElements(epoch: .j2000, frame: tree.root,
            centralBody: .earth, gravitationalParameter: mu, semiLatusRectumMetres: 7e6,
            eccentricity: 1e-7, inclination: 1e-7, longitudeOfAscendingNode: 1.0,
            argumentOfPeriapsis: 0.8, trueAnomaly: 0.3)
        let state = try original.cartesianState()
        let precise = try OsculatingOrbitalElements(state: state, centralBody: .earth,
            gravitationalParameter: mu,
            options: OrbitalElementConversionOptions(circularEccentricityTolerance: 0,
                                                       equatorialSineTolerance: 0))
        #expect(precise.angleConvention == .classical)
        #expect(norm(try precise.cartesianState().position - state.position) < 1e-6)
        let snapped = try OsculatingOrbitalElements(state: state, centralBody: .earth,
            gravitationalParameter: mu,
            options: OrbitalElementConversionOptions(circularEccentricityTolerance: 1e-6,
                                                       equatorialSineTolerance: 1e-6))
        #expect(snapped.angleConvention == .circularEquatorial)
        #expect(snapped.eccentricity == 0)
        #expect(snapped.inclination == 0)
        #expect(norm(try snapped.cartesianState().position - state.position) < 2)
        #expect(throws: OrbitalElementsError.invalidElements) {
            try OsculatingOrbitalElements(epoch: .j2000, frame: tree.root, centralBody: .earth,
                gravitationalParameter: mu, semiLatusRectumMetres: -1, eccentricity: 0,
                inclination: 0, longitudeOfAscendingNode: 0, argumentOfPeriapsis: 0, trueAnomaly: 0)
        }
    }

    @Test func rejectsDegenerateStatesAndUnreachableAnomalies() throws {
        let tree = FrameTree(), mu = try gm()
        for (r, v, error) in [
            (SIMD3<Double>.zero, SIMD3(1.0,0,0), OrbitalElementsError.zeroPosition),
            (SIMD3(1.0,0,0), SIMD3(2.0,0,0), .radialTrajectory),
            (SIMD3(1.0,0,0), SIMD3<Double>.zero, .radialTrajectory)
        ] {
            let state = try CartesianState(epoch: .j2000, frame: tree.root, position: r, velocity: v)
            #expect(throws: error) {
                try OsculatingOrbitalElements(state: state, centralBody: .earth, gravitationalParameter: mu)
            }
        }
        #expect(throws: OrbitalElementsError.anomalyOutsideConic) {
            try OsculatingOrbitalElements(epoch: .j2000, frame: tree.root, centralBody: .earth,
                gravitationalParameter: mu, semiLatusRectumMetres: 7e6, eccentricity: 2, inclination: 0,
                longitudeOfAscendingNode: 0, argumentOfPeriapsis: 0, trueAnomaly: .pi)
        }
        let rotating = FrameTree(pseudoInertial: false)
        let state = try CartesianState(epoch: .j2000, frame: rotating.root,
                                       position: SIMD3(1.0,0,0), velocity: SIMD3(0,1.0,0))
        #expect(throws: OrbitalElementsError.nonInertialFrame) {
            try OsculatingOrbitalElements(state: state, centralBody: .earth, gravitationalParameter: mu)
        }
        #expect(throws: OrbitalElementsError.invalidOptions) {
            try OrbitalElementConversionOptions(circularEccentricityTolerance: -1)
        }
    }
}
