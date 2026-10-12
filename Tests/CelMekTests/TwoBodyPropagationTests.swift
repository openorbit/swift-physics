import Foundation
import Testing
@testable import CelMek

@Suite struct TwoBodyPropagationTests {
    private func gm(_ value: Double, convention: GravitationalParameterConvention = .newtonian) throws -> GravitationalParameter {
        try GravitationalParameter(cubicMetresPerSecondSquared: value, convention: convention,
            provenance: PhysicalDataProvenance(source: "Test", version: "1", reference: "analytic model"))
    }
    private func norm(_ v: SIMD3<Double>) -> Double { hypot(hypot(v.x, v.y), v.z) }
    private func cross(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> SIMD3<Double> {
        SIMD3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x)
    }
    private func orbit(radius: Double = 7e6, speed: Double = 9000, mu: Double = 3.986004418e14,
                       options: TwoBodyPropagationOptions = .default) throws -> TwoBodyPropagator {
        let tree = FrameTree()
        let state = try CartesianState(epoch: .j2000, frame: tree.root,
                                       position: SIMD3(radius, 0, 0), velocity: SIMD3(0, speed, 0))
        return try TwoBodyPropagator(initialState: state, centralBody: .earth,
                                     gravitationalParameter: gm(mu), options: options)
    }

    @Test func independentConicReferences() throws {
        struct Case: Decodable {
            let name: String
            let radius, speed, gm, seconds: Double
            let position, velocity: [Double]
        }
        struct File: Decodable { let cases: [Case] }
        let url = Bundle.module.url(forResource: "two_body_reference", withExtension: "json", subdirectory: "Fixtures")!
        let fixtures = try JSONDecoder().decode(File.self, from: Data(contentsOf: url))
        let options = try TwoBodyPropagationOptions(absoluteTimeToleranceSeconds: 1e-13,
                                                   relativeTimeTolerance: 1e-14)
        for fixture in fixtures.cases {
            let model = try orbit(radius: fixture.radius, speed: fixture.speed, mu: fixture.gm, options: options)
            let epoch = try Epoch.j2000.advanced(by: Duration(seconds: fixture.seconds))
            let result = try model.propagation(at: epoch)
            let p = SIMD3(fixture.position[0], fixture.position[1], fixture.position[2])
            let v = SIMD3(fixture.velocity[0], fixture.velocity[1], fixture.velocity[2])
            #expect(norm(result.state.position - p) < max(1e-10, norm(p) * 3e-12), "\(fixture.name)")
            #expect(norm(result.state.velocity - v) < max(1e-11, norm(v) * 3e-12), "\(fixture.name)")
            #expect(result.state.epoch == epoch)
        }
    }

    @Test func nonPeriapsisAnalyticEllipseReference() throws {
        // mu=1, a=2, e=1/2 at eccentric anomaly pi/2. These
        // Cartesian values and elapsed times follow the conic equations directly.
        let tree = FrameTree()
        let initial = try CartesianState(epoch: .j2000, frame: tree.root,
            position: SIMD3(-1.0,sqrt(3),0), velocity: SIMD3(-1/sqrt(2.0),0,0))
        let model = try TwoBodyPropagator(initialState: initial, centralBody: .earth,
            gravitationalParameter: gm(1), options: TwoBodyPropagationOptions(absoluteTimeToleranceSeconds: 1e-13))
        let periapsisEpoch = try Epoch.j2000.advanced(by: Duration(seconds: -2*sqrt(2.0)*(.pi/2-0.5)))
        let apoapsisEpoch = try Epoch.j2000.advanced(by: Duration(seconds: 2*sqrt(2.0)*(.pi/2+0.5)))
        let periapsis = try model.state(at: periapsisEpoch)
        let apoapsis = try model.state(at: apoapsisEpoch)
        #expect(norm(periapsis.position-SIMD3(1.0,0,0)) < 1e-11)
        #expect(norm(periapsis.velocity-SIMD3(0.0,sqrt(1.5),0)) < 1e-11)
        #expect(norm(apoapsis.position-SIMD3(-3.0,0,0)) < 1e-11)
        #expect(norm(apoapsis.velocity-SIMD3(0.0,-sqrt(1.0/6),0)) < 1e-11)
    }

    @Test func conservationAndReversePropagation() throws {
        for speed in [5000.0, 9000, sqrt(2 * 3.986004418e14 / 7e6), 14000] {
            let model = try orbit(speed: speed), initial = model.initialState
            let mu = model.gravitationalParameter.cubicMetresPerSecondSquared
            let energy = (initial.velocity * initial.velocity).sum()/2 - mu/norm(initial.position)
            let h = cross(initial.position, initial.velocity)
            for seconds in [-30000.0, -1000, 0, 1000, 30000] {
                let epoch = try Epoch.j2000.advanced(by: Duration(seconds: seconds))
                let state = try model.state(at: epoch)
                let actualEnergy = (state.velocity * state.velocity).sum()/2 - mu/norm(state.position)
                #expect(abs(actualEnergy-energy) < 1e-4)
                #expect(norm(cross(state.position, state.velocity)-h) < norm(h)*2e-11)
                let reverse = try TwoBodyPropagator(initialState: state, centralBody: .earth,
                                                     gravitationalParameter: model.gravitationalParameter)
                let rebuilt = try reverse.state(at: .j2000)
                #expect(norm(rebuilt.position-initial.position) < 0.003)
                #expect(norm(rebuilt.velocity-initial.velocity) < 3e-6)
                let radius = norm(state.position)
                let acceleration = try #require(state.acceleration)
                #expect(norm(acceleration + (mu/radius/radius) * (state.position/radius)) < 1e-14)
            }
        }
    }

    @Test func arbitraryInclinedInitialStateAndQueryOrder() throws {
        let tree = FrameTree(), parameter = try gm(3.986004418e14)
        let elements = try OsculatingOrbitalElements(epoch: .j2000, frame: tree.root, centralBody: .earth,
            gravitationalParameter: parameter, semiLatusRectumMetres: 9e6, eccentricity: 0.6,
            inclination: 1.3, longitudeOfAscendingNode: 2, argumentOfPeriapsis: 0.7, trueAnomaly: -1)
        let initial = try elements.cartesianState()
        let model = try TwoBodyPropagator(initialState: initial, centralBody: .earth, gravitationalParameter: parameter)
        let epochs = try [5000.0, -2000, 0, 100, 5000].map { try Epoch.j2000.advanced(by: Duration(seconds: $0)) }
        let states = try epochs.map { try model.state(at: $0) }
        for (epoch, expected) in zip(epochs.reversed(), states.reversed()) {
            #expect(try model.state(at: epoch) == expected)
        }
        let end = states[0]
        let inverse = try TwoBodyPropagator(initialState: end, centralBody: .earth, gravitationalParameter: parameter)
        #expect(norm(try inverse.state(at: .j2000).position - initial.position) < 0.001)
    }

    @Test func multiRevolutionCircularAndZeroTime() throws {
        let model = try orbit(radius: 1, speed: 1, mu: 1)
        let epoch = try Epoch.j2000.advanced(by: Duration(seconds: 1000 * 2 * .pi + .pi/2))
        let state = try model.state(at: epoch)
        #expect(norm(state.position-SIMD3(0.0,1,0)) < 1e-10)
        #expect(norm(state.velocity-SIMD3(-1.0,0,0)) < 1e-10)
        let zero = try model.propagation(at: .j2000)
        #expect(zero.iterations == 0)
        #expect(zero.state.position == model.initialState.position)
        #expect(zero.state.velocity == model.initialState.velocity)
        #expect(zero.state.acceleration == SIMD3(-1.0,0,0))
    }

    @Test func explicitFailures() throws {
        #expect(throws: PropagationError.invalidOptions) { try TwoBodyPropagationOptions(maximumIterations: 0) }
        let model = try orbit(options: TwoBodyPropagationOptions(maximumIterations: 1))
        let epoch = try Epoch.j2000.advanced(by: Duration(seconds: 4000))
        #expect(throws: PropagationError.convergenceFailure(iterations: 1)) { try model.state(at: epoch) }
        let hyperbolic = try orbit(radius: 1, speed: 2, mu: 1,
                                  options: TwoBodyPropagationOptions(maximumBracketExpansions: 1))
        let distant = try Epoch.j2000.advanced(by: Duration(seconds: 1e8))
        #expect(throws: PropagationError.bracketFailure) { try hyperbolic.state(at: distant) }
        let tree = FrameTree()
        let radial = try CartesianState(epoch: .j2000, frame: tree.root,
                                       position: SIMD3(1.0,0,0), velocity: SIMD3(1.0,0,0))
        #expect(throws: PropagationError.radialTrajectory) {
            try TwoBodyPropagator(initialState: radial, centralBody: .earth, gravitationalParameter: gm(1))
        }
        let zero = try CartesianState(epoch: .j2000, frame: tree.root, position: .zero, velocity: SIMD3(0.0,1,0))
        #expect(throws: PropagationError.zeroPosition) {
            try TwoBodyPropagator(initialState: zero, centralBody: .earth, gravitationalParameter: gm(1))
        }
        let rotating = FrameTree(pseudoInertial: false)
        let rotatingState = try CartesianState(epoch: .j2000, frame: rotating.root,
            position: SIMD3(1.0,0,0), velocity: SIMD3(0.0,1,0))
        #expect(throws: PropagationError.nonInertialFrame) {
            try TwoBodyPropagator(initialState: rotatingState, centralBody: .earth, gravitationalParameter: gm(1))
        }
        let compatible = try TwoBodyPropagator(initialState: orbit().initialState, centralBody: .earth,
            gravitationalParameter: gm(3.986e14, convention: .ttCompatible))
        _ = try compatible.state(at: .j2000)
        let circular = try orbit()
        #expect(throws: PropagationError.incompatibleGMConvention) {
            try TwoBodyPropagator(initialState: circular.initialState, centralBody: .earth,
                                  gravitationalParameter: gm(3.986e14, convention: .tdbCompatible))
        }
    }
}
