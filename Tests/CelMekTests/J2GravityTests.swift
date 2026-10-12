// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct J2GravityTests {
    private let frame = Frame(name: "Synthetic inertial", pseudoInertial: true)
    private func model(_ j2: Double = 0.01, rotation: Quaternion = .identity, rate: Double = 0) throws -> J2Gravity {
        let provenance = try PhysicalDataProvenance(source: "Synthetic test parameters", version: "1")
        return try J2Gravity(body: .earth,
            gravitationalParameter: GravitationalParameter(cubicMetresPerSecondSquared: 100, convention: .newtonian, provenance: provenance),
            equatorialRadius: SourcedPhysicalValue(value: 2, provenance: provenance), j2: j2,
            coefficientProvenance: provenance,
            orientation: FixedAxisOrientationProvider(referenceFrame: frame, axis: SIMD3(0,0,1), rate: rate, initialRotation: rotation))
    }
    private func state(_ p: SIMD3<Double>, epoch: Epoch = .j2000) throws -> SpacecraftState {
        SpacecraftState(cartesian: try CartesianState(epoch: epoch, frame: frame, position: p, velocity: .zero),
            mass: try SpacecraftMassProperties(dryMassKilograms: 1))
    }
    private func acceleration(_ p: SIMD3<Double>, model: J2Gravity, epoch: Epoch = .j2000) throws -> SIMD3<Double> {
        try model.contribution(for: state(p, epoch: epoch), in: ForceEvaluationContext(frame: frame, origin: .earth)).accelerationMetresPerSecondSquared
    }
    @Test func equatorPoleAndSignedCoefficient() throws {
        let gravity = try model()
        let equator = try acceleration(SIMD3(10,0,0), model: gravity)
        let pole = try acceleration(SIMD3(0,0,10), model: gravity)
        #expect(abs(equator.x + 0.0006) < 1e-18)
        #expect(abs(pole.z - 0.0012) < 1e-18)
        #expect(equator.y == 0 && equator.z == 0 && pole.x == 0 && pole.y == 0)
        #expect(try acceleration(SIMD3(3,4,5), model: model(-0.01)) == -acceleration(SIMD3(3,4,5), model: gravity))
        #expect(try acceleration(SIMD3(3,4,5), model: model(0)) == .zero)
        #expect(throws: J2GravityError.invalidCoefficient) { try model(.nan) }
        #expect(throws: GravityModelError.singularPosition) { try acceleration(.zero, model: gravity) }
    }
    @Test func independentPotentialGradient() throws {
        // J2 potential per unit mass: mu*J2*R²/(2r³)*(3z²/r²-1).
        // Central differences of the potential independently check all force components.
        func potential(_ p: SIMD3<Double>) -> Double {
            let r = sqrt((p*p).sum())
            return 100 * 0.01 * 4 / (2*r*r*r) * (3*p.z*p.z/(r*r)-1)
        }
        let gravity = try model()
        for p in [SIMD3(3.0,4,5), SIMD3(-7.0,2,-3), SIMD3(0.0,0,8)] {
            let a = try acceleration(p, model: gravity)
            for axis in 0..<3 {
                var plus = p, minus = p
                plus[axis] += 1e-4; minus[axis] -= 1e-4
                let reference = -(potential(plus)-potential(minus))/2e-4
                #expect(abs(a[axis]-reference) < 1e-11)
            }
        }
    }
    @Test func tiltedPoleAndSpinInvariance() throws {
        let q = Quaternion(angle: 0.8, axis: SIMD3(1,2,3))
        let p = SIMD3(3.0,4,5)
        let expected = try q.invMatrix * acceleration(q.matrix*p, model: model())
        let actual = try acceleration(p, model: model(rotation: q))
        #expect(sqrt(((actual-expected)*(actual-expected)).sum()) < 1e-16)
        let later = try Epoch.j2000.advanced(by: Duration(seconds: 17))
        let spinning = try acceleration(p, model: model(rate: 2), epoch: later)
        let fixed = try acceleration(p, model: model())
        #expect(sqrt(((spinning-fixed)*(spinning-fixed)).sum()) < 1e-16)
        #expect(throws: GravityModelError.originMismatch) {
            try model().contribution(for: state(p), in: ForceEvaluationContext(frame: frame, origin: .moon))
        }
        #expect(throws: SpacecraftDynamicsError.frameMismatch) {
            try model().contribution(for: state(p), in: ForceEvaluationContext(frame: Frame(name: "Other", pseudoInertial: true), origin: .earth))
        }
    }

    @Test func orbitalNodalRegression() throws {
        // First-order secular prediction: ΔΩ = -3π J2 (R/p)² cos(i) per orbit.
        // https://ai-solutions.com/_freeflyeruniversityguide/j2_perturbation.htm
        // A weak synthetic coefficient limits omitted higher-order/short-period terms.
        let gravity = try model(1e-5)
        let inclination = Double.pi / 3
        let speed = sqrt(100.0 / 10)
        let initial = SpacecraftState(cartesian: try CartesianState(epoch: .j2000, frame: frame,
            position: SIMD3(10,0,0), velocity: SIMD3(0,speed*cos(inclination),speed*sin(inclination))),
            mass: try SpacecraftMassProperties(dryMassKilograms: 1))
        let central = try CentralGravity(body: .earth, gravitationalParameter: gravity.gravitationalParameter)
        let dynamics = try SpacecraftDynamics(models: [central,gravity],
            context: ForceEvaluationContext(frame: frame, origin: .earth))
        // Synthetic metre-scale orbit needs tighter absolute tolerances than game presets.
        let options = try NumericalPropagationOptions(absolutePositionToleranceMetres: 1e-11,
            absoluteVelocityToleranceMetresPerSecond: 1e-12, relativeTolerance: 1e-12,
            initialStepSeconds: 0.1, maximumStepSeconds: 1)
        let propagator = try NumericalSpacecraftPropagator(initialState: initial, dynamics: dynamics, options: options)
        let duration = 10 * 2 * Double.pi * sqrt(1000.0 / 100)
        let result = try propagator.state(at: Epoch.j2000.advanced(by: Duration(seconds: duration)))
        let r = result.position, v = result.velocity
        let hx = r.y*v.z-r.z*v.y, hy = r.z*v.x-r.x*v.z
        let measured = atan2(hx,-hy)
        let expected = -30 * Double.pi * 1e-5 * 0.04 * cos(inclination)
        #expect(measured < 0)
        #expect(abs((measured-expected)/expected) < 0.01)
        // Axisymmetry conserves angular momentum about the fixed pole.
        #expect(abs((r.x*v.y-r.y*v.x)-10*speed*cos(inclination)) < 1e-7)
    }
}
