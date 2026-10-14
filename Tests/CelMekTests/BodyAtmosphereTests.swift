// SPDX-License-Identifier: Apache-2.0
import Testing
import Math
import Atmospheres
@testable import CelMek
@Suite struct BodyAtmosphereTests {
    @Test func lookupRotationAndDrag() throws {
        let frame = Frame(name: "Atmosphere inertial", pseudoInertial: true)
        let earth = try #require(BodyIdentifier.earth.atmosphericModel)
        #expect(BodyIdentifier.jupiterBarycentre.atmosphericModel == nil)
        let atmosphere = try SphericalBodyAtmosphere(bodyIdentifier: .earth, profile: earth,
            referenceRadiusMetres: 6371000, frame: frame, angularVelocityRadiansPerSecond: SIMD3(0,0,1e-4))
        let context = try ForceEvaluationContext(frame: frame, origin: .earth)
        let state = SpacecraftState(cartesian: try CartesianState(epoch: .j2000,frame: frame,
            position: SIMD3(6371000,0,0), velocity: SIMD3(0,700,0)),
            mass: try SpacecraftMassProperties(dryMassKilograms: 100))
        let medium = try atmosphere.sample(for: state, in: context)
        #expect(abs(medium.densityKilogramsPerCubicMetre-1.225) < 1e-6)
        #expect(abs(medium.velocityMetresPerSecond.y-637.1) < 1e-10)
        let drag = try AtmosphericDrag(referenceAreaSquareMetres: 1, dragCoefficient: 2, atmosphere: atmosphere)
        let force = try drag.contribution(for: state,in: context)
        #expect(abs(force.accelerationMetresPerSecondSquared.y + medium.densityKilogramsPerCubicMetre*62.9*62.9/100) < 1e-8)
        #expect(throws: EnvironmentalForceError.originMismatch) {
            try atmosphere.sample(for: state,in: ForceEvaluationContext(frame: frame,origin: .mars))
        }
        #expect(throws: BodyModelError.bodyMismatch) {
            try SphericalBodyAtmosphere(bodyIdentifier: .mars,profile: earth,referenceRadiusMetres: 1,frame: frame)
        }
    }
}
