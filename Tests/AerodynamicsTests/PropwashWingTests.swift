// SPDX-License-Identifier: Apache-2.0
import Testing
import FluidDynamics
import Aerodynamics

@Suite struct PropwashWingTests {
    @Test func localWingResponseAndVelocityComposition() throws {
        struct Section: WingSectionCoefficientModel {
            func coefficients(for query: WingSectionQuery) throws -> AirfoilCoefficients {
                AirfoilCoefficients(lift: 1,drag: 0.1,moment: 0)
            }
        }
        let rotor = try ActuatorDiskSlipstream(centreMetres: .zero,downstreamAxis: SIMD3(-1,0,0),radiusMetres: 1,
            thrustNewtons: 200,fluidTorqueNewtonMetres: 10,densityKilogramsPerCubicMetre: 2,upstreamAxialSpeedMetresPerSecond: 10)
        let field = try SuperposedFluidVelocityField(backgroundMetresPerSecond: SIMD3(-10,0,0),increments: [rotor])
        let segment = try WingSegment(positionMetres: SIMD3(-2,0,0),axes: WingSectionAxes(),chordMetres: 1,
            areaSquareMetres: 1,model: Section())
        let base = try FluidFlow(relativeVelocityMetresPerSecond: SIMD3(10,0,0),densityKilogramsPerCubicMetre: 2)
        let local = try base.replacingRelativeVelocity(-field.velocityMetresPerSecond(at: segment.positionMetres))
        #expect(try segment.evaluate(bodyFlow: local).loads.forceNewtons.z < segment.evaluate(bodyFlow: base).loads.forceNewtons.z)
        #expect(try rotor.velocityMetresPerSecond(at: SIMD3(1,0,0)) == .zero)
        #expect(try rotor.velocityMetresPerSecond(at: SIMD3(-1,2,0)) == .zero)
    }
}
