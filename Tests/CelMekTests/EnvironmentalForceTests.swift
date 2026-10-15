// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct EnvironmentalForceTests {
    private let frame = Frame(name: "Environmental inertial", pseudoInertial: true)
    private func state(position: SIMD3<Double> = SIMD3(10,0,0), velocity: SIMD3<Double> = .zero, mass: Double = 10) throws -> SpacecraftState {
        SpacecraftState(cartesian: try CartesianState(epoch: .j2000, frame: frame, position: position, velocity: velocity),
            mass: try SpacecraftMassProperties(dryMassKilograms: mass))
    }
    private func context(_ origin: BodyIdentifier = .sun) throws -> ForceEvaluationContext {
        try ForceEvaluationContext(frame: frame, origin: origin)
    }
    private func radiation(_ fraction: Double = 1) throws -> SolarRadiationPressure {
        let provenance = try PhysicalDataProvenance(source: "Synthetic reference radiation", version: "1")
        return try SolarRadiationPressure(referencePressurePascals: SourcedPhysicalValue(value: 2, provenance: provenance),
            referenceDistanceMetres: SourcedPhysicalValue(value: 10, provenance: provenance),
            referenceAreaSquareMetres: 3, radiationCoefficient: 1.5,
            illumination: ConstantRadiationIllumination(illuminatedFraction: fraction))
    }
    @Test func dragRelativeVelocityAndAnalyticDecay() throws {
        let atmosphere = try UniformAtmosphere(frame: frame, origin: .sun, densityKilogramsPerCubicMetre: 0.1,
            velocityMetresPerSecond: SIMD3(1,0,0))
        let drag = try AtmosphericDrag(referenceAreaSquareMetres: 2, dragCoefficient: 2, atmosphere: atmosphere)
        let input = try state(velocity: SIMD3(4,4,0)) // relative (3,4,0), speed 5
        let force = try drag.contribution(for: input, in: context())
        #expect(abs(force.accelerationMetresPerSecondSquared.x + 0.3) < 1e-15)
        #expect(abs(force.accelerationMetresPerSecondSquared.y + 0.4) < 1e-15)
        #expect(force.propellantMassRateKilogramsPerSecond == 0)
        #expect(try drag.contribution(for: state(velocity: SIMD3(1,0,0)), in: context()).accelerationMetresPerSecondSquared == .zero)
        let dynamics = try SpacecraftDynamics(models: [drag], context: context())
        let seed = try state(position: .zero, velocity: SIMD3(6,0,0))
        let options = try NumericalPropagationOptions(absolutePositionToleranceMetres: 1e-9,
            absoluteVelocityToleranceMetresPerSecond: 1e-10, relativeTolerance: 1e-12)
        let end = try NumericalSpacecraftPropagator(initialState: seed, dynamics: dynamics, options: options)
            .spacecraftState(at: Epoch.j2000.advanced(by: Duration(seconds: 10)))
        // du/dt=-k*u², u0=5, k=.02; v=wind+u, x=wind*t+log(1+k*u0*t)/k.
        #expect(abs(end.cartesian.velocity.x-3.5) < 1e-8)
        #expect(abs(end.cartesian.position.x-(10+log(2)/0.02)) < 1e-7)
        #expect(end.mass == seed.mass)
    }
    @Test func radiationDirectionDistanceMassAndIllumination() throws {
        let model = try radiation()
        let full = try model.contribution(for: state(), in: context())
        #expect(abs(full.accelerationMetresPerSecondSquared.x-0.9) < 1e-15)
        let distant = try model.contribution(for: state(position: SIMD3(0,-20,0), mass: 20), in: context())
        #expect(abs(distant.accelerationMetresPerSecondSquared.y+0.1125) < 1e-15)
        #expect(distant.accelerationMetresPerSecondSquared.x == 0)
        #expect(try radiation(0).contribution(for: state(), in: context()).accelerationMetresPerSecondSquared == .zero)
        #expect(abs(try radiation(0.25).contribution(for: state(), in: context()).accelerationMetresPerSecondSquared.x-0.225) < 1e-15)
        #expect(full.propellantMassRateKilogramsPerSecond == 0)
    }
    private struct Source: EphemerisProvider {
        let metadata: EphemerisMetadata
        func state(at epoch: Epoch) throws -> CartesianState {
            try CartesianState(epoch: epoch, frame: frame, position: SIMD3(-10,0,0), velocity: .zero)
        }
    }
    @Test func ephemerisSourceAndCoverage() throws {
        let source = Source(metadata: try EphemerisMetadata(body: .sun, origin: .earth, frame: frame,
            coverage: Epoch.j2000...Epoch.j2000, coverageBasis: .callerSelected, modelTimeScale: .tt,
            velocityMethod: .analytic, accelerationMethod: .unavailable,
            provenance: PhysicalDataProvenance(source: "Synthetic Sun", version: "1")))
        let environment = try ForceEvaluationContext(frame: frame, origin: .earth, ephemerides: EphemerisCollection(providers: [source]))
        let force = try radiation().contribution(for: state(position: .zero), in: environment)
        #expect(abs(force.accelerationMetresPerSecondSquared.x-0.9) < 1e-15)
        #expect(throws: EnvironmentalForceError.originMismatch) {
            try radiation().contribution(for: state(), in: ForceEvaluationContext(frame: frame, origin: .moon,
                ephemerides: EphemerisCollection(providers: [source])))
        }
        #expect(throws: EnvironmentalForceError.missingEphemerides) { try radiation().contribution(for: state(), in: context(.earth)) }
    }
    private struct WrongSample: AtmosphereProvider {
        let sample: AtmosphereSample
        func sample(for state: SpacecraftState, in context: ForceEvaluationContext) throws -> AtmosphereSample { sample }
    }
    private struct BadIllumination: RadiationIlluminationProvider {
        func fraction(for state: SpacecraftState, sourcePosition: SIMD3<Double>, in context: ForceEvaluationContext) throws -> Double { .nan }
    }
    @Test func invalidEnvironmentAndContracts() throws {
        #expect(throws: EnvironmentalForceError.invalidDensity) {
            try UniformAtmosphere(frame: frame, origin: .earth, densityKilogramsPerCubicMetre: -1)
        }
        #expect(throws: EnvironmentalForceError.invalidIllumination) { try ConstantRadiationIllumination(illuminatedFraction: 1.1) }
        #expect(throws: EnvironmentalForceError.singularSource) { try radiation().contribution(for: state(position: .zero), in: context()) }
        let wrong = try AtmosphereSample(epoch: Epoch.j2000.advanced(by: Duration(seconds: 1)), frame: frame, origin: .sun,
            densityKilogramsPerCubicMetre: 1, velocityMetresPerSecond: .zero)
        #expect(throws: SpacecraftDynamicsError.epochMismatch) {
            try AtmosphericDrag(referenceAreaSquareMetres: 1, dragCoefficient: 1, atmosphere: WrongSample(sample: wrong))
                .contribution(for: state(), in: context())
        }
        let prototype = try radiation()
        let bad = try SolarRadiationPressure(referencePressurePascals: prototype.referencePressurePascals,
            referenceDistanceMetres: prototype.referenceDistanceMetres, referenceAreaSquareMetres: 1,
            radiationCoefficient: 1, illumination: BadIllumination())
        #expect(throws: EnvironmentalForceError.invalidIllumination) { try bad.contribution(for: state(), in: context()) }
    }
}
