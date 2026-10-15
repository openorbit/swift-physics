// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct GravityModelTests {
    private struct Body: EphemerisProvider {
        let metadata: EphemerisMetadata
        let position: SIMD3<Double>
        var velocity: SIMD3<Double> = .zero
        func state(at epoch: Epoch) throws -> CartesianState {
            try CartesianState(epoch: epoch,frame: frame,
                position: position+velocity*epoch.duration(since: .j2000).seconds,velocity: velocity)
        }
    }
    private func vector(_ values: [Double]) -> SIMD3<Double> { SIMD3(values[0],values[1],values[2]) }
    private func gm(_ value: Double) throws -> GravitationalParameter {
        try GravitationalParameter(cubicMetresPerSecondSquared: value,convention: .newtonian,
            provenance: PhysicalDataProvenance(source: "Explicit reference GM",version: "1"))
    }
    private func body(_ frame: Frame, position: SIMD3<Double>, origin: BodyIdentifier = .earth) throws -> Body {
        Body(metadata: try EphemerisMetadata(body: .moon,origin: origin,frame: frame,
            coverage: Epoch.j2000.advanced(by: Duration(seconds: -100))...Epoch.j2000.advanced(by: Duration(seconds: 100)),
            coverageBasis: .callerSelected,modelTimeScale: .tt,velocityMethod: .analytic,accelerationMethod: .unavailable,
            provenance: PhysicalDataProvenance(source: "Synthetic prescribed body",version: "1")),position: position)
    }
    private func state(_ frame: Frame,position: SIMD3<Double>,epoch: Epoch = .j2000, mass: Double = 1) throws -> SpacecraftState {
        SpacecraftState(cartesian: try CartesianState(epoch: epoch,frame: frame,position: position,velocity: .zero),
            mass: try SpacecraftMassProperties(dryMassKilograms: mass))
    }
    private struct Reference: Decodable {
        let position: [Double]; let bodyPosition: [Double]?; let gm: Double; let acceleration: [Double]
    }
    @Test func independentHighPrecisionReferences() throws {
        let url = try #require(Bundle.module.url(forResource: "gravity-reference",withExtension: "json",subdirectory: "Fixtures"))
        let rows = try JSONDecoder().decode([Reference].self,from: Data(contentsOf: url))
        let frame = Frame(name: "Earth inertial",pseudoInertial: true)
        for row in rows {
            let model: any ForceModel
            let collection: EphemerisCollection?
            if let position = row.bodyPosition {
                model = try ThirdBodyGravity(body: .moon,centralBody: .earth,gravitationalParameter: gm(row.gm))
                collection = try EphemerisCollection(providers: [body(frame,position: vector(position))])
            } else {
                model = try CentralGravity(body: .earth,gravitationalParameter: gm(row.gm)); collection = nil
            }
            let context = try ForceEvaluationContext(frame: frame,origin: .earth,ephemerides: collection)
            let input = try state(frame,position: vector(row.position))
            let actual = try model.contribution(for: input,in: context)
            for axis in 0..<3 {
                #expect(abs(actual.accelerationMetresPerSecondSquared[axis]-row.acceleration[axis]) <= max(abs(row.acceleration[axis])*1e-13,1e-32))
            }
            #expect(actual.propellantMassRateKilogramsPerSecond == 0)
            let heavier = try state(frame,position: vector(row.position),mass: 1e6)
            #expect(try model.contribution(for: heavier,in: context).accelerationMetresPerSecondSquared == actual.accelerationMetresPerSecondSquared)
        }
    }
    @Test func movingEphemerisAndForceComposition() throws {
        let frame = Frame(name: "Earth inertial",pseudoInertial: true)
        var provider = try body(frame,position: SIMD3(100,0,0)); provider.velocity = SIMD3(1,0,0)
        let context = try ForceEvaluationContext(frame: frame,origin: .earth,
            ephemerides: EphemerisCollection(providers: [provider]))
        let central = try CentralGravity(body: .earth,gravitationalParameter: gm(10))
        let third = try ThirdBodyGravity(body: .moon,centralBody: .earth,gravitationalParameter: gm(1))
        let dynamics = try SpacecraftDynamics(models: [central,third],context: context)
        let epoch = try Epoch.j2000.advanced(by: Duration(seconds: 20))
        let input = try state(frame,position: SIMD3(1,0,0),epoch: epoch)
        let result = try dynamics.derivative(at: input)
        let expected = -10 + 1/(119.0*119)-1/(120.0*120)
        #expect(abs(result.velocityRateMetresPerSecondSquared.x-expected) < 1e-14)
        #expect(result.velocityRateMetresPerSecondSquared.y == 0)
        #expect(try third.contribution(for: state(frame,position: .zero),in: context).accelerationMetresPerSecondSquared == .zero)
        #expect(throws: EphemerisError.outsideCoverage) {
            try third.contribution(for: state(frame,position: .zero,
                epoch: Epoch.j2000.advanced(by: Duration(seconds: 101))),in: context)
        }
    }
    @Test func originFrameBodyAndGMErrors() throws {
        let frame = Frame(name: "Earth inertial",pseudoInertial: true)
        let input = try state(frame,position: SIMD3(1,0,0))
        let third = try ThirdBodyGravity(body: .moon,centralBody: .earth,gravitationalParameter: gm(1))
        #expect(throws: GravityModelError.identicalBodies) { try ThirdBodyGravity(body: .earth,centralBody: .earth,gravitationalParameter: gm(1)) }
        let incompatible = try GravitationalParameter(cubicMetresPerSecondSquared: 1,convention: .tdbCompatible,
            provenance: PhysicalDataProvenance(source: "Test",version: "1"))
        #expect(throws: GravityModelError.incompatibleGMConvention) { try CentralGravity(body: .earth,gravitationalParameter: incompatible) }
        #expect(throws: GravityModelError.incompatibleGMConvention) { try ThirdBodyGravity(body: .moon,centralBody: .earth,gravitationalParameter: incompatible) }
        #expect(throws: GravityModelError.missingEphemerides) {
            try third.contribution(for: input,in: ForceEvaluationContext(frame: frame,origin: .earth))
        }
        #expect(throws: GravityModelError.originMismatch) {
            try third.contribution(for: input,in: ForceEvaluationContext(frame: frame,origin: .sun))
        }
        #expect(throws: GravityModelError.originMismatch) {
            try third.contribution(for: input,in: ForceEvaluationContext(frame: frame,origin: .earth,
                ephemerides: EphemerisCollection(providers: [body(frame,position: SIMD3(100,0,0),origin: .sun)])))
        }
        #expect(throws: SpacecraftDynamicsError.frameMismatch) {
            try third.contribution(for: input,in: ForceEvaluationContext(frame: frame,origin: .earth,
                ephemerides: EphemerisCollection(providers: [body(Frame(name: "other",pseudoInertial: true),position: SIMD3(100,0,0))])))
        }
        #expect(throws: EphemerisError.unsupportedBody) {
            try third.contribution(for: input,in: ForceEvaluationContext(frame: frame,origin: .earth,ephemerides: EphemerisCollection(providers: [])))
        }
        let central = try CentralGravity(body: .earth,gravitationalParameter: gm(1))
        #expect(throws: GravityModelError.originMismatch) {
            try central.contribution(for: input,in: ForceEvaluationContext(frame: frame,origin: .sun))
        }
        #expect(throws: GravityModelError.singularPosition) {
            try central.contribution(for: state(frame,position: .zero),in: ForceEvaluationContext(frame: frame,origin: .earth))
        }
        for position in [SIMD3<Double>.zero,SIMD3(1,0,0)] {
            #expect(throws: GravityModelError.singularPosition) {
                try third.contribution(for: input,in: ForceEvaluationContext(frame: frame,origin: .earth,
                    ephemerides: EphemerisCollection(providers: [body(frame,position: position)])))
            }
        }
    }
    @Test func combinedGravityIntegratesWithMovingProvider() throws {
        let frame = Frame(name: "relative inertial",pseudoInertial: true)
        let primaryGM = try gm(10)
        let seed = SpacecraftState(cartesian: try CartesianState(epoch: .j2000,frame: frame,
            position: SIMD3(1,0,0),velocity: SIMD3(0,sqrt(10),0)),mass: try SpacecraftMassProperties(dryMassKilograms: 1))
        var moon = try body(frame,position: SIMD3(100,0,0)); moon.velocity = SIMD3(1,0,0)
        let context = try ForceEvaluationContext(frame: frame,origin: .earth,ephemerides: EphemerisCollection(providers: [moon]))
        let primary = try CentralGravity(body: .earth,gravitationalParameter: primaryGM)
        let perturber = try ThirdBodyGravity(body: .moon,centralBody: .earth,gravitationalParameter: gm(1))
        let dynamics = try SpacecraftDynamics(models: [primary,perturber],context: context)
        let model = try NumericalSpacecraftPropagator(initialState: seed,dynamics: dynamics,accuracy: .full,
            options: NumericalPropagationOptions(absolutePositionToleranceMetres: 1e-12,
                absoluteVelocityToleranceMetresPerSecond: 1e-12,relativeTolerance: 1e-12,
                initialStepSeconds: 0.1,maximumStepSeconds: 0.1))
        let analytic = try TwoBodyPropagator(initialState: seed.cartesian,centralBody: .earth,gravitationalParameter: primaryGM)
        let date = try Epoch.j2000.advanced(by: Duration(seconds: 1))
        let actual = try model.spacecraftState(at: date), unperturbed = try analytic.state(at: date)
        let delta = actual.cartesian.position-unperturbed.position
        let difference = hypot(hypot(delta.x,delta.y),delta.z)
        #expect(difference > 1e-9 && difference < 1e-3)
        #expect(actual.mass == seed.mass)
        let derivative = try dynamics.derivative(at: actual)
        #expect(actual.cartesian.acceleration == derivative.velocityRateMetresPerSecondSquared)
        #expect(try model.state(at: date) == actual.cartesian)
    }

}
