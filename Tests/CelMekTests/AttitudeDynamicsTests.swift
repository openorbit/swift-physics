// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct AttitudeDynamicsTests {
    let reference = Frame(name: "Attitude inertial",pseudoInertial: true)
    func body() -> Frame { Frame(parent: reference,name: "Spacecraft body",pseudoInertial: false) }
    func epoch(_ seconds: Double) throws -> Epoch { try Epoch.j2000.advanced(by: CelMek.Duration(seconds: seconds)) }
    func context(_ body: Frame, inertia: RigidBodyInertia? = nil) throws -> RotationalEvaluationContext {
        try RotationalEvaluationContext(inertia: inertia ?? RigidBodyInertia(principalMomentsKilogramSquareMetres: SIMD3(2,3,4)),
            bodyFrame: body,referenceContext: ForceEvaluationContext(frame: reference,origin: .earth))
    }
    func attitude(_ body: Frame, q: Quaternion = .identity, omega: SIMD3<Double> = .zero) throws -> AttitudeState {
        try AttitudeState(epoch: .j2000,referenceFrame: reference,bodyFrame: body,bodyToReference: q,angularVelocityBodyRadiansPerSecond: omega)
    }
    func translation(_ p: SIMD3<Double> = .zero, _ v: SIMD3<Double> = .zero) throws -> SpacecraftState {
        try SpacecraftState(cartesian: CartesianState(epoch: .j2000,frame: reference,position: p,velocity: v),
            mass: SpacecraftMassProperties(dryMassKilograms: 5,propellantMassKilograms: 5))
    }
    func angleDifference(_ a: Quaternion, _ b: Quaternion) -> Double {
        let delta = Quaternion.dot(lhs: a,rhs: b) >= 0 ? a-b : a+b
        return 4*asin(min(1,sqrt(Quaternion.dot(lhs: delta,rhs: delta))/2))
    }
    struct Circle: StateProvider {
        let frame: Frame
        let rate: Double
        func state(at epoch: Epoch) throws -> CartesianState {
            let t = try epoch.duration(since: .j2000).seconds, angle = rate*t
            let p = SIMD3<Double>(10*cos(angle),10*sin(angle),0)
            return try CartesianState(epoch: epoch,frame: frame,position: p,
                velocity: SIMD3(-10*rate*sin(angle),10*rate*cos(angle),0),acceleration: -p*rate*rate)
        }
    }
    @Test func quaternionConventionAndFrameKinematics() throws {
        let b = body(), seed = try attitude(b,q: Quaternion(angle: .pi/2,axis: SIMD3(0,0,1)),omega: SIMD3(0,0,0.2))
        let spin = UniformSpinAttitudeProvider(initialState: seed), axes = try spin.attitude(at: .j2000)
        #expect(routeNorm(try axes.referenceVector(fromBody: SIMD3(1,0,0))-SIMD3(0,1,0)) < 1e-14)
        let orientation = try axes.orientation()
        #expect(orientation.angularVelocity == SIMD3(0,0,-0.2))
        let old = OrientationAttitudeProvider(orientation: AttitudeOrientationProvider(spin),bodyFrame: b)
        #expect(try angleDifference(old.attitude(at: epoch(2)).bodyToReference,spin.attitude(at: epoch(2)).bodyToReference) < 1e-14)
        let tree = try FrameTree(root: reference)
        let centre = LambertTransferPlannerTests.Linear(frame: reference,p: SIMD3(100,0,0),v: SIMD3(3,0,0))
        try tree.registerSpacecraftFrame(spacecraft: centre,attitude: spin,accelerationPolicy: .assumeZero)
        // A point fixed at body +X moves around the translating centre with ω×r.
        let local = try CartesianState(epoch: .j2000,frame: b,position: SIMD3(1,0,0),velocity: .zero,acceleration: .zero)
        let global = try tree.transform(from: b,to: reference,at: .j2000).applying(to: local,destination: reference)
        #expect(routeNorm(global.position-SIMD3(100,1,0)) < 1e-13)
        #expect(routeNorm(global.velocity-SIMD3(2.8,0,0)) < 1e-13)
        #expect(routeNorm(global.acceleration!-SIMD3(0,-0.04,0)) < 1e-13)
        #expect(throws: AttitudeError.angularAccelerationUnavailable) { try seed.orientation() }
    }
    @Test func prescribedTrackingAndDegeneracies() throws {
        let b = body(), circle = Circle(frame: reference,rate: 0.1)
        for direction in [AttitudeTrackingDirection.nadir,.velocity] {
            let provider = try TrackingAttitudeProvider(spacecraft: circle,bodyFrame: b,direction: direction,rollReference: SIMD3(0,0,1),derivativeStep: CelMek.Duration(seconds: 0.01))
            for t in [0.0,5,Double.pi/0.1] {
                let state = try provider.attitude(at: epoch(t)), orbital = try circle.state(at: epoch(t))
                let expected = direction == .nadir ? -orbital.position/10 : orbital.velocity
                #expect(routeNorm(try state.referenceVector(fromBody: SIMD3(1,0,0))-attitudeUnit(expected)) < 1e-12)
                #expect(routeNorm(state.angularVelocityBodyRadiansPerSecond-SIMD3(0,0,0.1)) < 1e-10)
                #expect(routeNorm(state.angularAccelerationBodyRadiansPerSecondSquared!) < 1e-7)
            }
        }
        let stationary = LambertTransferPlannerTests.Linear(frame: reference,p: .zero,v: .zero)
        let tracking = try TrackingAttitudeProvider(spacecraft: stationary,bodyFrame: b,direction: .target,target: circle,
            rollReference: SIMD3(0,0,1),derivativeStep: CelMek.Duration(seconds: 0.01))
        #expect(routeNorm(try tracking.attitude(at: .j2000).referenceVector(fromBody: SIMD3(1,0,0))-SIMD3(1,0,0)) < 1e-14)
        #expect(throws: AttitudeError.undefinedDirection) {
            try TrackingAttitudeProvider(spacecraft: stationary,bodyFrame: b,direction: .velocity,rollReference: SIMD3(0,0,1),derivativeStep: CelMek.Duration(seconds: 0.01)).attitude(at: .j2000)
        }
        #expect(throws: AttitudeError.undefinedRoll) {
            try TrackingAttitudeProvider(spacecraft: circle,bodyFrame: b,direction: .nadir,rollReference: SIMD3(1,0,0),derivativeStep: CelMek.Duration(seconds: 0.01)).attitude(at: .j2000)
        }
    }
    @Test func fullInertiaSolveAndPhysicalValidation() throws {
        // Rotation of principal moments (2,3,4) by 45 degrees around Z.
        let inertia = try RigidBodyInertia(xx: 2.5,yy: 2.5,zz: 4,xy: -0.5), v = SIMD3<Double>(1,2,3)
        #expect(routeNorm(try inertia.solve(inertia.matrix*v)-v) < 1e-14)
        #expect(abs(try inertia.kineticEnergyJoules(for: v)-23.25) < 1e-14)
        #expect(throws: AttitudeError.invalidInertia) { try RigidBodyInertia(xx: 1,yy: 1,zz: 10) }
        #expect(throws: AttitudeError.invalidInertia) { try RigidBodyInertia(xx: 1,yy: 1,zz: 1,xy: 2) }
        #expect(throws: AttitudeError.invalidInertia) { try RigidBodyInertia(xx: .nan,yy: 1,zz: 1) }
        let small = try RigidBodyInertia(xx: 2e-200,yy: 3e-200,zz: 4e-200)
        #expect(routeNorm(try small.solve(SIMD3(2e-200,6e-200,12e-200))-v) < 1e-14)
    }
    @Test func torqueFreeSymmetricTopAndConservation() throws {
        let b = body(), inertia = try RigidBodyInertia(xx: 2,yy: 2,zz: 3)
        let seed = try attitude(b,omega: SIMD3(0.2,0.3,0.4)), dynamics = RigidBodyDynamics(context: try context(b,inertia: inertia))
        let propagator = try RigidBodyAttitudePropagator(initialState: seed,dynamics: dynamics,options: .preset(for: .full))
        let end = try propagator.propagation(at: epoch(50)), w = end.attitude.angularVelocityBodyRadiansPerSecond
        let phase = 0.2*50, expected = SIMD3<Double>(0.2*cos(phase)-0.3*sin(phase),0.2*sin(phase)+0.3*cos(phase),0.4)
        #expect(routeNorm(w-expected) < 1e-9)
        #expect(abs(end.kineticEnergyJoules - (try inertia.kineticEnergyJoules(for: seed.angularVelocityBodyRadiansPerSecond))) < 1e-9)
        #expect(routeNorm(end.angularMomentumReferenceKilogramSquareMetresPerSecond - (try inertia.angularMomentumBody(for: seed.angularVelocityBodyRadiansPerSecond))) < 1e-9)
        #expect(abs(Quaternion.dot(lhs: end.attitude.bodyToReference,rhs: end.attitude.bodyToReference)-1) < 1e-14)
        let reversed = try RigidBodyAttitudePropagator(initialState: end.attitude,dynamics: dynamics,options: .preset(for: .full)).attitude(at: .j2000)
        #expect(angleDifference(reversed.bodyToReference,seed.bodyToReference) < 1e-8)
        #expect(routeNorm(reversed.angularVelocityBodyRadiansPerSecond-seed.angularVelocityBodyRadiansPerSecond) < 1e-9)
        let triaxial = try context(b), asymmetric = try RigidBodyAttitudePropagator(initialState: seed,dynamics: RigidBodyDynamics(context: triaxial),options: .preset(for: .full)).propagation(at: epoch(100))
        #expect(abs(asymmetric.kineticEnergyJoules - (try triaxial.inertia.kineticEnergyJoules(for: seed.angularVelocityBodyRadiansPerSecond))) < 1e-8)
        #expect(routeNorm(asymmetric.angularMomentumReferenceKilogramSquareMetresPerSecond - (try triaxial.inertia.angularMomentumBody(for: seed.angularVelocityBodyRadiansPerSecond))) < 1e-8)
    }
    @Test func constantTorqueAndController() throws {
        let b = body(), torque = try ConstantBodyTorque(bodyFrame: b,torqueBodyNewtonMetres: SIMD3(0,0,0.8))
        let seed = try attitude(b,omega: SIMD3(0,0,0.1))
        let end = try RigidBodyAttitudePropagator(initialState: seed,dynamics: RigidBodyDynamics(context: context(b),models: [torque]),options: .preset(for: .full)).attitude(at: epoch(3))
        #expect(routeNorm(end.angularVelocityBodyRadiansPerSecond-SIMD3(0,0,0.7)) < 1e-12)
        #expect(angleDifference(end.bodyToReference,Quaternion(angle: 1.2,axis: SIMD3(0,0,1))) < 1e-9)
        let goal = try InertialHoldAttitudeProvider(referenceFrame: reference,bodyFrame: b)
        let controller = try AttitudePDController(desired: goal,proportionalGainNewtonMetres: 2,derivativeGainNewtonMetreSeconds: 3,maximumTorqueNewtonMetres: 1)
        let turned = try attitude(b,q: Quaternion(angle: .pi/2,axis: SIMD3(0,0,1)))
        let command = try controller.contribution(for: turned,spacecraft: nil,in: context(b)).torqueBodyNewtonMetres
        #expect(abs(routeNorm(command)-1) < 1e-14 && command.z < 0)
        let equivalent = try attitude(b,q: -turned.bodyToReference)
        #expect(routeNorm(try controller.contribution(for: equivalent,spacecraft: nil,in: context(b)).torqueBodyNewtonMetres-command) < 1e-14)
        let halfTurn = try attitude(b,q: Quaternion(x: 0,y: 0,z: 1,w: 0))
        let negativeHalfTurn = try attitude(b,q: Quaternion(x: 0,y: 0,z: -1,w: 0))
        let a = try controller.contribution(for: halfTurn,spacecraft: nil,in: context(b)).torqueBodyNewtonMetres
        let c = try controller.contribution(for: negativeHalfTurn,spacecraft: nil,in: context(b)).torqueBodyNewtonMetres
        #expect(a == c)
        let controlled = try RigidBodyAttitudePropagator(initialState: turned,dynamics: RigidBodyDynamics(context: context(b),models: [controller])).attitude(at: epoch(40))
        #expect(angleDifference(controlled.bodyToReference,.identity) < 1e-5)
        #expect(routeNorm(controlled.angularVelocityBodyRadiansPerSecond) < 1e-5)
    }
    @Test func gravityGradientAndMagneticActuation() throws {
        let b = body(), ctx = try context(b), state = try attitude(b), spacecraft = try translation(SIMD3(1,1,0))
        let centre = try MovingWaypointRouteTarget(provider: LambertTransferPlannerTests.Linear(frame: reference,p: .zero,v: .zero),origin: .earth,coverage: epoch(-10)...epoch(10))
        let gm = try GravitationalParameter(cubicMetresPerSecondSquared: 1,convention: .newtonian,provenance: PhysicalDataProvenance(source: "Analytic",version: "1"))
        let gradient = try GravityGradientTorque(centre: centre,gravitationalParameter: gm)
        let torque = try gradient.contribution(for: state,spacecraft: spacecraft,in: ctx).torqueBodyNewtonMetres
        #expect(abs(torque.z-3/(4*sqrt(2))) < 1e-14)
        let aligned = try gradient.contribution(for: state,spacecraft: translation(SIMD3(1,0,0)),in: ctx)
        #expect(aligned.torqueBodyNewtonMetres == .zero)
        #expect(throws: GravityModelError.singularPosition) { try gradient.contribution(for: state,spacecraft: translation(),in: ctx) }
        let field = try UniformMagneticField(frame: reference,origin: .earth,fieldTeslas: SIMD3(0,0,2e-5))
        let magnetic = try MagnetorquerTorque(fieldProvider: field,maximumDipoleAmpereSquareMetres: 2) { _,_ in SIMD3(10,0,0) }
        let magneticTorque = try magnetic.contribution(for: state,spacecraft: spacecraft,in: ctx).torqueBodyNewtonMetres
        #expect(routeNorm(magneticTorque-SIMD3(0,-4e-5,0)) < 1e-15)
        #expect(abs((magneticTorque*field.fieldTeslas).sum()) < 1e-15)
        let rotated = try attitude(b,q: Quaternion(angle: .pi/2,axis: SIMD3(1,0,0)))
        let rotatedTorque = try magnetic.contribution(for: rotated,spacecraft: spacecraft,in: ctx).torqueBodyNewtonMetres
        #expect(routeNorm(rotatedTorque-SIMD3(0,0,4e-5)) < 1e-15)
    }
    @Test func prescribedThrustAndCoupledCentredEngine() throws {
        let b = body(), seed = try translation(), hold = try InertialHoldAttitudeProvider(referenceFrame: reference,bodyFrame: b,bodyToReference: Quaternion(angle: .pi/2,axis: SIMD3(0,0,1)))
        let engine = try BodyMountedThruster(bodyFrame: b,thrustNewtons: 1,directionBody: SIMD3(1,0,0),effectiveExhaustVelocityMetresPerSecond: 100)
        let prescribed = try NumericalSpacecraftPropagator(initialState: seed,dynamics: SpacecraftDynamics(models: [engine.forceModel(attitude: hold)],
            context: context(b).referenceContext),options: .preset(for: .full)).spacecraftState(at: epoch(2))
        let initial = try CoupledSpacecraftState(spacecraft: seed,attitude: hold.attitude(at: .j2000))
        let coupled = try CoupledSpacecraftPropagator(initialState: initial,dynamics: CoupledSpacecraftDynamics(context: context(b),wrenchModels: [engine]),
            attitudeOptions: .preset(for: .full),numericalOptions: .preset(for: .full)).state(at: epoch(2))
        #expect(abs(coupled.spacecraft.mass.totalMassKilograms-9.98) < 1e-12)
        #expect(routeNorm(coupled.spacecraft.cartesian.velocity-prescribed.cartesian.velocity) < 1e-9)
        #expect(abs(coupled.spacecraft.cartesian.velocity.y-100*log(10/9.98)) < 1e-9)
        #expect(abs(coupled.spacecraft.cartesian.velocity.x) < 1e-14)
        #expect(angleDifference(coupled.attitude.bodyToReference,initial.attitude.bodyToReference) < 1e-14)
    }
    @Test func offCentreEngineCouplesRotationAndTranslation() throws {
        let b = body(), seed = try CoupledSpacecraftState(spacecraft: translation(),attitude: attitude(b))
        let engine = try BodyMountedThruster(bodyFrame: b,thrustNewtons: 1,directionBody: SIMD3(1,0,0),positionBodyMetres: SIMD3(0,1,0),effectiveExhaustVelocityMetresPerSecond: 100)
        let propagator = try CoupledSpacecraftPropagator(initialState: seed,dynamics: CoupledSpacecraftDynamics(context: context(b),wrenchModels: [engine]),
            attitudeOptions: .preset(for: .full),numericalOptions: .preset(for: .full))
        let end = try propagator.state(at: epoch(2)), q = Quaternion(angle: -0.5,axis: SIMD3(0,0,1))
        #expect(angleDifference(end.attitude.bodyToReference,q) < 1e-9)
        #expect(routeNorm(end.attitude.angularVelocityBodyRadiansPerSecond-SIMD3(0,0,-0.5)) < 1e-12)
        #expect(end.spacecraft.cartesian.velocity.y < -0.01)
        // Independent composite Simpson integral of rotated F/m; Simpson is only
        // the test oracle, not the propagation implementation.
        let n = 10000, dt = 2.0/Double(n)
        var velocity = SIMD3<Double>.zero, position = SIMD3<Double>.zero
        for i in 0...n {
            let t = Double(i)*dt, angle = -0.125*t*t, mass = 10-0.01*t
            let a = SIMD3<Double>(cos(angle),sin(angle),0)/mass
            let weight = i == 0 || i == n ? 1.0 : i.isMultiple(of: 2) ? 2.0 : 4.0
            velocity += a*weight; position += a*((2-t)*weight)
        }
        velocity *= dt/3; position *= dt/3
        #expect(routeNorm(end.spacecraft.cartesian.velocity-velocity) < 1e-9)
        #expect(routeNorm(end.spacecraft.cartesian.position-position) < 1e-9)
        #expect(try angleDifference(propagator.attitudeProvider.attitude(at: epoch(1)).bodyToReference,propagator.state(at: epoch(1)).attitude.bodyToReference) < 1e-14)
        #expect(try propagator.spacecraftProvider.spacecraftState(at: epoch(2)).mass == end.spacecraft.mass)
    }
    struct StaleAttitude: AttitudeProvider {
        let seed: AttitudeState
        var referenceFrame: Frame { seed.referenceFrame }
        var bodyFrame: Frame { seed.bodyFrame }
        func attitude(at epoch: Epoch) throws -> AttitudeState { seed }
    }
    struct BadTorqueEpoch: TorqueModel {
        func contribution(for attitude: AttitudeState, spacecraft: SpacecraftState?, in context: RotationalEvaluationContext) throws -> TorqueContribution {
            try TorqueContribution(epoch: .j2000,bodyFrame: context.bodyFrame,torqueBodyNewtonMetres: .zero)
        }
    }
    @Test func invalidContractsAndWorkLimits() throws {
        let b = body(), seed = try attitude(b,omega: SIMD3(0,0,1)), dynamics = RigidBodyDynamics(context: try context(b))
        #expect(throws: AttitudeError.invalidRotation) { try attitude(b,q: Quaternion()) }
        let scaled = try attitude(b,q: Quaternion(x: 1e300,y: 0,z: 0,w: 1e300))
        #expect(abs(Quaternion.dot(lhs: scaled.bodyToReference,rhs: scaled.bodyToReference)-1) < 1e-14)
        #expect(throws: AttitudeError.epochMismatch) { try StaleAttitude(seed: seed).checkedAttitude(at: epoch(1)) }
        let future = try UniformSpinAttitudeProvider(initialState: seed).attitude(at: epoch(1))
        #expect(throws: AttitudeError.epochMismatch) { try RigidBodyDynamics(context: context(b),models: [BadTorqueEpoch()]).derivative(at: future) }
        #expect(throws: AttitudeError.frameMismatch) { try dynamics.derivative(at: attitude(body())) }
        #expect(throws: NumericalPropagationError.evaluationLimitExceeded) {
            try RigidBodyAttitudePropagator(initialState: seed,dynamics: dynamics,options: AttitudePropagationOptions(maximumDerivativeEvaluations: 1)).attitude(at: epoch(1))
        }
        #expect(throws: NumericalPropagationError.stepLimitExceeded) {
            try RigidBodyAttitudePropagator(initialState: seed,dynamics: dynamics,options: AttitudePropagationOptions(maximumStepAttempts: 1)).attitude(at: epoch(100))
        }
        #expect(throws: NumericalPropagationError.minimumStepReached) {
            try RigidBodyAttitudePropagator(initialState: attitude(b,omega: SIMD3(0,0,1e10)),dynamics: dynamics,
                options: AttitudePropagationOptions(initialStepSeconds: 0.001,minimumStepSeconds: 0.001,maximumStepSeconds: 1)).attitude(at: epoch(1))
        }
        #expect(throws: NumericalPropagationError.invalidOptions) { try AttitudePropagationOptions(absoluteAngleToleranceRadians: 0) }
        let zero = try RigidBodyAttitudePropagator(initialState: attitude(b),dynamics: dynamics).propagation(at: .j2000)
        #expect(zero.diagnostics.acceptedSteps == 0 && zero.diagnostics.derivativeEvaluations == 1)
    }
    @Test func tighterTolerancesImprovePhysicalRotation() throws {
        let b = body(), inertia = try RigidBodyInertia(xx: 1,yy: 1,zz: 1), seed = try attitude(b,omega: SIMD3(0.2,-0.3,0.4))
        let dynamics = RigidBodyDynamics(context: try context(b,inertia: inertia)), end = try epoch(50)
        let exact = try UniformSpinAttitudeProvider(initialState: seed).attitude(at: end)
        let low = try RigidBodyAttitudePropagator(initialState: seed,dynamics: dynamics,options: .preset(for: .low)).propagation(at: end)
        let full = try RigidBodyAttitudePropagator(initialState: seed,dynamics: dynamics,options: .preset(for: .full)).propagation(at: end)
        let lowError = angleDifference(low.attitude.bodyToReference,exact.bodyToReference), fullError = angleDifference(full.attitude.bodyToReference,exact.bodyToReference)
        #expect(fullError < lowError/10)
        #expect(fullError < 1e-9)
        #expect(full.diagnostics.derivativeEvaluations > low.diagnostics.derivativeEvaluations)
    }
    @Test func coupledBackwardReplayAndDryEngine() throws {
        let b = body(), initial = try CoupledSpacecraftState(spacecraft: translation(),attitude: attitude(b))
        let engine = try BodyMountedThruster(bodyFrame: b,thrustNewtons: 1,directionBody: SIMD3(1,0,0),positionBodyMetres: SIMD3(0,1,0),effectiveExhaustVelocityMetresPerSecond: 100)
        let dynamics = try CoupledSpacecraftDynamics(context: context(b),wrenchModels: [engine])
        let end = try CoupledSpacecraftPropagator(initialState: initial,dynamics: dynamics,attitudeOptions: .preset(for: .full),numericalOptions: .preset(for: .full)).state(at: epoch(2))
        let replay = try CoupledSpacecraftPropagator(initialState: end,dynamics: dynamics,attitudeOptions: .preset(for: .full),numericalOptions: .preset(for: .full)).state(at: .j2000)
        #expect(angleDifference(replay.attitude.bodyToReference,initial.attitude.bodyToReference) < 1e-9)
        #expect(routeNorm(replay.spacecraft.cartesian.position-initial.spacecraft.cartesian.position) < 1e-9)
        #expect(routeNorm(replay.spacecraft.cartesian.velocity-initial.spacecraft.cartesian.velocity) < 1e-9)
        #expect(abs(replay.spacecraft.mass.totalMassKilograms-10) < 1e-10)
        let dry = try CoupledSpacecraftState(spacecraft: SpacecraftState(cartesian: initial.spacecraft.cartesian,
            mass: SpacecraftMassProperties(dryMassKilograms: 5)),attitude: initial.attitude)
        #expect(throws: NumericalPropagationError.propellantDepleted) { try CoupledSpacecraftPropagator(initialState: dry,dynamics: dynamics).state(at: epoch(1)) }
        let wrongField = try UniformMagneticField(frame: reference,origin: .sun,fieldTeslas: SIMD3(0,0,1))
        let magnetic = try MagnetorquerTorque(fieldProvider: wrongField,maximumDipoleAmpereSquareMetres: 1) { _,_ in SIMD3(1,0,0) }
        #expect(throws: GravityModelError.originMismatch) { try magnetic.contribution(for: initial.attitude,spacecraft: initial.spacecraft,in: context(b)) }
    }

}
