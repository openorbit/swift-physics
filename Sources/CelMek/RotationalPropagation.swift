// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

/// Native Dormand–Prince local-error settings. Quaternion error is a physical
/// rotation angle, invariant under q/-q; no global accuracy bound is implied.
public struct AttitudePropagationOptions: Equatable {
    public let absoluteAngleToleranceRadians: Double
    public let absoluteAngularVelocityToleranceRadiansPerSecond: Double
    public let relativeTolerance: Double
    public let initialStepSeconds, minimumStepSeconds, maximumStepSeconds: Double
    public let maximumRotationPerStepRadians: Double
    public let maximumStepAttempts, maximumDerivativeEvaluations: Int
    public init(absoluteAngleToleranceRadians: Double = 1e-9, absoluteAngularVelocityToleranceRadiansPerSecond: Double = 1e-10,
        relativeTolerance: Double = 1e-9, initialStepSeconds: Double = 0.1, minimumStepSeconds: Double = 1e-8,
        maximumStepSeconds: Double = 10, maximumRotationPerStepRadians: Double = 0.25,
        maximumStepAttempts: Int = 100_000, maximumDerivativeEvaluations: Int = 700_001) throws {
        guard [absoluteAngleToleranceRadians,absoluteAngularVelocityToleranceRadiansPerSecond,initialStepSeconds,
               minimumStepSeconds,maximumStepSeconds,maximumRotationPerStepRadians].allSatisfy({ $0.isFinite && $0 > 0 }),
              relativeTolerance.isFinite, relativeTolerance >= 0, minimumStepSeconds <= initialStepSeconds,
              initialStepSeconds <= maximumStepSeconds, maximumRotationPerStepRadians <= .pi/2,
              maximumStepAttempts > 0, maximumDerivativeEvaluations > 0 else { throw NumericalPropagationError.invalidOptions }
        self.absoluteAngleToleranceRadians = absoluteAngleToleranceRadians
        self.absoluteAngularVelocityToleranceRadiansPerSecond = absoluteAngularVelocityToleranceRadiansPerSecond
        self.relativeTolerance = relativeTolerance; self.initialStepSeconds = initialStepSeconds
        self.minimumStepSeconds = minimumStepSeconds; self.maximumStepSeconds = maximumStepSeconds
        self.maximumRotationPerStepRadians = maximumRotationPerStepRadians
        self.maximumStepAttempts = maximumStepAttempts; self.maximumDerivativeEvaluations = maximumDerivativeEvaluations
    }
    public static func preset(for accuracy: AccuracyLevel) -> Self {
        let angle, velocity, relative: Double
        switch accuracy {
        case .low: (angle,velocity,relative) = (1e-4,1e-5,1e-5)
        case .balanced: (angle,velocity,relative) = (1e-7,1e-8,1e-7)
        case .high: (angle,velocity,relative) = (1e-9,1e-10,1e-9)
        case .full: (angle,velocity,relative) = (1e-11,1e-12,1e-11)
        }
        return try! Self(absoluteAngleToleranceRadians: angle,absoluteAngularVelocityToleranceRadiansPerSecond: velocity,relativeTolerance: relative)
    }
}
func attitudeValues(_ state: AttitudeState) -> [Double] {
    let q = state.bodyToReference, w = state.angularVelocityBodyRadiansPerSecond
    return [q.x,q.y,q.z,q.w,w.x,w.y,w.z]
}
func attitudeFromValues(_ values: [Double], offset: Int, epoch: Epoch, context: RotationalEvaluationContext,
    acceleration: SIMD3<Double>? = nil) throws -> AttitudeState {
    try AttitudeState(epoch: epoch,referenceFrame: context.referenceFrame,bodyFrame: context.bodyFrame,
        bodyToReference: Quaternion(x: values[offset],y: values[offset+1],z: values[offset+2],w: values[offset+3]),
        angularVelocityBodyRadiansPerSecond: SIMD3(values[offset+4],values[offset+5],values[offset+6]),
        angularAccelerationBodyRadiansPerSecondSquared: acceleration)
}
func rawQuaternionRate(_ values: [Double], offset: Int) -> [Double] {
    let q = Quaternion(x: values[offset],y: values[offset+1],z: values[offset+2],w: values[offset+3])
    let w = Quaternion(x: values[offset+4],y: values[offset+5],z: values[offset+6],w: 0)
    let d = 0.5*(q*w)
    return [d.x,d.y,d.z,d.w]
}
private struct RotationIntegrationResult {
    let values: [Double]
    let derivative: [Double]
    let diagnostics: NumericalPropagationDiagnostics
}
/// Shared native vector integrator for rotation alone or tightly coupled motion.
/// Accepted states normalize q; stage derivatives use raw q for qdot, unit q only
/// for physical vector transforms. No Euler-angle chart or quaternion sign reset.
private func integrateRotation(initial: [Double], from start: Epoch, to target: Epoch, quaternionOffset: Int,
    attitudeOptions: AttitudePropagationOptions, translationOptions: NumericalPropagationOptions?,
    derivative: ([Double],Epoch) throws -> [Double]) throws -> RotationIntegrationResult {
    var current = start, y = initial, accepted = 0, rejected = 0, evaluations = 0
    let forward = target >= start
    let minimumStep = max(attitudeOptions.minimumStepSeconds,translationOptions?.minimumStepSeconds ?? attitudeOptions.minimumStepSeconds)
    let maximumStep = min(attitudeOptions.maximumStepSeconds,translationOptions?.maximumStepSeconds ?? attitudeOptions.maximumStepSeconds)
    var nextStep = min(attitudeOptions.initialStepSeconds,translationOptions?.initialStepSeconds ?? attitudeOptions.initialStepSeconds)
    guard nextStep >= minimumStep, nextStep <= maximumStep else { throw NumericalPropagationError.invalidOptions }
    let attempts = min(attitudeOptions.maximumStepAttempts,translationOptions?.maximumStepAttempts ?? attitudeOptions.maximumStepAttempts)
    let evaluationLimit = min(attitudeOptions.maximumDerivativeEvaluations,translationOptions?.maximumDerivativeEvaluations ?? attitudeOptions.maximumDerivativeEvaluations)
    var smallest: Double?, largest: Double?, largestError = 0.0
    func query(_ values: [Double], _ epoch: Epoch) throws -> [Double] {
        guard evaluations < evaluationLimit else { throw NumericalPropagationError.evaluationLimitExceeded }
        evaluations += 1
        let result = try derivative(values,epoch)
        guard result.count == values.count, result.allSatisfy(\.isFinite) else { throw AttitudeError.numericalRangeExceeded }
        return result
    }
    func quaternion(_ values: [Double]) -> Quaternion {
        Quaternion(x: values[quaternionOffset],y: values[quaternionOffset+1],z: values[quaternionOffset+2],w: values[quaternionOffset+3])
    }
    let c = [0.0,1.0/5,3.0/10,4.0/5,8.0/9,1,1]
    let a: [[Double]] = [[],[1.0/5],[3.0/40,9.0/40],[44.0/45,-56.0/15,32.0/9],
        [19372.0/6561,-25360.0/2187,64448.0/6561,-212.0/729],
        [9017.0/3168,-355.0/33,46732.0/5247,49.0/176,-5103.0/18656],
        [35.0/384,0,500.0/1113,125.0/192,-2187.0/6784,11.0/84]]
    let weights = [71.0/57600,0,-71.0/16695,71.0/1920,-17253.0/339200,22.0/525,-1.0/40]
    while current != target {
        guard accepted+rejected < attempts else { throw NumericalPropagationError.stepLimitExceeded }
        let remaining = abs(try target.duration(since: current).seconds)
        let w = SIMD3(y[quaternionOffset+4],y[quaternionOffset+5],y[quaternionOffset+6]), speed = routeNorm(w)
        guard speed.isFinite else { throw AttitudeError.numericalRangeExceeded }
        let rotationStep = speed == 0 ? maximumStep : attitudeOptions.maximumRotationPerStepRadians/speed
        guard rotationStep >= minimumStep || remaining <= rotationStep else { throw NumericalPropagationError.minimumStepReached }
        let magnitude = min(nextStep,remaining,rotationStep)
        let end = magnitude >= remaining ? target : try current.advanced(by: Duration(seconds: forward ? magnitude : -magnitude))
        guard forward ? end > current : end < current else { throw NumericalPropagationError.epochResolutionExceeded }
        let h = try end.duration(since: current).seconds
        var stages: [[Double]] = [], candidate = y, invalidFuel = false
        for stage in 0..<7 {
            var values = y
            for i in y.indices { for j in 0..<stage { values[i] += h*a[stage][j]*stages[j][i] } }
            guard values.allSatisfy(\.isFinite) else { throw AttitudeError.numericalRangeExceeded }
            if translationOptions != nil && values[6] < 0 { invalidFuel = true; break }
            let date = stage == 0 ? current : c[stage] == 1 ? end : try current.advanced(by: Duration(seconds: c[stage]*h))
            if c[stage] > 0 && c[stage] < 1 && !(forward ? date > current && date < end : date < current && date > end) { throw NumericalPropagationError.epochResolutionExceeded }
            stages.append(try query(values,date))
            if stage == 6 { candidate = values }
        }
        var error = Double.infinity
        if !invalidFuel {
            var difference = y.map { _ in 0.0 }
            for i in y.indices { for stage in 0..<7 { difference[i] += h*weights[stage]*stages[stage][i] } }
            let low = zip(candidate,difference).map { $0.0-$0.1 }
            let highQ = try attitudeQuaternion(quaternion(candidate)), lowQ = try attitudeQuaternion(quaternion(low))
            let delta = Quaternion.dot(lhs: highQ,rhs: lowQ) >= 0 ? highQ-lowQ : highQ+lowQ
            let angleError = 4*asin(min(1,sqrt(Quaternion.dot(lhs: delta,rhs: delta))/2))
            let norm = sqrt(Quaternion.dot(lhs: quaternion(candidate),rhs: quaternion(candidate)))
            let angleScale = attitudeOptions.absoluteAngleToleranceRadians+attitudeOptions.relativeTolerance
            error = max(angleError,abs(norm-1))/angleScale
            for i in (quaternionOffset+4)..<(quaternionOffset+7) {
                let scale = attitudeOptions.absoluteAngularVelocityToleranceRadiansPerSecond+attitudeOptions.relativeTolerance*max(abs(y[i]),abs(candidate[i]))
                error = max(error,abs(difference[i])/scale)
            }
            if let options = translationOptions {
                for i in 0..<7 {
                    let absolute = i < 3 ? options.absolutePositionToleranceMetres : i < 6 ? options.absoluteVelocityToleranceMetresPerSecond : options.absolutePropellantToleranceKilograms
                    let scale = absolute+options.relativeTolerance*max(abs(y[i]),abs(candidate[i]))
                    error = max(error,abs(difference[i])/scale)
                }
            }
            var stageSpeed = speed
            for stage in stages.indices {
                var stageValues = y
                for component in (quaternionOffset+4)..<(quaternionOffset+7) {
                    for j in 0..<stage { stageValues[component] += h*a[stage][j]*stages[j][component] }
                }
                stageSpeed = max(stageSpeed,routeNorm(SIMD3(stageValues[quaternionOffset+4],stageValues[quaternionOffset+5],stageValues[quaternionOffset+6])))
            }
            error = max(error,abs(h)*stageSpeed/attitudeOptions.maximumRotationPerStepRadians)
            guard error.isFinite else { throw AttitudeError.numericalRangeExceeded }
        }
        if error <= 1 {
            let q = try attitudeQuaternion(quaternion(candidate))
            candidate[quaternionOffset] = q.x; candidate[quaternionOffset+1] = q.y
            candidate[quaternionOffset+2] = q.z; candidate[quaternionOffset+3] = q.w
            y = candidate; current = end; accepted += 1
            smallest = min(smallest ?? abs(h),abs(h)); largest = max(largest ?? abs(h),abs(h)); largestError = max(largestError,error)
            nextStep = max(minimumStep,min(maximumStep,abs(h)*(error == 0 ? 5 : min(5,max(0.2,0.9*pow(error,-0.2))))))
        } else {
            rejected += 1
            guard magnitude > minimumStep else { throw invalidFuel ? NumericalPropagationError.propellantDepleted : NumericalPropagationError.minimumStepReached }
            nextStep = max(minimumStep,magnitude*(invalidFuel ? 0.5 : max(0.1,min(0.5,0.9*pow(error,-0.2)))))
        }
    }
    let final = try query(y,current)
    return RotationIntegrationResult(values: y,derivative: final,diagnostics: NumericalPropagationDiagnostics(
        acceptedSteps: accepted,rejectedSteps: rejected,derivativeEvaluations: evaluations,
        smallestAcceptedStepSeconds: smallest,largestAcceptedStepSeconds: largest,maximumAcceptedNormalizedError: largestError))
}
public struct AttitudePropagationResult {
    public let attitude: AttitudeState
    public let diagnostics: NumericalPropagationDiagnostics
    public let kineticEnergyJoules: Double
    public let angularMomentumReferenceKilogramSquareMetresPerSecond: SIMD3<Double>
}
public struct RigidBodyAttitudePropagator: AttitudeProvider {
    public let initialState: AttitudeState
    public let dynamics: RigidBodyDynamics
    public let options: AttitudePropagationOptions
    private let spacecraft: (any SpacecraftStateProvider)?
    public var referenceFrame: Frame { initialState.referenceFrame }
    public var bodyFrame: Frame { initialState.bodyFrame }
    public init(initialState: AttitudeState, dynamics: RigidBodyDynamics,
        spacecraft: (any SpacecraftStateProvider)? = nil, options: AttitudePropagationOptions = .preset(for: .high)) throws {
        guard initialState.referenceFrame === dynamics.context.referenceFrame, initialState.bodyFrame === dynamics.context.bodyFrame,
              spacecraft == nil || spacecraft!.frame === initialState.referenceFrame else { throw AttitudeError.frameMismatch }
        self.initialState = initialState; self.dynamics = dynamics; self.spacecraft = spacecraft; self.options = options
    }
    public func attitude(at epoch: Epoch) throws -> AttitudeState { try propagation(at: epoch).attitude }
    public func propagation(at epoch: Epoch) throws -> AttitudePropagationResult {
        let result = try integrateRotation(initial: attitudeValues(initialState),from: initialState.epoch,to: epoch,quaternionOffset: 0,
            attitudeOptions: options,translationOptions: nil) { values,date in
            let attitude = try attitudeFromValues(values,offset: 0,epoch: date,context: dynamics.context)
            let translation = try spacecraft?.spacecraftState(at: date)
            let d = try dynamics.derivative(at: attitude,spacecraft: translation), alpha = d.angularAccelerationBodyRadiansPerSecondSquared
            return rawQuaternionRate(values,offset: 0)+[alpha.x,alpha.y,alpha.z]
        }
        let state = try attitudeFromValues(result.values,offset: 0,epoch: epoch,context: dynamics.context,
            acceleration: SIMD3(result.derivative[4],result.derivative[5],result.derivative[6]))
        return try AttitudePropagationResult(attitude: state,diagnostics: result.diagnostics,
            kineticEnergyJoules: dynamics.context.inertia.kineticEnergyJoules(for: state.angularVelocityBodyRadiansPerSecond),
            angularMomentumReferenceKilogramSquareMetresPerSecond: state.referenceVector(fromBody:
                dynamics.context.inertia.angularMomentumBody(for: state.angularVelocityBodyRadiansPerSecond)))
    }
}
public struct CoupledSpacecraftPropagationResult {
    public let state: CoupledSpacecraftState
    public let diagnostics: NumericalPropagationDiagnostics
}
/// Immutable replay, forward/backward, of translation, fuel and rotation at the
/// same RK stages. Continuous models only; split/reseed explicitly at command
/// discontinuities. Fuel exhaustion fails rather than silently firing while dry.
public struct CoupledSpacecraftPropagator {
    public let initialState: CoupledSpacecraftState
    public let dynamics: CoupledSpacecraftDynamics
    public let attitudeOptions: AttitudePropagationOptions
    public let numericalOptions: NumericalPropagationOptions
    public init(initialState: CoupledSpacecraftState, dynamics: CoupledSpacecraftDynamics,
        attitudeOptions: AttitudePropagationOptions = .preset(for: .high), numericalOptions: NumericalPropagationOptions = .preset(for: .high)) throws {
        guard initialState.attitude.bodyFrame === dynamics.context.bodyFrame, initialState.spacecraft.frame === dynamics.context.referenceFrame else { throw AttitudeError.frameMismatch }
        self.initialState = initialState; self.dynamics = dynamics; self.attitudeOptions = attitudeOptions; self.numericalOptions = numericalOptions
    }
    public func state(at epoch: Epoch) throws -> CoupledSpacecraftState { try propagation(at: epoch).state }
    public func propagation(at epoch: Epoch) throws -> CoupledSpacecraftPropagationResult {
        let s = initialState.spacecraft, p = s.cartesian.position, v = s.cartesian.velocity
        let initial = [p.x,p.y,p.z,v.x,v.y,v.z,s.mass.propellantMassKilograms]+attitudeValues(initialState.attitude)
        func state(_ values: [Double], at date: Epoch, derivative: [Double]? = nil) throws -> CoupledSpacecraftState {
            let translation = try SpacecraftState(cartesian: CartesianState(epoch: date,frame: dynamics.context.referenceFrame,
                position: SIMD3(values[0],values[1],values[2]),velocity: SIMD3(values[3],values[4],values[5]),
                acceleration: derivative.map { SIMD3($0[3],$0[4],$0[5]) }),
                mass: SpacecraftMassProperties(dryMassKilograms: s.mass.dryMassKilograms,propellantMassKilograms: values[6]))
            let attitude = try attitudeFromValues(values,offset: 7,epoch: date,context: dynamics.context,
                acceleration: derivative.map { SIMD3($0[11],$0[12],$0[13]) })
            return try CoupledSpacecraftState(spacecraft: translation,attitude: attitude)
        }
        let result = try integrateRotation(initial: initial,from: initialState.epoch,to: epoch,quaternionOffset: 7,
            attitudeOptions: attitudeOptions,translationOptions: numericalOptions) { values,date in
            let local = try state(values,at: date), d = try dynamics.derivative(at: local)
            let velocity = d.translation.positionRateMetresPerSecond, acceleration = d.translation.velocityRateMetresPerSecondSquared
            let alpha = d.rotation.angularAccelerationBodyRadiansPerSecondSquared
            return [velocity.x,velocity.y,velocity.z,acceleration.x,acceleration.y,acceleration.z,d.translation.propellantMassRateKilogramsPerSecond]
                + rawQuaternionRate(values,offset: 7)+[alpha.x,alpha.y,alpha.z]
        }
        return CoupledSpacecraftPropagationResult(state: try state(result.values,at: epoch,derivative: result.derivative),diagnostics: result.diagnostics)
    }
    /// Separate views avoid ambiguous StateProvider/AttitudeProvider overloads.
    public var attitudeProvider: CoupledAttitudeProvider { CoupledAttitudeProvider(propagator: self) }
    public var spacecraftProvider: CoupledTranslationProvider { CoupledTranslationProvider(propagator: self) }
}
public struct CoupledAttitudeProvider: AttitudeProvider {
    public let propagator: CoupledSpacecraftPropagator
    public var referenceFrame: Frame { propagator.initialState.attitude.referenceFrame }
    public var bodyFrame: Frame { propagator.initialState.attitude.bodyFrame }
    public func attitude(at epoch: Epoch) throws -> AttitudeState { try propagator.state(at: epoch).attitude }
}
public struct CoupledTranslationProvider: SpacecraftStateProvider {
    public let propagator: CoupledSpacecraftPropagator
    public var frame: Frame { propagator.initialState.spacecraft.frame }
    public func state(at epoch: Epoch) throws -> CartesianState { try spacecraftState(at: epoch).cartesian }
    public func spacecraftState(at epoch: Epoch) throws -> SpacecraftState { try propagator.state(at: epoch).spacecraft }
}
