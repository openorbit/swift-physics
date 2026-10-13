// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum AttitudeError: Error, Equatable {
    case invalidRotation, invalidVector, nonInertialReference, frameMismatch, epochMismatch
    case undefinedDirection, undefinedRoll, invalidConfiguration, angularAccelerationUnavailable
    case invalidInertia, numericalRangeExceeded
}
func attitudeCross(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> SIMD3<Double> {
    SIMD3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x)
}
func attitudeFinite(_ vector: SIMD3<Double>) -> Bool { [vector.x,vector.y,vector.z].allSatisfy(\.isFinite) }
func attitudeUnit(_ vector: SIMD3<Double>) throws -> SIMD3<Double> {
    guard attitudeFinite(vector) else { throw AttitudeError.invalidVector }
    let scale = max(abs(vector.x),abs(vector.y),abs(vector.z))
    guard scale > 0 else { throw AttitudeError.undefinedDirection }
    let scaled = vector/scale
    return scaled/routeNorm(scaled)
}
func attitudeQuaternion(_ quaternion: Quaternion) throws -> Quaternion {
    let scale = max(abs(quaternion.x),abs(quaternion.y),abs(quaternion.z),abs(quaternion.w))
    guard [quaternion.x,quaternion.y,quaternion.z,quaternion.w].allSatisfy(\.isFinite), scale > 0 else { throw AttitudeError.invalidRotation }
    let scaled = quaternion/scale
    return scaled/sqrt(Quaternion.dot(lhs: scaled,rhs: scaled))
}
/// Hamilton quaternion maps body vectors into inertial reference axes.
/// Angular velocity/acceleration are physical rates relative to that reference,
/// expressed in body axes (rad/s and rad/s²), not coordinate-map rates.
public struct AttitudeState {
    public let epoch: Epoch
    public let referenceFrame: Frame
    public let bodyFrame: Frame
    public let bodyToReference: Quaternion
    public let angularVelocityBodyRadiansPerSecond: SIMD3<Double>
    public let angularAccelerationBodyRadiansPerSecondSquared: SIMD3<Double>?
    public init(epoch: Epoch, referenceFrame: Frame, bodyFrame: Frame, bodyToReference: Quaternion,
        angularVelocityBodyRadiansPerSecond: SIMD3<Double> = .zero,
        angularAccelerationBodyRadiansPerSecondSquared: SIMD3<Double>? = nil) throws {
        guard referenceFrame.pseudoInertial else { throw AttitudeError.nonInertialReference }
        guard referenceFrame !== bodyFrame else { throw AttitudeError.frameMismatch }
        guard attitudeFinite(angularVelocityBodyRadiansPerSecond), angularAccelerationBodyRadiansPerSecondSquared.map(attitudeFinite) ?? true else { throw AttitudeError.invalidVector }
        self.epoch = epoch; self.referenceFrame = referenceFrame; self.bodyFrame = bodyFrame
        self.bodyToReference = try attitudeQuaternion(bodyToReference)
        self.angularVelocityBodyRadiansPerSecond = angularVelocityBodyRadiansPerSecond
        self.angularAccelerationBodyRadiansPerSecondSquared = angularAccelerationBodyRadiansPerSecondSquared
    }
    public func referenceVector(fromBody vector: SIMD3<Double>) throws -> SIMD3<Double> {
        guard attitudeFinite(vector) else { throw AttitudeError.invalidVector }
        let result = bodyToReference.matrix*vector
        guard attitudeFinite(result) else { throw AttitudeError.numericalRangeExceeded }
        return result
    }
    public func bodyVector(fromReference vector: SIMD3<Double>) throws -> SIMD3<Double> {
        guard attitudeFinite(vector) else { throw AttitudeError.invalidVector }
        let result = bodyToReference.conj.matrix*vector
        guard attitudeFinite(result) else { throw AttitudeError.numericalRangeExceeded }
        return result
    }
    public func orientation() throws -> BodyOrientation {
        guard let alpha = angularAccelerationBodyRadiansPerSecondSquared else { throw AttitudeError.angularAccelerationUnavailable }
        return try BodyOrientation(referenceFrame: referenceFrame,epoch: epoch,rotation: bodyToReference.conj,
            angularVelocity: -angularVelocityBodyRadiansPerSecond,angularAcceleration: -alpha)
    }
}
public protocol AttitudeProvider {
    var referenceFrame: Frame { get }
    var bodyFrame: Frame { get }
    func attitude(at epoch: Epoch) throws -> AttitudeState
}
extension AttitudeProvider {
    public func checkedAttitude(at epoch: Epoch) throws -> AttitudeState {
        let state = try attitude(at: epoch)
        guard state.epoch == epoch else { throw AttitudeError.epochMismatch }
        guard state.referenceFrame === referenceFrame, state.bodyFrame === bodyFrame else { throw AttitudeError.frameMismatch }
        return state
    }
}
public struct AttitudeOrientationProvider: OrientationProvider {
    public let attitudeProvider: any AttitudeProvider
    public var referenceFrame: Frame { attitudeProvider.referenceFrame }
    public init(_ provider: any AttitudeProvider) { attitudeProvider = provider }
    public func orientation(at epoch: Epoch) throws -> BodyOrientation { try attitudeProvider.checkedAttitude(at: epoch).orientation() }
}
/// Adapts existing planetary/local-orbital orientation providers without changing
/// their reference-to-body coordinate-map convention.
public struct OrientationAttitudeProvider: AttitudeProvider {
    public let orientationProvider: any OrientationProvider
    public let bodyFrame: Frame
    public var referenceFrame: Frame { orientationProvider.referenceFrame }
    public init(orientation: any OrientationProvider, bodyFrame: Frame) { orientationProvider = orientation; self.bodyFrame = bodyFrame }
    public func attitude(at epoch: Epoch) throws -> AttitudeState {
        let orientation = try orientationProvider.orientation(at: epoch)
        guard orientation.epoch == epoch else { throw AttitudeError.epochMismatch }
        guard orientation.referenceFrame === referenceFrame else { throw AttitudeError.frameMismatch }
        return try AttitudeState(epoch: epoch,referenceFrame: referenceFrame,bodyFrame: bodyFrame,bodyToReference: orientation.rotation.conj,
            angularVelocityBodyRadiansPerSecond: -orientation.angularVelocity,
            angularAccelerationBodyRadiansPerSecondSquared: -orientation.angularAcceleration)
    }
}
public struct InertialHoldAttitudeProvider: AttitudeProvider {
    public let referenceFrame: Frame
    public let bodyFrame: Frame
    public let bodyToReference: Quaternion
    public init(referenceFrame: Frame, bodyFrame: Frame, bodyToReference: Quaternion = .identity) throws {
        let seed = try AttitudeState(epoch: .j2000,referenceFrame: referenceFrame,bodyFrame: bodyFrame,bodyToReference: bodyToReference)
        self.referenceFrame = referenceFrame; self.bodyFrame = bodyFrame; self.bodyToReference = seed.bodyToReference
    }
    public func attitude(at epoch: Epoch) throws -> AttitudeState {
        try AttitudeState(epoch: epoch,referenceFrame: referenceFrame,bodyFrame: bodyFrame,bodyToReference: bodyToReference,
            angularAccelerationBodyRadiansPerSecondSquared: .zero)
    }
}
/// Prescribed constant physical spin in body axes; no torque model is inferred.
public struct UniformSpinAttitudeProvider: AttitudeProvider {
    public let initialState: AttitudeState
    public var referenceFrame: Frame { initialState.referenceFrame }
    public var bodyFrame: Frame { initialState.bodyFrame }
    public init(initialState: AttitudeState) { self.initialState = initialState }
    public func attitude(at epoch: Epoch) throws -> AttitudeState {
        let omega = initialState.angularVelocityBodyRadiansPerSecond, speed = routeNorm(omega)
        let angle = speed * (try epoch.duration(since: initialState.epoch).seconds)
        guard angle.isFinite else { throw AttitudeError.numericalRangeExceeded }
        let rotation = speed == 0 ? Quaternion.identity : Quaternion(angle: angle.truncatingRemainder(dividingBy: 2 * .pi),axis: omega)
        return try AttitudeState(epoch: epoch,referenceFrame: referenceFrame,bodyFrame: bodyFrame,
            bodyToReference: initialState.bodyToReference*rotation,angularVelocityBodyRadiansPerSecond: omega,
            angularAccelerationBodyRadiansPerSecondSquared: .zero)
    }
}
public enum AttitudeTrackingDirection { case target, velocity, nadir }
/// Body X follows the selected direction, Z follows the projected explicit roll
/// reference, Y = Z × X. Zero direction or parallel roll axes fail explicitly.
/// Five-point angular derivatives require provider coverage at epoch ±2*step.
public struct TrackingAttitudeProvider: AttitudeProvider {
    public let referenceFrame: Frame
    public let bodyFrame: Frame
    public let derivativeStep: Duration
    public let direction: AttitudeTrackingDirection
    public let rollReference: SIMD3<Double>
    private let spacecraft: any StateProvider
    private let target: (any StateProvider)?
    public init(spacecraft: any StateProvider, bodyFrame: Frame, direction: AttitudeTrackingDirection,
        target: (any StateProvider)? = nil, rollReference: SIMD3<Double>, derivativeStep: Duration) throws {
        guard spacecraft.frame.pseudoInertial else { throw AttitudeError.nonInertialReference }
        guard target == nil || target!.frame === spacecraft.frame, bodyFrame !== spacecraft.frame else { throw AttitudeError.frameMismatch }
        guard derivativeStep.seconds.isFinite, derivativeStep > .zero,
              (derivativeStep.seconds*derivativeStep.seconds).isFinite, derivativeStep.seconds*derivativeStep.seconds > 0 else { throw AttitudeError.invalidConfiguration }
        if case .target = direction, target == nil { throw AttitudeError.invalidConfiguration }
        referenceFrame = spacecraft.frame; self.bodyFrame = bodyFrame; self.direction = direction
        self.spacecraft = spacecraft; self.target = target; self.derivativeStep = derivativeStep
        self.rollReference = try attitudeUnit(rollReference)
    }
    private func matrix(at epoch: Epoch) throws -> [Double] {
        func checked(_ provider: any StateProvider) throws -> CartesianState {
            let state = try provider.state(at: epoch)
            guard state.epoch == epoch else { throw AttitudeError.epochMismatch }
            guard state.frame === referenceFrame else { throw AttitudeError.frameMismatch }
            return state
        }
        let state = try checked(spacecraft), centre = try target.map(checked)
        let look: SIMD3<Double>
        switch direction {
        case .target, .nadir: look = (centre?.position ?? .zero)-state.position
        case .velocity: look = state.velocity-(centre?.velocity ?? .zero)
        }
        let x = try attitudeUnit(look), projected = rollReference-x*(rollReference*x).sum()
        guard routeNorm(projected) > 1e-12 else { throw AttitudeError.undefinedRoll }
        let z = try attitudeUnit(projected), y = attitudeCross(z,x)
        return [x.x,x.y,x.z,y.x,y.y,y.z,z.x,z.y,z.z]
    }
    public func attitude(at epoch: Epoch) throws -> AttitudeState {
        let h = derivativeStep.seconds, r = try matrix(at: epoch)
        let samples = try [-2.0,-1,1,2].map { try matrix(at: epoch.advanced(by: Duration(seconds: $0*h))) }
        let d = (0..<9).map { (samples[0][$0]-8*samples[1][$0]+8*samples[2][$0]-samples[3][$0])/(12*h) }
        let dd = (0..<9).map { (-samples[3][$0]+16*samples[2][$0]-30*r[$0]+16*samples[1][$0]-samples[0][$0])/(12*h*h) }
        func axial(_ derivative: [Double]) -> SIMD3<Double> {
            func entry(_ i: Int, _ j: Int) -> Double { (0..<3).reduce(0) { $0+derivative[3*i+$1]*r[3*j+$1] } }
            return SIMD3((entry(2,1)-entry(1,2))/2,(entry(0,2)-entry(2,0))/2,(entry(1,0)-entry(0,1))/2)
        }
        return try AttitudeState(epoch: epoch,referenceFrame: referenceFrame,bodyFrame: bodyFrame,
            bodyToReference: Matrix3x3(elements: r).quaternion.conj,
            angularVelocityBodyRadiansPerSecond: -axial(d),angularAccelerationBodyRadiansPerSecondSquared: -axial(dd))
    }
}
extension FrameTree {
    /// Registers the provider's exact body-frame identity with its inertial parent.
    /// As with other frame providers, do not capture the owning tree strongly.
    public func registerSpacecraftFrame(spacecraft: any StateProvider, attitude: any AttitudeProvider,
        accelerationPolicy: OriginAccelerationPolicy = .requireProvided) throws {
        guard attitude.bodyFrame.parent === attitude.referenceFrame, spacecraft.frame === attitude.referenceFrame else { throw AttitudeError.frameMismatch }
        try registerFrame(attitude.bodyFrame,provider: BodyFixedTransformProvider(centre: spacecraft,
            orientation: AttitudeOrientationProvider(attitude),accelerationPolicy: accelerationPolicy))
    }
}
