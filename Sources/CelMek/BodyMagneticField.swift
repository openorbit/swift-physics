// SPDX-License-Identifier: Apache-2.0
import Math
import MagneticFields

public extension BodyIdentifier {
    /// Preset intrinsic field. Missing means no supplied model, never a zero total field.
    var magneticModel: BodyMagneticField? {
        guard case .naif(let code) = self else { return nil }
        return BodyMagneticField.model(naifBodyCode: code)
    }
}

public struct BodyMagneticSample {
    public let epoch: Epoch
    /// Body-fixed components in the provider's documented axes, in tesla.
    public let vectorTesla: SIMD3<Double>
    public init(epoch: Epoch, vectorTesla: SIMD3<Double>) throws {
        let validated = try MagneticFluxDensity(vectorTesla: vectorTesla)
        self.epoch = epoch; self.vectorTesla = validated.vectorTesla
    }
}

/// Body identity and epoch contract independent of the concrete field evaluator.
public protocol BodyMagneticProvider {
    var bodyIdentifier: BodyIdentifier { get }
    var metadata: MagneticFieldMetadata { get }
    func magneticField(at epoch: Epoch, displacementMetres: SIMD3<Double>) throws -> BodyMagneticSample
}

/// Labels a pluggable field model. Calendar conversion is injected explicitly:
/// future WMM/IGRF models require calendar decimal years, not Julian-year offsets.
public struct ModelBodyMagneticField: BodyMagneticProvider {
    public let bodyIdentifier: BodyIdentifier
    public let model: any MagneticFieldModel
    public var metadata: MagneticFieldMetadata { model.metadata }
    private let decimalYear: (Epoch) throws -> Double
    public init(bodyIdentifier: BodyIdentifier, model: any MagneticFieldModel,
                decimalYear: @escaping (Epoch) throws -> Double) throws {
        try bodyIdentifier.validate()
        if let identified = model as? BodyMagneticField,
           bodyIdentifier != .naif(identified.body.rawValue) { throw BodyModelError.bodyMismatch }
        self.bodyIdentifier = bodyIdentifier; self.model = model; self.decimalYear = decimalYear
    }
    public func magneticField(at epoch: Epoch, displacementMetres: SIMD3<Double>) throws -> BodyMagneticSample {
        let query = try MagneticFieldQuery(positionMetres: displacementMetres,decimalYear: decimalYear(epoch))
        let field = try model.checkedField(at: query)
        return try BodyMagneticSample(epoch: epoch, vectorTesla: field.vectorTesla)
    }
}

/// Rotates body-centred positions into model axes and B back to reference axes.
/// Caller must ensure the orientation axes agree with the magnetic model metadata.
/// Magnetic fields transform as vectors; no angular-velocity or origin-velocity terms.
public struct OrientedBodyMagneticField: MagneticFieldProvider {
    public let provider: any BodyMagneticProvider
    public let orientation: any BodyOrientationProvider
    public init(provider: any BodyMagneticProvider, orientation: any BodyOrientationProvider) throws {
        guard provider.bodyIdentifier == orientation.bodyIdentifier else { throw BodyModelError.bodyMismatch }
        self.provider = provider; self.orientation = orientation
    }
    /// Spacecraft state positions must be relative to this provider's body centre.
    public func field(at state: SpacecraftState) throws -> MagneticFieldSample {
        try magneticField(at: state.epoch,displacementMetres: state.cartesian.position,in: state.frame)
    }
    public func magneticField(at epoch: Epoch, displacementMetres: SIMD3<Double>,
                              in frame: Frame) throws -> MagneticFieldSample {
        guard frame === orientation.referenceFrame else { throw BodyFixedFrameError.frameMismatch }
        let axes = try orientation.orientation(at: epoch)
        guard axes.epoch == epoch else { throw FrameTransformError.epochMismatch }
        guard axes.referenceFrame === frame else { throw BodyFixedFrameError.frameMismatch }
        let field = try provider.magneticField(at: epoch,displacementMetres: axes.rotation.matrix*displacementMetres)
        guard field.epoch == epoch else { throw FrameTransformError.epochMismatch }
        let validated = try MagneticFluxDensity(vectorTesla: axes.rotation.conj.matrix*field.vectorTesla)
        return try MagneticFieldSample(epoch: epoch,frame: frame,origin: provider.bodyIdentifier,fieldTeslas: validated.vectorTesla)
    }
}
