// SPDX-License-Identifier: Apache-2.0
import Math

public enum IAUOrientationError: Error, Equatable {
  case unknownBody, legacyModelRequiresOptIn, outsideCoverage
}

/// IAU report orientation with ICRF reference axes, as in NAIF's generic PCK.
/// Coverage is an explicit caller-selected interval, not a library-certified
/// physical accuracy bound. Earth/Moon legacy approximations require opt-in.
public struct IAUOrientationProvider: OrientationProvider {
  public let body: String
  public let coverage: ClosedRange<Epoch>
  public let referenceFrame: Frame = CelestialFrames.icrf
  public let usesLegacyModel: Bool
  public let tdbModel: (any TDBTimeModel)?
  public static var availableBodies: [String] { iauRotationModels.keys.sorted() }

  public init(body: String, coverage: ClosedRange<Epoch>, allowLegacyModel: Bool = false, tdbModel: (any TDBTimeModel)? = nil) throws {
    guard let model = iauRotationModels[body] else { throw IAUOrientationError.unknownBody }
    let legacy = model.provenance == .legacy2009
    guard !legacy || allowLegacyModel else { throw IAUOrientationError.legacyModelRequiresOptIn }
    self.tdbModel = tdbModel
    self.body = body
    self.coverage = coverage
    usesLegacyModel = legacy
  }

  public func orientation(at epoch: Epoch) throws -> BodyOrientation {
    guard coverage.contains(epoch) else { throw IAUOrientationError.outsideCoverage }
    let sample = try tdbModel.map { try epoch.tdbCorrection(using: $0) }
    let date = try sample.map { try epoch.tdbDate(using: $0) } ?? epoch.julianDate(in: .tdb)
    let days = Double(date.wholeDay - 2_451_545) + date.fractionalDay
    let elements = try IAURotModel().elements(body: body, tdbDaysSinceJ2000: days)
    let rate = sample?.rate ?? epoch.coordinateTimeRate(in: .tdb)
    let acceleration = sample?.acceleration ?? epoch.coordinateTimeAcceleration(in: .tdb)
    func angleTransform(_ angle: Double, _ first: Double, _ second: Double,
                        axis: SIMD3<Double>) throws -> FrameTransform {
      try FrameTransform(epoch: epoch,
        rotation: Quaternion(angle: angle.truncatingRemainder(dividingBy: 2 * .pi), axis: axis),
        angularVelocity: axis * (first * rate),
        angularAcceleration: axis * (second * rate * rate + first * acceleration))
    }
    let ra = elements.rightAscension, dec = elements.declination, w = elements.primeMeridian
    let node = try angleTransform(-(ra.radians + .pi/2), -ra.rate, -ra.acceleration, axis: SIMD3(0, 0, 1))
    let tilt = try angleTransform(dec.radians - .pi/2, dec.rate, dec.acceleration, axis: SIMD3(1, 0, 0))
    let spin = try angleTransform(-w.radians, -w.rate, -w.acceleration, axis: SIMD3(0, 0, 1))
    let transform = try node.followed(by: tilt).followed(by: spin)
    return try BodyOrientation(referenceFrame: referenceFrame, epoch: epoch,
      rotation: transform.rotation, angularVelocity: transform.angularVelocity,
      angularAcceleration: transform.angularAcceleration)
  }
}
