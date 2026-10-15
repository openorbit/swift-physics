//
// SPDX-License-Identifier: Apache-2.0
//
// Copyright 2026 Mattias Holm
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
import Foundation
import Math

/// Coverage-aware VSOP87E coefficient selection, with an explicit TDB model.
public struct VSOP87Ephemeris: AccuracyAwareEphemerisProvider {
  public let metadata: EphemerisMetadata
  public let accuracy: AccuracyLevel
  public let effectiveAccuracy: AccuracyLevel
  public var availableAccuracyLevels: Set<AccuracyLevel> { Set(AccuracyLevel.allCases) }
  public let coefficientSelection: VSOP87CoefficientSelection
  private let provider: VSOP87StateProvider
  public init(body: vsop87_body, coverage: ClosedRange<Epoch>, tdbModel: any TDBTimeModel, accuracy: AccuracyLevel = .full) throws {
    self.accuracy = accuracy
    guard coverage.lowerBound >= tdbModel.coverage.lowerBound,
          coverage.upperBound <= tdbModel.coverage.upperBound else { throw EphemerisError.outsideCoverage }
    let id: BodyIdentifier
    switch body.body_id {
    case .vsop87_sun: id = .sun
    case .vsop87_mercury: id = .mercury
    case .vsop87_venus: id = .venus
    case .vsop87_earth: id = .earth
    case .vsop87_mars: id = .marsBarycentre
    case .vsop87_jupiter: id = .jupiterBarycentre
    case .vsop87_saturn: id = .saturnBarycentre
    case .vsop87_uranus: id = .uranusBarycentre
    case .vsop87_neptune: id = .neptuneBarycentre
    }
    let selection = try VSOP87CoefficientSelection(body: body,coverage: coverage,tdbModel: tdbModel,accuracy: accuracy)
    coefficientSelection = selection
    effectiveAccuracy = selection.retainedTermCount == selection.originalTermCount ? .full : accuracy
    provider = VSOP87StateProvider(body: selection.selectedBody,tdbModel: tdbModel)
    metadata = try EphemerisMetadata(body: id,origin: .solarSystemBarycentre,frame: body.frame,
      coverage: coverage,coverageBasis: .callerSelected,modelTimeScale: .tdb,
      velocityMethod: .analytic,accelerationMethod: .unavailable,
      provenance: PhysicalDataProvenance(source: "VSOP87E \(effectiveAccuracy) evaluation, \(selection.retainedTermCount)/\(selection.originalTermCount) terms; \(tdbModel.provenance)",
        version: "VSOP87",reference: "https://ftp.imcce.fr/pub/ephem/planets/vsop87/"),accuracy: accuracy,effectiveAccuracy: effectiveAccuracy)
  }
  public func state(at epoch: Epoch) throws -> CartesianState {
    guard metadata.coverage.contains(epoch) else { throw EphemerisError.outsideCoverage }
    return try provider.state(at: epoch)
  }
}

/// Native identities for the legacy lunar and Pluto approximations. Retain this
/// object and register its transforms to query these states in other frames.
public struct AnalyticEphemerisFrames {
  public let moonEclipticOfDate: Frame
  public let plutoEclipticJ2000: Frame
  internal init(moon: Frame, pluto: Frame) { moonEclipticOfDate = moon; plutoEclipticJ2000 = pluto }
}

extension FrameTree {
  /// Earth and Sun states must already be expressed in this tree's ICRF frame.
  public func makeAnalyticEphemerisFrames(earth: any StateProvider, sun: any StateProvider,
    accelerationPolicy: OriginAccelerationPolicy = .requireProvided,
    derivativeStep: Duration = try! Duration(seconds: 60)) throws -> AnalyticEphemerisFrames {
    guard earth.frame === CelestialFrames.icrf, sun.frame === CelestialFrames.icrf else {
      throw EphemerisError.frameMismatch
    }
    let dated = try MeanEclipticOfDateProvider(step: derivativeStep)
    let geo = try makeGeocentricFrames(earth: earth,accelerationPolicy: accelerationPolicy)
    let lunar = try makeFrame(parent: geo.eme2000,name: "Meeus Moon mean ecliptic/equinox of date",provider: dated)
    let ecliptic = Quaternion(angle: -84381.448 * .pi/(180*3600),axis: SIMD3(1,0,0)) * CelestialFrames.icrfToEME2000
    let orientation = try FixedAxisOrientationProvider(referenceFrame: CelestialFrames.icrf,axis: SIMD3(1,0,0),
      rate: 0,initialRotation: ecliptic)
    let heliocentric = try makeBodyFixedFrame(name: "Goffin heliocentric ecliptic J2000",centre: sun,
      orientation: orientation,accelerationPolicy: accelerationPolicy)
    return AnalyticEphemerisFrames(moon: lunar,pluto: heliocentric)
  }
}

/// Adapters for the existing truncated lunar and Pluto spherical series. TT is the
/// legacy ephemeris-time argument. Position coverage is caller selected; query
/// coverage is shortened by ±2*velocityStep for the five-point Cartesian derivative.
public struct SphericalAnalyticEphemeris: AccuracyAwareEphemerisProvider {
  public enum Model: Sendable { case meeusMoon, goffinPluto }
  public let metadata: EphemerisMetadata
  public let accuracy: AccuracyLevel
  public let effectiveAccuracy: AccuracyLevel = .full
  public var availableAccuracyLevels: Set<AccuracyLevel> { [.full] }
  public let model: Model
  public let positionCoverage: ClosedRange<Epoch>
  public let velocityStep: Duration
  public init(model: Model, frames: AnalyticEphemerisFrames, positionCoverage: ClosedRange<Epoch>,
              velocityStep: Duration, accuracy: AccuracyLevel = .full) throws {
    self.accuracy = accuracy
    guard velocityStep > .zero else { throw EphemerisError.invalidDerivativeMethod }
    let margin = try velocityStep.adding(velocityStep)
    let lower = try positionCoverage.lowerBound.advanced(by: margin)
    let upper = try positionCoverage.upperBound.advanced(by: margin.negated())
    guard lower <= upper else { throw EphemerisError.outsideCoverage }
    self.model = model; self.positionCoverage = positionCoverage; self.velocityStep = velocityStep
    let moon = model == .meeusMoon
    metadata = try EphemerisMetadata(body: moon ? .moon : .plutoBarycentre,origin: moon ? .earth : .sun,
      frame: moon ? frames.moonEclipticOfDate : frames.plutoEclipticJ2000,coverage: lower...upper,
      coverageBasis: .callerSelected,modelTimeScale: .tt,
      velocityMethod: .finiteDifference(step: velocityStep),accelerationMethod: .unavailable,
      provenance: PhysicalDataProvenance(source: moon ? "Meeus chapter 47 truncated lunar series, including constant light-time longitude term; legacy elp2000_82b" : "Meeus/Goffin Pluto spherical series; legacy goffin2000",
        version: "built-in legacy coefficients; TT adapter"),accuracy: accuracy,effectiveAccuracy: effectiveAccuracy)
  }
  private func position(at epoch: Epoch) throws -> SIMD3<Double> {
    guard positionCoverage.contains(epoch) else { throw EphemerisError.outsideCoverage }
    let t = try centuriesTT(epoch)
    let spherical = model == .meeusMoon ? meeusMoonPosition(centuries: t) : goffinPosition(centuries: t)
    let distance = spherical.z * (model == .meeusMoon ? 1 : CartesianUnits.metresPerAstronomicalUnit)
    return SIMD3(distance*cos(spherical.y)*cos(spherical.x),distance*cos(spherical.y)*sin(spherical.x),distance*sin(spherical.y))
  }
  public func state(at epoch: Epoch) throws -> CartesianState {
    guard metadata.coverage.contains(epoch) else { throw EphemerisError.outsideCoverage }
    let h = velocityStep.seconds
    let p = try position(at: epoch)
    let samples = try [-2.0,-1,1,2].map { try position(at: epoch.advanced(by: Duration(seconds: $0*h))) }
    let velocity = (samples[0]-8*samples[1]+8*samples[2]-samples[3])/(12*h)
    return try CartesianState(epoch: epoch,frame: frame,position: p,velocity: velocity)
  }
}

internal func centuriesTT(_ epoch: Epoch) throws -> Double {
  let date = try epoch.julianDate(in: .tt)
  return (Double(date.wholeDay-2451545)+date.fractionalDay)/36525
}

/// IAU 1976 precession and IAU 1980 mean obliquity, compatible with the legacy series.
/// Native Swift; rotation convention follows ERFA pmat76. No nutation is applied.
internal struct MeanEclipticOfDateProvider: TransformProvider {
  let step: Duration
  init(step: Duration) throws {
    guard step > .zero else { throw EphemerisError.invalidDerivativeMethod }
    self.step = step
  }
  func rotation(at epoch: Epoch) throws -> Quaternion {
    let t = try centuriesTT(epoch), arcsecond = Double.pi/(180*3600)
    let zeta = (2306.2181*t+0.30188*t*t+0.017998*t*t*t)*arcsecond
    let z = (2306.2181*t+1.09468*t*t+0.018203*t*t*t)*arcsecond
    let theta = (2004.3109*t-0.42665*t*t-0.041833*t*t*t)*arcsecond
    let eps = (84381.448-46.815*t-0.00059*t*t+0.001813*t*t*t)*arcsecond
    return Quaternion(angle: -eps,axis: SIMD3(1,0,0)) * Quaternion(angle: z,axis: SIMD3(0,0,1))
      * Quaternion(angle: -theta,axis: SIMD3(0,1,0)) * Quaternion(angle: zeta,axis: SIMD3(0,0,1))
  }
  func transform(at epoch: Epoch) throws -> FrameTransform {
    let q = try rotation(at: epoch), r = q.matrix, h = step.seconds
    let samples = try [-2.0,-1,1,2].map { try rotation(at: epoch.advanced(by: Duration(seconds: $0*h))).matrix }
    var d = [Double](), dd = [Double]()
    for i in 0..<3 { for j in 0..<3 {
      d.append((samples[0][i,j]-8*samples[1][i,j]+8*samples[2][i,j]-samples[3][i,j])/(12*h))
      dd.append((-samples[3][i,j]+16*samples[2][i,j]-30*r[i,j]+16*samples[1][i,j]-samples[0][i,j])/(12*h*h))
    }}
    func axial(_ v: [Double]) -> SIMD3<Double> {
      func e(_ i: Int,_ j: Int) -> Double { (0..<3).reduce(0) { $0+v[3*i+$1]*r[j,$1] } }
      return SIMD3((e(2,1)-e(1,2))/2,(e(0,2)-e(2,0))/2,(e(1,0)-e(0,1))/2)
    }
    return try FrameTransform(epoch: epoch,rotation: q,angularVelocity: axial(d),angularAcceleration: axial(dd))
  }
}
