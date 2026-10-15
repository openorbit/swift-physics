// Swift adaptation of ERFA 2.0.1: Copyright (C) 2013-2023 NumFOCUS.
// BSD terms and SOFA heritage: Sources/ERFASwift/LICENSE.txt.
import Foundation
import Math

/// Native Swift IAU 2006/2000A numerical model. All angles are radians.
package enum EarthOrientationModel {
  package enum Edge: Sendable {
    case celestialIntermediate, earthRotation, polarMotion, trueOfDate, apparentSidereal
  }
  package enum PoleReference: Sendable { case iau2000A, iau2006A }
  package enum EvaluationError: Error { case invalidParameters, correctionDidNotConverge }

  private static let arcsecond = Double.pi / (180 * 3600)
  private static let turn = 2 * Double.pi
  private static func polynomial(_ t: Double, _ coefficients: [Double]) -> Double {
    coefficients.reversed().reduce(0) { $0*t + $1 }
  }
  private static func argument(_ t: Double, _ coefficients: [Double]) -> Double {
    polynomial(t, coefficients).truncatingRemainder(dividingBy: 1296000) * arcsecond
  }
  private static func longitude(_ t: Double, _ a: Double, _ b: Double) -> Double {
    (a+b*t).truncatingRemainder(dividingBy: turn)
  }
  private static func fundamentalArguments(_ t: Double) -> [Double] {
    [argument(t, [485868.249036,1717915923.2178,31.8792,0.051635,-0.00024470]),
     argument(t, [1287104.793048,129596581.0481,-0.5532,0.000136,-0.00001149]),
     argument(t, [335779.526232,1739527262.8478,-12.7512,-0.001037,0.00000417]),
     argument(t, [1072260.703692,1602961601.2090,-6.3706,0.006593,-0.00003169]),
     argument(t, [450160.398036,-6962890.5431,7.4722,0.007702,-0.00005939]),
     longitude(t,3.176146697,1021.3285546211), longitude(t,1.753470314,628.3075849991),
     (0.024381750+0.00000538691*t)*t]
  }
  private static func phase(_ row: [Double], _ arguments: [Double]) -> Double {
    arguments.indices.reduce(0) { $0 + row[$1]*arguments[$1] }
  }
  private static func nutation(_ t: Double) -> (longitude: Double, obliquity: Double) {
    var args = Array(fundamentalArguments(t).prefix(5))
    // MHB2000 uses slightly different solar anomaly and elongation expressions.
    args[1] = argument(t,[1287104.79305,129596581.0481,-0.5532,0.000136,-0.00001149])
    args[3] = argument(t,[1072260.70369,1602961601.2090,-6.3706,0.006593,-0.00003169])
    var dp = 0.0, de = 0.0
    for row in EarthOrientationSeries.xls.reversed() {
      let a = phase(row,args).truncatingRemainder(dividingBy: turn)
      dp += (row[5]+row[6]*t)*sin(a) + row[7]*cos(a)
      de += (row[8]+row[9]*t)*cos(a) + row[10]*sin(a)
    }
    let lsLongitude = dp * (arcsecond/1e7), lsObliquity = de * (arcsecond/1e7)
    args = [longitude(t,2.35555598,8328.6914269554),longitude(t,1.627905234,8433.466158131),
      longitude(t,5.198466741,7771.3771468121),longitude(t,2.18243920,-33.757045),
      longitude(t,4.402608842,2608.7903141574),longitude(t,3.176146697,1021.3285546211),
      longitude(t,1.753470314,628.3075849991),longitude(t,6.203480913,334.0612426700),
      longitude(t,0.599546497,52.9690962641),longitude(t,0.874016757,21.3299104960),
      longitude(t,5.481293872,7.4781598567),longitude(t,5.321159000,3.8127774000),
      (0.024381750+0.00000538691*t)*t]
    dp = 0; de = 0
    for row in EarthOrientationSeries.xpl.reversed() {
      let a = phase(row,args).truncatingRemainder(dividingBy: turn)
      dp += row[13]*sin(a) + row[14]*cos(a)
      de += row[15]*sin(a) + row[16]*cos(a)
    }
    // IAU 2006 secular J2 correction and adjustment to 2000A longitude.
    let j2 = -2.7774e-6*t
    return ((lsLongitude+dp*(arcsecond/1e7))*(1+0.4697e-6+j2),
            (lsObliquity+de*(arcsecond/1e7))*(1+j2))
  }
  private static func cioLocator(_ t: Double, x: Double, y: Double) -> Double {
    let args = fundamentalArguments(t)
    var w = [94.00e-6,3808.65e-6,-122.68e-6,-72574.11e-6,27.98e-6,15.62e-6]
    let series = [EarthOrientationSeries.s0,EarthOrientationSeries.s1,EarthOrientationSeries.s2,
                  EarthOrientationSeries.s3,EarthOrientationSeries.s4]
    for order in series.indices {
      for row in series[order].reversed() {
        let a = phase(row,args)
        w[order] += row[8]*sin(a) + row[9]*cos(a)
      }
    }
    return polynomial(t,w)*arcsecond - x*y/2
  }
  /// Passive coordinate rotations, premultiplied in the published sequence.
  private static func rotation(_ angle: Double, axis: Int) -> Matrix3x3 {
    let c = cos(angle), s = sin(angle)
    switch axis {
    case 0: return Matrix3x3(elements: [1,0,0,0,c,s,0,-s,c])
    case 1: return Matrix3x3(elements: [c,0,-s,0,1,0,s,0,c])
    default: return Matrix3x3(elements: [c,s,0,-s,c,0,0,0,1])
    }
  }
  private static func earthRotationAngle(_ date1: Double, _ date2: Double) -> Double {
    let d1 = min(date1,date2), d2 = max(date1,date2)
    let t = d1 + (d2-2451545)
    let f = d1.truncatingRemainder(dividingBy: 1) + d2.truncatingRemainder(dividingBy: 1)
    let a = (turn*(f+0.7790572732640+0.00273781191135448*t)).truncatingRemainder(dividingBy: turn)
    return a < 0 ? a+turn : a
  }
  /// IAU 2000A bias/precession/nutation (ERFA bp00, obl80, pr00, numat).
  private static func precessionNutation2000(_ t: Double, n: (longitude: Double, obliquity: Double)) -> Matrix3x3 {
    let eps0 = 84381.448*arcsecond
    let dp = -0.29965*t*arcsecond, de = -0.02524*t*arcsecond
    let psi = polynomial(t,[0,5038.7784,-1.07259,-0.001147])*arcsecond+dp
    let omega = eps0+polynomial(t,[0,0,0.05127,-0.007726])*arcsecond+de
    let chi = polynomial(t,[0,10.5526,-2.38064,-0.001125])*arcsecond
    let bias = rotation(0.0068192*arcsecond,axis: 0)
      * rotation(-0.041775*arcsecond*sin(eps0),axis: 1) * rotation(-0.0146*arcsecond,axis: 2)
    let precession = rotation(chi,axis: 2) * rotation(-omega,axis: 0)
      * rotation(-psi,axis: 2) * rotation(eps0,axis: 0)
    let eps = polynomial(t,[84381.448,-46.8150,-0.00059,0.001813])*arcsecond+de
    let j2 = -2.7774e-6*t
    let nut = rotation(-(eps+n.obliquity/(1+j2)),axis: 0)
      * rotation(-n.longitude/(1+0.4697e-6+j2),axis: 2) * rotation(eps,axis: 0)
    return nut * precession * bias
  }
  package static func matrix(edge: Edge, tt1: Double, tt2: Double,
                     ut11: Double, ut12: Double, xp: Double, yp: Double,
                     dx: Double = 0, dy: Double = 0,
                     poleReference: PoleReference = .iau2006A) throws -> [Double] {
    let t = ((tt1-2451545)+tt2)/36525
    let r: Matrix3x3
    if edge == .polarMotion {
      r = rotation(-yp,axis: 0) * rotation(-xp,axis: 1) * rotation(-47e-6*t*arcsecond,axis: 2)
    } else if edge == .earthRotation {
      r = rotation(earthRotationAngle(ut11,ut12),axis: 2)
    } else {
      let n = nutation(t)
      let gamma = polynomial(t,[-0.052928,10.556378,0.4932044,-0.00031238,-0.000002788,0.0000000260])*arcsecond
      let phi = polynomial(t,[84381.412819,-46.811016,0.0511268,0.00053289,-0.000000440,-0.0000000176])*arcsecond
      let psi = polynomial(t,[-0.041775,5038.481484,1.5584175,-0.00018522,-0.000026452,-0.0000000148])*arcsecond+n.longitude
      let eps = polynomial(t,[84381.406,-46.836769,-0.0001831,0.00200340,-0.000000576,-0.0000000434])*arcsecond+n.obliquity
      func bpn(_ p: Double, _ e: Double) -> Matrix3x3 {
        rotation(-e,axis: 0) * rotation(-p,axis: 2) * rotation(phi,axis: 0) * rotation(gamma,axis: 2)
      }
      var b = bpn(psi, eps)
      // SOFA procedure: retain the theoretical locator s when applying observed CIP offsets.
      let s = cioLocator(t,x: b[2,0],y: b[2,1])
      var x = b[2,0] + dx, y = b[2,1] + dy
      if poleReference == .iau2000A {
        let old = precessionNutation2000(t, n: n)
        x = old[2,0] + dx; y = old[2,1] + dy
      }
      guard x.isFinite, y.isFinite, x*x+y*y < 1 else { throw EvaluationError.invalidParameters }
      if dx != 0 || dy != 0 || poleReference == .iau2000A {
        // Recover nutation angles for the observed pole, preserving the P03 equinox convention.
        var p = psi, e = eps
        var converged = false
        for _ in 0..<10 {
          b = bpn(p,e)
          let rx = x-b[2,0], ry = y-b[2,1]
          if max(abs(rx),abs(ry)) < 2e-16 { converged = true; break }
          let h = 1e-6
          let pp = bpn(p+h,e), pm = bpn(p-h,e), ep = bpn(p,e+h), em = bpn(p,e-h)
          let a = (pp[2,0]-pm[2,0])/(2*h), c = (pp[2,1]-pm[2,1])/(2*h)
          let bb = (ep[2,0]-em[2,0])/(2*h), d = (ep[2,1]-em[2,1])/(2*h)
          let determinant = a*d-bb*c
          guard abs(determinant) > 1e-12 else { throw EvaluationError.correctionDidNotConverge }
          p += (rx*d-bb*ry)/determinant; e += (a*ry-rx*c)/determinant
        }
        guard converged else { throw EvaluationError.correctionDidNotConverge }
      }
      if edge == .trueOfDate { r = b }
      else if edge == .celestialIntermediate {
        let radius = x*x+y*y
        let e = radius > 0 ? atan2(y,x) : 0
        let d = atan(sqrt(radius/(1-radius)))
        r = rotation(-(e+s),axis: 2) * rotation(d,axis: 1) * rotation(e,axis: 2)
      } else {
        let ax = x/(1+b[2,2])
        let v = SIMD3(1-ax*x,-ax*y,-x)
        let p = (b[0]*v).sum(), q = (b[1]*v).sum()
        let origins = p != 0 || q != 0 ? s-atan2(q,p) : s
        r = rotation(earthRotationAngle(ut11,ut12)-origins,axis: 2)
      }
    }
    return (0..<3).flatMap { i in (0..<3).map { r[i,$0] } }
  }
}
