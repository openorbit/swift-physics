// SPDX-License-Identifier: Apache-2.0
import Foundation

/// Internal-field potential V = a Σ (a/r)^(n+1) (g cos(mλ)+h sin(mλ)) P̅nm(cos θ).
/// Schmidt semi-normalization without the Condon–Shortley phase; B = −∇V.
/// All evaluation scratch storage is local. Coefficients remain immutable and reusable.
struct GeomagneticHarmonics: MagneticFieldModel {
    let coefficients: GeomagneticCoefficients
    let metadata: MagneticFieldMetadata
    private let normalization: [Double]
    init(coefficients: GeomagneticCoefficients, source: String) throws {
        self.coefficients = coefficients
        metadata = try MagneticFieldMetadata(name: coefficients.name,source: source,version: "File coefficients; epochs \(coefficients.epochs.first!)...\(coefficients.endYear)",
            axesDescription: "Earth-fixed geocentric XYZ: +X latitude 0 longitude 0, +Y longitude 90° east, +Z north",
            assumptions: "Schmidt semi-normalized internal harmonics; 6371200 m reference sphere. Linear epoch interpolation and final secular variation. No external/crustal field or geoid correction. Radius bound is a numerical evaluation bound, not model accuracy certification.",
            radiusRangeMetres: 6_000_000...7_500_000,
            decimalYearRange: coefficients.epochs.first!...coefficients.endYear)
        var factors = [Double](repeating: 0,count: GeomagneticCoefficients.index(coefficients.degree,coefficients.degree)+1)
        factors[0] = 1
        for n in 1...coefficients.degree {
            var s = 1.0
            factors[GeomagneticCoefficients.index(n,0)] = s
            for m in 1...n {
                s *= sqrt((m == 1 ? 2.0 : 1.0)/Double((n-m+1)*(n+m)))
                factors[GeomagneticCoefficients.index(n,m)] = s
            }
        }
        normalization = factors
    }
    func field(at query: MagneticFieldQuery) throws -> MagneticFluxDensity {
        try metadata.validate(query)
        let p = query.positionMetres, r = hypot(hypot(p.x,p.y),p.z)
        guard r > 0 else { throw MagneticFieldError.singularField }
        let horizontal = hypot(p.x,p.y), c = p.z/r, s = horizontal/r
        let longitude = horizontal == 0 ? 0 : atan2(p.y,p.x)
        let size = normalization.count
        var legendre = [Double](repeating: 0,count: size), derivative = legendre, overSin = legendre
        legendre[0] = 1
        let index = GeomagneticCoefficients.index
        for n in 1...coefficients.degree {
            for m in 0...n {
                let i = index(n,m)
                if n == m {
                    let prev = index(n-1,m-1), factor = Double(2*n-1)
                    legendre[i] = factor*s*legendre[prev]
                    derivative[i] = factor*(c*legendre[prev]+s*derivative[prev])
                    overSin[i] = factor*legendre[prev]
                } else {
                    let prev = index(n-1,m), denominator = Double(n-m)
                    let a = Double(2*n-1)
                    let old = n > m+1 ? index(n-2,m) : 0
                    let b = n > m+1 ? Double(n+m-1) : 0
                    legendre[i] = (a*c*legendre[prev]-b*legendre[old])/denominator
                    derivative[i] = (a*(-s*legendre[prev]+c*derivative[prev])-b*derivative[old])/denominator
                    if m > 0 { overSin[i] = (a*c*overSin[prev]-b*overSin[old])/denominator }
                }
            }
        }
        let values = coefficients.at(query.decimalYear), ratio = 6_371_200/r
        var radial = 0.0, colatitude = 0.0, azimuthal = 0.0, radialPower = ratio*ratio
        for n in 1...coefficients.degree {
            radialPower *= ratio
            for m in 0...n {
                let i = index(n,m), cosM = cos(Double(m)*longitude), sinM = sin(Double(m)*longitude)
                let g = values[2*i], h = values[2*i+1], term = g*cosM+h*sinM
                let scale = radialPower*normalization[i]
                radial += Double(n+1)*scale*term*legendre[i]
                colatitude -= scale*term*derivative[i]
                azimuthal += scale*Double(m)*(g*sinM-h*cosM)*overSin[i]
            }
        }
        let cl = cos(longitude), sl = sin(longitude)
        let vector = SIMD3(radial*s*cl+colatitude*c*cl-azimuthal*sl,
                           radial*s*sl+colatitude*c*sl+azimuthal*cl,
                           radial*c-colatitude*s)*1e-9
        return try MagneticFluxDensity(vectorTesla: vector)
    }
}
