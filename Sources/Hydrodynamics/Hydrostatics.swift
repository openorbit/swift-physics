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
import FluidDynamics

/// Flat hydrostatic free surface. Up points out of the liquid; gravity must point
/// in the opposite direction. Positions, geometry and gravity share caller axes.
public struct HydrostaticLiquid: Sendable {
    public let properties: LiquidProperties
    public let surfacePointMetres: SIMD3<Double>
    public let up: SIMD3<Double>
    public let gravityMetresPerSecondSquared: SIMD3<Double>
    public let surfaceAbsolutePressurePascals: Double
    public let maximumDepthMetres: Double
    public init(properties: LiquidProperties, surfacePointMetres: SIMD3<Double> = .zero, up: SIMD3<Double>,
                gravityMetresPerSecondSquared: SIMD3<Double>, surfaceAbsolutePressurePascals: Double,
                maximumDepthMetres: Double) throws {
        let length = fluidNorm(up), g = fluidNorm(gravityMetresPerSecondSquared)
        guard fluidFinite(surfacePointMetres), fluidFinite(up), length.isFinite, length > 0,
              fluidFinite(gravityMetresPerSecondSquared), g.isFinite, g > 0,
              surfaceAbsolutePressurePascals.isFinite, surfaceAbsolutePressurePascals > properties.vaporPressurePascals,
              maximumDepthMetres.isFinite, maximumDepthMetres > 0 else { throw FluidError.invalidParameter }
        guard fluidNorm(up/length+gravityMetresPerSecondSquared/g) < 1e-10 else { throw FluidError.unsupportedQuery }
        self.properties = properties; self.surfacePointMetres = surfacePointMetres; self.up = up/length
        self.gravityMetresPerSecondSquared = gravityMetresPerSecondSquared
        self.surfaceAbsolutePressurePascals = surfaceAbsolutePressurePascals; self.maximumDepthMetres = maximumDepthMetres
    }
    /// Signed distance above the liquid surface; negative inside the liquid.
    public func surfaceDistanceMetres(at positionMetres: SIMD3<Double>) throws -> Double {
        guard fluidFinite(positionMetres) else { throw FluidError.invalidGeometry }
        let distance = ((positionMetres-surfacePointMetres)*up).sum()
        guard distance.isFinite else { throw FluidError.nonFiniteResult }
        guard distance >= -maximumDepthMetres else { throw FluidError.outsideCoverage }
        return distance
    }
    /// Nil above the surface. Absolute liquid pressure includes atmosphere/surface pressure.
    public func absolutePressurePascals(at positionMetres: SIMD3<Double>) throws -> Double? {
        let distance = try surfaceDistanceMetres(at: positionMetres)
        guard distance <= 0 else { return nil }
        let pressure = surfaceAbsolutePressurePascals-properties.densityKilogramsPerCubicMetre*fluidNorm(gravityMetresPerSecondSquared)*distance
        guard pressure.isFinite else { throw FluidError.nonFiniteResult }; return pressure
    }
}

public struct DisplacedLiquidVolume: Sendable {
    public let volumeCubicMetres: Double
    public let centreOfBuoyancyMetres: SIMD3<Double>?
    /// Actual hull area below the surface; excludes the artificial clipping cap.
    public let wettedAreaSquareMetres: Double
    public init(volumeCubicMetres: Double, centreOfBuoyancyMetres: SIMD3<Double>?, wettedAreaSquareMetres: Double) throws {
        guard volumeCubicMetres.isFinite, volumeCubicMetres >= 0, wettedAreaSquareMetres.isFinite, wettedAreaSquareMetres >= 0,
              centreOfBuoyancyMetres.map(fluidFinite) ?? true,
              volumeCubicMetres == 0 || centreOfBuoyancyMetres != nil else { throw FluidError.invalidGeometry }
        self.volumeCubicMetres = volumeCubicMetres; self.centreOfBuoyancyMetres = volumeCubicMetres == 0 ? nil : centreOfBuoyancyMetres
        self.wettedAreaSquareMetres = wettedAreaSquareMetres
    }
    public func buoyancy(in liquid: HydrostaticLiquid, momentReferencePointMetres: SIMD3<Double> = .zero) throws -> FluidLoads {
        try Buoyancy(displacedVolumeCubicMetres: volumeCubicMetres,centreOfBuoyancyMetres: centreOfBuoyancyMetres ?? .zero)
            .evaluate(densityKilogramsPerCubicMetre: liquid.properties.densityKilogramsPerCubicMetre,
                field: .hydrostaticGravityMetresPerSecondSquared(liquid.gravityMetresPerSecondSquared),
                momentReferencePointMetres: momentReferencePointMetres)
    }
}

public protocol DisplacementGeometry: Sendable {
    func displacement(in liquid: HydrostaticLiquid) throws -> DisplacedLiquidVolume
}

/// Analytic spherical-cap displacement, useful for floats and buoyancy primitives.
public struct SphereDisplacement: DisplacementGeometry {
    public let centreMetres: SIMD3<Double>
    public let radiusMetres: Double
    public init(centreMetres: SIMD3<Double>, radiusMetres: Double) throws {
        guard fluidFinite(centreMetres), radiusMetres.isFinite, radiusMetres > 0,
              (4*Double.pi/3*radiusMetres*radiusMetres*radiusMetres).isFinite,
              4*Double.pi/3*radiusMetres*radiusMetres*radiusMetres > 0 else { throw FluidError.invalidGeometry }
        self.centreMetres = centreMetres; self.radiusMetres = radiusMetres
    }
    public func displacement(in liquid: HydrostaticLiquid) throws -> DisplacedLiquidVolume {
        let r = radiusMetres
        _ = try liquid.surfaceDistanceMetres(at: centreMetres-liquid.up*r)
        let d = try liquid.surfaceDistanceMetres(at: centreMetres)
        if d >= r { return try DisplacedLiquidVolume(volumeCubicMetres: 0,centreOfBuoyancyMetres: nil,wettedAreaSquareMetres: 0) }
        if d <= -r { return try DisplacedLiquidVolume(volumeCubicMetres: 4*Double.pi/3*r*r*r,centreOfBuoyancyMetres: centreMetres,wettedAreaSquareMetres: 4*Double.pi*r*r) }
        let h = r-d, volume = Double.pi*h*h*(r-h/3)
        let offset = -(2*r-h)*(2*r-h)/(4*(r-h/3))
        return try DisplacedLiquidVolume(volumeCubicMetres: volume,centreOfBuoyancyMetres: centreMetres+liquid.up*offset,wettedAreaSquareMetres: 2*Double.pi*r*h)
    }
}

/// Closed convex displacement geometry, with outward-wound planar polygon faces.
/// Clip against a supplied surface plane to obtain partial immersion and restoring
/// moments. Rotate the surface/gravity into body axes externally for a heeled hull.
public struct ConvexDisplacementHull: DisplacementGeometry {
    public let verticesMetres: [SIMD3<Double>]
    public let faces: [[Int]]
    public let volumeCubicMetres: Double
    public let volumeCentroidMetres: SIMD3<Double>
    private let tolerance: Double
    public init(verticesMetres: [SIMD3<Double>], faces: [[Int]]) throws {
        guard verticesMetres.count >= 4, verticesMetres.count <= 4096, faces.count >= 4, faces.count <= 8192,
              verticesMetres.allSatisfy(fluidFinite) else { throw FluidError.invalidGeometry }
        let reference = verticesMetres.reduce(SIMD3<Double>.zero,+)/Double(verticesMetres.count)
        let scale = verticesMetres.map { fluidNorm($0-reference) }.max()!
        guard scale.isFinite, scale > 0 else { throw FluidError.invalidGeometry }
        let epsilon = scale*1e-10
        var edges: [String:Int] = [:], polygons: [[SIMD3<Double>]] = []
        for face in faces {
            guard face.count >= 3, Set(face).count == face.count, face.allSatisfy({ verticesMetres.indices.contains($0) }) else { throw FluidError.invalidGeometry }
            let polygon = face.map { verticesMetres[$0] }, normal = fluidCross(polygon[1]-polygon[0],polygon[2]-polygon[0])
            let normalLength = fluidNorm(normal)
            guard normalLength.isFinite, normalLength > epsilon*epsilon else { throw FluidError.invalidGeometry }
            let n = normal/normalLength
            guard polygon.allSatisfy({ abs((($0-polygon[0])*n).sum()) <= epsilon }),
                  verticesMetres.allSatisfy({ (($0-polygon[0])*n).sum() <= epsilon }) else { throw FluidError.invalidGeometry }
            for i in polygon.indices {
                let a = polygon[i], edge = polygon[(i+1)%polygon.count]-a
                guard fluidNorm(edge) > epsilon, polygon.allSatisfy({ (fluidCross(edge,$0-a)*n).sum() >= -epsilon*fluidNorm(edge) }) else {
                    throw FluidError.invalidGeometry
                }
            }
            for i in face.indices {
                let a = face[i], b = face[(i+1)%face.count], key = "\(a):\(b)"
                edges[key,default: 0] += 1
            }
            polygons.append(polygon)
        }
        for (key,count) in edges {
            let parts = key.split(separator: ":")
            guard count == 1, edges["\(parts[1]):\(parts[0])"] == 1 else { throw FluidError.invalidGeometry }
        }
        let integral = try Self.integrate(polygons,reference: reference)
        guard integral.volume > epsilon*epsilon*epsilon else { throw FluidError.invalidGeometry }
        self.verticesMetres = verticesMetres; self.faces = faces; self.tolerance = epsilon
        self.volumeCubicMetres = integral.volume; self.volumeCentroidMetres = integral.centroid
    }
    public static func box(halfExtentsMetres e: SIMD3<Double>, centreMetres c: SIMD3<Double> = .zero) throws -> Self {
        guard fluidFinite(e), fluidFinite(c), e.x > 0, e.y > 0, e.z > 0 else { throw FluidError.invalidGeometry }
        let v = [SIMD3(-e.x,-e.y,-e.z),SIMD3(e.x,-e.y,-e.z),SIMD3(e.x,e.y,-e.z),SIMD3(-e.x,e.y,-e.z),
                 SIMD3(-e.x,-e.y,e.z),SIMD3(e.x,-e.y,e.z),SIMD3(e.x,e.y,e.z),SIMD3(-e.x,e.y,e.z)].map { $0+c }
        return try Self(verticesMetres: v,faces: [[0,3,2,1],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7]])
    }
    private static func integrate(_ polygons: [[SIMD3<Double>]], reference: SIMD3<Double>) throws -> (volume: Double, centroid: SIMD3<Double>) {
        var volume = 0.0, moment = SIMD3<Double>.zero
        for p in polygons where p.count >= 3 {
            for j in 1..<(p.count-1) {
                let a = p[0]-reference, b = p[j]-reference, c = p[j+1]-reference
                let v = (a*fluidCross(b,c)).sum()/6
                volume += v; moment += (a+b+c)*(v/4)
            }
        }
        guard volume.isFinite, fluidFinite(moment) else { throw FluidError.nonFiniteResult }
        return (volume,volume > 0 ? reference+moment/volume : reference)
    }
    public func displacement(in liquid: HydrostaticLiquid) throws -> DisplacedLiquidVolume {
        var polygons: [[SIMD3<Double>]] = [], cap: [SIMD3<Double>] = [], area = 0.0
        let distances = try verticesMetres.map { try liquid.surfaceDistanceMetres(at: $0) }
        if distances.allSatisfy({ $0 >= 0 }) { return try DisplacedLiquidVolume(volumeCubicMetres: 0,centreOfBuoyancyMetres: nil,wettedAreaSquareMetres: 0) }
        for face in faces {
            var clipped: [SIMD3<Double>] = []
            for i in face.indices {
                let ai = face[i], bi = face[(i+1)%face.count], a = verticesMetres[ai], b = verticesMetres[bi]
                let da = distances[ai], db = distances[bi]
                if da <= 0 { clipped.append(a) }
                if (da <= 0) != (db <= 0) {
                    let point = a+(b-a)*(da/(da-db)); clipped.append(point)
                    if !cap.contains(where: { fluidNorm($0-point) <= tolerance }) { cap.append(point) }
                }
            }
            // Remove repeated vertices caused by a clipping plane passing through a vertex.
            var unique: [SIMD3<Double>] = []
            for p in clipped where !unique.contains(where: { fluidNorm($0-p) <= tolerance }) { unique.append(p) }
            if unique.count >= 3 {
                polygons.append(unique)
                for j in 1..<(unique.count-1) { area += fluidNorm(fluidCross(unique[j]-unique[0],unique[j+1]-unique[0]))/2 }
            }
        }
        if cap.count >= 3 {
            let centre = cap.reduce(SIMD3<Double>.zero,+)/Double(cap.count)
            let candidate = abs(liquid.up.x) < 0.8 ? SIMD3(1.0,0,0) : SIMD3(0.0,1,0)
            let raw = fluidCross(candidate,liquid.up), u = raw/fluidNorm(raw), v = fluidCross(liquid.up,u)
            cap.sort { a,b in
                atan2(((a-centre)*v).sum(),((a-centre)*u).sum()) < atan2(((b-centre)*v).sum(),((b-centre)*u).sum())
            }
            polygons.append(cap)
        }
        let integral = try Self.integrate(polygons,reference: volumeCentroidMetres)
        guard integral.volume >= -tolerance*tolerance*tolerance, integral.volume <= volumeCubicMetres*(1+1e-9), area.isFinite else {
            throw FluidError.nonFiniteResult
        }
        let volume = max(0,min(volumeCubicMetres,integral.volume))
        return try DisplacedLiquidVolume(volumeCubicMetres: volume,centreOfBuoyancyMetres: volume > 0 ? integral.centroid : nil,wettedAreaSquareMetres: area)
    }
}
