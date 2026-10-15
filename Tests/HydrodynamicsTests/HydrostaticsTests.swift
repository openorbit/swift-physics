// SPDX-License-Identifier: Apache-2.0
import Foundation
import Testing
import FluidDynamics
import Hydrodynamics

@Suite struct HydrostaticsTests {
    func properties(density: Double = 1000) throws -> LiquidProperties {
        try LiquidProperties(identifier: "Synthetic test liquid",source: "Analytic reference",densityKilogramsPerCubicMetre: density,
            dynamicViscosityPascalSeconds: 0.001,speedOfSoundMetresPerSecond: 1500,vaporPressurePascals: 2000)
    }
    func liquid(surface: SIMD3<Double> = .zero, up: SIMD3<Double> = SIMD3(0,0,-1), density: Double = 1000, gravity: Double = 10) throws -> HydrostaticLiquid {
        try HydrostaticLiquid(properties: properties(density: density),surfacePointMetres: surface,up: up,
            gravityMetresPerSecondSquared: -up*gravity,surfaceAbsolutePressurePascals: 100000,maximumDepthMetres: 100)
    }
    @Test func halfImmersedBoxVolumeCentroidAndPressure() throws {
        let hull = try ConvexDisplacementHull.box(halfExtentsMetres: SIMD3(repeating: 1))
        let water = try liquid(), displacement = try hull.displacement(in: water)
        #expect(abs(hull.volumeCubicMetres-8) < 1e-12)
        #expect(abs(displacement.volumeCubicMetres-4) < 1e-12)
        #expect(abs(displacement.wettedAreaSquareMetres-12) < 1e-12)
        let centre = try #require(displacement.centreOfBuoyancyMetres)
        #expect(abs(centre.x)+abs(centre.y)+abs(centre.z-0.5) < 1e-12)
        let loads = try displacement.buoyancy(in: water,momentReferencePointMetres: SIMD3(1,0,0))
        #expect(abs(loads.forceNewtons.z+40000) < 1e-8)
        #expect(abs(loads.momentNewtonMetres.y+40000) < 1e-8)
        #expect(try water.absolutePressurePascals(at: SIMD3(0,0,2)) == 120000)
        #expect(try water.absolutePressurePascals(at: SIMD3(0,0,-1)) == nil)
        #expect(throws: FluidError.outsideCoverage) { try water.absolutePressurePascals(at: SIMD3(0,0,101)) }
    }
    @Test func tetrahedronClippingHasAnalyticVolumeAndCentroid() throws {
        let tetra = try ConvexDisplacementHull(verticesMetres: [.zero,SIMD3(1,0,0),SIMD3(0,1,0),SIMD3(0,0,1)],
            faces: [[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
        #expect(abs(tetra.volumeCubicMetres-1.0/6) < 1e-14)
        let small = try tetra.displacement(in: liquid(surface: SIMD3(0,0,0.5)))
        #expect(abs(small.volumeCubicMetres-1.0/48) < 1e-14)
        let c = try #require(small.centreOfBuoyancyMetres)
        #expect(abs(c.x-0.125)+abs(c.y-0.125)+abs(c.z-0.625) < 1e-12)
        let full = try tetra.displacement(in: liquid(surface: SIMD3(0,0,-1)))
        #expect(abs(full.volumeCubicMetres-tetra.volumeCubicMetres) < 1e-14)
        let dry = try tetra.displacement(in: liquid(surface: SIMD3(0,0,1)))
        #expect(dry.volumeCubicMetres == 0 && dry.centreOfBuoyancyMetres == nil)
    }
    @Test func heelRestoringSlopeAndTitanScaling() throws {
        let hull = try ConvexDisplacementHull.box(halfExtentsMetres: SIMD3(repeating: 1))
        for angle in [-1e-4,1e-4] {
            let up = SIMD3(0,sin(angle),-cos(angle)), water = try liquid(up: up)
            let displacement = try hull.displacement(in: water)
            #expect(abs(displacement.volumeCubicMetres-4) < 1e-12)
            let loads = try displacement.buoyancy(in: water,momentReferencePointMetres: SIMD3(0,0,0.5))
            // At G=B upright, GM=I_waterplane/V=(2*2³/12)/4=1/3 m.
            #expect(abs(loads.momentNewtonMetres.x/(40000*angle)-1.0/3) < 1e-7)
        }
        let titan = try liquid(density: 660,gravity: 1.352)
        let force = try hull.displacement(in: titan).buoyancy(in: titan).forceNewtons.z
        #expect(abs(force+660*1.352*4) < 1e-8)
        let translated = try ConvexDisplacementHull.box(halfExtentsMetres: SIMD3(repeating: 1),centreMetres: SIMD3(4,-3,7))
        let c = try #require(translated.displacement(in: liquid(surface: SIMD3(4,-3,7))).centreOfBuoyancyMetres)
        #expect(abs(c.x-4)+abs(c.y+3)+abs(c.z-7.5) < 1e-12)
    }
    @Test func sphereCapsAndPluggableDisplacement() throws {
        let shape: any DisplacementGeometry = try SphereDisplacement(centreMetres: .zero,radiusMetres: 2)
        let sample = try shape.displacement(in: liquid())
        #expect(abs(sample.volumeCubicMetres-16*Double.pi/3) < 1e-12)
        #expect(abs(sample.wettedAreaSquareMetres-8*Double.pi) < 1e-12)
        #expect(abs(try #require(sample.centreOfBuoyancyMetres).z-0.75) < 1e-12)
        let dry = try shape.displacement(in: liquid(surface: SIMD3(0,0,2)))
        #expect(dry.volumeCubicMetres == 0)
        #expect(throws: FluidError.invalidGeometry) { try DisplacedLiquidVolume(volumeCubicMetres: 1,centreOfBuoyancyMetres: nil,wettedAreaSquareMetres: 1) }
    }
    @Test func invalidMeshAndNonhydrostaticPlaneAreRejected() throws {
        let hull = try ConvexDisplacementHull.box(halfExtentsMetres: SIMD3(repeating: 1))
        #expect(throws: FluidError.invalidGeometry) { try ConvexDisplacementHull(verticesMetres: hull.verticesMetres,faces: Array(hull.faces.dropLast())) }
        var faces = hull.faces; faces[0].reverse()
        #expect(throws: FluidError.invalidGeometry) { try ConvexDisplacementHull(verticesMetres: hull.verticesMetres,faces: faces) }
        #expect(throws: FluidError.unsupportedQuery) {
            try HydrostaticLiquid(properties: properties(),up: SIMD3(0,0,-1),gravityMetresPerSecondSquared: SIMD3(0,1,10),
                surfaceAbsolutePressurePascals: 100000,maximumDepthMetres: 100)
        }
    }
}
