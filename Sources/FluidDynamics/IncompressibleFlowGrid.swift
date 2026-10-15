// SPDX-License-Identifier: Apache-2.0
import Foundation

/// Uniform Cartesian MAC-grid geometry and constant fluid properties, in SI.
/// The frame must translate at constant velocity without rotating; moving-frame
/// acceleration/rotation and moving solid geometry are not inferred.
public struct IncompressibleFlowConfiguration: Sendable {
    public let grid: CartesianFluidGrid
    public var dimensions: SIMD3<Int> { grid.dimensions }
    public var minimumCornerMetres: SIMD3<Double> { grid.minimumCornerMetres }
    public var cellSizeMetres: Double { grid.cellSizeMetres }
    public let densityKilogramsPerCubicMetre: Double
    public let speedOfSoundMetresPerSecond: Double
    public let kinematicViscositySquareMetresPerSecond: Double
    /// Prescribed normal velocity at all six outer faces; net box flux is zero.
    public let boundaryVelocityMetresPerSecond: SIMD3<Double>
    public var cellCount: Int { dimensions.x*dimensions.y*dimensions.z }
    public init(dimensions: SIMD3<Int> = SIMD3(repeating: 16), minimumCornerMetres: SIMD3<Double>,
                cellSizeMetres: Double, densityKilogramsPerCubicMetre: Double,
                speedOfSoundMetresPerSecond: Double, kinematicViscositySquareMetresPerSecond: Double = 0,
                boundaryVelocityMetresPerSecond: SIMD3<Double> = .zero) throws {
        guard [dimensions.x,dimensions.y,dimensions.z].allSatisfy({ (4...64).contains($0) }),
              fluidFinite(minimumCornerMetres), cellSizeMetres.isFinite, cellSizeMetres > 0,
              (cellSizeMetres*cellSizeMetres*cellSizeMetres).isFinite,
              cellSizeMetres*cellSizeMetres*cellSizeMetres > 0,
              densityKilogramsPerCubicMetre.isFinite, densityKilogramsPerCubicMetre > 0,
              speedOfSoundMetresPerSecond.isFinite, speedOfSoundMetresPerSecond > 0,
              kinematicViscositySquareMetresPerSecond.isFinite, kinematicViscositySquareMetresPerSecond >= 0,
              fluidFinite(boundaryVelocityMetresPerSecond) else { throw FluidError.invalidParameter }
        let extent = SIMD3(Double(dimensions.x),Double(dimensions.y),Double(dimensions.z))*cellSizeMetres
        guard fluidFinite(minimumCornerMetres+extent),
              (0..<3).allSatisfy({ minimumCornerMetres[$0]+cellSizeMetres > minimumCornerMetres[$0] }) else { throw FluidError.invalidGeometry }
        guard fluidNorm(boundaryVelocityMetresPerSecond)/speedOfSoundMetresPerSecond <= 0.3 else { throw FluidError.outsideCoverage }
        self.grid = try CartesianFluidGrid(dimensions: dimensions,minimumCornerMetres: minimumCornerMetres,cellSizeMetres: cellSizeMetres)
        self.densityKilogramsPerCubicMetre = densityKilogramsPerCubicMetre
        self.speedOfSoundMetresPerSecond = speedOfSoundMetresPerSecond
        self.kinematicViscositySquareMetresPerSecond = kinematicViscositySquareMetresPerSecond
        self.boundaryVelocityMetresPerSecond = boundaryVelocityMetresPerSecond
    }
}

public struct PressureProjectionReport: Sendable {
    public let iterations: Int
    public let relativePressureResidual: Double
    public let maximumDivergenceBeforePerSecond: Double
    public let maximumDivergenceAfterPerSecond: Double
}
public struct IncompressibleStepReport: Sendable {
    public let substeps: Int
    public let totalPressureIterations: Int
    public let maximumDivergenceAfterPerSecond: Double
    public let maximumSpeedBoundMetresPerSecond: Double
}

/// Optional CPU constant-density incompressible flow solver. Face-centred velocity,
/// cell pressure, first-order semi-Lagrangian advection, explicit viscosity and
/// Jacobi-preconditioned conjugate-gradient pressure projection. No turbulence,
/// boundary-layer, compressible shock or thermal model. Fixed stair-step solids
/// impose zero normal velocity (free-slip inviscid boundary, not resolved no-slip).
/// Owns fluid state, not scheduling, a vehicle assembly or rigid-body integration.
public struct IncompressibleFlowGrid: FluidVelocityField {
    public let configuration: IncompressibleFlowConfiguration
    public let solidCells: [Bool]
    /// Gauge pressure in Pa; zero mean over fluid cells after projection.
    public private(set) var pressurePascals: [Double]
    private var faces: [[Double]]
    private let neighbours: [Int]
    private let diagonal: [Double]
    private var nx: Int { configuration.dimensions.x }
    private var ny: Int { configuration.dimensions.y }
    private var nz: Int { configuration.dimensions.z }
    private var h: Double { configuration.cellSizeMetres }

    /// Initial velocity is total fluid velocity in grid axes, including background.
    /// Project after arbitrary initialization before treating it as divergence-free.
    /// Solids may not touch the outermost cell layer, preserving balanced box flux.
    public init(configuration: IncompressibleFlowConfiguration,
                solidAtCellCentre: (SIMD3<Double>) throws -> Bool = { _ in false },
                initialVelocity: ((SIMD3<Double>) throws -> SIMD3<Double>)? = nil) throws {
        self.configuration = configuration
        let n = configuration.dimensions, count = configuration.cellCount, origin = configuration.minimumCornerMetres
        let dx = configuration.cellSizeMetres
        var solids = [Bool](repeating: false,count: count)
        for k in 0..<n.z { for j in 0..<n.y { for i in 0..<n.x {
            let p = origin+SIMD3(Double(i)+0.5,Double(j)+0.5,Double(k)+0.5)*dx
            let solid = try solidAtCellCentre(p)
            guard !solid || (i > 0 && j > 0 && k > 0 && i < n.x-1 && j < n.y-1 && k < n.z-1) else {
                throw FluidError.invalidGeometry
            }
            solids[i+n.x*(j+n.y*k)] = solid
        } } }
        self.solidCells = solids; self.pressurePascals = [Double](repeating: 0,count: count)
        var adjacent = [Int](repeating: -1,count: count*6), diag = [Double](repeating: 0,count: count)
        let offsets = [SIMD3(-1,0,0),SIMD3(1,0,0),SIMD3(0,-1,0),SIMD3(0,1,0),SIMD3(0,0,-1),SIMD3(0,0,1)]
        for k in 0..<n.z { for j in 0..<n.y { for i in 0..<n.x {
            let index = i+n.x*(j+n.y*k)
            if solids[index] { continue }
            for d in 0..<6 {
                let q = SIMD3(i,j,k) &+ offsets[d]
                if q.x >= 0 && q.x < n.x && q.y >= 0 && q.y < n.y && q.z >= 0 && q.z < n.z {
                    let other = q.x+n.x*(q.y+n.y*q.z)
                    if !solids[other] { adjacent[6*index+d] = other; diag[index] += 1/(dx*dx) }
                }
            }
        } } }
        self.neighbours = adjacent; self.diagonal = diag
        self.faces = (0..<3).map { axis in
            var shape = n; shape[axis] += 1
            return [Double](repeating: 0,count: shape.x*shape.y*shape.z)
        }
        for axis in 0..<3 {
            let shape = faceShape(axis)
            for k in 0..<shape.z { for j in 0..<shape.y { for i in 0..<shape.x {
                let p = facePosition(axis,i,j,k)
                let v = try initialVelocity?(p) ?? configuration.boundaryVelocityMetresPerSecond
                guard fluidFinite(v) else { throw FluidError.invalidFlow }
                faces[axis][flat(i,j,k,shape)] = v[axis]
            } } }
        }
        enforceBoundaries()
        try validateSpeed()
    }
    private func flat(_ i: Int, _ j: Int, _ k: Int, _ shape: SIMD3<Int>) -> Int { i+shape.x*(j+shape.y*k) }
    private func faceShape(_ axis: Int) -> SIMD3<Int> {
        var shape = configuration.dimensions; shape[axis] += 1; return shape
    }
    private func facePosition(_ axis: Int, _ i: Int, _ j: Int, _ k: Int) -> SIMD3<Double> {
        var p = SIMD3(Double(i)+0.5,Double(j)+0.5,Double(k)+0.5); p[axis] -= 0.5
        return configuration.minimumCornerMetres+p*h
    }
    private func blocked(_ axis: Int, _ i: Int, _ j: Int, _ k: Int) -> Bool {
        let q = SIMD3(i,j,k), n = configuration.dimensions
        if q[axis] == 0 || q[axis] == n[axis] { return false }
        var a = q; a[axis] -= 1
        return solidCells[flat(a.x,a.y,a.z,n)] || solidCells[flat(q.x,q.y,q.z,n)]
    }
    private mutating func enforceBoundaries() {
        let n = configuration.dimensions
        for axis in 0..<3 {
            let shape = faceShape(axis)
            for k in 0..<shape.z { for j in 0..<shape.y { for i in 0..<shape.x {
                let q = SIMD3(i,j,k), index = flat(i,j,k,shape)
                if q[axis] == 0 || q[axis] == n[axis] { faces[axis][index] = configuration.boundaryVelocityMetresPerSecond[axis] }
                else if blocked(axis,i,j,k) { faces[axis][index] = 0 }
            } } }
        }
    }
    /// Conservative component-wise maximum speed bound, also used for low-Mach checks.
    public var maximumSpeedBoundMetresPerSecond: Double {
        fluidNorm(SIMD3(faces[0].reduce(0) { max($0,abs($1)) },faces[1].reduce(0) { max($0,abs($1)) },faces[2].reduce(0) { max($0,abs($1)) }))
    }
    private func validateSpeed() throws {
        guard faces.allSatisfy({ $0.allSatisfy(\.isFinite) }) else { throw FluidError.nonFiniteResult }
        let speed = maximumSpeedBoundMetresPerSecond
        guard speed.isFinite else { throw FluidError.nonFiniteResult }
        guard speed/configuration.speedOfSoundMetresPerSecond <= 0.3 else { throw FluidError.outsideCoverage }
    }
    private func gridCoordinate(_ p: SIMD3<Double>) -> SIMD3<Double> { (p-configuration.minimumCornerMetres)/h }
    private func cellIndex(at p: SIMD3<Double>) throws -> Int {
        guard fluidFinite(p) else { throw FluidError.invalidParameter }
        let q = gridCoordinate(p), n = configuration.dimensions
        guard fluidFinite(q), q.x >= 0, q.y >= 0, q.z >= 0,
              q.x <= Double(n.x), q.y <= Double(n.y), q.z <= Double(n.z) else { throw FluidError.outsideCoverage }
        return flat(min(n.x-1,Int(q.x)),min(n.y-1,Int(q.y)),min(n.z-1,Int(q.z)),n)
    }
    private func sampleComponent(_ field: [Double], axis: Int, coordinate: SIMD3<Double>) -> Double {
        let shape = faceShape(axis)
        var q = coordinate-SIMD3(repeating: 0.5); q[axis] += 0.5
        for d in 0..<3 { q[d] = max(0,min(Double(shape[d]-1),q[d])) }
        let lo = SIMD3(Int(q.x),Int(q.y),Int(q.z)), w = q-SIMD3(Double(lo.x),Double(lo.y),Double(lo.z))
        var value = 0.0
        for z in 0...1 { for y in 0...1 { for x in 0...1 {
            let weight = (x == 0 ? 1-w.x : w.x)*(y == 0 ? 1-w.y : w.y)*(z == 0 ? 1-w.z : w.z)
            let i = min(shape.x-1,lo.x+x), j = min(shape.y-1,lo.y+y), k = min(shape.z-1,lo.z+z)
            value += field[flat(i,j,k,shape)]*weight
        } } }
        return value
    }
    private func sample(_ field: [[Double]], at p: SIMD3<Double>) -> SIMD3<Double> {
        let coordinate = gridCoordinate(p)
        return SIMD3(sampleComponent(field[0],axis: 0,coordinate: coordinate),
                     sampleComponent(field[1],axis: 1,coordinate: coordinate),
                     sampleComponent(field[2],axis: 2,coordinate: coordinate))
    }
    public func velocityMetresPerSecond(at positionMetres: SIMD3<Double>) throws -> SIMD3<Double> {
        let index = try cellIndex(at: positionMetres)
        guard !solidCells[index] else { throw FluidError.unsupportedQuery }
        let velocity = sample(faces,at: positionMetres)
        guard fluidFinite(velocity) else { throw FluidError.nonFiniteResult }
        return velocity
    }
    /// Gauge pressure at the containing fluid cell. A pressure reference is needed
    /// before using this value for phase/cavitation or structural absolute pressure.
    public func gaugePressurePascals(at positionMetres: SIMD3<Double>) throws -> Double {
        let index = try cellIndex(at: positionMetres)
        guard !solidCells[index] else { throw FluidError.unsupportedQuery }
        return pressurePascals[index]
    }
    /// Cell-centred discrete divergence, 1/s. Solid entries are zero.
    public func divergencePerSecond() -> [Double] {
        var result = [Double](repeating: 0,count: configuration.cellCount)
        let us = faceShape(0), vs = faceShape(1), ws = faceShape(2), n = configuration.dimensions
        for k in 0..<nz { for j in 0..<ny { for i in 0..<nx {
            let index = flat(i,j,k,n)
            if solidCells[index] { continue }
            result[index] = (faces[0][flat(i+1,j,k,us)]-faces[0][flat(i,j,k,us)]
                + faces[1][flat(i,j+1,k,vs)]-faces[1][flat(i,j,k,vs)]
                + faces[2][flat(i,j,k+1,ws)]-faces[2][flat(i,j,k,ws)])/h
        } } }
        return result
    }
    private func matrix(_ p: [Double]) -> [Double] {
        var result = [Double](repeating: 0,count: p.count)
        let inverseCellSizeSquared = 1/(h*h)
        for i in p.indices where !solidCells[i] {
            var value = diagonal[i]*p[i]
            for d in 0..<6 { let j = neighbours[6*i+d]; if j >= 0 { value -= p[j]*inverseCellSizeSquared } }
            result[i] = value
        }
        return result
    }
    private func dot(_ a: [Double], _ b: [Double]) -> Double { zip(a,b).reduce(0) { $0+$1.0*$1.1 } }
    /// Atomic projection: failed solves leave the grid unchanged. Relative residual
    /// is L2; absolute divergence tolerance also supplies a stopping floor (1/s).
    public mutating func project(timeStepSeconds: Double, maximumIterations: Int = 300,
                                 relativeTolerance: Double = 1e-8, absoluteDivergenceTolerancePerSecond: Double = 1e-8) throws -> PressureProjectionReport {
        var next = self
        let report = try next.projectInPlace(dt: timeStepSeconds,maximumIterations: maximumIterations,
            relativeTolerance: relativeTolerance,absoluteTolerance: absoluteDivergenceTolerancePerSecond)
        self = next; return report
    }
    private mutating func projectInPlace(dt: Double, maximumIterations: Int, relativeTolerance: Double,
                                        absoluteTolerance: Double) throws -> PressureProjectionReport {
        guard dt.isFinite, dt > 0, (1...5000).contains(maximumIterations),
              relativeTolerance.isFinite, relativeTolerance > 0, relativeTolerance <= 1e-2,
              absoluteTolerance.isFinite, absoluteTolerance > 0 else { throw FluidError.invalidParameter }
        enforceBoundaries()
        let divergence = divergencePerSecond(), scale = configuration.densityKilogramsPerCubicMetre/dt
        var b = divergence.map { -scale*$0 }
        guard b.allSatisfy(\.isFinite) else { throw FluidError.nonFiniteResult }
        let fluidCount = solidCells.filter { !$0 }.count
        let mean = b.reduce(0,+)/Double(fluidCount)
        for i in b.indices where !solidCells[i] { b[i] -= mean }
        let normB = sqrt(dot(b,b)), threshold = max(relativeTolerance*normB,scale*absoluteTolerance)
        guard normB.isFinite, threshold.isFinite else { throw FluidError.nonFiniteResult }
        var p = [Double](repeating: 0,count: b.count), r = b
        var z = r.indices.map { diagonal[$0] > 0 ? r[$0]/diagonal[$0] : 0 }, direction = z
        var rz = dot(r,z), residual = normB, iterations = 0
        while residual > threshold && iterations < maximumIterations {
            let ad = matrix(direction), denominator = dot(direction,ad)
            guard denominator.isFinite, denominator > 0, rz.isFinite, rz > 0 else { throw FluidError.nonConvergence }
            let alpha = rz/denominator
            for i in p.indices { p[i] += alpha*direction[i]; r[i] -= alpha*ad[i] }
            iterations += 1; residual = sqrt(dot(r,r))
            if residual <= threshold { break }
            z = r.indices.map { diagonal[$0] > 0 ? r[$0]/diagonal[$0] : 0 }
            let nextRZ = dot(r,z), beta = nextRZ/rz
            for i in direction.indices { direction[i] = z[i]+beta*direction[i] }
            rz = nextRZ
        }
        guard residual.isFinite, residual <= threshold else { throw FluidError.nonConvergence }
        let pressureMean = p.reduce(0,+)/Double(fluidCount)
        pressurePascals = p.indices.map { solidCells[$0] ? 0 : p[$0]-pressureMean }
        guard pressurePascals.allSatisfy(\.isFinite) else { throw FluidError.nonFiniteResult }
        let n = configuration.dimensions
        for axis in 0..<3 {
            let shape = faceShape(axis)
            for k in 0..<shape.z { for j in 0..<shape.y { for i in 0..<shape.x {
                let q = SIMD3(i,j,k)
                if q[axis] == 0 || q[axis] == n[axis] || blocked(axis,i,j,k) { continue }
                var a = q; a[axis] -= 1
                let pressureDifference = pressurePascals[flat(q.x,q.y,q.z,n)]-pressurePascals[flat(a.x,a.y,a.z,n)]
                faces[axis][flat(i,j,k,shape)] -= (dt/configuration.densityKilogramsPerCubicMetre)*pressureDifference/h
            } } }
        }
        enforceBoundaries(); try validateSpeed()
        return PressureProjectionReport(iterations: iterations,relativePressureResidual: normB == 0 ? 0 : residual/normB,
            maximumDivergenceBeforePerSecond: divergence.map(abs).max() ?? 0,
            maximumDivergenceAfterPerSecond: divergencePerSecond().map(abs).max() ?? 0)
    }
    /// Caller supplies acceleration of the fluid (m/s²), evaluated at face positions.
    /// This is a source term, not a target velocity or an extra aircraft force.
    /// Advection is numerically dissipative; no sub-grid turbulence is invented.
    /// Atomic on failure, including substep-budget exhaustion and pressure failure.
    public mutating func advance(timeStepSeconds: Double, maximumSubsteps: Int = 64,
                                 maximumPressureIterations: Int = 300, relativePressureTolerance: Double = 1e-8,
                                 accelerationMetresPerSecondSquared: (SIMD3<Double>) throws -> SIMD3<Double> = { _ in .zero }) throws -> IncompressibleStepReport {
        guard timeStepSeconds.isFinite, timeStepSeconds > 0, (1...10000).contains(maximumSubsteps) else {
            throw FluidError.invalidParameter
        }
        var next = self, remaining = timeStepSeconds, substeps = 0, totalIterations = 0, maximumDivergence = 0.0
        var acceleration = next.faces
        for axis in 0..<3 {
            let shape = faceShape(axis)
            for k in 0..<shape.z { for j in 0..<shape.y { for i in 0..<shape.x {
                let value = try accelerationMetresPerSecondSquared(facePosition(axis,i,j,k))
                guard fluidFinite(value) else { throw FluidError.invalidFlow }
                acceleration[axis][flat(i,j,k,shape)] = value[axis]
            } } }
        }
        let accelerationBound = fluidNorm(SIMD3(acceleration[0].map(abs).max() ?? 0,
            acceleration[1].map(abs).max() ?? 0,acceleration[2].map(abs).max() ?? 0))
        guard accelerationBound.isFinite else { throw FluidError.nonFiniteResult }
        while remaining > 0 {
            guard substeps < maximumSubsteps else { throw FluidError.workBudgetExceeded }
            let speedBound = next.maximumSpeedBoundMetresPerSecond+accelerationBound*remaining
            guard speedBound.isFinite else { throw FluidError.nonFiniteResult }
            let advectiveLimit = speedBound > 0 ? 0.75*h/speedBound : remaining
            let nu = configuration.kinematicViscositySquareMetresPerSecond
            let viscousLimit = nu > 0 ? h*h/(6*nu) : remaining
            let dt = min(remaining,advectiveLimit,viscousLimit)
            guard dt.isFinite, dt > 0, remaining-dt < remaining else { throw FluidError.outsideCoverage }
            let old = next.faces
            var updated = old
            for axis in 0..<3 {
                let shape = faceShape(axis)
                for k in 0..<shape.z { for j in 0..<shape.y { for i in 0..<shape.x {
                    let q = SIMD3(i,j,k), index = flat(i,j,k,shape)
                    if q[axis] == 0 || q[axis] == configuration.dimensions[axis] || blocked(axis,i,j,k) { continue }
                    let position = facePosition(axis,i,j,k), departure = position-next.sample(old,at: position)*dt
                    // Avoid advecting values out of a solid cell. Coarse stair-step
                    // boundaries do not resolve near-wall flow or sub-cell geometry.
                    let departureCell = try? next.cellIndex(at: departure)
                    let advected = departureCell.map { next.solidCells[$0] } == true ? old[axis][index] : next.sampleComponent(old[axis],axis: axis,coordinate: next.gridCoordinate(departure))
                    var laplacian = 0.0
                    if nu > 0 {
                        for d in 0..<3 {
                            var lo = q, hi = q
                            lo[d] = max(0,q[d]-1); hi[d] = min(shape[d]-1,q[d]+1)
                            // Reflect tangential velocity across a solid: free slip.
                            let low = d != axis && blocked(axis,lo.x,lo.y,lo.z) ? old[axis][index] : old[axis][flat(lo.x,lo.y,lo.z,shape)]
                            let high = d != axis && blocked(axis,hi.x,hi.y,hi.z) ? old[axis][index] : old[axis][flat(hi.x,hi.y,hi.z,shape)]
                            laplacian += (low-2*old[axis][index]+high)/(h*h)
                        }
                    }
                    updated[axis][index] = advected+dt*(nu*laplacian+acceleration[axis][index])
                } } }
            }
            next.faces = updated; next.enforceBoundaries(); try next.validateSpeed()
            let report = try next.projectInPlace(dt: dt,maximumIterations: maximumPressureIterations,
                relativeTolerance: relativePressureTolerance,absoluteTolerance: 1e-8)
            totalIterations += report.iterations; maximumDivergence = max(maximumDivergence,report.maximumDivergenceAfterPerSecond)
            remaining -= dt; substeps += 1
        }
        self = next
        return IncompressibleStepReport(substeps: substeps,totalPressureIterations: totalIterations,
            maximumDivergenceAfterPerSecond: maximumDivergence,maximumSpeedBoundMetresPerSecond: maximumSpeedBoundMetresPerSecond)
    }
    /// Gauge-pressure force on the supplied fixed voxel solids, at grid-face centres.
    /// One-sided linear pressure reconstruction (cell value if no fluid neighbour).
    /// No viscous stress or sub-grid wing lift is included. This replaces the same
    /// pressure contribution in other models; it is not an extra force to add twice.
    public func pressureLoadsOnSolids(momentReferencePointMetres: SIMD3<Double> = .zero) throws -> FluidLoads {
        var force = SIMD3<Double>.zero, moment = SIMD3<Double>.zero
        let n = configuration.dimensions
        guard fluidFinite(momentReferencePointMetres) else { throw FluidError.invalidParameter }
        for axis in 0..<3 {
            let shape = faceShape(axis)
            for k in 0..<shape.z { for j in 0..<shape.y { for i in 0..<shape.x {
                let q = SIMD3(i,j,k)
                if q[axis] == 0 || q[axis] == n[axis] { continue }
                var a = q; a[axis] -= 1
                let left = flat(a.x,a.y,a.z,n), right = flat(q.x,q.y,q.z,n)
                if solidCells[left] == solidCells[right] { continue }
                let fluid = solidCells[right] ? left : right
                var away = solidCells[right] ? a : q
                away[axis] += solidCells[right] ? -1 : 1
                var facePressure = pressurePascals[fluid]
                if away[axis] >= 0 && away[axis] < n[axis] {
                    let other = flat(away.x,away.y,away.z,n)
                    if !solidCells[other] { facePressure += 0.5*(pressurePascals[fluid]-pressurePascals[other]) }
                }
                var f = SIMD3<Double>.zero
                f[axis] = (solidCells[right] ? facePressure : -facePressure)*h*h
                force += f; moment += fluidCross(facePosition(axis,i,j,k)-momentReferencePointMetres,f)
            } } }
        }
        return try FluidLoads(forceNewtons: force,momentNewtonMetres: moment,referencePointMetres: momentReferencePointMetres)
    }
}
