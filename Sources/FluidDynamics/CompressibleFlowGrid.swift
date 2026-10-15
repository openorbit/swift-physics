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

public enum CompressibleBoundary: Sendable {
    case periodic
    /// Zero-gradient ghost state; not a nonreflecting far-field boundary.
    case transmissive
    /// Stationary adiabatic slip wall, reversing normal momentum.
    case reflecting
    case fixed(ConservedFluidState)
}
public struct CompressibleFlowConfiguration: Sendable {
    public let grid: CartesianFluidGrid
    public let equationOfState: any CompressibleEquationOfState
    /// Low/high boundaries for X, then Y, then Z. Periodic boundaries must be paired.
    public let boundaries: [CompressibleBoundary]
    public let courantNumber: Double
    public init(grid: CartesianFluidGrid, equationOfState: any CompressibleEquationOfState,
                boundaries: [CompressibleBoundary] = Array(repeating: .transmissive,count: 6), courantNumber: Double = 0.4) throws {
        guard boundaries.count == 6, courantNumber.isFinite, courantNumber > 0, courantNumber <= 0.45 else { throw FluidError.invalidParameter }
        for axis in 0..<3 {
            let a: Bool, b: Bool
            if case .periodic = boundaries[2*axis] { a = true } else { a = false }
            if case .periodic = boundaries[2*axis+1] { b = true } else { b = false }
            guard a == b else { throw FluidError.invalidParameter }
        }
        for boundary in boundaries { if case .fixed(let state) = boundary { _ = try state.primitive(using: equationOfState) } }
        self.grid = grid; self.equationOfState = equationOfState; self.boundaries = boundaries; self.courantNumber = courantNumber
    }
}
public struct ConservedFluidTotals: Sendable {
    public let massKilograms: Double
    public let momentumKilogramMetresPerSecond: SIMD3<Double>
    public let energyJoules: Double
}
public struct CompressibleStepReport: Sendable {
    public let substeps: Int
    public let maximumCourantNumber: Double
    public let maximumSignalSpeedMetresPerSecond: Double
    /// Net outward domain flux integrated over this call, including wall impulse.
    public let boundaryMassOutflowKilograms: Double
    public let boundaryMomentumOutflowKilogramMetresPerSecond: SIMD3<Double>
    public let boundaryEnergyOutflowJoules: Double
    /// Step-average reaction from reflected solid-face momentum flux, not heat/shear.
    public let averageLoadsOnSolids: FluidLoads
}

/// First-order conservative finite-volume Euler solver with local Lax–Friedrichs
/// (Rusanov) fluxes. Captures shocks/expansion waves; evolves mass, three momentum
/// components and total energy. No viscous wall heat flux, turbulence or chemistry
/// solver. EOS and continuum applicability must be established separately.
public struct CompressibleFlowGrid: FluidVelocityField {
    public let configuration: CompressibleFlowConfiguration
    public let solidCells: [Bool]
    public private(set) var cells: [ConservedFluidState]
    public init(configuration: CompressibleFlowConfiguration,
                solidAtCellCentre: (SIMD3<Double>) throws -> Bool = { _ in false },
                initialState: (SIMD3<Double>) throws -> ConservedFluidState) throws {
        self.configuration = configuration
        let grid = configuration.grid, n = grid.dimensions
        var solids: [Bool] = [], states: [ConservedFluidState] = []
        for k in 0..<n.z { for j in 0..<n.y { for i in 0..<n.x {
            let p = grid.centre(SIMD3(i,j,k)), solid = try solidAtCellCentre(p)
            guard !solid || (i > 0 && j > 0 && k > 0 && i < n.x-1 && j < n.y-1 && k < n.z-1) else { throw FluidError.invalidGeometry }
            let state = try initialState(p); _ = try state.primitive(using: configuration.equationOfState)
            solids.append(solid); states.append(state)
        } } }
        self.solidCells = solids; self.cells = states
    }
    public func state(at positionMetres: SIMD3<Double>) throws -> CompressiblePrimitiveState {
        let index = try configuration.grid.cellIndex(at: positionMetres)
        guard !solidCells[index] else { throw FluidError.unsupportedQuery }
        return try cells[index].primitive(using: configuration.equationOfState)
    }
    /// Piecewise-constant sampling, matching the initial first-order spatial scheme.
    public func velocityMetresPerSecond(at positionMetres: SIMD3<Double>) throws -> SIMD3<Double> {
        try state(at: positionMetres).velocityMetresPerSecond
    }
    public func conservedTotals() throws -> ConservedFluidTotals {
        var total = SIMD8<Double>.zero
        for i in cells.indices where !solidCells[i] { total += cells[i].vector*configuration.grid.cellVolumeCubicMetres }
        guard (0..<5).allSatisfy({ total[$0].isFinite }) else { throw FluidError.nonFiniteResult }
        return ConservedFluidTotals(massKilograms: total[0],momentumKilogramMetresPerSecond: SIMD3(total[1],total[2],total[3]),energyJoules: total[4])
    }
    private func reflected(_ state: ConservedFluidState, axis: Int) throws -> ConservedFluidState {
        var vector = state.vector; vector[axis+1] = -vector[axis+1]; return try ConservedFluidState(vector: vector)
    }
    private func ghost(_ adjacent: ConservedFluidState, boundary: CompressibleBoundary, axis: Int) throws -> ConservedFluidState {
        switch boundary {
        case .transmissive: return adjacent
        case .fixed(let state): return state
        case .reflecting: return try reflected(adjacent,axis: axis)
        case .periodic: throw FluidError.invalidParameter // Periodic neighbours are explicitly wrapped.
        }
    }
    private func physicalFlux(_ state: ConservedFluidState, primitive: CompressiblePrimitiveState, axis: Int) -> SIMD8<Double> {
        let u = primitive.velocityMetresPerSecond[axis], pressure = primitive.thermodynamics.pressurePascals
        var flux = state.vector*u
        flux[axis+1] += pressure; flux[4] = (state.totalEnergyDensityJoulesPerCubicMetre+pressure)*u
        return flux
    }
    private func numericalFlux(_ left: ConservedFluidState, _ right: ConservedFluidState,
                               lp: CompressiblePrimitiveState, rp: CompressiblePrimitiveState, axis: Int) -> SIMD8<Double> {
        let signal = max(abs(lp.velocityMetresPerSecond[axis])+lp.thermodynamics.speedOfSoundMetresPerSecond,
                         abs(rp.velocityMetresPerSecond[axis])+rp.thermodynamics.speedOfSoundMetresPerSecond)
        return 0.5*(physicalFlux(left,primitive: lp,axis: axis)+physicalFlux(right,primitive: rp,axis: axis))
            - (0.5*signal)*(right.vector-left.vector)
    }
    /// All updates are transactional. A budget/positivity/EOS failure leaves the
    /// original grid unchanged. No density, pressure or energy floors repair the state.
    public mutating func advance(timeStepSeconds: Double, maximumSubsteps: Int = 128,
                                 momentReferencePointMetres: SIMD3<Double> = .zero) throws -> CompressibleStepReport {
        guard timeStepSeconds.isFinite, timeStepSeconds > 0, (1...100000).contains(maximumSubsteps),
              fluidFinite(momentReferencePointMetres) else { throw FluidError.invalidParameter }
        let grid = configuration.grid, n = grid.dimensions, h = grid.cellSizeMetres, area = h*h, eos = configuration.equationOfState
        var next = self, remaining = timeStepSeconds, substeps = 0, maximumCFL = 0.0, maxSignal = 0.0
        var boundaryIntegral = SIMD8<Double>.zero, solidImpulse = SIMD3<Double>.zero, solidAngularImpulse = SIMD3<Double>.zero
        while remaining > 0 {
            guard substeps < maximumSubsteps else { throw FluidError.workBudgetExceeded }
            let primitives = try next.cells.map { try $0.primitive(using: eos) }
            var signals = SIMD3<Double>.zero
            for p in primitives {
                for axis in 0..<3 { signals[axis] = max(signals[axis],abs(p.velocityMetresPerSecond[axis])+p.thermodynamics.speedOfSoundMetresPerSecond) }
            }
            for boundary in configuration.boundaries {
                if case .fixed(let state) = boundary {
                    let p = try state.primitive(using: eos)
                    for axis in 0..<3 { signals[axis] = max(signals[axis],abs(p.velocityMetresPerSecond[axis])+p.thermodynamics.speedOfSoundMetresPerSecond) }
                }
            }
            let sum = signals.sum(), dt = min(remaining,configuration.courantNumber*h/sum)
            guard sum.isFinite, dt.isFinite, dt > 0, remaining-dt < remaining else { throw FluidError.nonFiniteResult }
            maximumCFL = max(maximumCFL,dt*sum/h); maxSignal = max(maxSignal,max(signals.x,signals.y,signals.z))
            var change = [SIMD8<Double>](repeating: .zero,count: next.cells.count)
            for axis in 0..<3 {
                var shape = n; shape[axis] += 1
                for k in 0..<shape.z { for j in 0..<shape.y { for i in 0..<shape.x {
                    let face = SIMD3(i,j,k)
                    var lq = face; lq[axis] -= 1
                    let leftIndex: Int? = lq[axis] >= 0 ? grid.index(lq) : nil
                    let rightIndex: Int? = face[axis] < n[axis] ? grid.index(face) : nil
                    if let l = leftIndex, let r = rightIndex, solidCells[l] && solidCells[r] { continue }
                    var left: ConservedFluidState, right: ConservedFluidState
                    if let l = leftIndex { left = next.cells[l] }
                    else {
                        let adjacent = next.cells[rightIndex!]
                        if case .periodic = configuration.boundaries[2*axis] {
                            var wrapped = face; wrapped[axis] = n[axis]-1; left = next.cells[grid.index(wrapped)]
                        } else { left = try ghost(adjacent,boundary: configuration.boundaries[2*axis],axis: axis) }
                    }
                    if let r = rightIndex { right = next.cells[r] }
                    else {
                        let adjacent = next.cells[leftIndex!]
                        if case .periodic = configuration.boundaries[2*axis+1] {
                            var wrapped = face; wrapped[axis] = 0; right = next.cells[grid.index(wrapped)]
                        } else { right = try ghost(adjacent,boundary: configuration.boundaries[2*axis+1],axis: axis) }
                    }
                    let leftSolid = leftIndex.map { solidCells[$0] } ?? false, rightSolid = rightIndex.map { solidCells[$0] } ?? false
                    if leftSolid { left = try reflected(right,axis: axis) }
                    if rightSolid { right = try reflected(left,axis: axis) }
                    let lp = !leftSolid && leftIndex != nil ? primitives[leftIndex!] : try left.primitive(using: eos)
                    let rp = !rightSolid && rightIndex != nil ? primitives[rightIndex!] : try right.primitive(using: eos)
                    let flux = numericalFlux(left,right,lp: lp,rp: rp,axis: axis)
                    guard (0..<5).allSatisfy({ flux[$0].isFinite }) else { throw FluidError.nonFiniteResult }
                    if let l = leftIndex, !leftSolid { change[l] -= flux*(dt/h) }
                    if let r = rightIndex, !rightSolid { change[r] += flux*(dt/h) }
                    if leftIndex == nil { boundaryIntegral -= flux*(area*dt) }
                    if rightIndex == nil { boundaryIntegral += flux*(area*dt) }
                    if leftSolid != rightSolid {
                        let impulse = SIMD3(flux[1],flux[2],flux[3])*(area*dt*(rightSolid ? 1 : -1))
                        var position = grid.centre(face); position[axis] -= h/2
                        solidImpulse += impulse; solidAngularImpulse += fluidCross(position-momentReferencePointMetres,impulse)
                    }
                } } }
            }
            for index in next.cells.indices where !solidCells[index] {
                let updated = try ConservedFluidState(vector: next.cells[index].vector+change[index])
                _ = try updated.primitive(using: eos); next.cells[index] = updated
            }
            substeps += 1; remaining -= dt
        }
        let loads = try FluidLoads(forceNewtons: solidImpulse/timeStepSeconds,
            momentNewtonMetres: solidAngularImpulse/timeStepSeconds,referencePointMetres: momentReferencePointMetres)
        guard (0..<5).allSatisfy({ boundaryIntegral[$0].isFinite }) else { throw FluidError.nonFiniteResult }
        self = next
        return CompressibleStepReport(substeps: substeps,maximumCourantNumber: maximumCFL,maximumSignalSpeedMetresPerSecond: maxSignal,
            boundaryMassOutflowKilograms: boundaryIntegral[0],boundaryMomentumOutflowKilogramMetresPerSecond: SIMD3(boundaryIntegral[1],boundaryIntegral[2],boundaryIntegral[3]),
            boundaryEnergyOutflowJoules: boundaryIntegral[4],averageLoadsOnSolids: loads)
    }
}
