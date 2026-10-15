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

public enum SpiceSPKResolutionError: Error, Equatable, Sendable {
    case invalidOptions
    case invalidFrame
    case missingLink(bodyID: Int32, secondsTDB: Double)
    case coverageGap(bodyID: Int32, secondsTDB: Double)
    case cycle([Int32])
    case chainLimit(Int)
    case descriptorLimit(Int)
    case frameMismatch(bodyID: Int32, expected: Int32, actual: Int32)
}
public struct SpiceSPKResolutionOptions: Sendable {
    /// Bounds selected metadata edges explored across both branches, including unused ancestors.
    public let maximumChainLinks: Int
    /// Bounds descriptor inspections, including those that do not match the requested body.
    public let maximumDescriptorChecks: Int
    public init(maximumChainLinks: Int = 64, maximumDescriptorChecks: Int = 1_000_000) throws {
        guard maximumChainLinks > 0, maximumDescriptorChecks > 0 else { throw SpiceSPKResolutionError.invalidOptions }
        self.maximumChainLinks = maximumChainLinks; self.maximumDescriptorChecks = maximumDescriptorChecks
    }
}
public enum SpiceSPKChainSide: Equatable, Sendable { case target, centre }
public struct SpiceSPKChainContribution: Sendable {
    /// Add target-branch states and subtract centre-branch states.
    public let side: SpiceSPKChainSide
    /// Index in the snapshot's binaryKernels array, distinguishing repeated loads.
    public let kernelIndex: Int
    public let state: SpiceSPKState
}
/// Necessary primitive links, with independent native axes and no implicit frame relabelling.
public struct SpiceSPKChain: Sendable {
    public let targetID: Int32
    public let centreID: Int32
    public let secondsTDB: Double
    public let commonCentreID: Int32
    public let contributions: [SpiceSPKChainContribution]
    public let cacheScope: EphemerisCacheScope
    public let exploredChainLinks: Int
    public let inspectedDescriptors: Int
}
/// Geometric target-minus-centre state in explicitly requested, common native axes.
/// Kilometres and TDB argument-second derivatives; no light-time/aberration corrections.
public struct SpiceSPKResolvedState: Sendable {
    public let targetID: Int32
    public let centreID: Int32
    public let frameID: Int32
    public let secondsTDB: Double
    public let positionKilometres: SIMD3<Double>
    public let velocityKilometresPerTDBSecond: SIMD3<Double>
    public let commonCentreID: Int32
    public let contributions: [SpiceSPKChainContribution]
    public let cacheScope: EphemerisCacheScope
    public let exploredChainLinks: Int
    public let inspectedDescriptors: Int
}

struct SelectedSPKSegment {
    let file: SpiceDAFFile
    let kernelIndex: Int
    let arrayIndex: Int
    let descriptor: SpiceSPKDescriptor
}
enum SPKSelectionLookup { case selected(SelectedSPKSegment), missing, gap }
extension SpiceKernelSnapshot {
    func lookupSPK(targetID: Int32, seconds: Double, remainingChecks: inout Int, checkLimit: Int) throws -> SPKSelectionLookup {
        var hasBody = false
        for kernelIndex in binaryKernels.indices.reversed() {
            let file = binaryKernels[kernelIndex].file
            guard file.metadata.kind == "SPK" else { continue }
            for arrayIndex in file.metadata.arrays.indices.reversed() {
                guard remainingChecks > 0 else { throw SpiceSPKResolutionError.descriptorLimit(checkLimit) }
                remainingChecks -= 1
                guard let descriptor = file.metadata.arrays[arrayIndex].spk, descriptor.targetID == targetID else { continue }
                hasBody = true
                guard seconds >= descriptor.startSecondsTDB, seconds <= descriptor.endSecondsTDB else { continue }
                return .selected(SelectedSPKSegment(file: file, kernelIndex: kernelIndex, arrayIndex: arrayIndex, descriptor: descriptor))
            }
        }
        return hasBody ? .gap : .missing
    }

    /// Resolves two precedence-selected centre chains at the same epoch.
    /// Returned links retain independent native frame IDs for explicit downstream conversion.
    public func spkChain(targetID: Int32, centreID: Int32,
                         atSecondsTDB seconds: Double, options: SpiceSPKResolutionOptions? = nil) throws -> SpiceSPKChain {
        guard seconds.isFinite else { throw SpiceSPKError.invalidEpoch }
        let options = try options ?? SpiceSPKResolutionOptions()
        var target = SPKBranch(targetID), centre = SPKBranch(centreID)
        var links = 0, remainingChecks = options.maximumDescriptorChecks
        func advance(_ branch: inout SPKBranch) throws {
            guard branch.stop == nil else { return }
            let body = branch.nodes.last!
            if body == 0 { branch.stop = .root; return }
            guard links < options.maximumChainLinks else { throw SpiceSPKResolutionError.chainLimit(options.maximumChainLinks) }
            switch try lookupSPK(targetID: body, seconds: seconds, remainingChecks: &remainingChecks, checkLimit: options.maximumDescriptorChecks) {
            case .missing: branch.stop = .missing(body)
            case .gap: branch.stop = .gap(body)
            case .selected(let segment):
                links += 1
                let next = segment.descriptor.centreID
                if let index = branch.indices[next] {
                    branch.stop = .cycle(Array(branch.nodes[index...]) + [next])
                } else {
                    branch.edges.append(segment); branch.indices[next] = branch.nodes.count; branch.nodes.append(next)
                }
            }
        }
        var common: Int32? = targetID == centreID ? targetID : nil
        // Alternating metadata traversal can meet before an absolute SSB state is available.
        // Payloads on shared ancestors are never evaluated.
        while common == nil, target.stop == nil || centre.stop == nil {
            try advance(&target)
            if let node = target.nodes.last, centre.indices[node] != nil { common = node; break }
            try advance(&centre)
            if let node = centre.nodes.last, target.indices[node] != nil { common = node; break }
        }
        guard let common else {
            // Cycles take priority over a missing/disconnected opposite branch.
            for stop in [target.stop,centre.stop] { if case .cycle(let path)? = stop { throw SpiceSPKResolutionError.cycle(path) } }
            for stop in [target.stop,centre.stop] {
                if case .gap(let body)? = stop { throw SpiceSPKResolutionError.coverageGap(bodyID: body,secondsTDB: seconds) }
            }
            for stop in [target.stop,centre.stop] {
                if case .missing(let body)? = stop { throw SpiceSPKResolutionError.missingLink(bodyID: body,secondsTDB: seconds) }
            }
            throw SpiceSPKResolutionError.missingLink(bodyID: targetID,secondsTDB: seconds)
        }
        let targetEdges = target.edges.prefix(target.indices[common]!)
        let centreEdges = centre.edges.prefix(centre.indices[common]!)
        var contributions: [SpiceSPKChainContribution] = []
        for (side,edges) in [(SpiceSPKChainSide.target,targetEdges),(.centre,centreEdges)] {
            for edge in edges {
                let state = try edge.file.spkState(array: edge.arrayIndex,atSecondsTDB: seconds)
                contributions.append(SpiceSPKChainContribution(side: side,kernelIndex: edge.kernelIndex,state: state))
            }
        }
        return SpiceSPKChain(targetID: targetID,centreID: centreID,secondsTDB: seconds,commonCentreID: common,
            contributions: contributions,cacheScope: cacheScope,exploredChainLinks: links,
            inspectedDescriptors: options.maximumDescriptorChecks - remainingChecks)
    }

    /// Composes a chain only when every contributing native frame ID equals frameID.
    public func spkState(targetID: Int32, centreID: Int32, frameID: Int32 = 1,
                         atSecondsTDB seconds: Double, options: SpiceSPKResolutionOptions? = nil) throws -> SpiceSPKResolvedState {
        guard seconds.isFinite else { throw SpiceSPKError.invalidEpoch }
        guard frameID != 0 else { throw SpiceSPKResolutionError.invalidFrame }
        let chain = try spkChain(targetID: targetID,centreID: centreID,atSecondsTDB: seconds,options: options)
        var p = SPKCompensatedVector(), v = SPKCompensatedVector()
        for term in chain.contributions {
            let state = term.state
            guard state.descriptor.frameID == frameID else {
                throw SpiceSPKResolutionError.frameMismatch(bodyID: state.descriptor.targetID,expected: frameID,actual: state.descriptor.frameID)
            }
            let sign = term.side == .target ? 1.0 : -1.0
            try p.add(sign*state.positionKilometres); try v.add(sign*state.velocityKilometresPerTDBSecond)
        }
        return try SpiceSPKResolvedState(targetID: targetID,centreID: centreID,frameID: frameID,secondsTDB: seconds,
            positionKilometres: p.value(),velocityKilometresPerTDBSecond: v.value(),commonCentreID: chain.commonCentreID,
            contributions: chain.contributions,cacheScope: cacheScope,exploredChainLinks: chain.exploredChainLinks,
            inspectedDescriptors: chain.inspectedDescriptors)
    }
}
private enum SPKBranchStop: Equatable { case root, missing(Int32), gap(Int32), cycle([Int32]) }
private struct SPKBranch {
    var nodes: [Int32]
    var indices: [Int32:Int]
    var edges: [SelectedSPKSegment] = []
    var stop: SPKBranchStop?
    init(_ body: Int32) { nodes = [body]; indices = [body:0] }
}
struct SPKCompensatedVector {
    var sum = SIMD3<Double>.zero, correction = SIMD3<Double>.zero
    mutating func add(_ term: SIMD3<Double>) throws {
        for axis in 0..<3 {
            let next = sum[axis] + term[axis]
            correction[axis] += abs(sum[axis]) >= abs(term[axis]) ? (sum[axis]-next)+term[axis] : (term[axis]-next)+sum[axis]
            sum[axis] = next
            guard next.isFinite, correction[axis].isFinite else { throw SpiceSPKError.numericalRangeExceeded }
        }
    }
    func value() throws -> SIMD3<Double> {
        let value = sum + correction
        guard (0..<3).allSatisfy({ value[$0].isFinite }) else { throw SpiceSPKError.numericalRangeExceeded }
        return value
    }
}
