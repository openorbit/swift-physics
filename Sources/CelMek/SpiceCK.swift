// SPDX-License-Identifier: Apache-2.0
import Foundation
@preconcurrency import Math

public enum SpiceCKError: Error, Equatable {
    case numericalRangeExceeded
    case invalidQuery, malformedSegment, unsupportedType(Int32), resourceLimit
    case missingInstrument, coverageGap, angularVelocityUnavailable
    case invalidReferenceFrame, frameTreeUnavailable
}
/// Type-1 pointing is a discrete sample; types 2/3 are evaluated continuously. A tolerance never changes its timestamp.
public struct SpiceCKPointing: Sendable {
    public let requestedEncodedTicks: Double
    public let encodedTicks: Double
    /// Hamilton quaternion mapping base-reference coordinates into instrument coordinates.
    public let referenceToInstrument: Quaternion
    /// Physical instrument angular velocity, expressed in base-reference axes, radians/second.
    public let angularVelocityReferenceRadiansPerSecond: SIMD3<Double>?
    /// Type-2 model seconds per encoded tick; nil for types 1 and 3.
    public let secondsPerEncodedTick: Double?
    public let descriptor: SpiceCKDescriptor
    public let sourcePath: String
    public let recordIndex: Int
}

/// Public CK type-1 layout; no interpolation or extrapolation of pointing records.
public struct SpiceCKSegment {
    public let descriptor: SpiceCKDescriptor
    private let file: SpiceDAFFile
    private let array: Int
    private let times: [Double]
    private let width: Int
    public init(file: SpiceDAFFile, array: Int, maximumRecords: Int = 100_000) throws {
        guard file.metadata.arrays.indices.contains(array), let descriptor = file.metadata.arrays[array].ck else {
            throw SpiceCKError.malformedSegment
        }
        guard descriptor.dataType == 1 else { throw SpiceCKError.unsupportedType(descriptor.dataType) }
        guard maximumRecords > 0 else { throw SpiceCKError.resourceLimit }
        let a = file.metadata.arrays[array], length = a.finalAddress - a.initialAddress + 1
        let last = try file.readWords(array: array, offset: length-1, count: 1)[0]
        guard last.isFinite, last >= 1, last <= Double(maximumRecords), last < Double(Int.max), last.rounded() == last else {
            throw SpiceCKError.resourceLimit
        }
        let count = Int(last), width = descriptor.hasAngularVelocity ? 7 : 4
        guard count <= (length-1)/(width+1), length == count*(width+1)+(count-1)/100+1 else {
            throw SpiceCKError.malformedSegment
        }
        var times: [Double] = []
        var offset = 0
        while offset < count {
            let n = min(count-offset,file.limits.maximumReadWords)
            times += try file.readWords(array: array,offset: count*width+offset,count: n)
            offset += n
        }
        guard times.allSatisfy({ $0.isFinite && $0 >= descriptor.startEncodedTicks && $0 <= descriptor.endEncodedTicks }),
              zip(times,times.dropFirst()).allSatisfy({ $0 < $1 }) else { throw SpiceCKError.malformedSegment }
        for i in 1...max(1,(count-1)/100) where i <= (count-1)/100 {
            let stored = try file.readWords(array: array,offset: count*(width+1)+i-1,count: 1)[0]
            let midpoint = times[i*100-1] + (times[i*100]-times[i*100-1])/2
            guard stored.isFinite, abs(stored-midpoint) <= 8*max(midpoint.ulp,stored.ulp) else { throw SpiceCKError.malformedSegment }
        }
        self.file = file; self.array = array; self.descriptor = descriptor; self.times = times; self.width = width
    }
    /// Equal-distance ties select the later record. Tolerance is in encoded ticks.
    public func pointing(atEncodedTicks ticks: Double, toleranceTicks: Double = 0) throws -> SpiceCKPointing? {
        guard ticks.isFinite, toleranceTicks.isFinite, toleranceTicks >= 0 else { throw SpiceCKError.invalidQuery }
        guard ticks >= descriptor.startEncodedTicks, ticks <= descriptor.endEncodedTicks else { return nil }
        var lo = 0, hi = times.count
        while lo < hi { let mid = lo+(hi-lo)/2; if times[mid] < ticks { lo = mid+1 } else { hi = mid } }
        var index = min(lo,times.count-1)
        if lo > 0, abs(times[lo-1]-ticks) < abs(times[index]-ticks) { index = lo-1 }
        guard abs(times[index]-ticks) <= toleranceTicks else { return nil }
        let words = try file.readWords(array: array,offset: index*width,count: width)
        guard words.allSatisfy(\.isFinite) else { throw SpiceCKError.malformedSegment }
        let norm = words.prefix(4).reduce(0) { $0+$1*$1 }
        guard abs(norm-1) <= 1e-10 else { throw SpiceCKError.malformedSegment }
        let q = try attitudeQuaternion(Quaternion(x: words[1],y: words[2],z: words[3],w: words[0]))
        return SpiceCKPointing(requestedEncodedTicks: ticks,encodedTicks: times[index],referenceToInstrument: q,
            angularVelocityReferenceRadiansPerSecond: width == 7 ? SIMD3(words[4],words[5],words[6]) : nil, secondsPerEncodedTick: nil,
            descriptor: descriptor,sourcePath: file.url.path,recordIndex: index)
    }
}

extension SpiceKernelSnapshot {
    /// Latest-file/latest-segment priority; a discrete gap permits an older segment.
    /// Selected unsupported or malformed data fails explicitly.
    public func ckPointing(instrumentID: Int32, atEncodedTicks ticks: Double, toleranceTicks: Double = 0,
                           maximumDescriptorChecks: Int = 1_000_000) throws -> SpiceCKPointing {
        guard instrumentID != 0, ticks.isFinite, toleranceTicks.isFinite, toleranceTicks >= 0,
              maximumDescriptorChecks > 0 else { throw SpiceCKError.invalidQuery }
        var found = false, checks = 0
        for kernel in binaryKernels.reversed() where kernel.file.metadata.kind == "CK" {
            for i in kernel.file.metadata.arrays.indices.reversed() {
                guard checks < maximumDescriptorChecks else { throw SpiceCKError.resourceLimit }; checks += 1
                guard let d = kernel.file.metadata.arrays[i].ck, d.instrumentID == instrumentID else { continue }
                found = true
                guard ticks >= d.startEncodedTicks, ticks <= d.endEncodedTicks else { continue }
                let sample: SpiceCKPointing?
                if d.dataType == 2 {
                    sample = try SpiceCKType2Segment(file: kernel.file,array: i).pointing(atEncodedTicks: ticks)
                } else if d.dataType == 3 {
                    sample = try SpiceCKType3Segment(file: kernel.file,array: i).pointing(atEncodedTicks: ticks,toleranceTicks: toleranceTicks)
                } else {
                    sample = try SpiceCKSegment(file: kernel.file,array: i).pointing(atEncodedTicks: ticks,toleranceTicks: toleranceTicks)
                }
                if let sample { return sample }
            }
        }
        throw found ? SpiceCKError.coverageGap : SpiceCKError.missingInstrument
    }
}

/// CK rates have seconds units; the caller explicitly identifies their time scale.
public enum SpiceCKAngularRateTimeScale: Sendable { case atomic, tdb }
/// Type-1 exact samples and type-2/3 continuous attitude. Type-1/3 acceleration is unavailable.
/// Type-3 stored rates are independent of the quaternion interpolation; they cannot
/// supply a derivative-consistent orientation/frame mapping through this adapter.
/// For nearest-sample rendering use ckPointing and preserve its returned clock timestamp.
public final class SpiceCKAttitudeProvider: AttitudeProvider {
    public let referenceFrame: Frame
    public let bodyFrame: Frame
    public let snapshot: SpiceKernelSnapshot
    public let clock: SpiceSCLKClock
    public let instrumentID: Int32
    private let bindings: SpiceFrameBindings
    private let model: (any TDBTimeModel)?
    private let rateScale: SpiceCKAngularRateTimeScale
    private weak var tree: FrameTree?
    public init(snapshot: SpiceKernelSnapshot, instrumentID: Int32, clock: SpiceSCLKClock,
                referenceFrame: Frame, bodyFrame: Frame, frameBindings: SpiceFrameBindings,
                angularRateTimeScale: SpiceCKAngularRateTimeScale, tdbModel: (any TDBTimeModel)? = nil,
                frameTree: FrameTree? = nil) throws {
        guard instrumentID != 0 else { throw SpiceCKError.invalidQuery }
        guard referenceFrame.pseudoInertial, bodyFrame !== referenceFrame,
              frameBindings.bindings.allSatisfy({ $0.axesFrame.pseudoInertial }) else { throw SpiceCKError.invalidReferenceFrame }
        self.snapshot = snapshot; self.instrumentID = instrumentID; self.clock = clock
        self.referenceFrame = referenceFrame; self.bodyFrame = bodyFrame; self.bindings = frameBindings
        self.rateScale = angularRateTimeScale; self.model = tdbModel; self.tree = frameTree
    }
    public func attitude(at epoch: Epoch) throws -> AttitudeState {
        let clockSample = try clock.sample(at: epoch,using: model)
        let native = try snapshot.ckPointing(instrumentID: instrumentID,atEncodedTicks: clockSample.encodedTicks)
        guard var rate = native.angularVelocityReferenceRadiansPerSecond else { throw SpiceCKError.angularVelocityUnavailable }
        var acceleration: SIMD3<Double>?
        if let secondsPerTick = native.secondsPerEncodedTick {
            var parallelRate = 1.0, parallelAcceleration = 0.0
            if clock.parallelTimeScale == .tdb {
                guard let model else { throw SpiceSCLKError.tdbModelRequired }
                let correction = try epoch.tdbCorrection(using: model)
                parallelRate = correction.rate; parallelAcceleration = correction.acceleration
            }
            let factor = secondsPerTick/clockSample.parallelSecondsPerTick
            acceleration = rate*(factor*parallelAcceleration)
            rate *= factor*parallelRate
        } else if case .tdb = rateScale {
            guard let model else { throw SpiceSCLKError.tdbModelRequired }
            rate *= try epoch.tdbCorrection(using: model).rate
        }
        let binding = try bindings.binding(for: native.descriptor.frameID)
        var axes = FrameTransform.identity(at: epoch)
        if binding.axesFrame !== referenceFrame {
            guard let tree else { throw SpiceCKError.frameTreeUnavailable }
            axes = try tree.transform(from: referenceFrame,to: binding.axesFrame,at: epoch)
        }
        let c = native.referenceToInstrument
        if let acceleration {
            let nativeTransform = try FrameTransform(epoch: epoch,rotation: c,
                angularVelocity: -(c.matrix*rate),angularAcceleration: -(c.matrix*acceleration))
            let composed = try axes.followed(by: nativeTransform)
            return try AttitudeState(epoch: epoch,referenceFrame: referenceFrame,bodyFrame: bodyFrame,
                bodyToReference: composed.rotation.conj,angularVelocityBodyRadiansPerSecond: -composed.angularVelocity,
                angularAccelerationBodyRadiansPerSecondSquared: -composed.angularAcceleration)
        }
        return try AttitudeState(epoch: epoch,referenceFrame: referenceFrame,bodyFrame: bodyFrame,
            bodyToReference: (c*axes.rotation).conj,
            angularVelocityBodyRadiansPerSecond: c.matrix*(rate-axes.angularVelocity),
            angularAccelerationBodyRadiansPerSecondSquared: nil)
    }
}
