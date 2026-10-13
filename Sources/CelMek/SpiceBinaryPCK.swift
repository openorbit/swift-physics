// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum SpiceBinaryPCKError: Error, Equatable {
    case invalidEpoch, invalidSegment, outOfCoverage, numericalRangeExceeded
    case unsupportedType(Int32)
    case missingFrameClass(Int32)
    case coverageGap(Int32)
    case descriptorLimit
    case malformedSegment(path: String, array: Int, reason: String)
    case frameTreeRequired, frameTreeUnavailable
}

/// Native type-2 3-1-3 Euler angles, not text-PCK pole RA/declination values.
public struct SpiceBinaryPCKAngles: Sendable {
    public let secondsTDB: Double
    public let eulerAnglesRadians: SIMD3<Double>
    public let ratesRadiansPerTDBSecond: SIMD3<Double>
    public let accelerationsRadiansPerTDBSecondSquared: SIMD3<Double>
    public let descriptor: SpicePCKDescriptor
    public let sourcePath: String
    public let segmentName: String
    public let arrayIndex: Int
    public let recordIndex: Int
}
/// Validated type 2 directory and bounded on-demand record evaluation.
/// This value retains its reader; the underlying binary asset must remain unchanged.
public struct SpiceBinaryPCKSegment: Sendable {
    public let file: SpiceDAFFile
    public let arrayIndex: Int
    public let descriptor: SpicePCKDescriptor
    public let initialSecondsTDB: Double
    public let intervalSeconds: Double
    public let recordSizeWords: Int
    public let recordCount: Int
    public let polynomialDegree: Int

    public init(file: SpiceDAFFile, array index: Int) throws {
        guard file.metadata.kind == "PCK", file.metadata.arrays.indices.contains(index),
              let descriptor = file.metadata.arrays[index].pck else { throw SpiceBinaryPCKError.invalidSegment }
        guard descriptor.dataType == 2 else { throw SpiceBinaryPCKError.unsupportedType(descriptor.dataType) }
        let array = file.metadata.arrays[index], words = array.finalAddress - array.initialAddress + 1
        func malformed(_ reason: String) -> SpiceBinaryPCKError { .malformedSegment(path: file.url.path, array: index, reason: reason) }
        guard words >= 4 else { throw malformed("missing directory") }
        let directory = try file.readWords(array: index, offset: words - 4, count: 4)
        func integer(_ value: Double) throws -> Int {
            // DAF addresses already bound the possible count/size below Int32.max.
            guard value.isFinite, value >= 1, value <= Double(words - 4), value.rounded(.towardZero) == value else {
                throw malformed("invalid directory integer")
            }
            return Int(value)
        }
        let initial = directory[0], interval = directory[1]
        let size = try integer(directory[2]), count = try integer(directory[3])
        let components = 3
        guard size >= 2 + components, (size - 2) % components == 0,
              (words - 4) % size == 0, (words - 4) / size == count else { throw malformed("inconsistent record layout") }
        guard size <= file.limits.maximumReadWords else { throw SpiceDAFError.resourceLimit("PCK record words") }
        let end = initial + interval * Double(count)
        guard initial.isFinite, interval.isFinite, interval > 0, (interval/2) > 0,
              end.isFinite, end > initial, initial + interval > initial, end - interval < end else {
            throw malformed("invalid or unresolvable directory time grid")
        }
        guard descriptor.startSecondsTDB >= initial - pckTimeSlack(initial, interval),
              descriptor.endSecondsTDB <= end + pckTimeSlack(end, interval) else { throw malformed("descriptor exceeds record coverage") }
        self.file = file; self.arrayIndex = index; self.descriptor = descriptor
        self.initialSecondsTDB = initial; self.intervalSeconds = interval
        self.recordSizeWords = size; self.recordCount = count; self.polynomialDegree = (size - 2)/components - 1
    }

    /// Later record wins at an internal boundary; the final endpoint uses the last record.
    /// Only the selected record is read/validated. No extrapolation beyond descriptor coverage.
    public func angles(atSecondsTDB seconds: Double) throws -> SpiceBinaryPCKAngles {
        guard seconds.isFinite else { throw SpiceBinaryPCKError.invalidEpoch }
        guard seconds >= descriptor.startSecondsTDB, seconds <= descriptor.endSecondsTDB else { throw SpiceBinaryPCKError.outOfCoverage }
        func malformed(_ reason: String) -> SpiceBinaryPCKError { .malformedSegment(path: file.url.path, array: arrayIndex, reason: reason) }
        let ratio = (seconds - initialSecondsTDB)/intervalSeconds
        guard ratio.isFinite else { throw SpiceBinaryPCKError.numericalRangeExceeded }
        var record = Int(max(0, min(Double(recordCount - 1), floor(ratio))))
        // Correct a quotient rounded across one boundary using the directory's actual grid.
        if record + 1 < recordCount, seconds >= initialSecondsTDB + Double(record + 1)*intervalSeconds { record += 1 }
        if record > 0, seconds < initialSecondsTDB + Double(record)*intervalSeconds { record -= 1 }
        let values = try file.readWords(array: arrayIndex, offset: record * recordSizeWords, count: recordSizeWords)
        guard values.allSatisfy(\.isFinite), values[1] > 0 else { throw malformed("non-finite coefficients or invalid record radius") }
        let start = initialSecondsTDB + Double(record)*intervalSeconds
        let expectedMid = start + intervalSeconds/2
        guard abs(values[0] - expectedMid) <= pckTimeSlack(expectedMid, intervalSeconds),
              abs(values[1] - intervalSeconds/2) <= pckTimeSlack(intervalSeconds/2, intervalSeconds) else {
            throw malformed("record interval disagrees with directory")
        }
        let normalized = (seconds - values[0])/values[1]
        let slack = pckTimeSlack(values[0], intervalSeconds)/values[1] + 16*Double.ulpOfOne
        guard normalized.isFinite, normalized >= -1 - slack, normalized <= 1 + slack else { throw malformed("epoch outside selected coefficient interval") }
        let x = max(-1, min(1, normalized)), coefficientCount = polynomialDegree + 1
        var angles = SIMD3<Double>.zero, rates = SIMD3<Double>.zero, accelerations = SIMD3<Double>.zero
        for axis in 0..<3 {
            let polynomial = try pckPolynomial(values, offset: 2 + axis*coefficientCount, count: coefficientCount, x: x)
            angles[axis] = polynomial.value
            rates[axis] = polynomial.first/values[1]
            accelerations[axis] = polynomial.second/values[1]/values[1]
        }
        guard [angles,rates,accelerations].allSatisfy({ $0.x.isFinite && $0.y.isFinite && $0.z.isFinite }) else {
            throw SpiceBinaryPCKError.numericalRangeExceeded
        }
        return SpiceBinaryPCKAngles(secondsTDB: seconds, eulerAnglesRadians: angles,
            ratesRadiansPerTDBSecond: rates, accelerationsRadiansPerTDBSecondSquared: accelerations,
            descriptor: descriptor, sourcePath: file.url.path, segmentName: file.metadata.arrays[arrayIndex].name,
            arrayIndex: arrayIndex, recordIndex: record)
    }
}

private func pckTimeSlack(_ time: Double, _ interval: Double) -> Double { 8 * max(time.ulp, interval.ulp) }

// Clenshaw value and first/second differentiated recurrences, with full c0.
private func pckPolynomial(_ coefficients: [Double], offset: Int, count: Int, x: Double) throws -> (value: Double, first: Double, second: Double) {
    var b1 = 0.0, b2 = 0.0, d1 = 0.0, d2 = 0.0, e1 = 0.0, e2 = 0.0
    if count > 1 {
        for k in stride(from: count-1, through: 1, by: -1) {
            let b = coefficients[offset+k] + 2*x*b1 - b2
            let d = 2*b1 + 2*x*d1 - d2
            let e = 4*d1 + 2*x*e1 - e2
            guard b.isFinite, d.isFinite, e.isFinite else { throw SpiceBinaryPCKError.numericalRangeExceeded }
            b2 = b1; b1 = b; d2 = d1; d1 = d; e2 = e1; e1 = e
        }
    }
    let value = coefficients[offset] + x*b1 - b2, first = b1 + x*d1 - d2, second = 2*d1 + x*e1 - e2
    guard value.isFinite, first.isFinite, second.isFinite else { throw SpiceBinaryPCKError.numericalRangeExceeded }
    return (value,first,second)
}

extension SpiceDAFFile {
    public func pckAngles(array index: Int, atSecondsTDB seconds: Double) throws -> SpiceBinaryPCKAngles {
        try SpiceBinaryPCKSegment(file: self,array: index).angles(atSecondsTDB: seconds)
    }
}
extension SpiceKernelSnapshot {
    /// Later covering file/array wins; unsupported/corrupt selected data never falls back.
    public func pckAngles(frameClassID: Int32, atSecondsTDB seconds: Double,
                          maximumDescriptorChecks: Int = 1_000_000) throws -> SpiceBinaryPCKAngles {
        guard seconds.isFinite else { throw SpiceBinaryPCKError.invalidEpoch }
        guard frameClassID != 0 else { throw SpiceBinaryPCKError.invalidSegment }
        guard maximumDescriptorChecks > 0 else { throw SpiceBinaryPCKError.descriptorLimit }
        var checked = 0, found = false
        for kernel in binaryKernels.reversed() where kernel.file.metadata.kind == "PCK" {
            for index in kernel.file.metadata.arrays.indices.reversed() {
                guard checked < maximumDescriptorChecks else { throw SpiceBinaryPCKError.descriptorLimit }
                checked += 1
                guard let descriptor = kernel.file.metadata.arrays[index].pck, descriptor.frameClassID == frameClassID else { continue }
                found = true
                if seconds >= descriptor.startSecondsTDB, seconds <= descriptor.endSecondsTDB {
                    return try kernel.file.pckAngles(array: index,atSecondsTDB: seconds)
                }
            }
        }
        if found { throw SpiceBinaryPCKError.coverageGap(frameClassID) }
        throw SpiceBinaryPCKError.missingFrameClass(frameClassID)
    }
}

public struct SpiceBinaryPCKEvaluation {
    public let orientation: BodyOrientation
    public let nativeAngles: SpiceBinaryPCKAngles
    public let tdbCorrection: TDBCorrection
}

/// Binary type-2 orientation with explicit axes bindings and caller-selected coverage.
/// Retain frameTree externally when selected segments use different base axes.
public final class SpiceBinaryPCKOrientationProvider: OrientationProvider {
    public let referenceFrame: Frame
    public let frameClassID: Int32
    public let coverage: ClosedRange<Epoch>
    public let snapshot: SpiceKernelSnapshot
    public let provenance: PhysicalDataProvenance
    public let maximumDescriptorChecks: Int
    private let bindings: SpiceFrameBindings
    private let tdbModel: any TDBTimeModel
    private weak var tree: FrameTree?
    public init(snapshot: SpiceKernelSnapshot, frameClassID: Int32, referenceFrame: Frame,
                frameBindings: SpiceFrameBindings, frameTree: FrameTree? = nil,
                coverage: ClosedRange<Epoch>, tdbModel: any TDBTimeModel,
                provenance: PhysicalDataProvenance, maximumDescriptorChecks: Int = 1_000_000) throws {
        guard frameClassID != 0 else { throw SpiceBinaryPCKError.invalidSegment }
        guard maximumDescriptorChecks > 0 else { throw SpiceBinaryPCKError.descriptorLimit }
        guard referenceFrame.pseudoInertial else { throw SpiceTextProviderError.invalidReferenceFrame }
        guard tdbModel.coverage.contains(coverage.lowerBound), tdbModel.coverage.contains(coverage.upperBound) else { throw TDBModelError.outsideCoverage }
        for binding in frameBindings.bindings {
            guard binding.axesFrame.pseudoInertial else { throw SpiceTextProviderError.invalidReferenceFrame }
            if binding.axesFrame !== referenceFrame {
                guard let frameTree else { throw SpiceBinaryPCKError.frameTreeRequired }
                guard frameTree.isRegistered(binding.axesFrame), frameTree.isRegistered(referenceFrame) else {
                    throw FrameTransformError.disconnectedFrames
                }
            }
        }
        self.snapshot = snapshot; self.frameClassID = frameClassID; self.referenceFrame = referenceFrame
        self.bindings = frameBindings; self.tree = frameTree; self.coverage = coverage; self.tdbModel = tdbModel
        self.maximumDescriptorChecks = maximumDescriptorChecks
        self.provenance = try PhysicalDataProvenance(source: provenance.source + "; " + snapshot.binaryKernels.filter { $0.file.metadata.kind == "PCK" }.map { $0.source.path }.joined(separator: "; ") + "; time: " + tdbModel.provenance + "; axes: " + frameBindings.bindings.map { $0.provenance.source + " (" + $0.provenance.version + ")" }.joined(separator: "; "), version: provenance.version,reference: provenance.reference)
        _ = try evaluation(at: coverage.lowerBound)
        _ = try evaluation(at: coverage.upperBound)
    }
    public func orientation(at epoch: Epoch) throws -> BodyOrientation { try evaluation(at: epoch).orientation }
    public func evaluation(at epoch: Epoch) throws -> SpiceBinaryPCKEvaluation {
        guard coverage.contains(epoch) else { throw SpiceBinaryPCKError.outOfCoverage }
        let clock = try epoch.tdbCorrection(using: tdbModel)
        let seconds = try epoch.offsetFromJ2000.adding(clock.tdbMinusTT).seconds
        let native = try snapshot.pckAngles(frameClassID: frameClassID,atSecondsTDB: seconds,maximumDescriptorChecks: maximumDescriptorChecks)
        let binding = try bindings.binding(for: native.descriptor.frameID)
        var transform = FrameTransform.identity(at: epoch)
        for (index, axis) in [SIMD3<Double>(0,0,1),SIMD3(1,0,0),SIMD3(0,0,1)].enumerated() {
            let first = native.ratesRadiansPerTDBSecond[index], second = native.accelerationsRadiansPerTDBSecondSquared[index]
            let part = try FrameTransform(epoch: epoch,
                rotation: Quaternion(angle: -native.eulerAnglesRadians[index].truncatingRemainder(dividingBy: 2 * .pi),axis: axis),
                angularVelocity: axis * (-first * clock.rate),
                angularAcceleration: axis * (-second * clock.rate * clock.rate - first * clock.acceleration))
            transform = try transform.followed(by: part)
        }
        if binding.axesFrame !== referenceFrame {
            guard let tree else { throw SpiceBinaryPCKError.frameTreeUnavailable }
            let axes = try tree.transform(from: referenceFrame,to: binding.axesFrame,at: epoch)
            let rotationOnly = try FrameTransform(epoch: epoch,rotation: axes.rotation,
                angularVelocity: axes.angularVelocity,angularAcceleration: axes.angularAcceleration)
            transform = try rotationOnly.followed(by: transform)
        }
        let orientation = try BodyOrientation(referenceFrame: referenceFrame,epoch: epoch,rotation: transform.rotation,
            angularVelocity: transform.angularVelocity,angularAcceleration: transform.angularAcceleration)
        return SpiceBinaryPCKEvaluation(orientation: orientation,nativeAngles: native,tdbCorrection: clock)
    }
}
