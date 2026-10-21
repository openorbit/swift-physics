// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

/// Two fixed-ignition inertial burns corrected against a full terminal state.
/// Six duration-vector components are solved with finite differences and
/// pivoted elimination. `relativeControlStep` scales the duration-vector finite
/// differences (with a 1 ms floor); a smaller value can resolve long-burn
/// terminal encounters. Existing callers retain the 1e-4 default. Clearance and operational constraints remain caller checks.
public struct FiniteBurnStateTargetingReport {
    public let best: FiniteTransferEvaluation
    public let evaluations: Int
}

public enum FiniteBurnStateTargeting {
    public static func refine(initialState: SpacecraftState, initialSchedule: ScheduledThrust,
        target: CartesianState, context: ForceEvaluationContext, backgroundModels: [any ForceModel],
        positionToleranceMetres: Double = 10, velocityToleranceMetresPerSecond: Double = 0.01,
        minimumRemainingPropellantKilograms: Double = 0, maximumEvaluations: Int = 160,
        numericalOptions: NumericalPropagationOptions = .preset(for:.high),
        exhaustVelocityForDirection: ((SIMD3<Double>) throws -> Double)? = nil,
        relativeControlStep: Double = 1e-4) throws -> FiniteTransferEvaluation {
        try refineWithReport(initialState:initialState,initialSchedule:initialSchedule,target:target,context:context,
            backgroundModels:backgroundModels,positionToleranceMetres:positionToleranceMetres,
            velocityToleranceMetresPerSecond:velocityToleranceMetresPerSecond,minimumRemainingPropellantKilograms:minimumRemainingPropellantKilograms,
            maximumEvaluations:maximumEvaluations,numericalOptions:numericalOptions,exhaustVelocityForDirection:exhaustVelocityForDirection,relativeControlStep:relativeControlStep).best
    }
    public static func refineWithReport(initialState: SpacecraftState, initialSchedule: ScheduledThrust,
        target: CartesianState, context: ForceEvaluationContext, backgroundModels: [any ForceModel],
        positionToleranceMetres: Double = 10, velocityToleranceMetresPerSecond: Double = 0.01,
        minimumRemainingPropellantKilograms: Double = 0, maximumEvaluations: Int = 160,
        numericalOptions: NumericalPropagationOptions = .preset(for:.high),
        exhaustVelocityForDirection: ((SIMD3<Double>) throws -> Double)? = nil,
        relativeControlStep: Double = 1e-4) throws -> FiniteBurnStateTargetingReport {
        guard relativeControlStep.isFinite,relativeControlStep>0,relativeControlStep<=0.01,
            initialSchedule.burns.count == 2, (14...10_000).contains(maximumEvaluations),
            initialSchedule.frame === initialState.frame,
            positionToleranceMetres.isFinite, velocityToleranceMetresPerSecond.isFinite,
            minimumRemainingPropellantKilograms.isFinite,
            positionToleranceMetres > 0, velocityToleranceMetresPerSecond > 0,
            minimumRemainingPropellantKilograms >= 0,
            target.frame === initialState.frame, context.frame === initialState.frame,
            target.epoch > initialState.epoch else { throw FiniteTransferPlanningError.invalidInput }
        let seeds = initialSchedule.burns
        var controls = try seeds.flatMap { b -> [Double] in
            let v = b.thrust.direction * (try b.end.duration(since:b.start).seconds)
            return [v.x,v.y,v.z]
        }
        var count = 0
        func norm(_ v: SIMD3<Double>) -> Double { hypot(hypot(v.x,v.y),v.z) }
        func residual(_ e: FiniteTransferEvaluation) -> [Double] {
            let p = (e.arrivalState.cartesian.position-target.position)/positionToleranceMetres
            let v = (e.arrivalState.cartesian.velocity-target.velocity)/velocityToleranceMetresPerSecond
            return [p.x,p.y,p.z,v.x,v.y,v.z]
        }
        func score(_ e: FiniteTransferEvaluation) -> Double { residual(e).reduce(0) { $0+$1*$1 } }
        func evaluate(_ x: [Double]) throws -> FiniteTransferEvaluation? {
            guard count < maximumEvaluations else { return nil }
            var burns = [FiniteBurn](), fuel = 0.0
            for i in 0..<2 {
                let seed = seeds[i], vector = SIMD3(x[i*3],x[i*3+1],x[i*3+2]), duration = norm(vector)
                guard duration.isFinite, duration > 0.001 else { return nil }
                let end = try seed.start.advanced(by:Duration(seconds:duration))
                guard seed.start >= initialState.epoch, end <= (i == 0 ? seeds[1].start : target.epoch) else { return nil }
                let ve = try exhaustVelocityForDirection?(vector/duration) ?? seed.thrust.effectiveExhaustVelocityMetresPerSecond
                let thrust = try FiniteThrust(frame:initialState.frame,thrustNewtons:seed.thrust.thrustNewtons,
                    effectiveExhaustVelocityMetresPerSecond:ve,direction:vector)
                fuel += duration*thrust.propellantConsumptionKilogramsPerSecond
                burns.append(try FiniteBurn(start:seed.start,end:end,thrust:thrust))
            }
            guard fuel < initialState.mass.propellantMassKilograms,
                initialState.mass.propellantMassKilograms-fuel >= minimumRemainingPropellantKilograms else { return nil }
            count += 1
            let schedule = try ScheduledThrust(frame:initialState.frame,burns:burns)
            let state = try ScheduledThrustPropagator(initialState:initialState,schedule:schedule,backgroundModels:backgroundModels,
                context:context,options:numericalOptions).spacecraftState(at:target.epoch)
            return FiniteTransferEvaluation(schedule:schedule,arrivalState:state,positionErrorMetres:norm(state.cartesian.position-target.position),
                velocityErrorMetresPerSecond:norm(state.cartesian.velocity-target.velocity),propellantUsedKilograms:fuel)
        }
        guard var best = try evaluate(controls) else { throw FiniteTransferPlanningError.invalidInput }
        while best.positionErrorMetres > positionToleranceMetres || best.velocityErrorMetresPerSecond > velocityToleranceMetresPerSecond {
            guard count+13 < maximumEvaluations else { break }
            var matrix = Array(repeating:Array(repeating:0.0,count:7),count:6)
            let r = residual(best)
            for j in 0..<6 {
                let h = max(0.001,abs(controls[j])*relativeControlStep)
                var a = controls, b = controls; a[j] += h; b[j] -= h
                let ea = try evaluate(a), eb = try evaluate(b)
                guard ea != nil || eb != nil else { throw FiniteTransferPlanningError.invalidInput }
                if let ea, let eb {
                    let ra = residual(ea), rb = residual(eb)
                    for i in 0..<6 { matrix[i][j] = (ra[i]-rb[i])/(2*h) }
                } else if let ea {
                    let ra = residual(ea)
                    for i in 0..<6 { matrix[i][j] = (ra[i]-r[i])/h }
                } else if let eb {
                    let rb = residual(eb)
                    for i in 0..<6 { matrix[i][j] = (r[i]-rb[i])/h }
                }
            }
            for i in 0..<6 { matrix[i][6] = -r[i] }
            for j in 0..<6 {
                let pivot = (j..<6).max { abs(matrix[$0][j]) < abs(matrix[$1][j]) }!
                guard abs(matrix[pivot][j]) > 1e-12 else { throw FiniteTransferPlanningError.invalidInput }
                matrix.swapAt(j,pivot)
                let scale = matrix[j][j]
                for k in j...6 { matrix[j][k] /= scale }
                for i in 0..<6 where i != j {
                    let factor = matrix[i][j]
                    for k in j...6 { matrix[i][k] -= factor*matrix[j][k] }
                }
            }
            let correction = matrix.map { $0[6] }
            var fraction = 1.0, accepted = false
            for _ in 0..<12 {
                let x = zip(controls,correction).map { $0+$1*fraction }
                if let trial = try evaluate(x), score(trial) < score(best) {
                    controls = x; best = trial; accepted = true; break
                }
                fraction *= 0.5
            }
            if !accepted { break }
        }
        return FiniteBurnStateTargetingReport(best:best,evaluations:count)
    }
}
