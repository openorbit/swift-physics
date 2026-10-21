// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math

public enum FiniteBurnPositionTargetingStatus: Equatable {
    case positionConverged, evaluationLimitReached, stepLimitReached, singularSensitivity
}
public struct FiniteBurnPositionTargetingResult {
    public let status: FiniteBurnPositionTargetingStatus
    public let evaluations, rejectedTrials: Int
    public let initial, best: FiniteTransferEvaluation
}

/// Differential correction for one fixed-ignition, constant inertial burn.
/// Adjusts duration and direction using a numerically differentiated position
/// sensitivity and bounded backtracking. Velocity is diagnostic, NOT targeted.
/// Callers must independently verify clearance, attitude, resources and velocity.
public enum FiniteBurnPositionTargeting {
    public static func refine(initialState: SpacecraftState, initialBurn: FiniteBurn,
        target: CartesianState, context: ForceEvaluationContext, backgroundModels: [any ForceModel] = [],
        positionToleranceMetres: Double = 1, minimumRemainingPropellantKilograms: Double = 0,
        maximumEvaluations: Int = 120, burnWindow: ClosedRange<Epoch>? = nil,
        numericalOptions: NumericalPropagationOptions = .preset(for:.high)) throws -> FiniteBurnPositionTargetingResult {
        guard target.epoch > initialState.epoch, initialBurn.start >= initialState.epoch,
            initialBurn.end <= target.epoch, initialBurn.thrust.thrustNewtons > 0,
            positionToleranceMetres.isFinite, positionToleranceMetres > 0,
            minimumRemainingPropellantKilograms.isFinite, minimumRemainingPropellantKilograms >= 0,
            (1...10_000).contains(maximumEvaluations) else { throw FiniteTransferPlanningError.invalidInput }
        guard target.frame === initialState.frame, initialBurn.thrust.frame === initialState.frame,
            context.frame === initialState.frame else { throw SpacecraftDynamicsError.frameMismatch }
        let thrust = initialBurn.thrust
        var evaluations = 0, rejected = 0
        func norm(_ v: SIMD3<Double>) -> Double { hypot(hypot(v.x,v.y),v.z) }
        func cross(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> SIMD3<Double> {
            SIMD3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x)
        }
        func evaluate(_ control: SIMD3<Double>) throws -> FiniteTransferEvaluation? {
            let seconds = norm(control)
            guard seconds.isFinite, seconds > 1e-8 else { rejected += 1; return nil }
            let end = try initialBurn.start.advanced(by:Duration(seconds:seconds))
            let used = seconds*thrust.propellantConsumptionKilogramsPerSecond
            guard end <= target.epoch, initialState.mass.propellantMassKilograms-used > 0,
                initialState.mass.propellantMassKilograms-used >= minimumRemainingPropellantKilograms,
                burnWindow.map({ $0.contains(initialBurn.start) && $0.contains(end) }) ?? true else { rejected += 1; return nil }
            let burn = try FiniteBurn(start:initialBurn.start,end:end,thrust:FiniteThrust(frame:thrust.frame,
                thrustNewtons:thrust.thrustNewtons,effectiveExhaustVelocityMetresPerSecond:thrust.effectiveExhaustVelocityMetresPerSecond,direction:control))
            let schedule = try ScheduledThrust(frame:initialState.frame,burns:[burn])
            evaluations += 1
            do {
                let state = try ScheduledThrustPropagator(initialState:initialState,schedule:schedule,
                    backgroundModels:backgroundModels,context:context,options:numericalOptions).spacecraftState(at:target.epoch)
                return FiniteTransferEvaluation(schedule:schedule,arrivalState:state,
                    positionErrorMetres:norm(state.cartesian.position-target.position),
                    velocityErrorMetresPerSecond:norm(state.cartesian.velocity-target.velocity),
                    propellantUsedKilograms:initialState.mass.propellantMassKilograms-state.mass.propellantMassKilograms)
            } catch is NumericalPropagationError { rejected += 1; return nil }
        }
        var control = thrust.direction * (try initialBurn.end.duration(since:initialBurn.start).seconds)
        guard let initial = try evaluate(control) else { throw FiniteTransferPlanningError.invalidInput }
        var best = initial, status = FiniteBurnPositionTargetingStatus.evaluationLimitReached
        while best.positionErrorMetres > positionToleranceMetres {
            guard evaluations+7 <= maximumEvaluations else { break }
            let step = max(1e-4,norm(control)*1e-4)
            var columns: [SIMD3<Double>] = []
            for axis in 0..<3 {
                var plus = control, minus = control
                plus[axis] += step; minus[axis] -= step
                let a = try evaluate(plus), b = try evaluate(minus)
                if let a, let b { columns.append((a.arrivalState.cartesian.position-b.arrivalState.cartesian.position)/(2*step)) }
                else if let a { columns.append((a.arrivalState.cartesian.position-best.arrivalState.cartesian.position)/step) }
                else if let b { columns.append((best.arrivalState.cartesian.position-b.arrivalState.cartesian.position)/step) }
                else { break }
            }
            guard columns.count == 3 else { status = .stepLimitReached; break }
            let lengths = columns.map(norm)
            guard lengths.allSatisfy({$0.isFinite && $0 > 1e-15}) else { status = .singularSensitivity; break }
            let a = columns[0]/lengths[0], b = columns[1]/lengths[1], c = columns[2]/lengths[2]
            let determinant = (a*cross(b,c)).sum()
            guard determinant.isFinite, abs(determinant) > 1e-10 else { status = .singularSensitivity; break }
            let residual = target.position-best.arrivalState.cartesian.position
            var correction = SIMD3((residual*cross(b,c)).sum()/lengths[0],
                (a*cross(residual,c)).sum()/lengths[1],(a*cross(b,residual)).sum()/lengths[2])/determinant
            let length = norm(correction)
            guard length.isFinite else { status = .singularSensitivity; break }
            correction *= min(1,max(0.01,norm(control)*0.2)/max(length,1e-30))
            var accepted = false, fraction = 1.0
            for _ in 0..<10 where evaluations < maximumEvaluations {
                let candidate = control+correction*fraction
                if let trial = try evaluate(candidate), trial.positionErrorMetres < best.positionErrorMetres {
                    control = candidate; best = trial; accepted = true; break
                }
                fraction *= 0.5
            }
            guard accepted else { status = evaluations >= maximumEvaluations ? .evaluationLimitReached : .stepLimitReached; break }
        }
        if best.positionErrorMetres <= positionToleranceMetres { status = .positionConverged }
        return FiniteBurnPositionTargetingResult(status:status,evaluations:evaluations,rejectedTrials:rejected,initial:initial,best:best)
    }
}

