// SPDX-License-Identifier: Apache-2.0
import Foundation
import Math
import Testing
@testable import CelMek

@Suite struct UnpoweredFlybyTests {
    let frame = Frame(name: "Body-centred flyby",pseudoInertial: true)
    func request(incoming: SIMD3<Double> = SIMD3(1,0,0), outgoing: SIMD3<Double> = SIMD3(0,1,0),
        altitude: Double = 0, tolerance: Double = 1e-12) throws -> FlybyEncounterRequest {
        try FlybyEncounterRequest(body: .earth,periapsisEpoch: .j2000,frame: frame,
            gravitationalParameter: GravitationalParameter(cubicMetresPerSecondSquared: 1,convention: .newtonian,
                provenance: PhysicalDataProvenance(source: "Analytic hyperbola",version: "1")),bodyRadiusMetres: 0.1,
            minimumPeriapsisAltitudeMetres: altitude,incomingExcessVelocityMetresPerSecond: incoming,
            outgoingExcessVelocityMetresPerSecond: outgoing,relativeSpeedTolerance: tolerance)
    }
    @Test func independentRightAngleHyperbolaAndFiniteStates() throws {
        let result = try UnpoweredFlyby.solve(request())
        // mu=vInf=p=1, e=sqrt(2), rp=sqrt(2)-1, vp=sqrt(2)+1.
        #expect(abs(result.turningAngleRadians - .pi/2) < 1e-14)
        #expect(abs(result.periapsisRadiusMetres-(sqrt(2)-1)) < 1e-14)
        #expect(abs(result.eccentricity-sqrt(2)) < 1e-14)
        #expect(abs(result.periapsisSpeedMetresPerSecond-(sqrt(2)+1)) < 1e-14)
        #expect(routeNorm(result.orbitalNormal-SIMD3(0,0,1)) < 1e-14)
        let r = SIMD3<Double>(1,-1,0)/sqrt(2), t = SIMD3<Double>(1,1,0)/sqrt(2)
        #expect(routeNorm(result.periapsisRadialDirection-r) < 1e-14)
        #expect(routeNorm(result.periapsisProgradeDirection-t) < 1e-14)
        // Hyperbolic anomaly at true anomaly pi/3: aMagnitude=1.
        let f = Double.pi/3, e = sqrt(2.0)
        let h = 2*atanh(sqrt((e-1)/(e+1))*tan(f/2)), time = e*sinh(h)-h
        for sign in [-1.0,1.0] {
            let epoch = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: sign*time))
            let actual = try result.trajectory.state(at: epoch)
            let angle = sign*f, radius = 1/(1+e*cos(angle))
            let expectedPosition = (r*cos(angle)+t*sin(angle))*radius
            let expectedVelocity = -r*sin(angle)+t*(e+cos(angle))
            #expect(routeNorm(actual.position-expectedPosition) < 1e-9)
            #expect(routeNorm(actual.velocity-expectedVelocity) < 1e-9)
        }
    }
    @Test func incomingOutgoingAsymptotesAndCentralVelocityChange() throws {
        for outgoing in [SIMD3<Double>(0,1,0),SIMD3(0,0,-1)] {
            let result = try UnpoweredFlyby.solve(request(outgoing: outgoing))
            for (sign,expected) in [(-1.0,SIMD3<Double>(1,0,0)),(1.0,outgoing)] {
                let epoch = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: sign*1e6))
                let actual = try result.trajectory.state(at: epoch)
                #expect(routeNorm(actual.velocity-expected) < 2e-6)
                #expect(routeNorm(actual.position) > 1e5)
            }
            #expect(abs(routeNorm(result.outgoingExcessVelocityMetresPerSecond)-1) < 1e-14)
        }
        let bodyVelocity = SIMD3<Double>(0,2,0)
        let incomingCentral = bodyVelocity+SIMD3<Double>(1,0,0), outgoingCentral = bodyVelocity+SIMD3<Double>(0,1,0)
        #expect(abs(routeNorm(incomingCentral)-sqrt(5)) < 1e-14)
        #expect(routeNorm(outgoingCentral) == 3)
    }
    @Test func safetySpeedConservationAndDegenerateRequests() throws {
        let base = try UnpoweredFlyby.solve(request())
        let boundary = try UnpoweredFlyby.solve(request(altitude: base.periapsisAltitudeMetres))
        #expect(boundary.periapsisRadiusMetres == base.periapsisRadiusMetres)
        #expect(throws: FlybyError.unsafePeriapsis) { try UnpoweredFlyby.solve(request(altitude: 0.4)) }
        #expect(throws: FlybyError.excessSpeedMismatch) { try UnpoweredFlyby.solve(request(outgoing: SIMD3(0,2,0))) }
        #expect(throws: FlybyError.degenerateTurn) { try UnpoweredFlyby.solve(request(outgoing: SIMD3(1,0,0))) }
        #expect(throws: FlybyError.degenerateTurn) { try UnpoweredFlyby.solve(request(outgoing: SIMD3(-1,0,0))) }
        #expect(throws: FlybyError.invalidVelocity) { try request(incoming: .zero) }
        #expect(throws: FlybyError.invalidAltitude) { try request(altitude: -1) }
        let accepted = try UnpoweredFlyby.solve(request(outgoing: SIMD3(0,1+1e-10,0),tolerance: 1e-9))
        #expect(accepted.outgoingExcessVelocityMetresPerSecond == SIMD3(0,1,0))
        #expect(accepted.outgoingVelocityResidualMetresPerSecond > 0)
        #expect(abs(accepted.requestedOutgoingSpeedDifferenceMetresPerSecond-1e-10) < 1e-15)
    }
    @Test func smallTurnAndReversalGeometry() throws {
        let angle = 1e-6
        let small = try UnpoweredFlyby.solve(request(outgoing: SIMD3(cos(angle),sin(angle),0)))
        #expect(abs((small.periapsisRadialDirection*small.periapsisProgradeDirection).sum()) < 1e-14)
        #expect(abs(small.turningAngleRadians-angle) < 1e-14)
        #expect(abs(small.periapsisRadiusMetres/(1/sin(angle/2)-1)-1) < 1e-12)
        let first = try UnpoweredFlyby.solve(request())
        let reversed = try UnpoweredFlyby.solve(request(incoming: SIMD3(0,-1,0),outgoing: SIMD3(-1,0,0)))
        #expect(routeNorm(reversed.periapsisRadialDirection-first.periapsisRadialDirection) < 1e-14)
        #expect(routeNorm(reversed.periapsisProgradeDirection+first.periapsisProgradeDirection) < 1e-14)
        #expect(routeNorm(reversed.orbitalNormal+first.orbitalNormal) < 1e-14)
    }
}
