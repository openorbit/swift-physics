import Foundation
import Testing
@testable import CelMek

@Suite struct ConicAnomalyTests {
    @Test func anomalyRoundTripsAcrossUnity() throws {
        for e in [0.0, 0.3, 0.9999999999999999, 1, 1.0000000000000002, 1.5, 10] {
            for nu in [-0.9, -1e-6, 0, 1e-6, 0.9] {
                let anomaly = try ConicAnomaly.from(trueAnomaly: nu, eccentricity: e)
                #expect(abs(try anomaly.trueAnomaly(eccentricity: e) - nu) < 2e-14)
                let m = try anomaly.meanAnomaly(eccentricity: e)
                // Tight absolute anomaly tolerance is necessary because E/H is
                // tiny near e=1; the default tolerance does not bound nu error.
                let options = try KeplerSolverOptions(absoluteTolerance: 1e-25, relativeTolerance: 1e-14,
                                                       maximumIterations: 256)
                let solved = try ConicAnomaly.from(meanAnomaly: m, eccentricity: e, options: options)
                #expect(abs(try solved.trueAnomaly(eccentricity: e) - nu) < 2e-12)
            }
        }
    }

    @Test func knownAnomaliesAndRevolutions() throws {
        let elliptic = try ConicAnomaly.from(trueAnomaly: .pi/2, eccentricity: 0.5)
        guard case .eccentric(let eccentric) = elliptic else { Issue.record("Wrong conic"); return }
        #expect(abs(eccentric - .pi/3) < 1e-14)
        let hyperbolic = try ConicAnomaly.from(trueAnomaly: .pi/2, eccentricity: 2)
        guard case .hyperbolic(let h) = hyperbolic else { Issue.record("Wrong conic"); return }
        #expect(abs(h - acosh(2)) < 1e-14)
        #expect(abs(try hyperbolic.meanAnomaly(eccentricity: 2) - (2*sqrt(3)-acosh(2))) < 1e-14)
        let parabolic = try ConicAnomaly.from(trueAnomaly: .pi/2, eccentricity: 1)
        guard case .barker(let d) = parabolic else { Issue.record("Wrong conic"); return }
        #expect(abs(d - 1) < 1e-14)
        #expect(abs(try parabolic.meanAnomaly(eccentricity: 1) - 4.0/3) < 1e-14)
        for b in [-1e100, -4.0/3, -1e-100, 0, 1e-100, 4.0/3, 1e100, Double.greatestFiniteMagnitude] {
            let solution = try ConicAnomaly.from(meanAnomaly: b, eccentricity: 1)
            if b == .greatestFiniteMagnitude {
                // Cubing the rounded root can overflow at this boundary;
                // inversion still gives a finite parameter.
                if case .barker(let x) = solution { #expect(x.isFinite) }
            } else {
                let reconstructed = try solution.meanAnomaly(eccentricity: 1)
                #expect(abs(reconstructed-b) <= abs(b)*1e-13)
            }
        }
        let mean = 8 * Double.pi + 0.1
        let phase = try ConicAnomaly.from(meanAnomaly: mean, eccentricity: 0.9)
        #expect(abs(try phase.meanAnomaly(eccentricity: 0.9) - mean) < 1e-12)
        #expect(abs(try phase.trueAnomaly(eccentricity: 0.9)) < .pi)
    }

    @Test func nearUnityMeanAnomaliesMatchIndependentDecimalReferences() throws {
        // Inverse of the 80-digit fixtures from generate_kepler_fixtures.py.
        let elliptic = ConicAnomaly.eccentric(1.8169983969771785e-6)
        let hyperbolic = ConicAnomaly.hyperbolic(1.8168762011230227e-6)
        #expect(abs(try elliptic.meanAnomaly(eccentricity: 1.0.nextDown) - 1e-18) < 1e-32)
        #expect(abs(try hyperbolic.meanAnomaly(eccentricity: 1.0.nextUp) - 1e-18) < 1e-32)
    }

    @Test func nearParabolicApoapsisRemainsFinite() throws {
        let tree = FrameTree()
        let mu = try GravitationalParameter(cubicMetresPerSecondSquared: 1, convention: .newtonian,
            provenance: PhysicalDataProvenance(source: "Test", version: "1", reference: "analytic"))
        let e = 1.0.nextDown
        let elements = try OsculatingOrbitalElements(epoch: .j2000, frame: tree.root, centralBody: .earth,
            gravitationalParameter: mu, semiLatusRectumMetres: 1, eccentricity: e, inclination: 0,
            longitudeOfAscendingNode: 0, argumentOfPeriapsis: 0, trueAnomaly: .pi)
        let state = try elements.cartesianState()
        #expect(abs(state.position.x + 1/(1-e)) < 4)
        let almostPi = Double.pi.nextDown
        let denominator = try OsculatingOrbitalElements.radialDenominator(eccentricity: 1, trueAnomaly: almostPi)
        #expect(denominator > 0)
        #expect(1 + cos(almostPi) == 0) // The naive formula loses this finite branch.
    }

    @Test func rejectsWrongKindsAndAsymptotes() throws {
        #expect(throws: OrbitalElementsError.invalidElements) {
            try ConicAnomaly.barker(1).meanAnomaly(eccentricity: 0.5)
        }
        #expect(throws: OrbitalElementsError.invalidElements) {
            try ConicAnomaly.eccentric(.nan).trueAnomaly(eccentricity: 0.5)
        }
        for (e, nu) in [(1.0, Double.pi), (2.0, acos(-0.5)), (2.0, Double.pi)] {
            #expect(throws: OrbitalElementsError.anomalyOutsideConic) {
                try ConicAnomaly.from(trueAnomaly: nu, eccentricity: e)
            }
        }
        #expect(throws: OrbitalElementsError.anomalyOutsideConic) {
            try ConicAnomaly.barker(1e100).trueAnomaly(eccentricity: 1)
        }
        #expect(throws: OrbitalElementsError.numericalRangeExceeded) {
            try ConicAnomaly.hyperbolic(1000).meanAnomaly(eccentricity: 2)
        }
    }
}
