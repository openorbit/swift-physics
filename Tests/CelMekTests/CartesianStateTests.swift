// SPDX-License-Identifier: Apache-2.0
import Testing
import CelMek
import Math

@Suite
struct CartesianStateTests {
  let frame = Frame(name: "Test", pseudoInertial: true)

  @Test
  func normalizesUnitsAndPreservesMissingAcceleration() throws {
    let vector = SIMD3<Double>(1, -2, 3)
    let si = try CartesianState(epoch: .j2000, frame: frame, position: vector, velocity: vector)
    #expect(si.position == vector && si.velocity == vector)
    #expect(si.acceleration == nil)
    let km = try CartesianState(epoch: .j2000, frame: frame, position: vector, velocity: vector,
                                acceleration: .zero, units: .kilometresAndSeconds)
    #expect(km.position == vector * 1000 && km.velocity == vector * 1000)
    #expect(km.acceleration == .zero)
    let au = try CartesianState(epoch: .j2000, frame: frame, position: vector, velocity: vector,
                                acceleration: vector, units: .astronomicalUnitsAndDays)
    for axis in 0..<3 {
      #expect(au.position[axis] == vector[axis] * 149_597_870_700)
      #expect(abs(au.velocity[axis] - vector[axis] * 149_597_870_700 / 86_400) < 1e-9)
      #expect(abs(au.acceleration![axis] - vector[axis] * 149_597_870_700 / (86_400 * 86_400)) < 1e-13)
    }
  }

  @Test
  func equalityIncludesEpochAndFrameIdentity() throws {
    let state = try CartesianState(epoch: .j2000, frame: frame, position: .zero, velocity: .zero)
    #expect(state == state)
    let otherFrame = Frame(name: frame.name, pseudoInertial: true)
    #expect(state != (try CartesianState(epoch: .j2000, frame: otherFrame, position: .zero, velocity: .zero)))
    let later = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: 1))
    #expect(state != (try CartesianState(epoch: later, frame: frame, position: .zero, velocity: .zero)))
  }

  @Test(arguments: [Double.nan, .infinity, -.infinity], [0, 1, 2])
  func rejectsNonfiniteComponents(value: Double, axis: Int) {
    var vector = SIMD3<Double>.zero
    vector[axis] = value
    #expect(throws: CartesianStateError.nonFinitePosition) {
      try CartesianState(epoch: .j2000, frame: frame, position: vector, velocity: .zero)
    }
    #expect(throws: CartesianStateError.nonFiniteVelocity) {
      try CartesianState(epoch: .j2000, frame: frame, position: .zero, velocity: vector)
    }
    #expect(throws: CartesianStateError.nonFiniteAcceleration) {
      try CartesianState(epoch: .j2000, frame: frame, position: .zero, velocity: .zero, acceleration: vector)
    }
  }

  @Test
  func rejectsConversionOverflow() {
    #expect(throws: CartesianStateError.nonFinitePosition) {
      try CartesianState(epoch: .j2000, frame: frame,
                         position: SIMD3(repeating: Double.greatestFiniteMagnitude), velocity: .zero,
                         units: .astronomicalUnitsAndDays)
    }
  }

  @Test
  func publicProviderReturnsSIReferenceState() throws {
    let epoch = try Epoch(julianDate: 2451545, scale: .tdb)
    let provider: any StateProvider = earth
    let state = try provider.state(at: epoch)
    #expect(state.epoch == epoch)
    #expect(state.frame === provider.frame && state.frame === CelestialFrames.vsop87E)
    #expect(state.acceleration == nil)
    let position = SIMD3<Double>(-0.1842769826, 0.9644534530, 0.0002022407) * 149_597_870_700
    let velocity = SIMD3<Double>(-0.0172022466, -0.0031661954, 0.0000000125) * (149_597_870_700 / 86_400)
    for axis in 0..<3 {
      #expect(abs(state.position[axis] - position[axis]) < 15)
      #expect(abs(state.velocity[axis] - velocity[axis]) < 0.0002)
    }
  }

  @Test(arguments: [0.0, -3_155_760_000, 3_155_760_000])
  func providerVelocitiesFollowElapsedSeconds(offset: Double) throws {
    let epoch = try Epoch.j2000.advanced(by: CelMek.Duration(seconds: offset))
    let before = try epoch.advanced(by: CelMek.Duration(seconds: -100))
    let after = try epoch.advanced(by: CelMek.Duration(seconds: 100))
    let before2 = try epoch.advanced(by: CelMek.Duration(seconds: -200))
    let after2 = try epoch.advanced(by: CelMek.Duration(seconds: 200))
    // A five-point difference reduces orbital-curvature error while a longer
    // interval limits cancellation in distant-epoch series evaluation.
    for body in [sun, mercury, venus, earth, mars, jupiter, saturn, uranus, neptune] {
      let state = try body.state(at: epoch)
      let difference = try (8 * (body.state(at: after).position - body.state(at: before).position)
        - (body.state(at: after2).position - body.state(at: before2).position)) / 1200
      for axis in 0..<3 {
        #expect(abs(state.velocity[axis] - difference[axis]) < 0.001)
      }
    }
  }
}
