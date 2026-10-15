// SPDX-License-Identifier: Apache-2.0
import Testing
import CelMek

@Suite
struct CalendarAndStateAPIContractTests {
  private func requireSendable<T: Sendable>(_ value: T) {}

  @Test
  func renamedCalendarCasesPreserveRawValuesAndDescriptions() {
    #expect(Month.january.rawValue == 1)
    #expect(Month.december.rawValue == 12)
    #expect(Month.march.description == "March")
    #expect(Weekday.sunday.rawValue == 0)
    #expect(Weekday.monday.rawValue == 1)
    #expect(JewishMonth.tishri.rawValue == 1)
    #expect(MoslemMonth.muharram.rawValue == 1)
    requireSendable(Month.january)
    requireSendable(Weekday.monday)
  }

  @Test
  func cartesianStateCanCrossIsolationBoundaries() throws {
    let state = try CartesianState(epoch: .j2000, frame: CelestialFrames.icrf,
      position: .zero, velocity: .zero)
    requireSendable(state)
    #expect(state.frame === CelestialFrames.icrf)
  }
}
