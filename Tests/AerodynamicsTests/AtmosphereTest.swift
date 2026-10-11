//
//  AtmosphereTest.swift
//  
//
//  Created by Mattias Holm on 2023-11-05.
//

import Testing
import Aerodynamics
@Suite
struct AtmosphereTest {

    @Test
    func testExample() throws {
      let atm = EarthAtmosphere()
      #expect(1.789 == atm.dynamicViscosity(at: 0))
    }
}
