import Testing
import Atmospheres
@Suite struct AtmosphereTest {
    @Test func legacyFacadeUsesSI() {
        let atm = EarthAtmosphere()
        #expect(abs(atm.dynamicViscosity(at: 0)-1.7894e-5) < 1e-9)
        #expect(abs(atm.pressure(at: 0)-101325) < 1e-8)
        #expect(abs(atm.density(at: 0)-1.225) < 1e-6)
        #expect(atm.temperature(at: 0) == 288.15)
        #expect(atm.dynamicViscosity(at: -10000).isNaN)
    }
}
