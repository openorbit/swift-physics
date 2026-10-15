# Atmospheres

Swift target with no external dependencies. Import `Atmospheres` via
swift-celmek's `Atmospheres` or `Physics` product.

```swift
import Atmospheres
let earth = ISAAtmosphere()
let seaLevel = try earth.properties(at: 0) // Pa, K, kg/m³, m/s, Pa s
let jupiter = BodyAtmosphere.model(for: .jupiter)
let belowClouds = try jupiter.properties(at: -50_000) // below the 1-bar datum
```

Earth uses dry-air ISA/US Standard Atmosphere 1976 lower layers, through 84.852 km
geopotential height, converting geometric input altitude internally. All profiles
throw outside their explicit coverage; coverage is not an atmospheric edge or a
promise of vacuum above it. Custom `HydrostaticAtmosphere` profiles accept increasing
temperature levels and a reference pressure at one of those levels.

Built-ins cover Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Titan, Triton,
Io and Callisto. Non-Earth solid-body profiles are **nominal isothermal approximations**;
they reproduce a reference state, not a measured vertical temperature profile.
Gas giants fit the published 1-bar and 0.1-bar temperatures with a hydrostatic
linear-temperature layer, a dry adiabat below 1 bar (limited to 100 bar), and an
isothermal upper extension. Neither cloud microphysics nor deep interiors are modeled.

Io is highly variable; its pressure is an order-of-magnitude nominal value, with a
chosen 150 K gas temperature. Callisto represents the historically observed CO2
component as an exponential exosphere proxy. Its pressure is equivalent kinetic
pressure; continuum sound speed and viscosity are absent. Non-Earth viscosity is
unavailable rather than guessed. These fixed profiles omit latitude, seasons, solar
activity, winds, plumes and atmospheric escape. Additional moons can be represented
with custom profiles; no poorly constrained defaults are silently supplied.

Every model exposes its source URL, assumptions, altitude datum and coverage.
Reference data: [NASA planetary factsheets](https://nssdc.gsfc.nasa.gov/planetary/factsheet/),
[US Standard Atmosphere 1976](https://ntrs.nasa.gov/citations/19770009539),
[Titan nominal reference](https://www.nasa.gov/wp-content/uploads/2015/04/niac_sunspiral_superballbot_phaseii_finalreport_tagged.pdf),
[Triton reference](https://ntrs.nasa.gov/citations/19950012719),
[Io review](https://doi.org/10.1006/icar.1996.0186),
[Callisto CO2 observation](https://doi.org/10.1126/science.283.5403.820).

For CelMek, `BodyIdentifier.earth.atmosphericModel` returns the built-in profile.
`SphericalBodyAtmosphere` adapts a profile to `AtmosphericDrag`. Supply a compatible
body-centred inertial frame, the radius of the profile's datum (1 bar for giants),
and optionally angular velocity in those axes. It computes density and co-rotation;
it neither transforms frames nor obtains ephemerides. Exosphere density should be
paired with a suitable rarefied-gas force model rather than assumed continuum drag.

Migration: `EarthAtmosphere` remains accessible through `import Aerodynamics`.
Its formerly incorrect scalar results now use SI: Kelvin rather than Celsius,
pressure in Pa, density in kg/m³ and viscosity in Pa s (about 1.789e-5 at sea level).
Legacy scalar queries return NaN outside coverage; prefer the throwing profile API.

Run `swift test --filter ProfileTests` from the repository root.
