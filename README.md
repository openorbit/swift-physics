# swift-celmek

A port of the celmek library to swift.

The celmek library implements several astronomical algorithms.

| Algorithm | Status
| --- | --- 
| VSOP87 | Done
| ELP2000 82b | Done
| Goffin2000 | Done

## Status

CelMek is very much in development, and the scope has grown to implement more astronomical algorithms.
It is currently NOT in a 1.0 state, so expect that API breakage happens.

A companion app for interactivelly working with this library is developed here:
[https://github.com/lorrden/CelMek-Viewer]

## Evolution

In addition to the celestial mechanics algorithms,
various functionality related to Aerodynamics have been added.
It is probable that we will rename the package `swift-physics` at some point.

## Modules

Use the `CelMek`, `Math`, or `Aerodynamics` library product for the modules your
application needs. The existing `Physics` product bundles all three. Products select
dependencies; Swift imports remain `import CelMek`, `import Math`, and
`import Aerodynamics`. See [MODULES.md](MODULES.md) for responsibilities and dependency rules.

## Time API

Precise epochs, durations, UTC leap seconds, and calendar readings with explicit
TAI/TT/TDB/UTC scales are documented in [TIME.md](TIME.md).

Coordinate handedness, quaternion and angle conventions, transform direction, and
error handling are documented in [CONVENTIONS.md](CONVENTIONS.md).

## Testing

The tests use Swift Testing and require a Swift 6 or newer toolchain. Run all
test targets with:

```sh
swift test
```

## License

Apache 2.0

Shared Cartesian states use metres and metres per second, with explicit unit conversion and a VSOP87 provider adapter. See [STATES.md](STATES.md).

Time-dependent frame trees, rotating transforms, and their conventions are documented in [FRAMES.md](FRAMES.md).

The IAU coefficient review, updated data, and independent SPICE reference tests are documented in [IAU_ROTATION_AUDIT.md](IAU_ROTATION_AUDIT.md).
