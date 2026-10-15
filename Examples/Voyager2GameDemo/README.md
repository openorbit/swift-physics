# Voyager 2 inspired game planning example

Run from the repository root:

```sh
swift run Voyager2GameDemo > /tmp/voyager2-report.json
```

In a restricted development sandbox use the same module-cache configuration as the
other demos:

```sh
env CLANG_MODULE_CACHE_PATH=/tmp/celmek-clang-cache \
    SWIFTPM_MODULECACHE_OVERRIDE=/tmp/celmek-swift-cache \
    swift run --disable-sandbox Voyager2GameDemo > /tmp/voyager2-report.json
```

No DE/SPICE assets, network service, or CSPICE installation are required. The example
uses full VSOP87E and the injected Fairhead TDB model, plus an explicit built-in
leap-second table for historical UTC input/output. It subtracts the Sun's barycentric
position and velocity to form Sun-relative states in the same ecliptic J2000 axes.
No metre/AU or barycentric/heliocentric relabelling is implicit.

## What it plans

The [NASA mission fact sheet](https://science.nasa.gov/mission/voyager/fact-sheet/)
records launch on 1977-08-20 and Jupiter/Saturn/Uranus/Neptune encounters on
1979-07-09, 1981-08-25, 1986-01-24 and 1989-08-25. The
[NASA Voyager 2 mission page](https://science.nasa.gov/mission/voyager/voyager-2/)
supplies launch time 14:29:44 UT. The example uses that launch timestamp and **noon
UTC** on the encounter reference days, rather than exact historical closest-approach
instants.

It searches all 81 combinations of -30, 0 and +30 days around the four encounter
dates. Thirty precomputed zero-revolution, prograde Lambert connectors join Earth,
Jupiter, Saturn, Uranus and Neptune in the Sun's two-body field. Candidate ranking
minimizes illustrative parking-orbit injection delta-v plus unresolved encounter
velocity mismatches. The dates are constraints with a small search window, so matching
the historical flight time is not an independent prediction of the mission schedule.

A 200 km circular Earth parking orbit is converted into a local hyperbolic departure
with the Lambert leg's outgoing excess vector. The example samples the actual
Earth-centred parking/burn/escape propagator and translates those states into the
heliocentric plotting frame. Its departure calculation includes injection delta-v
and C3; it does not simulate the launch vehicle, ascent, launch-site plane or staging.

At Jupiter, Saturn and Uranus, incoming/outgoing Lambert velocities are measured
relative to the moving body. A local unpowered flyby preserves incoming excess speed,
turns toward the next connector, and respects the configured spherical periapsis
limit. If the requested turn exceeds that limit, the model clamps it. The report
shows requested/achieved turn, signed excess-speed mismatch, inferred periapsis altitude
and the remaining outgoing vector mismatch. That mismatch is **not an implemented
maneuver or a periapsis burn cost**, and is not silently assigned to spacecraft fuel.

Neptune has no following Lambert endpoint. A final local flyby targets a chosen
out-of-ecliptic departure at approximately the demo's 5,000 km minimum altitude.
This demonstrates a fourth flyby without inventing a historical interstellar target.

## Expected result and practical limits

On the development run the chosen dates are the reference encounter days: total
flight time is about **4,387.90 days (12.01 Julian years)**. The illustrative Earth
injection is about **7.19 km/s**, with **C3 103.01 km²/s²**. Jupiter/Saturn/Uranus
connector velocity residuals are approximately **102 / 134 / 13 m/s**. Propagating
each solar Lambert leg independently reaches its requested endpoint within about
0.1 m on that run; that is solver agreement, not a bound on mission-model error.
These are diagnostics, not historical reference fixtures or cross-platform exact values.

The output declares `isContinuousMissionTrajectory: false`. Local flyby, parking/escape
and solar coast samples are separate views. Finite sphere-of-influence patch times,
positions and velocities are not solved, and transitions may have discontinuities.
The solar coasts reach system barycentres; local flybys approximate those as planetary
centres. Satellite motion, atmospheric/ring constraints, planetary perturbations during
coasts, navigation corrections and historical spacecraft mass/propellant are omitted.
Rounded Newtonian GM and spherical equatorial radii are explicit illustrative constants.
The arbitrary parking plane, terminal turn and safety altitudes are not recovered mission
parameters. Do not interpret this example as a continuous Voyager reconstruction.

The library's finite encounter-sequence and numerical verification APIs can be used
for the next refinement: optimize B-plane geometry and dates, connect finite local arcs
with explicit residual/correction constraints, and propagate the resulting itinerary
under chosen perturbing forces. The current example intentionally exposes the residuals
that refinement must resolve.

## Using the samples in a game

The JSON includes legs, flyby diagnostics, flight-time comparison, clock playback and sampled SI
position/velocity states. Solar coasts contain 49 uniformly spaced samples each;
parking/departure and local encounters have their own labelled sample groups. Separate
`phase` values before drawing: joining all rows as one polyline would connect unrelated
local views and discontinuous patches. Encounter dates are not rest stops or rendezvous.

For playback, retain each leg's `TwoBodyPropagator` and query `state(at:)` at the atomic
simulation epoch advanced by your time-warp scale. Use the same epoch for body ephemeris
queries. For prediction query a future covered epoch without mutating the current state;
choose the appropriate coast segment and display transition residuals explicitly. The
returned solar propagator itself is not coverage-limited, so the example/application
must enforce each leg's departure/arrival bounds. Avoid linearly interpolating distant
render samples as a replacement for orbital propagation.

Rendering can subtract the camera's **Double** position before converting to `Float`,
with an explicit scene scale. The sample coordinates are Sun-relative metres, not
camera-relative graphics units; a separate renderer adapter remains roadmap work.
The `playback` rows exercise the reusable `SimulationClock`: one real second advances
one 365-day simulation year, followed by pause, reverse and a seek to the Jupiter patch.
Each row queries the spacecraft and all five prescribed planets at the same epoch.
A 30-day future query leaves the clock untouched; predictions beyond Neptune are
omitted instead of extrapolating the itinerary. Segment selection uses the following
coast at an exact patch epoch, so velocity jumps remain explicit. Orbital display
samples use `TrajectorySampler` rather than hidden interpolation.

See [time guide](../../Sources/CelMek/CelMek.docc/TimeGuide.md) for clock semantics, bounds and floating-point limits.
Renderer adapters, persistence and the remaining release-readiness work are separate
Milestone 9 items.

Playback also exports `spacecraftCentredPlanetPositionsRenderUnits`: five planet
positions relative to the spacecraft, converted to `Float` after subtracting in
`Double`, at 1e9 metres per render unit. Axes remain right-handed heliocentric
axes; the spacecraft renders at zero. See [rendering contracts](../../Sources/CelMek/CelMek.docc/GAME_INTEGRATION.md).
