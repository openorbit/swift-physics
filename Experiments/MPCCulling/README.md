# MPC culling experiment

Run from the repository root:

```sh
swift run -c release MPCCullingExperiment
swift run -c release MPCCullingExperiment --limit 10000
```

Options: `--data Data/MPC`, `--half-window-days 1`, `--cell-au 0.25`.
`--limit` limits asteroids only (the text reader still scans the file); all comets are included.
Requires the caller-supplied MPCORB.DAT and cometels.json, which are not packaged as library resources.

This standalone experiment streams fixed-width asteroid records and decodes comet JSON.
It prepares Newtonian Sun-centered two-body propagators in ecliptic J2000 axes.
Asteroids start at their packed TT element epoch; comets start at perihelion, retaining
optional element epochs separately. Dates are interpreted as TT for this experiment.
It does not include planetary perturbations or comet nongravitational forces and is
not a precision replacement for MPC/JPL ephemerides. No source GM/time-scale conversion
is asserted. Names, photometry and full source metadata are not yet imported.

The index groups positions at a center epoch into Cartesian grid cells. Bounds are
expanded by perihelion speed times the half-window duration, plus a numerical guard.
This global two-body speed bound covers motion throughout the window, including
hyperbolic and parabolic orbits, without relying on endpoint samples. Per-body radial
bounds intersect swept-radius bounds with whole-orbit perihelion/aphelion bounds.
Queries scan cell radial bounds, optionally reject sphere/AABB nonintersection, then
filter individual radial intervals before evaluating selected states. Spatial indexing
therefore avoids propagation per query but requires one full pass when rebuilt.
Cell scans are linear in occupied cells; this is a baseline for a future BVH.

The benchmark reports load/preparation, build, selection and evaluation separately.
Four sphere queries at both window endpoints and its center are checked against
brute-force exact queries. It also checks stale-index rejection. All geometric
coordinates and query radii use metres in the same heliocentric axes.
Bounds are conservative for the chosen ideal two-body model; the floating-point
padding is pragmatic and not a formally certified error bound. Wide windows or
small-perihelion comets loosen bounds. Bodies are not truncated or silently skipped
on loader/evaluator failures. Rebuild after jumps outside the window.

Source data credit: Minor Planet Center. Preserve the complete MPCORB header,
per-record references and orbit-computer attribution when redistributing source or
derived datasets. Update information is in Data/MPC/MPCORB.DAT's header. This
experiment does not write or redistribute a compiled dataset.

## Preliminary full-catalog measurement

2026-10-13, local arm64 release executable, after compilation completed.
These are indicative single-run measurements, not repeated benchmark medians.
The full transcript is in `results/full-catalog.txt` (timings recorded to microseconds).

- 1,574,501 records; asteroid load/preparation 6.014 s, comet load/preparation 0.008 s.
- One-day half-window, 0.25 AU cells: 25,832 occupied cells; build 2.977 s.
- Full evaluation oracle: 2.403–2.643 s per sphere query (includes distance checks/output collection).
- Belt probe at (2.7, 0, 0) AU with radius 0.1 AU: 294,308 radial candidates,
  5,276 spatial candidates, 535 exact hits at center epoch. Spatial candidate
  evaluation took 0.021 s at center epoch and 0.009–0.021 s across the three epochs.
- All 12 full-catalog queries matched brute force; synthetic cases and stale-index rejection passed.

This baseline deliberately uses prepared existing propagators per record, rather
than the eventual compact structure-of-arrays catalog. It measures the existing
solver and prioritizes verifying culling correctness before memory optimization.
