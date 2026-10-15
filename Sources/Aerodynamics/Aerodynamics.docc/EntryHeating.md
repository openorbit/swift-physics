# Entry surface and leading-edge heating

Use explicit recovery/chemistry data for local surface heating and a bounded swept-cylinder approximation at the leading-edge stagnation line.

## Enthalpy-based surface data

``EnthalpySurfaceHeatingModel`` is separate from low-speed temperature-difference
convection. ``EnthalpySurfaceHeatingQuery`` requires local edge flow, streamline
distance, gas identity and wall/adiabatic-wall enthalpies sharing the dataset's
reference. Adiabatic-wall enthalpy includes recovery; it is not automatically the
edge static or total enthalpy. Real-gas chemistry and catalytic-wall assumptions
must be supplied with the data rather than replaced by a constant cp times temperature.

``TabulatedEnthalpySurfaceHeating`` evaluates
`q = rho_edge V_edge St (h_adiabaticWall - h_wall)` in W/m². The sign preserves
convective cooling. Its dimensionless Stanton table uses, in order, axes named
`mach`, `log10Re` and `wallEnthalpyRatio`, each with units `1`. Reynolds uses edge
viscosity/density and the supplied streamline distance. The last coordinate is
wall enthalpy divided by adiabatic-wall enthalpy, so the chosen enthalpy reference
must remain consistent. Positive adiabatic-wall enthalpy is required.

A dataset explicitly declares laminar or turbulent boundary layer; selection and
transition detection remain caller responsibilities. Each table needs gas identity,
source, version and assumptions covering geometry, pressure gradient, chemistry
and catalytic state. Its rectangle must be wholly supported by the data, including
intermediate interpolated conditions. There is no extrapolation, automatic
laminar/turbulent switch or built-in calibrated high-enthalpy surface dataset.
``EnthalpySurfaceHeatingCoverage`` also bounds edge temperature, edge pressure,
absolute adiabatic-wall enthalpy and streamline distance. These checks prevent
a Mach/Reynolds match from silently admitting an unsupported thermochemical state.
Kn must not exceed the declared continuum limit (at most 0.01). Missing edge temperature/pressure, local sound
speed, viscosity or mean free path throws. Rest/zero Reynolds is outside coverage.

## Swept leading edge

``LeadingEdgeConvectiveHeatingModel`` accepts a ``LeadingEdgeHeatingQuery`` so
higher-fidelity geometry/chemistry models can replace the initial provider.
``SweptCylinderStagnationHeating`` implements the idealized laminar **Earth-air**
geometry approximation `q_line = q_sphere / sqrt(2) * sqrt(Rsphere/RleadingEdge)
* cos(effectiveSweep)^1.1`. Its reference sphere evaluator and radius are explicit;
the query radius denotes the actual leading-edge radius. The sphere provider
retains its own gas, wall and flow restrictions. Effective sweep is relative to
incident flow. Declared speed/density/radius/Mach/Kn limits are checked using the
actual leading-edge radius; implementation sweep is restricted to at most ±60°.
This restriction is not an experimental accuracy guarantee. It does not calculate
sweep from vehicle attitude or discover the stagnation line. Other gases require
separately justified models or data.

The source geometry approximation is [NASA CR-181922, equation 1](https://ntrs.nasa.gov/api/citations/19900007739/downloads/19900007739.pdf).
Its implementation is verified against unswept and 60° reference calculations.
It predicts a stagnation **line**, with no circumferential distribution, shock
interaction, transition or acreage distribution. Do not use a small leading-edge
radius with a sphere model's Kn check alone: the smaller radius can be rarefied
while the reference sphere is still in continuum flow.

## Incoming radiation

``TabulatedRadiativeHeating`` implements ``RadiativeHeatingModel`` using a supplied
cold-wall shock-layer dataset. Its axes are `speed` (`m/s`), `density` (`kg/m^3`)
and spherical nose `radius` (`m`), in that order, with values in `W/m^2`. Gas identity,
provenance, version, assumptions, freestream-temperature range and declared
continuum coverage are mandatory.
It accepts cold-wall queries without additional stagnation-gas state. It does not
infer radiative equilibrium, chemistry, optical thickness or wall emissivity, and
no universal planetary radiation correlation is bundled. Missing radiation data
still means unknown, not zero. Convective and radiative contributions can be
combined externally only when they describe the same state and position.

All these evaluations are stateless. Temperature evolution, radiative emission,
conduction, ablation and changes in geometry/mass belong to external thermal and
systems packages. The supplied-data tests use analytic synthetic fixtures to
verify interpolation, units, signs, identity and coverage; they do not validate a
real vehicle database. See <doc:ModelCoverage> for the runnable entry harness and
its distinction between integration error and physical-model error.

The enthalpy/Stanton formulation is described in
[NPTEL compressible boundary-layer heat transfer](https://archive.nptel.ac.in/content/storage2/courses/101103003/module5/lec26/4.html).
