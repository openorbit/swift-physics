# Model coverage and validation

Select the physical approximation and data coverage explicitly, then assess numerical convergence separately.

## Selection and cost

| Provider | Supported approximation / required data | Evaluation cost |
| --- | --- | --- |
| Checked airfoil/control/body tables | Supplied rectangular incidence, Mach, Reynolds and optional control grid; no implied validity between different regimes | Up to 16 body-grid corners; section/angle lookup also searches samples |
| SegmentedWing | Independent strips and sampled local flow; no automatic wake | O(N) |
| FiniteWingAirfoil | Uniform equivalent attached low-speed wing, AR ≥ 4 | Iterative scalar root |
| LiftingLineWing | Straight planar unswept wing, linear sections, uniform Mach ≤ 0.3 flow | O(N³) setup, O(N²) evaluation |
| NonlinearLiftingLine | Supplied nonlinear sections, nonuniform flow, unswept fixed wakes and mutual surface coupling, Mach ≤ 0.3 | O(N²) setup/memory, O(I N²) solve |
| ModifiedNewtonianPressure | Windward calorically perfect gas, declared Mach ≥ 5 and Kn ≤ 0.01 | O(1) per panel |
| MaxwellianFreeMolecularTraction | Single-species Maxwellian, specified reflection/wall temperature, Kn ≥ 10 | O(1) per panel |
| TabulatedRarefiedTraction | Supplied consistent gas/surface data across declared Kn/incidence/speed/wall-temperature range | At most 16 corners per pressure/shear lookup |
| LaminarFlatPlateHeating / TurbulentFlatPlateHeating | Low-Mach isothermal attached zero-pressure-gradient plate; explicit film properties and boundary-layer choice | O(1) per point |
| TabulatedEnthalpySurfaceHeating | Supplied gas/geometry/recovery/catalytic assumptions, explicit laminar/turbulent Stanton data, Kn ≤ 0.01 | At most 8 corners |
| SuttonGravesHeating / mixture / H₂-He | Distinct gas-specific stagnation correlations with their own bounds and wall assumptions | O(1) per point |
| SweptCylinderStagnationHeating | Idealized laminar Earth-air stagnation line, supplied sphere model, bounded effective sweep | One sphere evaluation plus O(1) |
| TabulatedRadiativeHeating | Supplied cold-wall shock-layer radiation data and gas-specific continuum coverage | At most 8 corners |
| DragAreaAerodynamics / buoyancy | Supplied current drag area/coefficient and displaced volume/local pressure gradient | Provider cost plus O(1) |

Tables carry source, version and assumptions. ``AerodynamicScalarTable`` supports
1–8 named/unit-labelled axes, strict ordering, finite values and complete grids,
with last-axis-fastest storage, multilinear interpolation and no extrapolation.
JSON decoding revalidates structure. Physical wrappers additionally check axis
names, units, signs and gas/geometry conventions. A rectangular interpolation grid
must not span an unsupported hole in the physical dataset. Single-node dimensions
are exact constraints, not unconstrained coordinates. Floating-point query values
must remain within the declared endpoints.

``CoveredBodyCoefficients`` binds a provider to reference lengths and explicit
Mach/Reynolds/incidence/sideslip/Kn coverage. ``MachRegimeBodyCoefficients``
dispatches by caller-chosen thresholds to bounded providers. A threshold belongs
to the upper interval; it does not blend, extrapolate or hide discontinuities.
Transonic/high-incidence support therefore comes from appropriate supplied body
data, rather than a low-speed or hypersonic approximation stretched across all
Mach numbers. No validated generic capsule/spaceplane database is bundled.

## Gas and transport

``PerfectGasMixture`` stores explicit mass fractions of ``PerfectGasSpecies`` and
mass-weights cp and R. Each species declares a constant-cp temperature range; the
mixture uses their intersection. It supplies ideal-gas density and sound speed and
can construct film transport properties with **caller-supplied** viscosity and
conductivity. There is no viscosity/conductivity averaging rule or high-temperature
chemistry inferred from composition. Species identifiers must be unique, positive
mass fractions must sum to one, and property queries outside the shared temperature
interval throw. Use a richer provider for dissociation, ionization, variable cp,
non-equilibrium or multicomponent diffusion. A mixture-average R is not a resolved
multispecies free-molecular model.

## Failure and zero-flow behavior

Checked constructors reject invalid geometry/nonfinite coefficients and incomplete
data. Unsupported physical regimes throw `outsideCoverage`; missing state throws
`missingFlowProperty`; incompatible identities/conventions throw `unsupportedQuery`.
The nonlinear wake solver additionally reports `nonConvergence`.

Whole-body/strip/drag-area evaluators return zero aerodynamic load at zero dynamic
pressure without querying undefined incidence. The steady lifting-line solvers,
forced-convection plate models and supplied enthalpy/radiation datasets reject
rest when it lies outside their declared domains. Free-molecular traction retains
thermal molecular impacts at zero bulk speed. Sutton–Graves cold-wall stagnation
models retain their documented zero-speed convention. There is no universal
zero-density-to-vacuum shortcut, fallback model or assumed zero heat flux when a
provider is unavailable. Buoyancy is evaluated independently of airspeed.

## Runnable entry comparison

Run `swift run AerodynamicsDemo`. In addition to wing/control/parachute/balloon
examples it integrates 25 seconds of planar prescribed-attitude entry for a
lifting capsule and a pressure-panel wing. The harness owns RK4, gravity, mass,
attitude prescription and heat-load integration. Aerodynamics only returns loads
and flux. It uses ISA below 86 km, explicitly approximate mean free path,
synthetic capsule coefficients and an idealized pressure-only flat wing. Heating
uses Earth-air Sutton–Graves and a swept-cylinder stagnation-line estimate.

It reports peak deceleration, dynamic pressure, local heat flux and integrated
local heat load, plus final-speed/heat-load changes when halving the step from
0.2 to 0.1 seconds. The starting altitude is 70 km: the small leading-edge radius
would fail continuum coverage at 80 km. Coverage checks remain active throughout
the run. Heat load is J/m² at the chosen local point, not total vehicle energy or
a thermal-protection temperature/survival prediction. Attitude is prescribed;
this harness does not assess trim, structural response or closed-loop guidance.

Step refinement measures integration sensitivity, not physical accuracy. Other
tests separately compare spanwise convergence, classical lifting-line references,
independent molecular quadrature, published stagnation-heating CFD and similarity
heat-transfer values. Supplied-data tests establish interpolation and contracts;
real vehicle data still needs independent calibration/validation.

## Deliberate limits and later fidelity

The initial physics/API workstream includes replaceable providers for the remaining
regimes. It does not include a general CFD/DSMC/radiation solver or universal
validated vehicle datasets. Free-wake rollup, swept/separated-wing fidelity,
chemistry-resolved transport and high-enthalpy/transitional/radiative databases
remain specialist extensions. Atmosphere coverage above ISA, thermal response and
vehicle systems are separate workstreams, not missing state machines to place in
Aerodynamics. See <doc:NonlinearWakes>, <doc:RarefiedFlow> and <doc:EntryHeating>.
