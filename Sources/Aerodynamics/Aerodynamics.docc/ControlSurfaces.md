# Wing and tail control-surface data

``TabulatedControlSurface`` implements the existing section coefficient contract
for caller-supplied geometric control deflection. It uses a set of
``ControlSurfaceSetting`` datasets, each containing a checked incidence/Reynolds/Mach
``TabulatedAirfoil``. Lift, drag and pitching moment interpolate linearly across
deflection after evaluating the neighboring polars. No incidence offset or assumed
control-effectiveness factor is substituted for measured or computed data.

Deflection coordinates are finite radians within ±π, sorted on construction, with
duplicates rejected. All settings share the same Reynolds and Mach coordinate grid;
incidence samples may differ. Intermediate deflections require coverage in both
neighboring polars. Exact deflection settings use only their own dataset, avoiding
irrelevant neighboring coverage restrictions. A single setting represents one
supported deflection, not an assumed constant response at all deflections.

The deflection boundary policy defaults to rejection. Explicit clamping affects
only the deflection coordinate; each underlying airfoil retains its own declared
incidence/Re/Mach policy. Clamping is a caller-selected numerical approximation,
not a calibrated extrapolation. No transition or smoothing beyond piecewise linear
interpolation is inferred.

## Conventions and integration

Source/version, deflection sign and moment-reference descriptions are mandatory.
For example, positive deflection may mean trailing edge down and CM may refer to
quarter chord. All datasets must share reference chord/area conventions, section
axes and moment origin. Metadata records those assumptions but cannot verify the
physical meaning of externally supplied data.

Pass the model directly to ``WingSegment`` or ``SegmentedWing``. Segments already
accept their own `controlDeflectionRadians`, so an external caller can supply
differential aileron, elevator or rudder geometry without adding actuator state here.
Horizontal and vertical tail surfaces use the same section/strip contracts with
the appropriate supplied section axes and positions. Do not apply another incidence
or camber correction for the same deflection when the dataset already includes it.

```swift
let controlled = try TabulatedControlSurface(settings: [
    ControlSurfaceSetting(deflectionRadians: -0.1, airfoil: negativeDeflectionData),
    ControlSurfaceSetting(deflectionRadians: 0, airfoil: neutralData),
    ControlSurfaceSetting(deflectionRadians: 0.1, airfoil: positiveDeflectionData)
], source: "Caller-supplied control polar database", version: "1",
   deflectionConvention: "Positive trailing edge down",
   momentReferenceDescription: "Quarter chord in section axes")
let segment = try WingSegment(positionMetres: SIMD3(-3, 1, 0),
    axes: WingSectionAxes(incidenceRadians: 0), chordMetres: 1, areaSquareMetres: 1,
    controlDeflectionRadians: 0.05, model: controlled)
let loads = try segment.evaluate(bodyFlow: flow).loads
```

The example datasets and local flow are provided by the caller. No measured flap,
rudder or tail database is bundled. Section moments include the dataset CM plus
the existing force lever arm, with no separate centre-of-pressure addition.
Independent strips still do not compute spanwise circulation, wake interaction or
finite-wing control effectiveness automatically.

## Stall and data loading

Nonlinear lift curves, rising drag and post-stall pitching moments can be represented
by supplied polar samples. They are preserved rather than replaced by a lift-slope
formula. This is a quasi-steady interpolation model, not a dynamic-stall, hysteresis,
hinge-moment, flutter or unsteady separation model. Sparse sampling can miss extrema;
coverage does not itself prove accuracy. The low-speed attached-flow limits of
``FiniteWingAirfoil`` must not be extended simply because a control polar has
post-stall samples.

``ControlSurfaceTableLoader`` reads JSON. Codable decoding revalidates the table,
deflection settings and nested airfoil grids. Tests cover four-coordinate affine
interpolation, exact settings, different incidence coverage, rejection/clamping,
duplicate coordinates, invalid decoded data and a supplied post-stall curve.
