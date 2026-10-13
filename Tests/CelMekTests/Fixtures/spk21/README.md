# Type 21 binary oracle fixtures

Generated with binary SpiceyPy 8.3.0 and CSPICE_N0067. No toolkit implementation
source was inspected or translated. Scripts/generate_spk21_oracle.py calls the
published SPKW21 writer contract via its Fortran entry point (no C wrapper), then
queries geometric J2000 states through spkez. Synthetic difference tables have
independent random coefficients, per-axis orders, interleaved reference states,
and variable steps. Dimension 15 has 201 records to exercise both directory blocks;
dimension 25 has two records to exercise maximum supported integration orders.
Reference epochs are interior to their records; queries lie on both sides.

The real fixture is NASA/JPL Horizons' Ceres export, obtained through
Tools/SmallBodies/extract.py kernels with requested 2000-01-01 through 2000-01-03
TDB. Its original kernel comments are preserved. It is target 20000001 relative to
Sun (10), frame J2000 (1), type 21. It is test data, not a production coverage choice.
The real reference states were evaluated against the binary oracle, including both
sides of every internal boundary. References pin file SHA-256, native centre/axes,
record selections, positions in km and velocities in km/TDB second. Synthetic
position tolerance is 6e-8 km; real tolerance is 1e-7 km; velocity tolerances are
2e-13 and 1e-12 km/s respectively. These are implementation comparison tolerances.

Regenerate synthetic fixtures with an external optional binary installation:

```sh
PYTHONPATH=/path/to/oracle python3 Scripts/generate_spk21_oracle.py
```

Add or replace the real references with:

```sh
PYTHONPATH=/path/to/oracle python3 Scripts/generate_spk21_oracle.py --horizons /path/to/20000001.bsp
```

The public-format synthetic Swift fixtures separately cover big-endian input and
malformed directories/records. There is no toolkit runtime dependency in CelMek.
