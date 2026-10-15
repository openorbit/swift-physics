# Pinned SPICE black-box references

Generated on 2026-10-13 with **CSPICE_N0067**, binary SpiceyPy **8.3.0**, NumPy
2.5.4, on macOS arm64. The input states, attitude samples, polynomial coefficients,
clock correlations and abbreviated leap-second history are synthetic test data.
They do not represent a mission or a physical accuracy claim.

`Scripts/generate_spice_oracle.py` calls public writer/query APIs against the binary
oracle. No NASA/NAIF toolkit implementation source was inspected, translated or
used to derive the native Swift algorithms. Writer output supplies independent
small standard DAF files; query output supplies reference values for those files.
The script uses `m2q` to create attitude input quaternions independently of Swift.

| Asset | Coverage and purpose |
| --- | --- |
| states.bsp | SPK 2/3, targets 399/499 relative to centre 3, plus 3 relative to 0; -20 to +20 TDB seconds |
| attitudes.bc | CK 1/2/3 with rates, CK 3 without rates; encoded ticks 0 to 100, multiple intervals/gaps |
| rotation.bpc | Binary PCK 2, frame class 3000 in J2000; -20 to +20 TDB seconds |
| frames.tf | Oracle-only class-2 frame registration for querying binary PCK via pxform |
| time.tls | Synthetic post-1971 leap history, including the 1972 positive leap labels |
| clock.tsc | Two partitions and two TDB correlation slopes; ticks 0 to 2000 |

`reference.json` pins the toolkit/wrapper identity, inputs, expected outputs and
SHA-256 for all six kernel files. Its own SHA-256 is
`666a60aa902cd7bfffb479f8cf9ecc136575b144e54b3c29caf74fbdc17850be`.
The Swift test pins that hash and verifies asset hashes when CryptoKit is available;
all numeric checks run regardless. The fixture folder is copied by the existing
test-resource declaration, and needs no network or oracle installation to test.

Reference units/conventions:

- SPK: geometric J2000 Cartesian states, km and km/TDB second, aberration `NONE`;
  primitive centre composition and SI/atomic-epoch adapter are both compared.
- CK: reference-to-instrument row-major matrices, encoded requested/selected times,
  stored physical angular rates in reference axes (rad/model second).
- Binary PCK: J2000-to-body row-major orientation matrices from `pxform`.
- UTC: UTC labels converted by `str2et` then `unitim(ET, TT)` into TT seconds past
  J2000. Absolute-double precision at these distant epochs sets the 0.3 microsecond
  comparison tolerance; Swift UTC round trips remain exact.
- SCLK: fractional ticks converted by `sct2e` into TDB seconds, including correlation
  boundaries. This is a synthetic kernel-specific convention, not a universal clock.

Typical numeric comparison tolerances are 2e-12 km, 2e-13 km/s, and 2e-13 matrix/rate
components. Boundary and timestamp selection is exact. Synthetic public-format tests
separately cover both byte orders, malformed inputs, priority and unload/reload.
This oracle set contains the platform's little-endian writer output. Large real
mission/planetary assets remain optional; these fixtures do not validate deferred
SPK/CK/PCK types, all toolkit frame classes, or aberration corrections.

## Optional regeneration

Install binary wheels only into an external environment (not a library dependency),
then run:

```sh
PYTHONPATH=/path/to/optional/oracle python3 Scripts/generate_spice_oracle.py
```

The script replaces its generated files. Review numerical changes and update the
pinned reference hash deliberately. The required binary oracle and Python wrapper
are not distributed here. The bundled kernel data and reference outputs are generated
test artifacts from project-authored synthetic inputs, rather than redistributed JPL
mission datasets. Their file structure follows public SPICE documentation.

Public format/behavior references (not implementation source):

- [DAF Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html)
- [SPK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html)
- [CK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/ck.html), including type-3 independent rate interpolation
- [PCK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/pck.html)
- [Time Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/time.html)
- [SCLK Required Reading](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/sclk.html)

Public API contracts used by the optional generator: spkopn/spkw02/spkw03/spkcls,
spkez, ckopn/ckw01/ckw02/ckw03/ckcls, ckgp/ckgpav, rotate/m2q,
pckopn/pckw02/pckcls, furnsh/pxform, str2et/unitim and sct2e.
