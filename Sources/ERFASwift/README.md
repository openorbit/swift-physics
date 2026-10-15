# ERFASwift

This internal Swift module contains the numerical code and coefficient tables
adapted from ERFA 2.0.1:

- IAU Earth orientation, precession, nutation, CIO locator, and Earth rotation.
- The full geocentric Fairhead–Bretagnon TT–TDB series and analytic derivatives.

The ERFA copyright, BSD terms, and SOFA heritage notice are retained verbatim in
[LICENSE.txt](LICENSE.txt). This is a modified Swift implementation derived from
ERFA, not the official SOFA library. Retain the notice with distributions as
specified in that file. A documentation copy remains at
`Documentation/licenses/ERFA.txt`.

The module depends on the package's Apache-2.0 `Math` module and does not depend
on `CelMek`. Its entry points have Swift `package` access: this target is an
implementation detail, with no standalone library product or public API.

`CelMek` owns epochs, coverage checks, frame providers, and its public errors;
it adapts those to the numerical module without changing the public API.
Regenerate `FairheadTDBSeries.swift` using `Tools/generate_tdb_series.py` with a
caller-supplied ERFA 2.0.1 `dtdb.c` file. No external ERFA binary is required.
