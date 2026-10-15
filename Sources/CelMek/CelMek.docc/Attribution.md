# Licensing and attribution

CelMek is copyright Mattias Holm and distributed under Apache License 2.0.

The repository's `LICENSE` is the authoritative package license. Third-party data
and reference algorithms carry their own provenance and attribution in the relevant
sources and fixture notes. ERFA-derived numerical routines and coefficient tables
are isolated in the internal `ERFASwift` module under their BSD terms. The full
copyright, license, and SOFA heritage notice is colocated at
`Sources/ERFASwift/LICENSE.txt`, with a documentation copy at
`Documentation/licenses/ERFA.txt`. CelMek depends on this module for Earth
orientation and the full Fairhead TDB model; their public APIs remain in CelMek.

Aerodynamics uses Yams, and the documentation build uses the Swift-DocC plugin.
WMM2025 coefficients are distributed by NOAA; IGRF-14 coefficients are attributed
to IAGA in the MagneticFields guide. Optional JPL/SPICE datasets are supplied
separately and documented in <doc:SPICE_KERNELS> and <doc:DE441>.
