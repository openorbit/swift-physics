WMM2025_TEST_VALUES.txt is NOAA's published WMM2025 test table from the user's
WMM2025 Linux download. Height is kilometres above the WGS84 ellipsoid.
Source: https://www.ncei.noaa.gov/products/world-magnetic-model

IGRF14_REFERENCE.txt contains 48 vectors from NOAA's igrf14syn Fortran reference
routine (https://www.ngdc.noaa.gov/IAGA/vmod/igrf14.f), fetched 2026-10-14.
It covers historical epoch interpolation, final secular variation, WGS84 conversion,
altitude and both poles. Regenerate with Tools/MagneticFields/generate_igrf_fixtures.py and gfortran.
The small differences allowed by the test reflect rounded constants in that routine.
