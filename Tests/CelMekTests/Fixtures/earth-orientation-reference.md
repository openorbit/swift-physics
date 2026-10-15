Earth-orientation matrices generated using ERFA 2.0.1 in an external temporary
shared library. Generator: Tools/generate_earth_orientation_fixtures.py.
No C code or ERFA library is required to build or test this Swift package.

150 matrices cover all five edges at 30 split TT dates spanning J2000 ±1000
Julian years, including fractional dates near integer-day boundaries. These
are numerical port regression cases, not a claim of physical accuracy over
that entire interval. UT1−TT = −64.184 s; xp = 2.55060238e−7 rad;
yp = 1.860359247e−6 rad. Adaptation license: Documentation/licenses/ERFA.txt.
