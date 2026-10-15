"""Generate observed-CIP matrices with an external ERFA 2.0.1 shared library.

Usage: python3 Tools/generate_observed_pole_fixtures.py /path/to/liberfa.dylib
No C runtime or build target is added to CelMek.
"""
import argparse
import ctypes
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('library')
args = parser.parse_args()
lib = ctypes.CDLL(args.library)
D = ctypes.c_double
Matrix = (D * 3) * 3
M = ctypes.POINTER(D * 3)
for name, parameters in {
    'eraPnm06a': [D,D,M], 'eraBp00': [D,D,M,M,M],
    'eraNut00a': [D,D,ctypes.POINTER(D),ctypes.POINTER(D)],
    'eraBpn2xy': [M,ctypes.POINTER(D),ctypes.POINTER(D)],
    'eraS06': [D,D,D,D], 'eraObl80': [D,D],
    'eraIr': [M], 'eraRx': [D,M], 'eraRz': [D,M],
    'eraRxr': [M,M,M], 'eraC2ixys': [D,D,D,M],
}.items():
    function = getattr(lib,name)
    function.argtypes = parameters
    function.restype = D if name in ['eraS06','eraObl80'] else None

def reference(day, dx, dy, old):
    b, r = Matrix(), Matrix()
    x,y = D(),D()
    lib.eraPnm06a(2400000.5,day,b)
    lib.eraBpn2xy(b,ctypes.byref(x),ctypes.byref(y))
    s = lib.eraS06(2400000.5,day,x.value,y.value)
    if old:
        rb,rp,bp,n = Matrix(),Matrix(),Matrix(),Matrix()
        dp,de = D(),D()
        lib.eraBp00(2400000.5,day,rb,rp,bp)
        lib.eraNut00a(2400000.5,day,ctypes.byref(dp),ctypes.byref(de))
        t = ((2400000.5-2451545)+day)/36525
        eps = lib.eraObl80(2400000.5,day)-0.02524*t*3.141592653589793/(180*3600)
        lib.eraIr(n); lib.eraRx(eps,n); lib.eraRz(-dp.value,n); lib.eraRx(-(eps+de.value),n)
        lib.eraRxr(n,bp,b)
        lib.eraBpn2xy(b,ctypes.byref(x),ctypes.byref(y))
    lib.eraC2ixys(x.value+dx,y.value+dy,s,r)
    return [r[i][j] for i in range(3) for j in range(3)]

rows = []
for day in [37665.,51544.5,53736.,57754.,60000.,62502.]:
    for old in [False,True]:
        for dx,dy in [(0.,0.),(1.2e-9,-2.3e-9),(-8e-9,5e-9)]:
            rows.append(dict(mjd=day,old=old,dx=dx,dy=dy,matrix=reference(day,dx,dy,old)))
Path('Tests/CelMekTests/Fixtures/observed-pole-erfa.json').write_text(json.dumps(rows,indent=2)+'\n')
