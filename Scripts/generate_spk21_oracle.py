#!/usr/bin/env python3
"""Optional binary CSPICE oracle; public SPKW21/SPKE21 contracts only, no toolkit source.
PYTHONPATH=/path/to/binary/spiceypy python3 Scripts/generate_spk21_oracle.py
"""
import ctypes
import hashlib
import json
from pathlib import Path
import numpy as np
import spiceypy as spice
from spiceypy.utils.libspicehelper import libspice

out = Path(__file__).resolve().parents[1] / 'Tests/CelMekTests/Fixtures/spk21'
out.mkdir(parents=True, exist_ok=True)
refs = dict(schema=1, toolkit=spice.tkvrsn('TOOLKIT'), wrapper=spice.__version__,
            basis='Published API contracts; binary oracle only', files=[], states=[])
for dimension, count in [(15, 201), (25, 2)]:
    size = 4 * dimension + 11
    epochs = np.arange(count, dtype=np.float64) * 10 + 10
    lines = np.zeros((count, size), dtype=np.float64)
    rng = np.random.default_rng(21000 + dimension)
    for i, end in enumerate(epochs):
        line = lines[i]
        line[0] = end - 4
        line[1:dimension+1] = np.cumsum(rng.uniform(8, 12, dimension))
        line[1+dimension:7+dimension] = [1e8+i*100, 2, -2e8+i*200, -3, 3e7-i*300, 4]
        line[7+dimension:7+4*dimension] = rng.uniform(-0.02, 0.02, 3*dimension)
        line[7+4*dimension] = dimension+1
        line[8+4*dimension:11+4*dimension] = [dimension, dimension-4, dimension-7]
    path = out / f'dimension{dimension}.bsp'
    path.unlink(missing_ok=True)
    h = spice.spkopn(str(path), 'CELMEK TYPE21 BINARY ORACLE', 0)
    # No CSPICE C wrapper exists; call the documented Fortran writer entry point.
    ints = [ctypes.c_int(x) for x in [h, 900000+dimension, 0, count, size]]
    first, last = ctypes.c_double(0), ctypes.c_double(epochs[-1])
    libspice.spkw21_(ctypes.byref(ints[0]), ctypes.byref(ints[1]), ctypes.byref(ints[2]),
        ctypes.c_char_p(b'J2000'), ctypes.byref(first), ctypes.byref(last), ctypes.c_char_p(b'TYPE21 ORACLE'),
        ctypes.byref(ints[3]), ctypes.byref(ints[4]), lines.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        epochs.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), ctypes.c_int(5), ctypes.c_int(13))
    if spice.failed():
        message = spice.getmsg('LONG'); spice.reset(); raise RuntimeError(message)
    spice.spkcls(h)
    spice.furnsh(str(path))
    queries = sorted(set([0., 1.25, 6., 9.9, 10., np.nextafter(10., np.inf), float(epochs[-1])] +
        ([990., 1000., np.nextafter(1000., np.inf), 1999., 2000., np.nextafter(2000., np.inf)] if count > 100 else [])))
    for t in queries:
        state, _ = spice.spkez(900000+dimension, t, 'J2000', 'NONE', 0)
        refs['states'].append(dict(file=path.name, t=float(t), record=int(np.searchsorted(epochs,t)), state=state.tolist()))
    refs['files'].append(dict(name=path.name, dimension=dimension, count=count,
                             sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    spice.unload(str(path))
(out/'reference.json').write_text(json.dumps(refs, indent=2)+'\n')
print('Generated', len(refs['states']), 'type 21 oracle states')

# Optionally pin a real small-body kernel obtained by the extraction script.
# Existing real fixtures are retained when this argument is omitted.
if len(__import__('sys').argv) == 3 and __import__('sys').argv[1] == '--horizons':
    source = Path(__import__('sys').argv[2])
    destination = out / 'horizons-ceres.bsp'
    destination.write_bytes(source.read_bytes())
    h = spice.dafopr(str(destination)); spice.dafbfs(h)
    if not spice.daffna(): raise RuntimeError('Empty real kernel')
    descriptor = spice.dafgs(5)
    body, center, frame, kind, start, end, first, last = spice.spkuds(descriptor)
    if kind != 21 or body != 20000001: raise RuntimeError('Expected Ceres type 21')
    dimension, count = map(int, spice.dafgda(h, last-1, last))
    size = 4*dimension+11
    epochs = spice.dafgda(h, first+size*count, first+size*count+count-1)
    spice.dafcls(h); spice.furnsh(str(destination))
    queries = [start, start+(end-start)/7, (start+end)/2, end]
    for boundary in epochs:
        if start < boundary < end:
            queries += [float(np.nextafter(boundary, -np.inf)), float(boundary), float(np.nextafter(boundary, np.inf))]
    rows = []
    for t in sorted(set(queries)):
        state, _ = spice.spkez(body,t,'J2000','NONE',center)
        rows.append(dict(t=t, record=int(np.searchsorted(epochs,t)), state=state.tolist()))
    reference = dict(schema=1, toolkit=spice.tkvrsn('TOOLKIT'), wrapper=spice.__version__,
        file=destination.name, body=body, center=center, frame=frame,
        sha256=hashlib.sha256(destination.read_bytes()).hexdigest(),
        origin='NASA/JPL Horizons API; Ceres; requested 2000-01-01 through 2000-01-03 TDB',
        states=rows)
    (out/'horizons-reference.json').write_text(json.dumps(reference,indent=2)+'\n')
    spice.unload(str(destination))
    print('Generated',len(rows),'real Horizons oracle states')
