#!/usr/bin/env python3
"""Independent 80-digit conic propagation fixtures (Python standard library).

Uses eccentric/hyperbolic anomalies, not universal variables or Swift code.
Exact binary64 initial state, GM and elapsed-time inputs are lifted to Decimal.
Run from the repository root to regenerate Tests/CelMekTests/Fixtures/two_body_reference.json.
"""
import json
from decimal import Decimal, localcontext
from pathlib import Path

PI = Decimal('3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628')
D = Decimal.from_float


def sincos(x, hyperbolic=False):
    sin, cos = x, Decimal(1)
    st, ct = sin, cos
    for n in range(1, 400):
        sign = 1 if hyperbolic else -1
        st *= sign*x*x/Decimal(2*n*(2*n+1))
        ct *= sign*x*x/Decimal((2*n-1)*2*n)
        ns, nc = sin+st, cos+ct
        if ns == sin and nc == cos:
            return sin, cos
        sin, cos = ns, nc
    raise RuntimeError('Series did not converge')


def solve(f, target):
    sign = -1 if target < 0 else 1
    target = abs(target)
    lo, hi = Decimal(0), Decimal(1)
    while f(hi) < target:
        hi *= 2
    for _ in range(320):
        mid = (lo+hi)/2
        if f(mid) < target:
            lo = mid
        else:
            hi = mid
    return sign*(lo+hi)/2


def state(r0, v0, mu, dt):
    r0, v0, mu, dt = map(D, (r0, v0, mu, dt))
    e = v0*v0*r0/mu-1
    if e == 1:
        p = 2*r0
        b = 2*(mu/(p*p*p)).sqrt()*dt
        d = solve(lambda x: x+x*x*x/3, b)
        x, y = r0*(1-d*d), 2*r0*d
        scale = (mu/p).sqrt()
        vx, vy = -2*scale*d/(1+d*d), 2*scale/(1+d*d)
    elif e < 1:
        a = r0/(1-e)
        mean = (mu/(a*a*a)).sqrt()*dt
        mean %= 2*PI
        if mean > PI:
            mean -= 2*PI
        if mean < -PI:
            mean += 2*PI
        anomaly = solve(lambda x: x-e*sincos(x)[0], mean)
        s, c = sincos(anomaly)
        q = (1-e*e).sqrt()
        radius = a*(1-e*c)
        x, y = a*(c-e), a*q*s
        scale = (mu*a).sqrt()/radius
        vx, vy = -scale*s, scale*q*c
    else:
        a = r0/(e-1)
        mean = (mu/(a*a*a)).sqrt()*dt
        anomaly = solve(lambda x: e*sincos(x, True)[0]-x, mean)
        s, c = sincos(anomaly, True)
        q = (e*e-1).sqrt()
        radius = a*(e*c-1)
        x, y = a*(e-c), a*q*s
        scale = (mu*a).sqrt()/radius
        vx, vy = -scale*s, scale*q*c
    return [float(x), float(y), 0], [float(vx), float(vy), 0]


with localcontext() as context:
    context.prec = 80
    fixtures = []
    for name, r0, v0, mu in [
        ('circular', 1., 1., 1.),
        ('elliptic', 1., 1.2, 1.),
        ('near-parabolic-bound', 1., 1.414213562373095, 1.),
        ('parabolic', 1., 1., 0.5),
        ('near-parabolic-unbound', 1., 1.4142135623730951, 1.),
        ('hyperbolic', 1., 2., 1.),
        ('earth-elliptic', 7e6, 9000., 3.986004418e14),
    ]:
        for dt in ([0.1, -1., 10., -100.] if r0 == 1 else [100., -1000., 10000.]):
            position, velocity = state(r0, v0, mu, dt)
            fixtures.append(dict(name=name, radius=r0, speed=v0, gm=mu, seconds=dt,
                                 position=position, velocity=velocity))
    Path('Tests/CelMekTests/Fixtures/two_body_reference.json').write_text(
        json.dumps(dict(provenance='80-digit Decimal conic equations; generate_two_body_fixtures.py',
                        cases=fixtures), indent=2)+'\n')
