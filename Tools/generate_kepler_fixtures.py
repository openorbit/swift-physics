#!/usr/bin/env python3
"""Independent 80-digit Decimal bisection of Kepler's equations.

Print Swift fixture tuples using exact binary64 inputs. No production solver or
external numerical package is used. Small anomalies exercise cancellation.
"""
from decimal import Decimal, localcontext


def trig(x, hyperbolic):
    term = total = x
    for n in range(1, 400):
        term *= x * x / Decimal((2*n)*(2*n+1))
        if not hyperbolic:
            term = -term
        updated = total + term
        if updated == total:
            return total
        total = updated
    raise RuntimeError('Series did not converge')


with localcontext() as context:
    context.prec = 80
    for hyperbolic, eccentricity, mean in [
        (False, 0.9, 0.1), (False, 0.9999999999999999, 1e-18),
        (False, 0.999999, 1e-8), (False, 0.5, 2.0),
        (True, 1.0000000000000002, 1e-18),
        (True, 1.000001, 1e-8), (True, 2.0, 10.0), (True, 1.1, 1000.0),
    ]:
        e, m = map(Decimal.from_float, (eccentricity, mean))
        low, high = Decimal(0), Decimal(20 if hyperbolic else 4)
        for _ in range(300):
            x = (low + high)/2
            f = e*trig(x, True)-x-m if hyperbolic else x-e*trig(x, False)-m
            if f < 0:
                low = x
            else:
                high = x
        print(f'({str(hyperbolic).lower()}, {eccentricity!r}, {mean!r}, {float((low+high)/2)!r}),')
