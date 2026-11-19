#!/usr/bin/env python3
"""The interior Schwarzschild metric solves Einstein's equations (Appendix A).

Purpose  Verify symbolically each step of Appendix A: B(r) from the tt equation,
         the pressure p(r) solving the Tolman-Oppenheimer-Volkoff equation,
         hydrostatic equilibrium, matching to the exterior metric, the rr
         equation, and the divergence of the central pressure at mu = 8/9.
Paper    Eqs. (9)-(10) and Appendix A, Eqs. (A1)-(A5).
Inputs   none.  Units G = c = 1, constant density rho, radius R.
Outputs  "PASS | ..." / "FAIL | ..." lines and a summary; exit status 0 if all pass.
Runtime  a few seconds (SymPy).
"""
import sys

import sympy as sp

r, R, rho, x = sp.symbols("r R rho x", positive=True)
results = []


def check(label, ok):
    results.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'} | {label}", flush=True)


m = sp.Rational(4, 3) * sp.pi * rho * r**3            # mass inside radius r
M = m.subs(r, R)
rs = 2 * M                                            # Schwarzschild radius
mu = sp.simplify(rs / R)                              # compactness
a = sp.sqrt(1 - mu)
beta = sp.sqrt(1 - rs * r**2 / R**3)

B = 1 / (1 - 2 * m / r)                               # (A1): tt equation
p = rho * (beta - a) / (3 * a - beta)                 # (A5): pressure
A = ((3 * a - beta) / 2) ** 2                         # Eq. (9)

tt = sp.diff(r * (1 - 1 / B), r) - 8 * sp.pi * r**2 * rho     # tt Einstein equation for a static star
check("(A1) B(r) solves the tt equation and equals Eq. (10)",
      sp.simplify(tt) == 0 and sp.simplify(B - 1 / (1 - rs * r**2 / R**3)) == 0)
tov = sp.diff(p, r) + (rho + p) * (m + 4 * sp.pi * r**3 * p) / (r * (r - 2 * m))
check("(A4)-(A5) p(r) solves the Tolman-Oppenheimer-Volkoff equation", sp.simplify(tov) == 0)
check("(A5) the pressure vanishes at the surface, p(R) = 0", sp.simplify(p.subs(r, R)) == 0)
check("(A3) hydrostatic equilibrium: (rho + p) sqrt(A) = rho sqrt(1 - mu)", sp.simplify((rho + p) * (3 * a - beta) / 2 - rho * a) == 0)
check("matching to the exterior metric: A(R) = 1 - mu", sp.simplify(A.subs(r, R) - (1 - mu)) == 0)
rr = sp.diff(A, r) / A - 2 * (m + 4 * sp.pi * r**3 * p) / (r * (r - 2 * m))
check("the rr Einstein equation is satisfied by Eq. (9)", sp.simplify(rr) == 0)
p0 = (p.subs(r, 0) / rho).subs(rho, 3 * x / (8 * sp.pi * R**2))   # central pressure as a function of mu = x
check("the central pressure p(0) diverges as mu -> 8/9 (Buchdahl bound)",
      sp.limit(p0, x, sp.Rational(8, 9), dir="-") == sp.oo)

print(f"SUMMARY | {sum(results)}/{len(results)} PASS")
sys.exit(0 if all(results) else 1)
