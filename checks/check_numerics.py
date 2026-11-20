#!/usr/bin/env python3
"""Numerical results of the paper, recomputed from the exact Snell law.

Purpose  Recompute the numerical results listed below from the exact Snell law
         and compare them with expected_values.json.
Paper    Table I and Sec. VI F (accuracy of NLO/NNLO), Sec. IV (proper-time
         example), Sec. VII B (mu_c, mu_s, antipodal tunnels, mu = 0.888).
Inputs   none (reads ../expected_values.json).
Outputs  one line per check, "PASS | ..." or "FAIL | ...", and a summary line;
         exit status 0 if everything passes.
Runtime  a few seconds.
"""
import json
import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import gravity_train as g  # noqa: E402

EXP = json.loads((pathlib.Path(__file__).resolve().parents[1] / "expected_values.json").read_text())
results = []


def check(label, ok, detail=""):
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'} | {label}{(' | ' + detail) if detail else ''}", flush=True)


# Table I
e = EXP["table_I_delta90"]
for m, ref in e["values"].items():
    q = g.find_q(math.pi / 2, float(m))
    check(f"Table I, mu = {m}: exact q", abs(q - ref) < e["tolerance"], f"computed {q:.6f}, paper {ref}")

# Newtonian limit
check("Newtonian limit q = 1 - Delta/pi (Delta = 90 deg)", abs(g.find_q(math.pi / 2, 0.0) - 0.5) < 1e-10)

# Accuracy of the NNLO formula over 10..170 deg
e = EXP["nnlo_max_rel_error"]
for m, bound in e["upper_bounds"].items():
    err = max(abs(g.q_series(math.radians(d), float(m))[2] - g.find_q(math.radians(d), float(m)))
              / g.find_q(math.radians(d), float(m)) for d in range(10, 171, 10))
    check(f"Sec. VI F, mu = {m}: max NNLO relative error < {bound:g}", err < bound, f"computed {err:.2e}")

# Proper-time versus coordinate-time optimum (Sec. IV)
e = EXP["proper_time_example"]
q_tau = g.find_q(math.pi / 2, 0.3, clock="proper")
q_t = g.find_q(math.pi / 2, 0.3, clock="coordinate")
check("Sec. IV: proper-time optimum q_tau (Delta = 90 deg, mu = 0.3)", abs(q_tau - e["q_tau"]) < e["tolerance"], f"{q_tau:.5f}")
check("Sec. IV: coordinate-time optimum q_t", abs(q_t - e["q_t"]) < e["tolerance"], f"{q_t:.5f}")

# Thresholds (Sec. VII B)
for key, fn in (("mu_c", g.mu_c), ("mu_s", g.mu_s)):
    v = fn()
    check(f"{EXP[key]['where']}: {key}", abs(v - EXP[key]["value"]) < EXP[key]["tolerance"], f"{v:.8f}")
check("Newtonian limit of K(mu) is -pi/2", abs(g.K(1e-7) + math.pi / 2) < 1e-4)

# Antipodal tunnels at mu = 0.8
e = EXP["antipodal_mu_0p8"]
sols = g.fastest_tunnel(math.pi, 0.8)
best, diam = sols[0], [s for s in sols if s.get("note") == "diameter"][0]
check("Sec. VII B, mu = 0.8: curved tunnel is the fastest", best.get("note") != "diameter")
check("Sec. VII B, mu = 0.8: q of curved tunnel", abs(best["q"] - e["q_curved"]) < e["tolerance"], f"{best['q']:.6f}")
check("Sec. VII B, mu = 0.8: time of curved tunnel", abs(best["time"] - e["t_curved"]) < e["tolerance"], f"{best['time']:.4f} R/c")
check("Sec. VII B, mu = 0.8: time of diameter", abs(diam["time"] - e["t_diameter"]) < e["tolerance"], f"{diam['time']:.4f} R/c")

# Threshold mu_s from direct quadrature, independent of the formula for K (Sec. VII B)
from scipy.optimize import brentq  # noqa: E402


def slope(mu, q=0.002):
    """(Delta(q) - pi)/q for a tunnel passing close to the centre; tends to 2K(mu) as q -> 0."""
    return (g.opening(q, mu) - math.pi) / q


for m in (0.74, 0.78):
    s, k2 = slope(m), 2 * g.K(m)
    check(f"Sec. VII B, mu = {m}: direct (Delta - pi)/q at q = 0.002 equals 2K", abs(s - k2) < 1e-4,
          f"{s:+.5f} vs {k2:+.5f}")
mu_num = brentq(lambda m: slope(m, 0.001), 0.765, 0.772, xtol=1e-9)
check("Sec. VII B: sign change of the direct slope reproduces mu_s", abs(mu_num - g.mu_s()) < 1e-5,
      f"{mu_num:.6f} vs {g.mu_s():.6f}")


def antipodal_curved(mu):
    """Curved antipodal tunnels: roots of Delta(q) = pi on a grid of turning radii."""
    qs = [i / 100 for i in range(1, 96)]
    v = [g.opening(q, mu) - math.pi for q in qs]
    return [brentq(lambda q: g.opening(q, mu) - math.pi, qs[i], qs[i + 1])
            for i in range(len(qs) - 1) if v[i] * v[i + 1] < 0]


check("Sec. VII B, mu = 0.76 < mu_s: no curved antipodal tunnel", not antipodal_curved(0.76))
roots = antipodal_curved(0.775)
ok = len(roots) == 1 and g.travel_time(roots[0], 0.775) < g.diameter_time(0.775)
check("Sec. VII B, mu = 0.775 > mu_s: one curved antipodal tunnel, faster than the diameter", ok,
      f"q = {roots[0]:.4f}, t = {g.travel_time(roots[0], 0.775):.4f} vs {g.diameter_time(0.775):.4f} R/c" if roots else "none")

# Near the Buchdahl bound
e = EXP["near_buchdahl_mu_0p888"]
for d, key in ((90, "q_delta90"), (180, "q_delta180")):
    q = g.fastest_tunnel(math.radians(d), 0.888)[0]["q"]
    check(f"Sec. VII B, mu = 0.888: fastest q for Delta = {d} deg", abs(q - e[key]) < e["tolerance"], f"{q:.4f}")

print(f"SUMMARY | {sum(results)}/{len(results)} PASS")
sys.exit(0 if all(results) else 1)
