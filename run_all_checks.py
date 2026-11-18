#!/usr/bin/env python3
"""Run every check of the paper in one go and print a summary table.

Purpose  Reproduce, from scratch, every analytic and numerical result of the paper.
Inputs   optional: --with-mathematica  also run checks/verify_all.wl through
         wolframscript (requires a Wolfram Engine or Mathematica installation).
Outputs  the output of each check, then a table with one line per check script
         (PASS/FAIL, number of checks, time); exit status 0 if everything passes.
Runtime  about 20 s without Mathematica.
"""
import argparse
import pathlib
import re
import shutil
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent
CHECKS = [
    ("Appendix A: interior Schwarzschild metric", [sys.executable, "checks/check_interior_metric.py"]),
    ("Secs. VI-VII: analytic results (SymPy)", [sys.executable, "checks/check_derivations.py"]),
    ("Table I, Secs. IV, VI F, VII B: numerical results", [sys.executable, "checks/check_numerics.py"]),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--with-mathematica", action="store_true")
    args = ap.parse_args()
    checks = list(CHECKS)
    if args.with_mathematica:
        if shutil.which("wolframscript") is None:
            sys.exit("wolframscript not found: install the Wolfram Engine or drop --with-mathematica")
        checks.append(("Secs. VI-VII: analytic results (Mathematica)", ["wolframscript", "-file", "checks/verify_all.wl"]))
    rows = []
    for name, cmd in checks:
        print(f"\n=== {name} ===", flush=True)
        t0 = time.time()
        proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
        out = proc.stdout + proc.stderr
        print(out.rstrip())
        n_pass = len(re.findall(r"^PASS \||^True$", out, re.M))
        n_fail = len(re.findall(r"^FAIL \||^False$", out, re.M))
        ok = proc.returncode == 0 and n_fail == 0 and n_pass > 0
        rows.append((name, ok, n_pass, n_fail, time.time() - t0))
    print("\n" + "=" * 78)
    print(f"{'check':52s} {'result':6s} {'passed':>7s} {'failed':>7s} {'time':>5s}")
    for name, ok, p, f, dt in rows:
        print(f"{name:52s} {'PASS' if ok else 'FAIL':6s} {p:7d} {f:7d} {dt:4.0f}s")
    print("=" * 78)
    all_ok = all(r[1] for r in rows)
    print("ALL CHECKS PASSED" if all_ok else "SOME CHECKS FAILED")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
