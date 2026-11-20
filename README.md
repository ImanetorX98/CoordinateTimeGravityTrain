# The gravitational train in the interior Schwarzschild metric

Companion code for

> I. Rosignoli, *The Gravitational Train in the Interior Schwarzschild Metric: NLO and NNLO
> Perturbative Expansions*, submitted to European Journal of Physics (2026).

A particle is released from rest at one point of the surface of a homogeneous star and slides
without friction through a tunnel to another surface point. Which tunnel minimises the travel time?
In general relativity the answer is a geodesic of an optical metric, and the fastest tunnel
follows a radial Snell law. This repository reproduces from scratch **the results of the paper listed below**
and contains the **computational exercise for students** described in Sec. VI F.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24378845.svg)](https://doi.org/10.5281/zenodo.24378845)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ImanetorX98/CoordinateTimeGravityTrain/blob/main/notebooks/gravity_train_exercise.ipynb)

---

## Quick start

```bash
python3 -m pip install -r requirements.txt     # numpy, scipy, sympy, mpmath, matplotlib
python3 run_all_checks.py                      # about 20 s
```

Expected final output:

```
check                                                result  passed  failed  time
Appendix A: interior Schwarzschild metric            PASS         7       0    0s
Secs. VI-VII: analytic results (SymPy)               PASS        40       0    7s
Table I, Secs. IV, VI F, VII B: numerical results    PASS        26       0    5s
ALL CHECKS PASSED
```

If Mathematica or the free Wolfram Engine is installed, `python3 run_all_checks.py --with-mathematica`
also runs an independent implementation of the symbolic checks in Wolfram Language.

To regenerate the curves of Figs. 1 and 2:

```bash
python3 figures/fig1_exact_tunnels.py           # about a minute
python3 figures/fig2_high_compactness.py       # tens of seconds
```

They are written to `figures/output/` (PDF and PNG).

---

## Where the results are checked

| Result in the paper | Script | What is checked |
|---|---|---|
| Interior Schwarzschild metric, Eqs. (9)-(10), Appendix A | `checks/check_interior_metric.py` | Einstein's equations, TOV equation, matching at the surface, Buchdahl bound |
| Expansions of A, B and n, Sec. VI A-B | `checks/check_derivations.py` | exact series in the compactness |
| Snell integrand I0, I1, I2, Eqs. (47)-(49) | `checks/check_derivations.py` | exact series |
| Closed-form integrals and Eq. (54), Sec. VI D | `checks/check_derivations.py` | derivatives of the primitives, end-point values |
| q1, q2 and Eq. (59), Sec. VI E | `checks/check_derivations.py` | perturbative inversion; proof that q1, q2 > 0 |
| Table I (exact values) and the accuracy of NNLO, Sec. VI F | `checks/check_numerics.py` | exact Snell law versus the series, at 10° steps from 10° to 170° |
| Earth example, Sec. VI G | `checks/check_derivations.py` | depth shift and travel time |
| Proper-time versus coordinate-time optimum, Sec. IV | `checks/check_numerics.py` | q = 0.527 versus 0.559 (Δ = 90°, μ = 0.3) |
| Buchdahl asymptotics, Sec. VII A | `checks/check_derivations.py` | A(0) and n(0) as μ → 8/9 |
| Circular optical orbits, μ_c, Eq. (70) | `checks/check_derivations.py`, `checks/check_numerics.py` | exact cubic and its root |
| Loss of minimality of the diameter, μ_s, Eqs. (71)-(73) | `checks/check_derivations.py`, `checks/check_numerics.py` | Jacobi equation, K(μ_s) = 0; independent check by direct quadrature: the sign change of (Δ − π)/q reproduces μ_s, no curved antipodal tunnel at μ = 0.76, a faster one at μ = 0.775 and 0.8 |
| Fastest tunnels near the Buchdahl bound, Sec. VII B | `checks/check_numerics.py` | turning radii at μ = 0.888 |
| Fig. 1 | `figures/fig1_exact_tunnels.py` | exact tunnels, 70° apart, for 15 values of µ |
| Fig. 2 | `figures/fig2_high_compactness.py` | fastest antipodal tunnels, and turning radius versus compactness for 90° and 180° |

The reference values and tolerances used by `check_numerics.py` are listed, with the place where
they appear in the paper, in `expected_values.json`.

---

## Repository layout

| Path | Contents |
|---|---|
| `gravity_train.py` | Core physics shared by all scripts: metric functions, optical index, Snell integrals, travel times, perturbative formulas, the high-compactness analysis. Documented function by function. |
| `run_all_checks.py` | Runs every check and prints a summary table. |
| `checks/` | `check_interior_metric.py` and `check_derivations.py` (symbolic, SymPy), `check_numerics.py` (numerical), `verify_all.wl` (symbolic, Wolfram Language). |
| `figures/` | Scripts for Figs. 1 and 2, and the shared plotting style. |
| `notebooks/gravity_train_exercise.ipynb` | The computational exercise for students (supplementary material of the paper). |
| `expected_values.json` | Values quoted in the paper, with tolerances. |

Every Python script in `checks/` and `figures/`, and `run_all_checks.py`, starts with a short header
stating its **purpose**, the **part of the paper** it reproduces, its **inputs** (command-line options,
all with defaults), its **outputs** and its **runtime**. The Python check scripts print one
`PASS | ...` or `FAIL | ...` line per check and exit with status 0 only if everything passes, so they
can be used in automated tests; the Wolfram Language check prints a label and `True`/`False` per check.

---

## The notebook for students

`notebooks/gravity_train_exercise.ipynb` guides a student, in six parts, from the optical index to
the exact relativistic tunnels: Newtonian limit, the two numerical difficulties of the problem and how
to handle them, the exact values of Table I, comparison of Newtonian, NLO, NNLO and exact turning radii, and two
open explorations (proper time; antipodal tunnels at high compactness). Each part ends with a check
value and a conceptual question.

* **In the browser, no installation:** click the *Open in Colab* badge above (a Google account is needed).
* **Locally:** `python3 -m pip install jupyterlab numpy scipy matplotlib`, then
  `jupyter lab notebooks/gravity_train_exercise.ipynb`.

The notebook is self-contained: it does not import `gravity_train.py`, so that students see every line
of the computation.

---

## Units and conventions

The code uses R = c = 1 (the paper keeps c explicit), x = r/R, compactness μ = r_s/R with 0 < μ < 8/9, turning radius
q = r_*/R, station separation Δ. The optical index is used through w = 1/n², which stays finite at
the surface, where n diverges because the particle starts at rest. Times are in units of R/c.

---

## Versions

* The scripts released with the original submission (June 2026) are preserved under the tag
  `v1.0-as-submitted`.
* Between July and November 2026 this repository temporarily hosted the code of a different project,
  now maintained at [ImanetorX98/NonStationaryBrachistochrone](https://github.com/ImanetorX98/NonStationaryBrachistochrone).
  That material is still available in the history of this repository.

## How to cite

Please cite the paper and, for the code, the Zenodo archive:

* Software, all versions: [doi:10.5281/zenodo.24378845](https://doi.org/10.5281/zenodo.24378845);
  version 2.0.2, the one used in the revised manuscript: [doi:10.5281/zenodo.24396991](https://doi.org/10.5281/zenodo.24396991)
  (version 2.0.1:
  [doi:10.5281/zenodo.24383706](https://doi.org/10.5281/zenodo.24383706); version 2.0.0:
  [doi:10.5281/zenodo.24378846](https://doi.org/10.5281/zenodo.24378846)).
* Paper (preprint): [doi:10.5281/zenodo.20617415](https://doi.org/10.5281/zenodo.20617415).

`CITATION.cff` records the all-versions software DOI and the preprint DOI.

## License

MIT, see `LICENSE`.
