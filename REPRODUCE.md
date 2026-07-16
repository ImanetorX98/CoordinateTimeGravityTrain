# Reproducing the paper

Constrained-worldline brachistochrones in non-stationary spacetimes
(Kodama rails, closed forms, plunge inversion). This repository is self-contained
for the **non-stationary** metrics of the paper: FLRW, Vaidya, and Thakurta--Kerr,
together with the stationary Kerr machinery they build on.

## Layout
- `paper/` — LaTeX sources (`main.tex` = CQG/iopart, `main_prd_revtex.tex` = PRD),
  bibliographies, class/style files. (Figures are regenerated from the scripts
  below; per repository policy binaries are not tracked.)
- `FLRWmetric/`, `VaidyaMetric/`, `ThakurtaMetric/` — analysis + figure scripts for
  the three non-stationary metrics.
- `KerrScripts/` — the stationary-Kerr scripts the paper uses (genus-2 machinery,
  quasi-periods, third-kind theta, separatrices, penetration, colormaps, and the
  figure generators for the Kerr panels). Copied from the companion
  `KerrBrachistochrone` repository so this repo reproduces the paper on its own.
- `reproduce_reductions.py` — one-file check of the algebraic heart: solves the
  11x11 reduction identity and prints the reference coefficients `c_k`
  (App. "Reproducing the reductions").
- `paper_style.py` — shared matplotlib style helper imported by the figure scripts.

## Install
```
pip install -r requirements.txt
```
SageMath (+ the `abelfunctions` library) is needed **only** for the two `.sage`
scripts in `KerrScripts/` (genus-2 Riemann surface / theta). Everything else is
pure Python.

## Reproduce the closed-form results
```
python reproduce_reductions.py                       # algebraic c_k, residual 0
python KerrScripts/pipeline_completa_deltaphi.py      # full (M,a,E,J) -> delta phi
python KerrScripts/kerr_tbranch_psi_assembly.py       # t-branch on-curve assembly
python VaidyaMetric/vaidya_fully_explicit.py          # Vaidya explicit delta phi
python VaidyaMetric/vaidya_horizon_dilog.py           # horizon dilogarithm D_k
python VaidyaMetric/vaidya_ell_dilog_match.py         # separatrix Bloch-Wigner (sigma,zeta)
python VaidyaMetric/vaidya_asymmetry.py               # accretion/evaporation split
python FLRWmetric/flrw_perturbed_genus_jump.py        # perturbed-FLRW genus jump (proof of concept)
```
Table "Result-to-script map" in the paper appendix lists every closed form with the
script that generates and verifies it, and its achieved numerical residual.

## Regenerate the figures, then compile
Each `fig_*.py` / `genera_*.py` / `colormap_*.py` script writes its figure(s); run
them, place the output in `paper/Immagini/`, then
```
cd paper && pdflatex main.tex && pdflatex main.tex
```
(likewise `main_prd_revtex.tex`).
