# Constrained-worldline brachistochrones in non-stationary spacetimes

Source code and manuscript for the paper *Constrained-worldline brachistochrones in
non-stationary spacetimes: Kodama rails, closed forms, and plunge inversion*
(CQG/iopart and PRD/revtex versions in `paper/`).

The fastest constrained "rail" worldline between two events is extended from
stationary to **non-stationary** spacetimes — FLRW, Vaidya, and Thakurta--Kerr —
built on the stationary Kerr optical/genus-2 machinery. The repository is
self-contained for these results: every closed form has a script that generates and
numerically verifies it, and every figure in the paper has a generating script.

> Per repository policy, generated binaries (`*.png`, `*.pdf`, videos) are **not**
> tracked. Figures are regenerated from the scripts below.

---

## Layout

| Path | Contents |
|------|----------|
| `paper/` | LaTeX sources: `main.tex` (CQG/iopart), `main_prd_revtex.tex` (PRD), bibliographies, class/style files. |
| `FLRWmetric/` | FLRW (degenerate base) analysis + figure scripts, incl. the perturbed-FLRW genus-jump proof of concept. |
| `VaidyaMetric/` | Vaidya (dynamical black hole) analysis + figure scripts: explicit $\delta\varphi$, horizon dilogarithm, separatrix Bloch–Wigner, accretion/evaporation asymmetry. |
| `ThakurtaMetric/` | Thakurta--Kerr (conformal rotating) analysis + figure scripts: adiabatic hybrid, breathing families, cusp, inversion. |
| `KerrScripts/` | Stationary-Kerr machinery the paper builds on: genus-2 reduction, quasi-periods, third-kind theta, separatrices, penetration, colormaps, and the Kerr figure generators. Copied from the companion `KerrBrachistochrone` repo so this repo reproduces on its own. |
| `reproduce_reductions.py` | One-file check of the algebraic heart (App. "Reproducing the reductions"). |
| `paper_style.py` | Shared matplotlib style helper imported by the figure scripts. |

---

## Install

```bash
python3 -m pip install -r requirements.txt      # numpy scipy sympy matplotlib mpmath Pillow
```

**SageMath** (with the `abelfunctions` library) is required **only** for the four
`.sage` scripts in `KerrScripts/` (genus-2 Riemann surface, quasi-periods,
third-kind theta, Abel map) — these belong to the stationary-Kerr machinery.
**Everything else — all FLRW, Vaidya and Thakurta–Kerr analysis, every closed-form
check, and all figures — runs on the pure-Python stack** (sympy for the algebra,
`scipy.integrate.quad` for the abelian building blocks, `mpmath` for the elliptic
$\sigma,\zeta$ on the separatrix). No Sage needed there.

---

## How to use the scripts

All scripts are **standalone** — run with `python3 <script>.py` (or `sage <script>.sage`).
There are two kinds:

- **Verification scripts** print numerical residuals to stdout (a closed form vs. its
  direct integral); a residual near machine precision means the identity holds.
- **Figure scripts** (`fig_*.py`, `genera_*.py`, `colormap_*.py`) compute and save a
  figure (into `paper/Immagini/` or beside the script); they need `paper_style.py`
  on the path (run them from the repo root, e.g. `python3 KerrScripts/fig_....py`).

### Closed-form verification (reproduce the paper's results)

| Command | What it does | Expected residual |
|---------|--------------|-------------------|
| `python3 reproduce_reductions.py` | Solves the 11×11 reduction identity; prints the reference coefficients $c_k$ (rational, exact) for the TK $\tau$/$t$ branches and Vaidya. | `0` exactly |
| `python3 KerrScripts/pipeline_completa_deltaphi.py` | Full chain $(M,a,E,J)\to\delta\varphi$: solve → $c_k,b,Q_{kj}$ → building blocks $U_k,W_{kj},\mathcal I_j$ → assembled $\delta\varphi$ vs direct. | $4\times10^{-15}$ |
| `python3 KerrScripts/kerr_tbranch_psi_assembly.py` | $t$-branch on-curve clock and full $\psi$ assembly. | $9\times10^{-13}$ |
| `python3 KerrScripts/kerr_psi_explicit_verified.py` | TK $\tau$ reduction + $\psi=\tfrac12\hat E\sum Q_{kj}W_{kj}$. | $5\times10^{-14}$ |
| `sage KerrScripts/kerr_quasiperiods_bel.sage` | Genus-2 quasi-periods $\eta$ (Legendre relation). | $10^{-12}$ |
| `sage KerrScripts/kerr_thirdkind_theta_closed.sage` | Third-kind integral $L$ closed as a $\theta$-ratio (Fay). | $5\times10^{-8}$ |
| `python3 VaidyaMetric/vaidya_fully_explicit.py` | Vaidya fully explicit $\delta\varphi$ (all building blocks). | $2\times10^{-13}$ |
| `python3 VaidyaMetric/vaidya_horizon_dilog.py` | Horizon dilogarithm $\mathcal D_k=U_k\ln-\mathcal G_k$. | $10^{-12}$ |
| `python3 VaidyaMetric/vaidya_ell_dilog_match.py` | Separatrix elliptic closed form in Weierstrass $\sigma,\zeta$; Bloch–Wigner. | $10^{-14}$ |
| `python3 VaidyaMetric/vaidya_asymmetry.py` | Accretion/evaporation split $\mathcal A_\infty,\mathcal B_{\rm hor}$. | $10^{-15}$ |
| `python3 FLRWmetric/flrw_perturbed_genus_jump.py` | Perturbed-FLRW proof of concept: genus 0→1→2 with the tidal multipole. | — |

The paper appendix table *Result-to-script map* lists every closed form with its
script and achieved residual.

### Figures

Each figure script writes the panels of one paper figure. The mapping is by name:
`fig_<name>.py` (or a `genera_*`/`colormap_*` script) produces `Immagini/fig_<name>.*`.
For example:

```bash
python3 FLRWmetric/genera_figure_flrw.py            # FLRW kinematics / variational panels
python3 VaidyaMetric/genera_figure_vaidya.py        # Vaidya trajectories, timing, bounce
python3 ThakurtaMetric/genera_figure_thakurta.py    # Thakurta-Kerr branches, breathing, cusp
python3 KerrScripts/fig_penetrate_bounce_orbit.py   # ergosphere penetrate-and-bounce orbit
python3 KerrScripts/kerr_penetration_threshold.py   # penetration phase diagram
python3 KerrScripts/fig_atlas_tau_t.py              # t- and tau-branch atlases
python3 KerrScripts/colormap_conforme_estremi_fissi.py  # fixed-endpoint conformal colormap
```

### Compile the paper

Regenerate the figures into `paper/Immagini/`, then:

```bash
cd paper
pdflatex main.tex && pdflatex main.tex                    # CQG version
pdflatex main_prd_revtex.tex && pdflatex main_prd_revtex.tex   # PRD version
```

---

## Notes

- The stationary Kerr study lives in its own repository, `KerrBrachistochrone`; the
  subset used by this paper is mirrored in `KerrScripts/`.
- Earlier gravity-train / interior-Schwarzschild scripts from the initial release
  belong to a separate study and are not part of the non-stationary paper.
