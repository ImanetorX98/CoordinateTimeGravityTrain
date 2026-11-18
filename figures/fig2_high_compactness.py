#!/usr/bin/env python3
"""Fig. 2 of the paper: fastest coordinate-time tunnels at high compactness.

Purpose  (a) fastest antipodal tunnels for several mu: the diameter below mu_s,
         one of two mirror-image curved tunnels above; (b) turning radius q of
         the fastest tunnel versus mu for Delta = 90 and 180 deg, with the NNLO
         formula and the thresholds mu_s, mu_c and the Buchdahl bound.
         For every mu all Snell solutions joining the stations (winding ones and
         the diameter included) are compared and the fastest is kept.
Paper    Fig. 2, Sec. VII B.
Inputs   optional: --quick (coarser mu grid for panel b).
Outputs  figures/output/fig2_high_compactness.pdf and .png
Runtime  about 15 s (a few seconds with --quick).
"""
import argparse
import math
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "figures"))
import gravity_train as g  # noqa: E402
import style  # noqa: E402


def fastest_q(delta, mu):
    return 1 - delta / math.pi if mu == 0 else g.fastest_tunnel(delta, mu)[0]["q"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    mus_s, mus_c = g.mu_s(), g.mu_c()
    plt = style.apply()
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.2, 3.6), gridspec_kw=dict(width_ratios=[1, 1.25]))

    th = np.linspace(0, 2 * math.pi, 600)
    ax.plot(np.cos(th), np.sin(th), color=style.INK, lw=1.2)
    curve_mus = [0.70, 0.78, 0.80, 0.83, 0.86, 0.88, 0.888]
    for mu, c in zip(curve_mus, np.linspace(0.0, 0.85, len(curve_mus))):
        col = style.PLASMA(c)
        q = fastest_q(math.pi, mu)
        if q == 0.0:
            ax.plot([0, 0], [-1, 1], color=col, lw=2, label=rf"$\mu={mu:g}$ (diameter)")
            continue
        X, Y = g.tunnel_curve(q, mu)
        ax.plot(Y, X, color=col, lw=2, label=rf"$\mu={mu:g}$")          # stations at the poles
        ax.plot(-Y, X, color=col, lw=0.7, ls=(0, (1.5, 2)), alpha=0.55)   # mirror image
    ax.plot([0, 0], [1, -1], "o", color=style.INK, ms=4)
    ax.set_aspect("equal")
    ax.set_xlim(-1.08, 1.08)
    ax.set_ylim(-1.08, 1.08)
    ax.axis("off")
    ax.set_title("(a) fastest antipodal tunnels", fontsize=10, loc="left")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=2, fontsize=7.5, frameon=False)

    if args.quick:
        mus = np.concatenate([np.linspace(0, 0.74, 6), np.linspace(0.75, 0.885, 8), [mus_s + 1e-3, 0.888]])
    else:
        mus = np.concatenate([np.linspace(0, 0.74, 15), np.linspace(0.75, 0.885, 28), [mus_s + 1e-3, 0.887, 0.888]])
    mus = np.unique(np.round(mus, 6))
    for d, col in ((90, style.PLASMA(0.15)), (180, style.PLASMA(0.68))):
        qs = [fastest_q(math.radians(d), float(m)) for m in mus]
        bx.plot(mus, qs, color=col, lw=2, label=rf"$\Delta={d}^\circ$, exact")
    m = np.linspace(0, 0.6, 100)
    bx.plot(m, [g.q_series(math.pi / 2, v)[2] for v in m], color=style.PLASMA(0.15), lw=1, ls="--",
            label=r"$\Delta=90^\circ$, NNLO")
    for v, lab, ha in ((mus_s, r"$\mu_s$", "center"), (mus_c, r"$\mu_c$", "right"), (g.MU_BUCHDAHL, r"$\mu_B$", "left")):
        bx.axvline(v, color=style.MUTED, lw=0.8, ls=":")
        bx.text(v, 1.02, lab, ha=ha, va="bottom", fontsize=9, color=style.MUTED, transform=bx.get_xaxis_transform())
    bx.set_xlim(0, 0.92)
    bx.set_ylim(-0.02, 1.0)
    bx.set_xlabel(r"compactness $\mu = r_s/R$")
    bx.set_ylabel(r"turning radius $q = r_*/R$")
    bx.grid(True, color="#e6e6e6", lw=0.6)
    bx.spines[["top", "right"]].set_visible(False)
    bx.legend(loc="upper left", fontsize=8, frameon=False)
    bx.set_title("(b) fastest tunnel versus compactness", fontsize=10, loc="left", pad=14)
    out = ROOT / "figures" / "output"
    out.mkdir(exist_ok=True)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(out / f"fig2_high_compactness.{ext}", dpi=300, bbox_inches="tight")
    print("saved", out / "fig2_high_compactness.pdf")


if __name__ == "__main__":
    main()
