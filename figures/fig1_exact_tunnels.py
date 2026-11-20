#!/usr/bin/env python3
"""Fig. 1 of the paper: exact tunnels for stations 70 deg apart, for a family of compactness values.

Purpose  Draw the Newtonian hypocycloid and the exact relativistic tunnels for 15 values of mu,
         in the plane of the tunnel, rotated so that the symmetry axis is the bisector y = x.
Paper    Fig. 1; Eq. (34) for the exact turning radius, Eqs. (4)-(5) for the Newtonian curve.
Inputs   optional: --delta 70 (degrees), --mu-min 0.05, --mu-max 0.75, --n 15.
Outputs  figures/output/fig1_exact_tunnels.pdf and .png
Runtime  about a minute.
"""
import argparse
import math
import pathlib
import sys

import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "figures"))
import gravity_train as g  # noqa: E402
import style  # noqa: E402


def rotate(X, Y):
    """From stations symmetric about the Y axis to stations symmetric about y = x."""
    return (Y - X) / math.sqrt(2), (Y + X) / math.sqrt(2)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--delta", type=float, default=70.0, help="station separation in degrees")
    ap.add_argument("--mu-min", type=float, default=0.05)
    ap.add_argument("--mu-max", type=float, default=0.75)
    ap.add_argument("--n", type=int, default=15, help="number of compactness values")
    args = ap.parse_args()
    delta = math.radians(args.delta)
    mus = np.linspace(args.mu_min, args.mu_max, args.n)
    # plasma palette stopped at orange-yellow, so that the last curve stays visible on white
    cmap = LinearSegmentedColormap.from_list("plasma_cut", style.PLASMA(np.linspace(0.0, 0.85, 256)))
    norm = Normalize(mus[0], mus[-1])

    plt = style.apply()
    fig, ax = plt.subplots(figsize=(3.6, 3.1))
    th = np.linspace(0, math.pi / 2, 400)
    ax.plot(np.cos(th), np.sin(th), color=style.INK, lw=1.2)
    ax.plot([0, 1.05], [0, 1.05], color=style.MUTED, lw=0.6, ls=(0, (4, 3)))
    for mu in mus:
        ax.plot(*rotate(*g.tunnel_curve(g.find_q(delta, mu), mu)), color=cmap(norm(mu)), lw=1.3)
    q0 = 1 - delta / math.pi
    u = np.linspace(0, math.pi / 2, 600)
    r = np.sqrt(q0**2 + (1 - q0**2) * np.sin(u) ** 2)                 # Eq. (4)
    phi = np.arctan2(np.sin(u), q0 * np.cos(u)) - q0 * u              # Eq. (5)
    rr, pp = np.concatenate([r[::-1], r[1:]]), np.concatenate([-phi[::-1], phi[1:]])
    ax.plot(*rotate(rr * np.sin(pp), rr * np.cos(pp)), color=style.NEWTON, lw=2.2, label=r"Newtonian, $\mu=0$")
    for s in (1, -1):
        a = math.pi / 4 + s * delta / 2
        ax.plot(math.cos(a), math.sin(a), "o", color=style.INK, ms=4, zorder=5)
    lo = math.sin(math.pi / 4 - delta / 2) - 0.08
    ax.set_aspect("equal")
    ax.set_xlim(lo, 1.04)
    ax.set_ylim(lo, 1.04)
    ax.set_xlabel(r"$x/R$")
    ax.set_ylabel(r"$y/R$")
    ax.grid(True, color="#e6e6e6", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="lower left", fontsize=7.5, frameon=False)
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(r"$\mu = r_s/R$")
    cb.outline.set_edgecolor(style.MUTED)
    out = ROOT / "figures" / "output"
    out.mkdir(exist_ok=True)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(out / f"fig1_exact_tunnels.{ext}", dpi=300, bbox_inches="tight")
    print("saved", out / "fig1_exact_tunnels.pdf")


if __name__ == "__main__":
    main()
