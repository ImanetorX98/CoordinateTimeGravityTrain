#!/usr/bin/env python3
"""Fig. 1 of the paper: exact versus NNLO tunnels for stations 90 deg apart.

Purpose  Draw the Newtonian tunnel, the exact relativistic tunnels (solid) and the
         NNLO tunnels from the closed-form angular functions (dashed).
Paper    Fig. 1; Eqs. (51)-(53) for the NNLO shape, Eq. (59) for q.
Inputs   optional: --mu 0.1 0.3 0.5 0.7 (compactness values), --delta 90 (degrees).
Outputs  figures/output/fig1_exact_vs_nnlo.pdf and .png
Runtime  a few seconds.
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


def phi_nnlo(x, q, mu):
    """phi0 + mu phi1 + mu^2 phi2, Eqs. (51)-(53)."""
    U = np.arcsin(np.sqrt(np.clip((x**2 - q**2) / (1 - q**2), 0, 1)))
    S = np.sqrt(np.clip((x**2 - q**2) * (1 - x**2), 0, None))
    P, Q = -144 * q**4 + 12 * q**2 + 132, 147 * q**4 - 354 * q**2 + 399
    p0 = np.arctan2(np.sqrt(np.clip(x**2 - q**2, 0, None)), q * np.sqrt(np.clip(1 - x**2, 0, None))) - q * U
    p1 = q * (15 - 7 * q**2) / 32 * (S + (1 - q**2) * U)
    p2 = q / 512 * (P / 2 * (S + (1 - q**2) * U)
                    + Q / 8 * ((2 * x**2 + 3 * q**2 - 1) * S + (1 - q**2) * (1 + 3 * q**2) * U))
    return p0 + mu * p1 + mu**2 * p2


def whole(x, phi):
    r = np.concatenate([x[::-1], x[1:]])
    p = np.concatenate([-phi[::-1], phi[1:]])
    return r * np.sin(p), r * np.cos(p)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mu", type=float, nargs="+", default=[0.1, 0.3, 0.5, 0.7])
    ap.add_argument("--delta", type=float, default=90.0, help="station separation in degrees")
    args = ap.parse_args()
    delta = math.radians(args.delta)
    plt = style.apply()
    fig, ax = plt.subplots(figsize=(3.4, 2.9))
    th = np.linspace(-math.pi / 2, math.pi / 2, 400)
    ax.plot(np.sin(th), np.cos(th), color=style.INK, lw=1.2)
    u = np.linspace(0, math.pi / 2, 400)
    q0 = 1 - delta / math.pi
    x0 = q0 + (1 - q0) * np.sin(u) ** 2
    ax.plot(*whole(x0, phi_nnlo(x0, q0, 0.0)), color=style.NEWTON, lw=2, label=r"Newtonian, $\mu=0$")
    cols = np.linspace(0.2, 0.82, len(args.mu))
    for mu, c in zip(args.mu, cols):
        col = style.PLASMA(c)
        ax.plot(*g.tunnel_curve(g.find_q(delta, mu), mu), color=col, lw=2, label=rf"$\mu={mu:g}$")
        qn = g.q_series(delta, mu)[2]
        xn = qn + (1 - qn) * np.sin(u) ** 2
        ax.plot(*whole(xn, phi_nnlo(xn, qn, mu)), color=col, lw=1, ls=(0, (3, 2)))
    for s in (1, -1):
        ax.plot(s * math.sin(delta / 2), math.cos(delta / 2), "o", color=style.INK, ms=4)
    ax.plot([], [], color=style.INK, lw=1, ls=(0, (3, 2)), label="NNLO")
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x/R$")
    ax.set_ylabel(r"$y/R$")
    ax.grid(True, color="#e6e6e6", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.42), ncol=3, fontsize=7.5, frameon=False)
    out = ROOT / "figures" / "output"
    out.mkdir(exist_ok=True)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(out / f"fig1_exact_vs_nnlo.{ext}", dpi=300, bbox_inches="tight")
    print("saved", out / "fig1_exact_vs_nnlo.pdf")


if __name__ == "__main__":
    main()
