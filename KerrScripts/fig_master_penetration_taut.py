# -*- coding: utf-8 -*-
"""
Master plot: classificazione del comportamento all'ergosfera delle brachistocrone
tau e t in funzione del momento J, per la forma GENERALE di phi(r).

Ramo tau (TRICOTOMIA):  penetra solo J=+Jc (misura nulla); |J|>Jc rimbalzo
  liscio (fuori r_e); |J|<Jc e J=-Jc rimbalzo con CUSPIDE (a r_e).
Ramo t (DICOTOMIA):  penetra l'intervallo (Jc-, Jc+); fuori rimbalzo liscio.

Forme generali (BL, senza shift Doran):
  dphi/dr|_tau = J r sqrt(wf)/(Δ sqrt(Δ - J^2 w))
  dphi/dr|_t   = K(r) r sqrt(wf)/(Δ sqrt(Δ - K^2 w)),  K(r)=(fJ+2Ma/r)/E
Muro sqrt(wf): f<0 dentro r_e => tau confinato a r>=r_e (t no: auto-tuning
K(r_e)=a/E lo porta a criticita' per ogni J).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from paper_style import COL, set_style, savefig
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

set_style()
HERE = os.path.dirname(os.path.abspath(__file__))
M, a, E = 1.0, 0.9, 1.2
re = 2*M
Jc = a/E
f = lambda r: 1 - 2*M/r
w = lambda r: E**2 - f(r)
Dl = lambda r: r**2 - 2*M*r + a**2
Kt = lambda r, J: (f(r)*J + 2*M*a/r)/E
rr = np.linspace(re + 1e-6, 20, 500000)

# soglie t (bordi finestra penetrante)
pen_t = np.array([np.sum((Dl(rr) - Kt(rr, J)**2*w(rr)) < -1e-9) == 0
                  for J in np.linspace(-12, 6, 36001)])
Jscan = np.linspace(-12, 6, 36001)
edges = Jscan[np.where(np.diff(pen_t.astype(int)))[0]]
Jcm, Jcp = edges.min(), edges.max()
print(f"tau: Jc=+{Jc:.3f} (unico penetrante)")
print(f"t:   finestra ({Jcm:.3f}, {Jcp:.3f})")

# ---------------------------------------------------------------- figura
fig, ax = plt.subplots(figsize=(COL, 2.9))
Jmin, Jmax = -10, 6

# colori
c_pen = '#2ca02c'      # penetra
c_smooth = '#1f77b4'   # rimbalzo liscio
c_cusp = '#d62728'     # rimbalzo cuspide

# --- barra ramo t (y=1) ---
yt = 1.0
ax.add_patch(plt.Rectangle((Jmin, yt-0.18), Jcm-Jmin, 0.36, color=c_smooth))
ax.add_patch(plt.Rectangle((Jcm, yt-0.18), Jcp-Jcm, 0.36, color=c_pen))
ax.add_patch(plt.Rectangle((Jcp, yt-0.18), Jmax-Jcp, 0.36, color=c_smooth))
ax.plot([Jcm, Jcm], [yt-0.18, yt+0.18], 'k', lw=1.0)
ax.plot([Jcp, Jcp], [yt-0.18, yt+0.18], 'k', lw=1.0)
ax.text(Jcm, yt+0.28, r'$J_c^-=%.2f$' % Jcm, ha='center', fontsize=6)
ax.text(Jcp, yt+0.28, r'$J_c^+=%.2f$' % Jcp, ha='center', fontsize=6)
ax.text(Jmin-0.3, yt, r'$t$', ha='right', va='center', fontsize=9)

# --- barra ramo tau (y=0) ---
ytau = 0.0
ax.add_patch(plt.Rectangle((Jmin, ytau-0.18), (-Jc)-Jmin, 0.36, color=c_smooth))
ax.add_patch(plt.Rectangle((-Jc, ytau-0.18), 2*Jc, 0.36, color=c_cusp))
ax.add_patch(plt.Rectangle((Jc, ytau-0.18), Jmax-Jc, 0.36, color=c_smooth))
# singolo punto penetrante +Jc
ax.plot([Jc], [ytau], marker='*', color=c_pen, ms=13, mec='k', mew=0.5, zorder=5)
ax.plot([-Jc, -Jc], [ytau-0.18, ytau+0.18], 'k', lw=1.0)
ax.plot([Jc, Jc], [ytau-0.18, ytau+0.18], 'k', lw=1.0)
ax.text(-Jc, ytau-0.32, r'$-J_c$', ha='center', fontsize=6)
ax.text(Jc, ytau-0.32, r'$+J_c=%.2f$' % Jc, ha='center', fontsize=6)
ax.text(Jmin-0.3, ytau, r'$\tau$', ha='right', va='center', fontsize=9)

ax.set_xlim(Jmin-1.2, Jmax+0.5)
ax.set_ylim(-0.7, 1.7)
ax.set_yticks([])
ax.set_xlabel(r'angular momentum $J$')
ax.set_title('ergosphere behaviour vs $J$: $\\tau$ trichotomy, $t$ dichotomy')

legend = [Patch(color=c_pen, label='penetrates $r_e$'),
          Patch(color=c_smooth, label='smooth bounce (outside $r_e$)'),
          Patch(color=c_cusp, label='cusp bounce (at $r_e$)'),
          plt.Line2D([], [], marker='*', color=c_pen, mec='k', ms=9, ls='',
                     label=r'$\tau$: single penetrating $J=+J_c$')]
ax.legend(handles=legend, fontsize=5.5, loc='upper left',
          bbox_to_anchor=(0.0, -0.28), ncol=2, framealpha=0.9)
savefig(fig, HERE, 'fig_master_penetration_taut')
print('FATTO.')
