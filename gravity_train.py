"""Core physics of the relativistic gravity train in the interior Schwarzschild metric.

Shared by the checks in ``checks/`` and the figure scripts in ``figures/``.

Paper: I. Rosignoli, "The Gravitational Train in the Interior Schwarzschild Metric:
NLO and NNLO Perturbative Expansions", submitted to European Journal of Physics (2026).

Units and notation: R = c = 1 (the paper keeps c explicit), x = r/R, compactness mu = r_s/R
with 0 < mu < 8/9, turning radius q = r_*/R, station separation Delta (radians).
A particle is released from rest at the surface and slides without friction
through a planar tunnel; the tunnel minimising the coordinate travel time is a
geodesic of the optical metric n^2 (B dr^2 + r^2 dphi^2), Eqs. (26)-(27).

Numerically it is convenient to use w(x) = 1/n(x)^2, which is finite everywhere
(n diverges at the surface, where the particle starts at rest):
    w(x) = A(x) [1 - A(x)/A(1)]                     (Eq. 26)
Non-radial Snell integrals use the substitution x = q + (1-q) sin^2 u, which removes the
integrable singularities at the turning point (x = q) and at the surface (x = 1).
"""
import contextlib
import math
import warnings

import numpy as np
from scipy.integrate import IntegrationWarning, quad
from scipy.optimize import brentq

MU_BUCHDAHL = 8.0 / 9.0
_QUAD = dict(limit=800, epsabs=1e-12, epsrel=1e-11)
_U_NODES = [0.0, math.pi / 8, math.pi / 4, 3 * math.pi / 8, math.pi / 2]


@contextlib.contextmanager
def _scanning():
    """Silence SciPy IntegrationWarning while scanning or bracketing turning radii.

    Root bracketing and branch scans evaluate the Snell integrals at extreme q (very
    close to 0, or at the edge of the forbidden band of Sec. VII B), where the integrand
    is steep and QUADPACK may warn. Those values only decide signs; every quoted result
    is checked against high-precision references in checks/. The single-q functions
    (half_angle, travel_time, ...) never silence warnings.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", IntegrationWarning)
        yield


# ---------------------------------------------------------------- metric (Eqs. 9-10)
def A(x, mu):
    """Lapse g_tt/(-c^2) of the interior Schwarzschild metric, Eq. (9)."""
    return 0.25 * (3.0 * math.sqrt(1.0 - mu) - math.sqrt(1.0 - mu * x * x)) ** 2


def B(x, mu):
    """Radial metric function g_rr, Eq. (10)."""
    return 1.0 / (1.0 - mu * x * x)


def gap(x, mu):
    """A(1) - A(x), written without cancellation: mu (1 - x^2) (5a - beta) / [4 (beta + a)]."""
    a, beta = math.sqrt(1.0 - mu), math.sqrt(1.0 - mu * x * x)
    return mu * (1.0 - x * x) * (5.0 * a - beta) / (4.0 * (beta + a))


def w(x, mu, clock="coordinate"):
    """1/n^2, the inverse squared optical index.

    clock="coordinate": Schwarzschild coordinate time, Eq. (26) (the stations' clocks read a constant multiple of it); mu = 0 gives the Newtonian 1 - x^2.
    clock="proper":     the passenger's proper time, from d tau in Eq. (22): A(1)/A(x) - 1.
    """
    if mu == 0.0:
        return 1.0 - x * x
    if clock == "proper":
        return gap(x, mu) / A(x, mu)
    return A(x, mu) * gap(x, mu) / A(1.0, mu)


def optical_index(x, mu):
    """Optical index n(x) of Eq. (26), up to the irrelevant constant factor for mu = 0."""
    return 1.0 / math.sqrt(w(x, mu))


# ---------------------------------------------------------------- Snell integrals
def _parts(u, q, mu, clock="coordinate"):
    s, c = math.sin(u), math.cos(u)
    x = q + (1.0 - q) * s * s
    wx, wq = w(x, mu, clock), w(q, mu, clock)
    rad = x * x - wx * q * q / wq          # proportional to n^2 x^2 - C^2
    jac = 2.0 * (1.0 - q) * s * c          # dx/du
    return x, wx, wq, max(rad, 0.0), jac


def half_angle(q, mu, clock="coordinate"):
    """Half opening Delta/2 of the Snell tunnel with admissible turning radius q (Eq. 34)."""
    def f(u):
        x, wx, wq, rad, jac = _parts(u, q, mu, clock)
        if rad == 0.0:
            return 0.0
        return q * math.sqrt(B(x, mu) * wx / wq) / (x * math.sqrt(rad)) * jac
    return sum(quad(f, a, b, **_QUAD)[0] for a, b in zip(_U_NODES, _U_NODES[1:]))


def opening(q, mu):
    """Total opening Delta(q) in radians."""
    return 2.0 * half_angle(q, mu)


def travel_time(q, mu):
    """Coordinate travel time t = int n dl of the Snell tunnel with admissible turning radius q, in units R/c."""
    def f(u):
        x, wx, wq, rad, jac = _parts(u, q, mu)
        if rad == 0.0 or wx == 0.0:
            return 0.0
        return math.sqrt(B(x, mu) / wx) * x / math.sqrt(rad) * jac
    return 2.0 * sum(quad(f, a, b, **_QUAD)[0] for a, b in zip(_U_NODES, _U_NODES[1:]))


def diameter_time(mu):
    """Coordinate travel time along the diameter (antipodal stations), in units R/c."""
    def f(s):                              # x = 1 - s^2 removes the surface singularity
        x = 1.0 - s * s
        return math.sqrt(B(x, mu) / w(x, mu)) * 2.0 * s if s > 0 else 0.0
    return 2.0 * quad(f, 0.0, 1.0, **_QUAD)[0]


# ---------------------------------------------------------------- perturbative result (Eqs. 57-59)
def q_series(delta, mu):
    """Newtonian, NLO and NNLO turning radius for separation delta (radians), Eq. (59)."""
    q0 = 1.0 - delta / math.pi
    q1 = q0 * (1 - q0**2) * (15 - 7 * q0**2) / 32
    q2 = q0 * (1 - q0**2) * (1827 - 3489 * q0**2 + 2457 * q0**4 - 539 * q0**6) / 4096
    return q0, q0 + mu * q1, q0 + mu * q1 + mu**2 * q2


def find_q(delta, mu, clock="coordinate"):
    """Turning radius for a monotonic opening curve (mu <~ 0.75), searched in 0.0001 <= q <= 0.999.

    For larger mu the Snell law can have several solutions; use fastest_tunnel().
    """
    with _scanning():
        return brentq(lambda q: half_angle(q, mu, clock) - delta / 2.0, 1e-4, 0.999, xtol=1e-13)


# ---------------------------------------------------------------- high compactness (Sec. VII)
def admissible(qs, mu, n=20001):
    """Boolean mask: turning radius q reachable from the surface, Eq. (69): n(x) x > n(q) q for x > q."""
    xs = np.linspace(1e-4, 1 - 1e-6, n)
    G = np.array([x * x / w(x, mu) for x in xs])
    run_min = np.minimum.accumulate(G[::-1])[::-1]
    out = []
    for q in qs:
        i = np.searchsorted(xs, q, side="right")
        out.append(i < n and q * q / w(q, mu) < run_min[i] * (1 - 1e-12))
    return np.array(out)


def fastest_tunnel(delta, mu, qs=None):
    """Snell tunnels joining two stations delta apart found by a finite grid scan, sorted by travel time.

    A tunnel may wind: its total opening can be delta + 2 pi k or 2 pi - delta + 2 pi k; the scan tries k = 0, 1, 2, 3.
    For delta = pi the radial diameter is included as well. Returns a list of dicts
    with keys q, total_opening_deg, time (and note for the diameter).
    """
    with _scanning():
        return _fastest_tunnel(delta, mu, qs)


def _fastest_tunnel(delta, mu, qs):
    if qs is None:
        qs = np.concatenate([np.geomspace(1e-4, 0.05, 40), np.linspace(0.05, 0.999, 500)[1:]])
    mask = admissible(qs, mu)
    Dq = np.array([opening(q, mu) if m else np.nan for q, m in zip(qs, mask)])
    targets = sorted({round(v, 12) for k in range(4)
                      for v in (delta + 2 * math.pi * k, 2 * math.pi - delta + 2 * math.pi * k)})
    sols = []
    for T in targets:
        for i in range(len(qs) - 1):
            if mask[i] and mask[i + 1] and (Dq[i] - T) * (Dq[i + 1] - T) < 0:
                q = brentq(lambda s: opening(s, mu) - T, qs[i], qs[i + 1], xtol=1e-13)
                sols.append(dict(q=q, total_opening_deg=math.degrees(T), time=travel_time(q, mu)))
    if abs(delta - math.pi) < 1e-12:
        sols.append(dict(q=0.0, total_opening_deg=180.0, time=diameter_time(mu), note="diameter"))
    return sorted(sols, key=lambda s: s["time"])


def mu_c():
    """Onset of circular optical orbits, Eq. (70): mu_c = 1 - z_c, 3645 z^3 + 3888 z^2 - 1080 z + 64 = 0."""
    roots = np.roots([3645, 3888, -1080, 64])
    z = [r.real for r in roots if abs(r.imag) < 1e-12 and 1 / 9 < r.real < 1]
    assert len(z) == 1
    return 1.0 - z[0]


def K(mu):
    """Coefficient K(mu) of Delta(q) = pi + 2 K q + O(q^3) for near-radial tunnels, Eq. (72)."""
    W0 = A(0.0, mu) * gap(0.0, mu)

    def f(theta):                          # x = sin(theta) removes the sqrt(1 - x) at the surface
        x = math.sin(theta)
        if x == 0.0:
            return 0.0
        h = math.sqrt(B(x, mu) * A(x, mu) * gap(x, mu) / W0)
        return (h - 1.0) / (x * x) * math.cos(theta)
    return quad(f, 0.0, math.pi / 2, limit=400)[0] - 1.0


def mu_s():
    """Loss of local minimality of the diameter (conjugate point at the antipode): K(mu_s) = 0."""
    return brentq(K, 0.7, 0.8, xtol=1e-12)


def tunnel_curve(q, mu, n=400):
    """Exact tunnel shape for admissible q: returns (X, Y) of the whole tunnel, stations symmetric about the Y axis."""
    us = np.linspace(0.0, math.pi / 2, n)

    def f(u):
        x, wx, wq, rad, jac = _parts(u, q, mu)
        if rad == 0.0:
            return 0.0
        return q * math.sqrt(B(x, mu) * wx / wq) / (x * math.sqrt(rad)) * jac
    phi = np.concatenate([[0.0], np.cumsum([quad(f, a, b, limit=200)[0] for a, b in zip(us[:-1], us[1:])])])
    x = q + (1.0 - q) * np.sin(us) ** 2
    r = np.concatenate([x[::-1], x[1:]])
    p = np.concatenate([-phi[::-1], phi[1:]])
    return r * np.sin(p), r * np.cos(p)
