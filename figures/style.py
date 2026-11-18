"""Shared figure style: matplotlib "plasma" palette, as in the paper."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PLASMA = matplotlib.colormaps["plasma"]
NEWTON = "#1f77b4"          # Newtonian curve
INK, MUTED = "#222222", "#6b6b6b"


def apply():
    matplotlib.rcParams.update({"font.family": "serif", "font.size": 10, "axes.labelsize": 11,
                                "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                                "axes.labelcolor": INK, "text.color": INK})
    return plt
