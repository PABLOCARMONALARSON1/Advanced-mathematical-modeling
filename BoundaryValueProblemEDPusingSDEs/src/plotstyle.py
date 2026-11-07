"""Shared matplotlib style for the notebooks (colour-blind-safe categorical colours)."""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Categorical slots, used in this fixed order (never cycled)
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#ffffff"

# Diverging (signed fields): two hues and a neutral midpoint
DIVERGING = LinearSegmentedColormap.from_list("blue_orange", [BLUE, "#f3f3f1", ORANGE])
# Sequential (magnitudes): one hue, light -> dark
SEQUENTIAL = LinearSegmentedColormap.from_list("blues", ["#e8f0fb", BLUE, "#0f3f7a"])


def apply():
    """Set the global matplotlib defaults used by every figure of the project."""
    mpl.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": MUTED, "axes.labelcolor": MUTED, "axes.titlecolor": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "axes.axisbelow": True, "axes.titlesize": 11, "axes.titleweight": "bold",
        "axes.labelsize": 10, "legend.frameon": False, "legend.fontsize": 9,
        "lines.linewidth": 2.0, "lines.markersize": 6,
        "figure.dpi": 110, "savefig.dpi": 160, "savefig.bbox": "tight",
        "font.size": 10,
    })
