"""Shared matplotlib style: restrained, print-ready, colour-blind-safe (Okabe-Ito)."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OI = {"black": "#000000", "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73", "yellow": "#F0E442",
      "blue": "#0072B2", "vermillion": "#D55E00", "purple": "#CC79A7"}
CYCLE = [OI["blue"], OI["vermillion"], OI["green"], OI["purple"], OI["orange"], OI["sky"]]


def apply_style():
    plt.rcParams.update({
        "font.family": "serif", "mathtext.fontset": "stix", "font.size": 9.5,
        "axes.labelsize": 10, "axes.titlesize": 10, "legend.fontsize": 8, "legend.frameon": False,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "axes.prop_cycle": matplotlib.cycler(color=CYCLE), "axes.linewidth": 0.8,
        "lines.linewidth": 1.5, "grid.color": "#CCCCCC", "grid.linewidth": 0.5, "axes.grid": True,
        "grid.alpha": 0.6, "figure.dpi": 120, "savefig.dpi": 220, "savefig.bbox": "tight",
    })


def panel_label(ax, s):
    ax.text(-0.01, 1.03, s, transform=ax.transAxes, va="bottom", ha="left", fontweight="bold", fontsize=10)


def save(fig, outdir, name):
    import os
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, name + ".png"))
    fig.savefig(os.path.join(outdir, name + ".pdf"), metadata={"CreationDate": None})
    plt.close(fig)
