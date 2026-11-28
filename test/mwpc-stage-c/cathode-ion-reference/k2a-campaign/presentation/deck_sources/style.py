"""Shared plot style: ATLAS-like axes, Public Sans, deck palette. S0 = black, K2A = red everywhere."""
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "/home/claude/wt-val/test/mwpc-stage-c/cathode-ion-reference/timing-ablation")
import t_shaped_compare as TS  # noqa: E402  (the validation tool's own shaper, unchanged)

INK, INK2, MUTED, GRID = "#15191F", "#4A525C", "#8A9099", "#C9CED4"
S0C, K2C = "#15191F", "#C8102E"          # S0 reference, K2A prediction
UC, VC = "#1F4E9E", "#D9822B"            # U plane, V plane
GREEN = "#2E7D38"
FIGDIR = "/home/claude/pres/fig"

from matplotlib import font_manager as _fm
for _f in ("/root/.fonts/PublicSans-Regular.ttf", "/root/.fonts/PublicSans-Bold.ttf"):
    _fm.fontManager.addfont(_f)
import logging as _lg
_lg.getLogger("matplotlib.font_manager").setLevel(_lg.ERROR)

mpl.rcParams.update({
    "font.family": ["Public Sans", "DejaVu Sans"], "font.size": 11, "axes.labelsize": 12, "axes.titlesize": 12,
    "mathtext.fontset": "dejavusans", "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
    "ytick.right": True, "xtick.minor.visible": True, "ytick.minor.visible": True, "xtick.major.size": 6,
    "ytick.major.size": 6, "xtick.minor.size": 3, "ytick.minor.size": 3, "axes.linewidth": 1.0,
    "legend.frameon": False, "legend.fontsize": 9.5, "figure.dpi": 100, "savefig.dpi": 300,
    "axes.unicode_minus": True, "figure.facecolor": "white", "axes.facecolor": "white",
    "lines.linewidth": 1.6,
})


def label(ax, x=0.04, y=0.95, extra=None, ha="left", size=11):
    """'NA60+ Simulation' label in the ATLAS manner (bold italic experiment name)."""
    if ha == "right":
        t2 = ax.text(x, y, " Simulation", transform=ax.transAxes, ha="right", va="top", fontsize=size)
        t = ax.annotate("NA60+", xy=(0, 0), xycoords=t2, ha="right", va="bottom", fontsize=size, fontweight="bold", fontstyle="italic", family="DejaVu Sans")
    else:
        t = ax.text(x, y, "NA60+", transform=ax.transAxes, ha="left", va="top", fontsize=size, fontweight="bold", fontstyle="italic", family="DejaVu Sans")
        ax.annotate(" Simulation", xy=(1, 0), xycoords=t, ha="left", va="bottom", fontsize=size)
    if extra:
        for i, e in enumerate(extra if isinstance(extra, list) else [extra]):
            ax.text(x, y - 0.075 * (i + 1), e, transform=ax.transAxes, ha=ha, va="top", fontsize=size - 1.5, color=INK)
    return t


def ms(x, unit="", f="{:.3g}"):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return f"mean {f.format(x.mean())}{unit}, RMS {f.format(x.std(ddof=0))}{unit}"


def save(fig, name):
    fig.savefig(f"{FIGDIR}/{name}.png", bbox_inches="tight", pad_inches=0.04, facecolor="white")
    plt.close(fig)


def load_variant(d, prompt=True):
    lab, t, cur = TS.load_wave(f"{d}/replay/waveforms.csv")
    if prompt:
        l2, t2, pc = TS.load_wave(f"{d}/prompt/prompt_waveforms.csv")
        assert l2 == lab and np.array_equal(t, t2)
        cur = cur + pc
    return lab, t, cur


def shaped(cur, dt_in, tp, window_ns=6000.0, dt=1.0):
    h, _ = TS.impulse(tp, 3, dt, max(12 * tp, 2000.0))
    nb = int(round(window_ns / dt_in))
    return TS.shape(cur[:nb], dt_in, h, dt)
