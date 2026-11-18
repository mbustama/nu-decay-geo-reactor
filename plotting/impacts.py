# The paper's plotting code, its functions that read the analysis's results replaced by ones that read
# the release's data/ files (at the end of this file).
"""Figures E.1 and E.2 of the paper: how much each parameter moves the 90% limit from data on tau_1/m_1.

A row per parameter, its bars the relative change of the limit when the parameter is held at its
best fit plus or minus its fitted width, every other parameter fitted.  Rows are ranked by the
larger of the two changes; those under :data:`SHOWN_ABOVE` are left out, and :func:`shown` counts
them.  A side that would leave the parameter's range has no bar.  The figure of all four
experiments combined has a second panel, on the left, with each parameter's prior and its fit to
data, as shifts from the prior's center in widths of the prior.

The notebook draws the figures with :func:`figure_combined_all` and :func:`figure_experiments`.
"""

from __future__ import annotations

import glob
import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon, Rectangle
from matplotlib.ticker import MultipleLocator


from . import style as figure_style
from . import names as pn
from .analysis import ss

#: The free-power fit's own parameter, the Earth's radiogenic power, which ``names`` does not name.
pn.NAMES.setdefault("bse_power_tw", "Radiogenic power")

#: The records drawn: the study's, with the mixing the paper is drawn from (``dl.PAPER_MIXING``).
MIXING = ss.dl.PAPER_MIXING
ss.use_mixing(MIXING)

#: Each fit's title: the experiments with the years of their data, and the combinations by every name.
TITLES = {"borexino": "Borexino (2020)", "kamland": "KamLAND (2022)", "snoplus": "SNO+ (2026)",
          "juno": "JUNO (2026)", "combined_bx_kl_snoplus": "Borexino, KamLAND, and SNO+ combined",
          "combined_all": "Borexino, KamLAND, SNO+, and JUNO combined"}
#: The two sides, each with its legend entry and color (Okabe and Ito: vermilion and sky blue, which
#: no experiment's color is).
SIDES = (("+1", r"Best fit $+1\sigma$", "#D55E00"), ("-1", r"Best fit $-1\sigma$", "#56B4E9"))
#: Rows whose larger change is below this share are left out; the captions count them.
SHOWN_ABOVE = 0.01
#: Bar height, in rows.
BAR = 0.7
#: Where all four's figure cuts its impact axis, in percent: a longer bar stops there in a triangle with
#: its value written in it, as the money plot's bounds beyond its axis.
CUT_ALL = 100.0
#: The triangle's length, as a share of the axis's span, and the gap between its tip and the edge.
TRIANGLE, TIP_GAP = 0.025, 0.004
#: Major ticks, with a vertical grid line each, every this many percent.
TICK = 10.0
#: The most tick labels an axis prints, on a full-width panel and on a narrow one of the appendix;
#: beyond, only every second, fifth or tenth tick is labelled.
MOST_LABELS, MOST_LABELS_NARROW = 12, 6
#: The legend's font size.
LEGEND_SIZE = 25
#: The pull panel's colors: each prior's one-sigma range in a light tint (45%) of the green of the Brazilian plot's 68% band,
#: the range a parameter is expected in, and each fit to data in black, as the fits to data are there.
PRIOR_COLOR, FIT_COLOR = "#8CD7B1", "black"
#: The prior bars' height, in rows, and the pull axis's reach either side of the prior's center, in widths.
PRIOR_BAR, PULL_REACH = 0.45, 2.6
#: The empty rows under the last parameter that hold the pull panel's legend, inside the panel.
LEGEND_ROWS = 1.7
#: What the pull panel says in the row of a parameter whose prior is flat.
FLAT = "Flat prior"


def shown(site: str) -> tuple[list[tuple[str, dict[str, float]]], int]:
    """The rows a fit's figure shows, and how many it leaves out.

    Args:
        site: An experiment or a combination, such as ``"combined_all"``.

    Returns:
        ``(rows, left out)``: the parameters whose larger change reaches
        :data:`SHOWN_ABOVE`, largest first, and the number of the others.
    """
    every = impacts(site)
    rows = [(n, d) for n, d in every if max(abs(v) for v in d.values()) >= SHOWN_ABOVE]
    return rows, len(every) - len(rows)


def draw(ax, site: str, label_size: float = 21.0, title_size: float = 23.0, most_labels: int = MOST_LABELS,
         right_edge: float | None = None) -> int:
    """One fit's bars on an axes.

    Where both sides move the limit the same way, the longer bar is drawn first,
    so that the shorter one shows in front of it.

    Args:
        ax: The axes.
        site: An experiment or a combination, such as ``"combined_all"``.
        label_size: The parameters' names' font size.
        title_size: The title's font size.
        most_labels: The most tick labels the axis prints (:func:`grid`).
        right_edge: Where to cut the axis, in percent, or None to fit every bar.  A bar beyond it
            stops in a triangle just inside the edge, its value written in it (:data:`CUT_ALL`).

    Returns:
        How many parameters are left out.
    """
    rows, left_out = shown(site)
    beyond = []
    for y, (_name, d) in enumerate(rows):
        bars = [(100.0 * d[side], color) for side, _text, color in SIDES if side in d]
        for value, color in sorted(bars, key=lambda b: -abs(b[0])):
            if right_edge is not None and value > right_edge:
                beyond.append((y, value, color))
            else:
                ax.barh(y, value, height=BAR, color=color, lw=0, alpha=0.9)
    if right_edge is not None:
        low = min(ax.get_xlim()[0], 0.0)
        ax.set_xlim(low, right_edge)
        span = right_edge - low
        tip = right_edge - TIP_GAP * span
        base = tip - TRIANGLE * span
        for y, value, color in beyond:
            # Drawn behind any shorter bar of the same row, as the longer bar always is.
            ax.add_patch(Rectangle((0.0, y - BAR / 2), base, BAR, color=color, lw=0, alpha=0.9, zorder=0.9))
            ax.add_patch(Polygon([(base, y - BAR / 2), (base, y + BAR / 2), (tip, y)], closed=True, color=color,
                                 lw=0, alpha=0.9, zorder=0.9))
            ax.text(base - 0.012 * span, y, rf"${value:+.0f}\%$", ha="right", va="center_baseline",
                    fontsize=label_size, color="black", zorder=3)
    ax.axvline(0.0, color="0.25", lw=1.0)
    grid(ax, most_labels)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([pn.label(n, site in ss.dl.COMBINATIONS) for n, _d in rows], fontsize=label_size)
    ax.tick_params(axis="y", which="both", left=False, right=False)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_title(TITLES[site], fontsize=title_size, pad=12)
    return left_out


def draw_pulls(ax, site: str, label_size: float = 23.0) -> None:
    """One fit's pull panel on an axes, its rows those of :func:`draw`, with its legend inside it.

    Each row has its prior's one-sigma range as a green bar from -1 to 1 and
    the fit to data as a black dot and line; a flat prior's row says so.  The
    legend goes in :data:`LEGEND_ROWS` of empty rows under the last parameter,
    which the axes' limits leave.

    Args:
        ax: The axes, sharing its rows with the impact panel's.
        site: An experiment or a combination, such as ``"combined_all"``.
        label_size: The parameters' names' font size, which a flat prior's row also takes.
    """
    for y, pull in enumerate(pulls(site)):
        if pull is None:
            ax.text(0.0, y, FLAT, ha="center", va="center", fontsize=label_size, color="0.3")
            continue
        best, low, high = pull
        ax.barh(y, 2.0, left=-1.0, height=PRIOR_BAR, color=PRIOR_COLOR, lw=0, zorder=2)
        ax.errorbar(best, y, xerr=[[best - low], [high - best]], fmt="o", color=FIT_COLOR, ms=8, lw=2.4,
                    capsize=0, zorder=3)
    ax.axvline(0.0, color="0.6", lw=0.9, zorder=1)
    ax.set_xlim(-PULL_REACH, PULL_REACH)
    ax.set_xticks(np.arange(-2, 3))
    ax.set_xticklabels([f"${t}$" for t in range(-2, 3)])
    ax.set_xlabel("Shift from the prior's center [prior widths]")
    ax.tick_params(axis="y", which="both", left=False, right=False)
    handles = [Patch(color=PRIOR_COLOR, lw=0, label=r"Prior, $\pm1\sigma$"),
               Line2D([], [], color=FIT_COLOR, marker="o", ms=8, lw=2.4, label=r"Fit to data, $\pm1\sigma$")]
    box = ax.legend(handles=handles, loc="lower center", ncol=2, fontsize=LEGEND_SIZE, framealpha=1.0,
                    fancybox=False, handlelength=1.6, columnspacing=1.2)
    box.get_frame().set_edgecolor("black")
    box.get_frame().set_linewidth(1.0)


def grid(ax, most_labels: int = MOST_LABELS) -> None:
    """Major ticks every :data:`TICK` percent, each with a vertical grid line behind the bars.

    Where that would print more than ``most_labels`` numbers (JUNO's axis
    spans almost 300), only every second, fifth or tenth tick is labelled; every
    tick keeps its line.

    Args:
        ax: The axes, its bars drawn.
        most_labels: The most numbers the axis prints.
    """
    ax.xaxis.set_major_locator(MultipleLocator(TICK))
    ax.set_axisbelow(True)
    ax.grid(axis="x", which="major", color="0.85", lw=0.9)
    low, high = ax.get_xlim()
    ticks = np.arange(np.ceil(low / TICK) * TICK, high + 1e-9, TICK)
    step = next(k for k in (1, 2, 5, 10, 20) if len(ticks) / k <= most_labels)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"${t:.0f}$" if round(t / TICK) % step == 0 else "" for t in ticks])


def legend(ax, loc: str = "lower right", anchor=None) -> None:
    """The two sides' legend, in a box with a black border.

    Args:
        ax: The axes it goes in, or the figure, with ``anchor`` in figure coordinates.
        loc: Where in them, or which point of the box sits at ``anchor``.
        anchor: Where the box goes, for a legend outside every panel; None inside ``ax``.
    """
    box = ax.legend(handles=[Patch(color=c, alpha=0.9, lw=0, label=t) for _s, t, c in SIDES], loc=loc,
                    bbox_to_anchor=anchor, title="Parameter held at", fontsize=LEGEND_SIZE,
                    title_fontsize=LEGEND_SIZE, framealpha=1.0, fancybox=False)
    box.get_frame().set_edgecolor("black")
    box.get_frame().set_linewidth(1.0)


def xlabel(ax) -> None:
    """The axis label every panel shares.

    Args:
        ax: The axes.
    """
    ax.set_xlabel(r"Relative change of the $90\%$ C.L. limit [\%]")


def figure_combined_all():
    """All four's figure, for a full page of the paper: the pulls on the left, the impacts on the right.

    Returns:
        The matplotlib figure.
    """
    figure_style.use()
    site = "combined_all"
    fig, (left, right) = plt.subplots(1, 2, figsize=(22.0, 26.0), sharey=True,
                                      gridspec_kw={"width_ratios": [0.9, 1.25], "wspace": 0.06})
    draw(right, site, label_size=23.0, right_edge=CUT_ALL)
    right.set_title("")
    right.tick_params(labelleft=False)
    xlabel(right)
    legend(right)
    rows = shown(site)[0]
    left.set_yticks(range(len(rows)))
    left.set_yticklabels([pn.label(n, True) for n, _d in rows], fontsize=23.0)
    draw_pulls(left, site)
    left.set_ylim(len(rows) - 0.4 + LEGEND_ROWS, -0.6)
    # One title over both panels.
    a, b = left.get_position(), right.get_position()
    fig.text(0.5 * (a.x0 + b.x1), a.y1 + 0.006, TITLES[site], ha="center", va="bottom", fontsize=27.0)
    return fig


def figure_experiments():
    """The four experiments' figure, for the appendix: Borexino and KamLAND on top, SNO+ and JUNO below.

    Every row has the same height in all four panels, so their bars compare: a panel with fewer rows than
    its neighbour leaves the rest empty below, and Borexino's holds the legend there.  JUNO's axis stops
    at :data:`CUT_ALL`, as all four's does: its rate of geoneutrinos and sin2 th12 move its limit beyond
    +100%, which would squeeze every other bar.  BX+KL+SNO+ has no panel: its impacts are those of
    all four.

    Returns:
        The matplotlib figure.
    """
    figure_style.use()
    n = {site: len(shown(site)[0]) for site in ("borexino", "kamland", "snoplus", "juno")}
    top, bottom = max(n["borexino"], n["kamland"]), max(n["snoplus"], n["juno"])
    fig = plt.figure(figsize=(20.0, 21.0))
    grid = fig.add_gridspec(2, 2, height_ratios=[top + 0.2, bottom + 0.2], wspace=1.05, hspace=0.16)
    axes = {}
    for (row, col), site in (((0, 0), "borexino"), ((0, 1), "kamland"), ((1, 0), "snoplus"), ((1, 1), "juno")):
        axes[site] = fig.add_subplot(grid[row, col])
        draw(axes[site], site, label_size=19.0, most_labels=MOST_LABELS_NARROW,
             right_edge=CUT_ALL if site == "juno" else None)
        rows = top if row == 0 else bottom
        axes[site].set_ylim(rows - 0.4, -0.6)
    legend(axes["borexino"], loc="lower center")
    # One label for the four panels.
    fig.supxlabel(r"Relative change of the $90\%$ C.L. limit [\%]", fontsize=25, y=0.072)
    return fig


# The data, read from the release's data/ files.
from . import data as _data  # noqa: E402


def impacts(*args, **kwargs):
    """The paper's ``impacts``, read from the data file of the figure in use."""
    return getattr(_data.current(), 'impacts')(*args, **kwargs)


def pulls(*args, **kwargs):
    """The paper's ``pulls``, read from the data file of the figure in use."""
    return getattr(_data.current(), 'pulls')(*args, **kwargs)
