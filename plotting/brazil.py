# The paper's plotting code, its functions that read the analysis's results replaced by ones that read
# the release's data/ files (at the end of this file).
"""Figure C.1 of the paper: the limits from data against the bands of simulated limits.

For the decay nu_1 -> nu_3 in the inverted ordering, per experiment and combination: the limits of
simulated experiments without decay on the geochemical Earth (the bands), the noise-free expected
limit, and the limit from data.  One panel has the Earth's radiogenic power assumed, the other has
it free.  The green band holds the middle 68% of the simulated limits and the yellow one the middle
95%; the dashed line is their median.

The notebook draws the figure with :func:`figure_paper`.
"""

from __future__ import annotations

import json
import math
import os
import sys

import matplotlib.pyplot as plt


from .analysis import dl
from . import style as figure_style
from .analysis import BSE_PROFILE_RANGE_TW

SITES = (("borexino", "Borexino"), ("kamland", "KamLAND"), ("snoplus", "SNO+"), ("juno", "JUNO"))
#: The combined fits, which share the lifetime and the mixing across experiments; their columns
#: stay empty until their records exist.
COMBINED = (("combined_bx_kl_snoplus", "BX+KL+SNO+"), ("combined_all", "All four"))
COLUMNS = SITES + COMBINED
#: The channels with a Brazilian plot: the case, the lifetime limited and the label. Row 1 first,
#: then the six where some fit's data set a 90% limit.
#: Every fit has the mixing the paper is drawn from and the Earth's power assumed.
CHANNELS = (("io_n1_n3", "log10_tau1_over_m1", r"$\nu_1 \to \nu_3$, IO"),
            ("io_n2_n3", "log10_tau2_over_m2", r"$\nu_2 \to \nu_3$, IO"),
            ("io_n2_n1", "log10_tau2_over_m2", r"$\nu_2 \to \nu_1$, IO"),
            ("io_n2_n1_n3", "log10_tau1_over_m1", r"$\nu_2 \to \nu_1 \to \nu_3$, IO: $\tau_1/m_1$"),
            ("io_n2_n1_n3", "log10_tau2_over_m2", r"$\nu_2 \to \nu_1 \to \nu_3$, IO: $\tau_2/m_2$"),
            ("no_n2_n1", "log10_tau2_over_m2", r"$\nu_2 \to \nu_1$, NO"),
            ("no_n3_n2_n1", "log10_tau2_over_m2", r"$\nu_3 \to \nu_2 \to \nu_1$, NO: $\tau_2/m_2$"))


#: The colors of the bands, as in every Brazilian plot.
YELLOW, GREEN = "#FFDF00", "#00A651"
#: The fits to data: per Earth model, the marker, its fill, and its offset from the band's center.  Every
#: Earth model is open: a filled marker is the Earth's power free, in the overview.
DATA_MARKERS = (("cosmochemical", "Cosmochemical", "v", "white", -0.12),
                ("geochemical", "Geochemical", "o", "white", 0.0),
                ("geodynamical", "Geodynamical", "^", "white", 0.12))
#: The fit to data with the Earth's power free, started from the geochemical Earth: the overview's
#: filled diamond, alone.
FREE_MARKERS = (("geochemical", r"Power free, ${:g}$ to ${:g}$ TW".format(*BSE_PROFILE_RANGE_TW), "D", "k", 0.0),)
#: Where the simulations' legend goes with the power free: Borexino's band runs to the floor there, so
#: the legend goes to the top, beside the channel's label, in two columns.
FREE_SIMULATIONS_LEGEND = {"loc": "upper left", "bbox_to_anchor": (0.245, 1.0), "ncol": 2}


def markers(treatment: str) -> tuple:
    """The fits to data a plot shows.

    Args:
        treatment: ``"assumed"``, or ``"profiled"`` for the Earth's power free.

    Returns:
        :data:`DATA_MARKERS` or :data:`FREE_MARKERS`.
    """
    return DATA_MARKERS if treatment == "assumed" else FREE_MARKERS


def span(channel, treatment: str = "assumed", models=None) -> tuple[float, float]:
    """A panel's range: the whole decades around everything it draws, two at least.

    The lowest value keeps a margin of half its decade from the floor, and a
    range of one decade is opened downward to two, where matplotlib would label
    the minor ticks.

    Args:
        channel: An entry of :data:`CHANNELS`.
        treatment: ``"assumed"``, or ``"profiled"`` for the Earth's power free.
        models: The Earth models whose limits from data are drawn; every one when None.

    Returns:
        ``(floor, top)`` in s/eV; a band edge without a limit runs down to the floor.
    """
    values = []
    for site, _ in COLUMNS:
        toys = simulated(site, channel, treatment)
        if toys:
            b = dl.band(toys, "90")
            values += [v for v in (*b["95"], *b["68"], b["median"]) if v]
        shown = tuple(m[0] for m in markers(treatment) if models is None or m[0] in models)
        for data, bses in (("asimov", ("geochemical",)), ("real", shown)):
            for bse in bses:
                r = record(site, data, bse, channel, treatment)
                if r and r["limits_s_per_ev"]["90"]:
                    values.append(r["limits_s_per_ev"]["90"])
    if not values:
        return 1e-11, 1e-7
    floor = 10.0 ** math.floor(math.log10(min(values) / 1.5))
    top = 10.0 ** math.ceil(math.log10(1.5 * max(values)))
    return min(floor, top / 100.0), top


def draw(ax, channel=CHANNELS[0], legends: bool = True, treatment: str = "assumed", models=None,
         columns=None, label=None, label_size=None, ylabel: bool = True, limits=None,
         label_xy=(0.025, 0.982)) -> None:
    """Draw one channel's bands, expected limits and limits from data.

    Args:
        ax: The axes to draw on.
        channel: An entry of :data:`CHANNELS`.
        legends: Whether to draw the legends in the panel.
        treatment: ``"assumed"``, or ``"profiled"`` for the Earth's power free.
        models: The Earth models whose limits from data are marked; every one when None.
        columns: The columns' labels; :data:`COLUMNS`' names when None.
        label: The panel's label; the channel's, as :data:`CHANNELS` names it, when None.
        label_size: The label's font size; the rc's when None.
        ylabel: Whether to label the y-axis.
        limits: The panel's range ``(floor, top)``; its own (:func:`span`) when None.  A band edge
            without a limit runs down to the floor.
        label_xy: Where the panel's label's upper left corner goes, in axes fractions.
    """
    floor, top = limits or span(channel, treatment, models)
    width = 0.32
    for k, (site, _name) in enumerate(COLUMNS):
        # A fit's limits are drawn before its bands exist; the bands follow.
        toys = simulated(site, channel, treatment)
        if toys:
            b = dl.band(toys, "90")
            lo95, hi95 = (floor if v is None else v for v in b["95"])
            lo68, hi68 = (floor if v is None else v for v in b["68"])
            ax.fill_between([k - width, k + width], lo95, hi95, color=YELLOW, lw=0, zorder=1)
            ax.fill_between([k - width, k + width], lo68, hi68, color=GREEN, lw=0, zorder=2)
            if b["median"]:
                ax.plot([k - width, k + width], [b["median"]] * 2, "k--", lw=1.5, zorder=3)
            for lo in (lo95, lo68):
                if lo == floor:
                    ax.annotate("", (k, floor * 1.05), (k, floor * 4),
                                arrowprops=dict(arrowstyle="->", color="0.3", lw=1))
        elif treatment == "profiled":
            # Just above the fit to data, or below it near the top; midway up the panel without one.
            decades = math.log10(top / floor)
            fit = record(site, "real", FREE_MARKERS[0][0], channel, treatment)
            v = fit["limits_s_per_ev"]["90"] if fit else None
            if v is None:
                y, va = math.sqrt(floor * top), "center"
            elif math.log10(top / v) > 0.25 * decades:
                y, va = v * 10 ** (0.06 * decades), "bottom"
            else:
                y, va = v / 10 ** (0.06 * decades), "top"
            ax.text(k, y, "running", ha="center", va=va, fontsize=13, color="0.45", style="italic", zorder=3,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none"))
        asimov = record(site, "asimov", "geochemical", channel, treatment)
        if asimov and asimov["limits_s_per_ev"]["90"]:
            ax.plot(k + width + 0.08, asimov["limits_s_per_ev"]["90"], marker="<", ms=8, mfc="white",
                    mec="k", ls="none", zorder=4)
        for bse, _, marker, fill, offset in markers(treatment):
            if models is not None and bse not in models:
                continue
            offset = offset if models is None or len(models) > 1 else 0.0
            r = record(site, "real", bse, channel, treatment)
            if r and r["limits_s_per_ev"]["90"]:
                ax.plot(k + offset, r["limits_s_per_ev"]["90"], marker=marker, ms=9 if bse == "geochemical" else 8,
                        mfc=fill, mec="k", mew=1.4, ls="none", zorder=6 if bse == "geochemical" else 5)
    ax.set_yscale("log")
    ax.set_ylim(floor, top)
    ax.set_xlim(-0.55, len(COLUMNS) - 0.45)
    ax.set_xticks(range(len(COLUMNS)))
    ax.set_xticklabels(columns or [name for _, name in COLUMNS])
    ax.axvline(len(SITES) - 0.5, color="0.5", lw=1, ls=":", zorder=0)
    ax.tick_params(axis="x", which="both", top=False, bottom=False)
    lifetime = channel[1].replace("log10_tau", "").split("_over_m")[0]
    if ylabel:
        ax.set_ylabel(rf"$\tau_{lifetime}/m_{lifetime}$, 90\% lower limit [s~eV$^{{-1}}$]")
    ax.text(*label_xy, label or channel[2] + ("" if treatment == "assumed" else ", Earth's power free"),
            transform=ax.transAxes, ha="left", va="top", fontsize=label_size, zorder=10,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="black"))
    if legends:
        draw_legends(ax, ("lower left" if treatment == "assumed" else FREE_SIMULATIONS_LEGEND, "lower right"),
                     treatment)


def draw_legends(ax, places, treatment: str = "assumed") -> None:
    """The two legends: the simulations, and the fits to data on each Earth model or with the power free.

    Args:
        ax: The axes to draw them on.
        places: Where each goes: matplotlib's ``loc``, or a mapping of the legend's placement options.
        treatment: ``"assumed"``, or ``"profiled"`` for the Earth's power free.
    """
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    simulations = ax.legend(
        [Patch(color=YELLOW, lw=0), Patch(color=GREEN, lw=0), Line2D([], [], color="k", ls="--", lw=1.5),
         Line2D([], [], marker="<", ms=8, mfc="white", mec="k", ls="none")],
        [r"95\% band", r"68\% band", "Median", "Expected (noise-free)"],
        title="Simulations", edgecolor="black", framealpha=1, **placed(places[0]))
    ax.add_artist(simulations)
    shown = markers(treatment)
    ax.legend([Line2D([], [], marker=marker, ms=9 if bse == "geochemical" else 8, mfc=fill, mec="k", mew=1.4,
                      ls="none") for bse, _, marker, fill, _ in shown],
              [name for _, name, _, _, _ in shown],
              title="Fits to data" if len(shown) > 1 else "Fit to data", edgecolor="black", framealpha=1,
              **placed(places[1]))


def placed(place) -> dict:
    """A legend's placement options.

    Args:
        place: matplotlib's ``loc``, or a mapping of placement options.

    Returns:
        The options, as keyword arguments to ``legend``.
    """
    return dict(place) if isinstance(place, dict) else {"loc": place}


#: The paper's figure: each channel named as the other figures name a channel, a chain's panel with
#: its lifetime, and each experiment's column with the year of its analysis.
PAPER_LABELS = {("io_n1_n3", "log10_tau1_over_m1"): r"$\nu_1\to\nu_3$ (IO)",
                ("io_n2_n3", "log10_tau2_over_m2"): r"$\nu_2\to\nu_3$ (IO)",
                ("io_n2_n1", "log10_tau2_over_m2"): r"$\nu_2\to\nu_1$ (IO)",
                ("io_n2_n1_n3", "log10_tau1_over_m1"): r"$\nu_2\to\nu_1\to\nu_3$ (IO), $\tau_1/m_1$",
                ("io_n2_n1_n3", "log10_tau2_over_m2"): r"$\nu_2\to\nu_1\to\nu_3$ (IO), $\tau_2/m_2$",
                ("no_n2_n1", "log10_tau2_over_m2"): r"$\nu_2\to\nu_1$ (NO)",
                ("no_n3_n2_n1", "log10_tau2_over_m2"): r"$\nu_3\to\nu_2\to\nu_1$ (NO), $\tau_2/m_2$"}
#: The room above the top decade, in decades, that each panel's label sits in, and where the label's
#: upper left corner goes, in axes fractions.
PAPER_HEADROOM_DECADES = 0.7
PAPER_LABEL_XY = (0.04, 0.95)
PAPER_COLUMNS = ("Borexino\n(2020)", "KamLAND\n(2022)", "SNO+\n(2026)", "JUNO\n(2026)", "BX+KL+\nSNO+", "All\nfour")


#: The paper's panels: row 1 with the Earth assumed and with its power free.
PAPER_PANELS = ((CHANNELS[0], "assumed"), (CHANNELS[0], "profiled"))


def figure_paper(panels=PAPER_PANELS):
    """The Brazilian plots of the paper, two to a row, on one range.

    Only the geochemical Earth's limit from data is marked, the Earth model the bands are made on,
    and the power-free fit's in its panel.  The experiments are named under
    the bottom row only; two legends side by side across the figure name the simulated experiments
    and the limits from data; one label serves every y-axis.

    Args:
        panels: ``(channel, treatment)`` per panel, a channel of :data:`CHANNELS` and ``"assumed"``
            or ``"profiled"`` (the Earth's power free).

    Returns:
        The matplotlib figure.
    """
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    figure_style.use()
    size = plt.rcParams["axes.labelsize"]
    spans = [span(channel, treatment, ("geochemical",)) for channel, treatment in panels]
    limits = (min(lo for lo, _ in spans), max(hi for _, hi in spans))
    rows = (len(panels) + 1) // 2
    # Each row tall enough for the one y-axis label, which a single row must hold alone.
    fig, axes = plt.subplots(rows, 2, figsize=(20.0, 6.5 * rows + 1.6), squeeze=False,
                             gridspec_kw={"hspace": 0.06, "wspace": 0.05})
    cells = list(axes.flat)
    for ax, (channel, treatment) in zip(cells, panels):
        free = treatment == "profiled"
        draw(ax, channel, legends=False, treatment=treatment, models=None if free else ("geochemical",),
             columns=PAPER_COLUMNS, label=(r"$\nu_1\to\nu_3$ (IO), free radiogenic power" if free
                                           else PAPER_LABELS[channel[:2]]),
             label_size=size, ylabel=False, limits=limits, label_xy=PAPER_LABEL_XY)
        # Room above the top decade for the panel's label, which would otherwise hide the bands'
        # tops and Borexino's marker in row 1 and in the chain's tau_1 panel.
        ax.set_ylim(limits[0], limits[1] * 10 ** PAPER_HEADROOM_DECADES)
    for ax in cells[:-2]:
        ax.tick_params(axis="x", labelbottom=False)
    for ax in cells[-2:]:
        ax.tick_params(axis="x", labelsize=19)
    for ax in cells[1::2]:
        ax.tick_params(axis="y", labelleft=False)          # the left column's labels serve both
    # The legends take a fixed height above the panels, the experiments' names one below them.
    height = fig.get_size_inches()[1]
    fig.subplots_adjust(top=1.0 - 1.35 / height, left=0.08, right=0.98, bottom=0.95 / height)
    # One label for every y-axis: on the middle row's left panel, next to its ticks, moved to the
    # middle of the stack.
    top, bottom = cells[0].get_position().y1, cells[-1].get_position().y0
    anchor = cells[2 * ((rows - 1) // 2)]
    anchor.set_ylabel(r"Lower limit (90\% C.L.) on $\tau_i/m_i$ [s~eV$^{-1}$]")
    anchor.yaxis.label.set_y((0.5 * (top + bottom) - anchor.get_position().y0) / anchor.get_position().height)
    # Two legends side by side across the figure: the simulated experiments, and the limits from
    # data.  Split at the gap between the columns, each read as its column's own.  The second hugs its
    # two entries instead, at matplotlib's spacing and with the
    # money plot's short handle for markers, flush with the right column; the first takes the rest
    # of the row, the columns' gap apart.
    common = dict(borderaxespad=0.0, fontsize=20, title_fontsize=size, edgecolor="black", framealpha=1,
                  fancybox=False)
    simulated = ([Patch(color=YELLOW, lw=0), Patch(color=GREEN, lw=0), Line2D([], [], color="k", ls="--", lw=1.5),
                  Line2D([], [], marker="<", ms=10, mfc="white", mec="k", ls="none")],
                 [r"95\% band", r"68\% band", "Median", "Expected"])
    data = ([Line2D([], [], marker="o", ms=11, mfc="white", mec="k", mew=1.4, ls="none"),
             Line2D([], [], marker="D", ms=10, mfc="k", mec="k", ls="none")],
            ["Geochemical Earth", "Free radiogenic power"])
    left, right = cells[0].get_position().x0, cells[1].get_position().x1
    gap = cells[1].get_position().x0 - cells[0].get_position().x1
    legend = fig.legend(*data, title="Limits from data", ncol=2, loc="lower right", handlelength=1.2,
                        bbox_to_anchor=(right, top + 0.008), **common)
    width = legend.get_window_extent(fig.canvas.get_renderer()).width / fig.bbox.width
    fig.legend(*simulated, title="Simulated experiments, geochemical Earth", ncol=4, loc="lower left",
               mode="expand", bbox_to_anchor=(left, top + 0.008, right - width - gap - left, 0.0), **common)
    return fig


# The data, read from the release's data/ files.
from . import data as _data  # noqa: E402


def record(*args, **kwargs):
    """The paper's ``record``, read from the data file of the figure in use."""
    return getattr(_data.current(), 'record')(*args, **kwargs)


def simulated(*args, **kwargs):
    """The paper's ``simulated``, read from the data file of the figure in use."""
    return getattr(_data.current(), 'simulated')(*args, **kwargs)
