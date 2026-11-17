# The paper's plotting code, its functions that read the analysis's results replaced by ones that read
# the release's data/ files (at the end of this file).
"""Figure 1 of the paper: lower limits on tau_i/m_i from visible decay, ours and published ones.

Three panels on one log axis: nu_1 on top, nu_2 in the middle, nu_3 at the bottom.  In each, a bar
per constraint runs from the left edge, the excluded side, to its limit.  Ours come first, one row
per channel: the combined fit of all four experiments with the Earth's radiogenic power free, and
each experiment's own limit from the same fit as a marker (:func:`our_limits`).  A bar is hatched
from the limit that the reactor antineutrinos alone set to its end: the part that the geoneutrinos
add.  The limits on tau_3/m_3 are those of JUNO's reactor antineutrinos.  Published bounds on
visible decay follow, hatched, dark to the weakest of the cases their paper quotes and light to the
strongest; an astrophysical bound, beyond the axis, fills its row and ends in a triangle with its
value.  Bounds on invisible decay are left out.

The notebook draws the figure with :func:`figure`.
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon, Rectangle


from . import style as figure_style
from . import limits_overview as ov

#: The color of this work's bars, as in the spectrum figures.
OURS = "#C1272D"
#: Published bounds: a neutral gray, hatched.
PUBLISHED = "#8C8C8C"
#: Each experiment's tick: its color from the overview and a shape of its own, a second cue
#: where two colors are close for color-blind readers.
EXPERIMENTS = (("borexino", "Borexino (2020)", "#0072B2", "o"), ("kamland", "KamLAND (2022)", "#E69F00", "s"),
               ("snoplus", "SNO+ (2026)", "#009E73", "D"), ("juno", "JUNO (2026)", "#CC79A7", "P"))
#: The combined fit whose bars these are: all four experiments.
COMBINED = "combined_all"
BSES = ("cosmochemical", "geochemical", "geodynamical")
#: The axis range, s/eV: this work's limits and the terrestrial ones; astrophysical bounds
#: lie beyond it.
LEFT, RIGHT = 1e-12, 1e-7
#: hbar in eV s, to turn alpha = m/tau into tau/m.
HBAR_EV_S = 6.582119569e-16

#: The panels: the decaying state and its rows.  A row of ours
#: is ("ours", case, lifetime, label); a published one is ("published", plot_row, label,
#: source), where plot_row names its rows in the literature table.  A row ("none", cases,
#: lifetime, label) says that the combination limits none of the cases; no row is one now.
PANELS = (
    (r"$\nu_1$", (
        ("ours", "io_n1_n3", "log10_tau1_over_m1", r"$\nu_1\to\nu_3$ (IO)"),
        ("ours", "io_n2_n1_n3", "log10_tau1_over_m1", r"$\nu_2\to\nu_1\to\nu_3$ (IO)"),
        ("published", "glashow_nu1", r"$\nu_1\to\nu_3$ (IO)", "High-energy astro., IceCube, Glashow"),
        ("published", "sn1987a_nu1", r"$\nu_1\to\nu_3$ (IO)", "SN1987A"),
    )),
    (r"$\nu_2$", (
        ("ours", "io_n2_n3", "log10_tau2_over_m2", r"$\nu_2\to\nu_3$ (IO)"),
        ("ours", "io_n2_n1", "log10_tau2_over_m2", r"$\nu_2\to\nu_1$ (IO)"),
        ("ours", "io_n2_n1_n3", "log10_tau2_over_m2", r"$\nu_2\to\nu_1\to\nu_3$ (IO)"),
        ("ours", "no_n2_n1", "log10_tau2_over_m2", r"$\nu_2\to\nu_1$ (NO)"),
        ("ours", "no_n3_n2_n1", "log10_tau2_over_m2", r"$\nu_3\to\nu_2\to\nu_1$ (NO)"),
        ("published", "kamland_reactor_nu2", r"$\nu_2,\nu_1\to\nu_3$ (IO)", "Reactors, KamLAND"),
        ("published", "solar_antinu_nu2", r"$\nu_2\to\bar\nu_1$" " (IO or NO, Majorana)", r"Sun, $\bar\nu_e$ search"),
        ("published", "solar_nu2_dirac", r"$\nu_2\to\nu_1$" " (IO or NO, Dirac)", r"Sun, $\nu$ spectra"),
        ("published", "glashow_nu2", r"$\nu_2\to\nu_3$ (IO)", "High-energy astro., IceCube, Glashow"),
        ("published", "sn1987a_nu2", r"$\nu_2\to\nu_1,\nu_3$ (IO)", "SN1987A"),
    )),
    (r"$\nu_3$", (
        ("ours", "no_n3_n1", "log10_tau3_over_m3", r"$\nu_3\to\nu_1$ (NO)"),
        ("ours", "no_n3_n2", "log10_tau3_over_m3", r"$\nu_3\to\nu_2$ (NO)"),
        ("ours", "no_n3_n2_n1", "log10_tau3_over_m3", r"$\nu_3\to\nu_2\to\nu_1$ (NO)"),
        ("published", "t2k_minos_nu3", r"$\nu_3\to\nu_1$ (NO)", "T2K, MINOS"),
        ("published", "solar_antinu_nu3", r"$\nu_3\to\bar\nu_1$" " (NO, Majorana)", r"Sun, $\bar\nu_e$ spectra"),
    )),
)
#: Each published bound's paper, which the short labels in the panels leave to the caption.
CITATIONS = {
    "glashow_nu1": "IceCube's Glashow-resonance events (Bustamante 2020)",
    "glashow_nu2": "IceCube's Glashow-resonance events (Bustamante 2020)",
    "sn1987a_nu1": "SN1987A, with tau_1/m_1 = tau_2/m_2 (Ivanez-Ballesteros and Volpe 2023)",
    "sn1987a_nu2": "SN1987A, with tau_1/m_1 = tau_2/m_2 (Ivanez-Ballesteros and Volpe 2023)",
    "kamland_reactor_nu2": "KamLAND's reactor antineutrinos (Porto-Silva et al. 2020)",
    "solar_antinu_nu2": "KamLAND's search for solar antineutrinos (KamLAND 2004)",
    "solar_nu2_dirac": "solar neutrinos, for Dirac neutrinos (our reading of de Gouvea et al. 2024's bound on the coupling)",
    "t2k_minos_nu3": "T2K and MINOS (Gago et al. 2017)",
    "solar_antinu_nu3": "solar antineutrino spectra (Picoreti et al. 2022)",
}
#: The text sizes: the channels, the bounds' labels, the notes and the legend.
SIZES = {"channel": 23, "source": 23, "note": 23, "legend": 23, "inside": 23}
#: Half the thickness of a bar, in rows.
HALF = 0.36
#: The unit of tau/m as every label writes it.
UNIT = r"s~eV$^{-1}$"
#: Where the in-panel notes start, s eV^-1: the shaded band's label and the row without a limit.
NOTE_X = 1.3e-12
#: The fit our rows show: the radiogenic power free, ``(bse, treatment)`` as the records name it.
MODEL = ("geochemical", "profiled")
#: Every published bound's bar, in one gray, dark enough to carry a white label.
PUBLISHED_STYLE = dict(fc="#8C8C8C", hatch="///", ec="#A8A8A8", lw=0)
#: The decades a bar must span past the labels' start to hold its label inside.
INSIDE_DECADES = 2.0
#: The hatch of the part of a bar that the geoneutrinos add, and its name in the legend.
ADDED_HATCH, ADDED_LABEL = "\\\\\\", "Added by geoneutrinos"
#: The legend's title.
LEGEND_TITLE = r"Disfavored by geoneutrinos (free radiogenic power) + reactors (90\% C.L., this work)"
#: Where the short bar of a bound beyond the axis starts, in decades below the right edge.
BEYOND_DECADES = 0.7


def tex_number(x: float) -> str:
    """A number as LaTeX, in powers of ten.

    Args:
        x: The number.

    Returns:
        For example ``1.2\\times10^{5}``.
    """
    exponent = math.floor(math.log10(x))
    mantissa = f"{x / 10 ** exponent:.2f}".rstrip("0").rstrip(".")
    return rf"{mantissa}\times10^{{{exponent}}}"


def our_limits(case: str, param: str) -> tuple:
    """The limits a row of ours shows: all four's, and each experiment's own, with the power free.

    Args:
        case: The channel's case.
        param: The lifetime limited.

    Returns:
        ``(combined limit, the experiments' limits by site)``, s/eV; an experiment that sets
        none has ``0.0``.

    Raises:
        ValueError: If a fit has not been run, or the combination sets no limit.
    """
    bse, treatment = MODEL
    combined = ov.limit(COMBINED, case, param, bse, treatment)
    singles = {site: ov.limit(site, case, param, bse, treatment) for site, *_ in EXPERIMENTS}
    if combined is None or any(v is None for v in singles.values()):
        raise ValueError(f"{case} {param} with the power free: not every fit has been run")
    if not combined:
        raise ValueError(f"{case} {param}: the combination sets no limit with the power free")
    return combined, singles


def ours(ax, y: float, case: str, param: str) -> None:
    """One row of ours: the combined fit with the power free, what the geoneutrinos add, and each experiment.

    Args:
        ax: The panel.
        y: The row.
        case: The channel's case.
        param: The lifetime limited.
    """
    combined, limits = our_limits(case, param)
    ax.add_patch(Rectangle((LEFT, y - HALF), combined - LEFT, 2 * HALF, color=OURS, alpha=0.85, lw=0, zorder=2))
    reactors = reactors_alone(case, param)
    if reactors is not None and reactors < combined:
        start = max(reactors, LEFT)
        ax.add_patch(Rectangle((start, y - HALF), combined - start, 2 * HALF, fc="none", ec="white", hatch=ADDED_HATCH,
                               lw=0, zorder=3))
    for site, _, color, marker in EXPERIMENTS:
        if limits[site]:
            ax.plot(limits[site], y, marker=marker, ms=12, mfc=color, mec="white", mew=0.9, ls="none", zorder=4)


def none_limited(ax, y: float, cases: tuple, param: str) -> None:
    """One row for channels that the combination does not limit on any Earth model, said so.

    Args:
        ax: The panel.
        y: The row.
        cases: The channels.
        param: The lifetime.

    Raises:
        ValueError: If any of them is limited on any Earth model, or not run.
    """
    for case in cases:
        for bse in BSES:
            v = ov.limit(COMBINED, case, param, bse)
            if v is None or v:
                raise ValueError(f"{case} on the {bse} Earth: {v}; a limited channel needs a row of its own")
    ax.text(NOTE_X, y, r"No 90\% C.L.~limit from geoneutrinos", va="center_baseline", ha="left", fontsize=SIZES["source"],
            color="k", zorder=3)


def published(ax, y: float, bound: dict, source: str) -> None:
    """One row of a published bound: hatched, dark to its weakest case and light to its strongest.

    Every published bar is one gray.  A bound beyond the axis fills its row from the left edge
    and ends in a right-pointing triangle just inside the right edge.  The label goes inside the
    bar, in white, where the bar is long enough for it, and beside it otherwise.

    Args:
        ax: The panel.
        y: The row.
        bound: From ``published_bounds``.
        source: The data and the paper.
    """
    low, high = bound["low"], bound["high"]
    if low >= RIGHT:
        tip = RIGHT / 10 ** 0.03                    # the triangle's tip, just inside the edge
        base = tip / 10 ** 0.16
        ax.add_patch(Rectangle((LEFT, y - HALF), base - LEFT, 2 * HALF, zorder=2, **PUBLISHED_STYLE))
        ax.add_patch(Polygon([(base, y - HALF), (base, y + HALF), (tip, y)], closed=True, fc=PUBLISHED_STYLE["fc"],
                             lw=0, zorder=2))
        value = rf"$\geq {tex_number(low)}$" + (rf"--${tex_number(high)}$" if high > low else "") + " " + UNIT
        ax.text(NOTE_X, y, f"{source}: {value}", va="center_baseline", ha="left", fontsize=SIZES["inside"], color="white",
                zorder=3)
        return
    ax.add_patch(Rectangle((LEFT, y - HALF), low - LEFT, 2 * HALF, zorder=2, **PUBLISHED_STYLE))
    if high > low:
        ax.add_patch(Rectangle((low, y - HALF), high - low, 2 * HALF, alpha=0.45, zorder=2, **PUBLISHED_STYLE))
    if math.log10(low / NOTE_X) > INSIDE_DECADES:
        ax.text(NOTE_X, y, source, va="center_baseline", ha="left", fontsize=SIZES["inside"], color="white", zorder=3)
    else:
        ax.text(high * 1.6, y, source, va="center_baseline", ha="left", fontsize=SIZES["source"], color="0.2", zorder=3)


def figure():
    """The money plot, in the repository's style, two columns wide.

    Returns:
        The matplotlib figure.
    """
    figure_style.use()
    bounds = published_bounds()
    heights = [len(rows) for _, rows in PANELS]
    fig, axes = plt.subplots(3, 1, figsize=(17.0, 15.0), sharex=True,
                             gridspec_kw={"height_ratios": heights, "hspace": 0.06})
    for ax, (state, rows) in zip(axes, PANELS):
        labels, kinds = [], []
        for y, row in enumerate(rows):
            if row[0] == "ours":
                _, case, param, label = row
                ours(ax, y, case, param)
                kinds.append("ours")
            elif row[0] == "none":
                _, cases, param, label = row
                none_limited(ax, y, cases, param)
                kinds.append("ours")
            else:
                _, key, label, source = row
                published(ax, y, bounds[key], source)
                # Terrestrial and astrophysical bounds form one group: a line separates ours from them only.
                kinds.append("published")
            labels.append(label)
        for i in range(len(rows) - 1):
            if kinds[i] != kinds[i + 1]:
                ax.axhline(i + 0.5, color="0.6", lw=0.8, zorder=1)
        ax.set_ylim(len(rows) - 0.5, -0.5)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels(labels, fontsize=SIZES["channel"])
        for tick in ax.get_yticklabels():
            tick.set_multialignment("right")              # a two-line label keeps both lines at the axis
        ax.tick_params(axis="y", which="both", left=False, right=False)
    axes[-1].set_xscale("log")
    axes[-1].set_xlim(LEFT, RIGHT)
    axes[-1].set_xlabel(rf"Neutrino lifetime, $\tau_i/m_i$ [{UNIT}]")
    # The decaying states, left-aligned with one another, left of the widest channel label.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    left = min(t.get_window_extent(renderer).x0 for ax in axes for t in ax.get_yticklabels())
    x = axes[0].transAxes.inverted().transform((left, 0.0))[0]
    for ax, (state, _rows) in zip(axes, PANELS):
        ax.text(x - 0.02, 0.5, state, transform=ax.transAxes, ha="right", va="center", fontsize=30)
    handles = ([Line2D([], [], marker=m, ms=12, mfc=c, mec="white", mew=0.9, ls="none") for _, _, c, m in EXPERIMENTS]
               + [Patch(color=OURS, alpha=0.85, lw=0),
                  Patch(fc=OURS, alpha=0.85, ec="white", hatch=ADDED_HATCH, lw=0)])
    names = [name.split(" (")[0] for _, name, _, _ in EXPERIMENTS] + ["All four", ADDED_LABEL]
    # borderaxespad=0: the frame spans the panel's width exactly, with no inset at either side.
    legend = axes[0].legend(handles, names, title=LEGEND_TITLE, loc="lower left", bbox_to_anchor=(0.0, 1.035, 1.0, 0.0),
                            mode="expand", borderaxespad=0.0, ncol=len(names), edgecolor="black", framealpha=1,
                            fontsize=SIZES["legend"], title_fontsize=plt.rcParams["axes.labelsize"],
                            handlelength=1.2, handletextpad=0.35, borderpad=0.35)
    legend.set_clip_on(False)
    return fig


# The data, read from the release's data/ files.
from . import data as _data  # noqa: E402


def published_bounds(*args, **kwargs):
    """The paper's ``published_bounds``, read from the data file of the figure in use."""
    return getattr(_data.current(), 'published_bounds')(*args, **kwargs)


def reactors_alone(*args, **kwargs):
    """The paper's ``reactors_alone``, read from the data file of the figure in use."""
    return getattr(_data.current(), 'reactors_alone')(*args, **kwargs)
