# The paper's plotting code, its functions that read the analysis's results replaced by ones that read
# the release's data/ files (at the end of this file).
"""Figure D.1 of the paper: the coverage of the 90% C.L. limits, tested with simulated experiments.

Above, an example: one fit's simulated limits, made with the true lifetime at its expected limit,
and the shares on either side of the truth.  Below, for every fit and every true lifetime the test
was run at, the share of the simulated experiments whose limit lies above the truth, with its
binomial 68% interval (Wilson's); one panel has the Earth's radiogenic power assumed, the other has
it free.

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
from .analysis import rct

#: Each fit's color, as the overview draws it (Okabe-Ito for the experiments, grays for the combinations).
COLORS = {"borexino": "#0072B2", "kamland": "#E69F00", "snoplus": "#009E73", "juno": "#CC79A7",
          "combined_bx_kl_snoplus": "0.55", "combined_all": "0.15"}
#: The truths' markers, longest lifetime first: at the expected limit open, deeper filled.
OFFSETS = (-0.14, 0.14)
#: The paper's fits: every fit of the test, all four last; one not on file is marked as running.
PAPER_FITS = rct.FITS
#: The paper's example, its upper panel: a fit and its treatment of the Earth's power, at the
#: truth at its expected limit.  KamLAND, since every one of its simulated experiments sets a
#: limit there: the two shares are then the bars on either side of the truth.
EXAMPLE = ("kamland", "assumed")
#: The simulated limits that exclude the true lifetime, and those that do not; the shares that
#: exclude it are drawn in the first color below as well.
EXCLUDES, ALLOWS = "#D55E00", "0.72"
#: The width of the example's bins, in decades of tau/m, and the room above its
#: tallest bar, as a factor: the panel's top is that, rounded up to ten simulated experiments.
EXAMPLE_BIN = 0.1
EXAMPLE_HEADROOM = 1.3
#: Where a fit's two points sit about its column in the paper: the truth at the expected limit,
#: then half a decade shorter.
PAPER_OFFSETS = {True: -0.12, False: 0.12}
#: The simulated experiments a test holds when complete; one with fewer, still running or a test
#: on one of them, is not drawn.
PAPER_SIMULATED = 200


def wilson(k: int, n: int, z: float = 1.0) -> tuple[float, float]:
    """Wilson's binomial interval for a share.

    Args:
        k: Successes.
        n: Trials.
        z: Its half-width in standard deviations, one for $68\\%$.

    Returns:
        ``(low, high)`` as shares.
    """
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return max(centre - half, 0.0), min(centre + half, 1.0)


def paper_points() -> list:
    """Every point of the paper's lower panel.

    Returns:
        ``(column, site, treatment, at_expected, excluded, n)`` per fit, treatment and truth; the
        column is the fit's place in :data:`PAPER_FITS`.
    """
    out = []
    for k, (site, _name) in enumerate(PAPER_FITS):
        for treatment in ("assumed", "profiled"):
            for truth, records, _ in rct.points(site, treatment):
                excluded, n = dl.coverage(records, truth, "90")
                out.append((k, site, treatment, at_expected_limit(site, treatment, truth), excluded, n))
    return out


def corner_label(ax, text: str, size: float, gap: float = 0.2) -> None:
    """A rounded-box label in a panel's upper left corner, its box as far from both spines.

    Args:
        ax: The panel, laid out.
        text: The label; a line break starts a second line, left-aligned.
        size: Its font size, in points.
        gap: The space between the box and each spine, in inches.
    """
    width, height = ax.figure.get_size_inches()
    pos = ax.get_position()
    inset = gap + 0.3 * size / 72.0                     # the gap, then the box's own padding
    ax.text(inset / (pos.width * width), 1.0 - inset / (pos.height * height), text, transform=ax.transAxes,
            ha="left", va="top", multialignment="left", fontsize=size, zorder=10,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="k", lw=0.8))


def figure_paper():
    """The coverage test as the paper shows it, at 90% C.L.: an example, then every fit.

    Only the 90% limits, which are the paper's; an example first, so that the share is seen before
    it is plotted; one panel per treatment of the Earth's power, so that each point's marker says
    only its true lifetime; which side of the target is which, written beside it.  Top: the simulated
    90% limits of :data:`EXAMPLE` with the true lifetime at its expected limit, and the shares on
    either side of it.  Middle and bottom: the share that excludes the true lifetime at every fit,
    with the Earth's power assumed and free, with its binomial 68% interval (Wilson's), against
    the 10% a correct 90% limit gives.

    Returns:
        The matplotlib figure.
    """
    import numpy as np
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MultipleLocator

    from . import brazil as plot_brazil
    figure_style.use()
    size = plt.rcParams["axes.labelsize"]
    names = dict(zip((site for site, _ in PAPER_FITS), plot_brazil.PAPER_COLUMNS))
    fig = plt.figure(figsize=(10.0, 18.0))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.75], hspace=0.2)
    upper = fig.add_subplot(outer[0])
    inner = outer[1].subgridspec(2, 1, hspace=0.06)
    middle = fig.add_subplot(inner[0])
    lower = fig.add_subplot(inner[1], sharex=middle, sharey=middle)

    # Top: one fit's simulated limits at one truth, cut at the truth.
    site, treatment = EXAMPLE
    truth, records, _ = next(p for p in rct.points(site, treatment) if at_expected_limit(site, treatment, p[0]))
    excluded, n = dl.coverage(records, truth, "90")
    values = np.log10([r["limits_s_per_ev"]["90"] for r in records if r["limits_s_per_ev"]["90"] is not None])
    assert values.size == n, "the example's shares are its bars: every simulated experiment must set a limit"
    # Bins with an edge at the truth, so that each lies wholly on one side of it.
    edges = truth + EXAMPLE_BIN * np.arange(-math.ceil((truth - values.min()) / EXAMPLE_BIN),
                                            math.ceil((values.max() - truth) / EXAMPLE_BIN) + 2)
    counts, _ = np.histogram(values, edges)
    beyond = edges[:-1] > truth - 1e-9
    assert counts[beyond].sum() == excluded, (counts[beyond].sum(), excluded)
    upper.stairs(np.where(beyond, 0, counts), 10.0 ** edges, fill=True, color=ALLOWS, lw=0)
    upper.stairs(np.where(beyond, counts, 0), 10.0 ** edges, fill=True, color=EXCLUDES, lw=0)
    upper.axvline(10.0 ** truth, color="k", ls="--", lw=1.6, zorder=3)
    upper.set_xscale("log")
    upper.set_xlim(10.0 ** edges[0], 10.0 ** edges[-1])
    top = 10.0 * math.ceil(EXAMPLE_HEADROOM * counts.max() / 10.0)
    upper.set_ylim(0.0, top)
    upper.yaxis.set_minor_locator(MultipleLocator(1.0))
    upper.set_xlabel(r"Simulated $90\%$ C.L. lower limit on $\tau_1/m_1$ [s~eV$^{-1}$]")
    upper.set_ylabel("Simulated experiments")
    # The line's name along it, above the bars right of it; the shares on either side of it, above
    # every bar.
    upper.text(10.0 ** (truth + 0.025), counts[beyond].max() + 0.12 * top, "True lifetime", rotation=90,
               ha="left", va="bottom", fontsize=20)
    height = counts.max() + 0.04 * top
    upper.text(10.0 ** (truth - 0.03), height, rf"${100 * (n - excluded) / n:g}\%$", ha="right", va="bottom",
               fontsize=size, color="0.35")
    upper.text(10.0 ** (truth + 0.03), height, rf"${100 * excluded / n:g}\%$", ha="left", va="bottom",
               fontsize=size, color=EXCLUDES)
    corner_label(upper, names[site].replace("\n", " ") + "\n" + r"$\nu_1\to\nu_3$ (IO)", size)

    # Middle and bottom: the share that excludes the truth at every fit, the Earth's power assumed
    # (open markers) and free (filled), as in the paper's other figures.
    points = paper_points()
    for ax, treatment, name in ((middle, "assumed", "Geochemical Earth"), (lower, "profiled", "Free radiogenic power")):
        shown = set()
        for k, site, which, expected, excluded, n in points:
            if which != treatment or n < PAPER_SIMULATED:
                continue
            low, high = wilson(excluded, n)
            share = excluded / n
            ax.errorbar(k + PAPER_OFFSETS[expected], 100 * share, yerr=[[100 * (share - low)], [100 * (high - share)]],
                        fmt="o" if expected else "s", ms=10, mec=EXCLUDES,
                        mfc=EXCLUDES if treatment == "profiled" else "white", ecolor=EXCLUDES, elinewidth=1.6,
                        capsize=3, mew=1.6, zorder=4)
            shown.add(k)
        for k in set(range(len(PAPER_FITS))) - shown:
            ax.text(k, 5.0, "running", ha="center", va="center", fontsize=17, color="0.45", style="italic",
                    zorder=3, bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none"))
        ax.axhline(10.0, color="0.35", lw=1.4, ls="--", zorder=1)
        ax.axvline(3.5, color="0.5", lw=1, ls=":", zorder=0)
        # Which side of the target is which, outside the right spine, where no point can reach.
        side = ax.get_yaxis_transform()
        ax.text(1.015, 14.0, r"Too strong $\rightarrow$", rotation=90, transform=side, ha="left", va="center",
                fontsize=20, color="0.3")
        ax.text(1.015, 5.0, r"$\leftarrow$ Conservative", rotation=90, transform=side, ha="left", va="center",
                fontsize=20, color="0.3")
        ax.yaxis.set_minor_locator(MultipleLocator(1.0))
        ax.tick_params(axis="x", which="both", top=False, bottom=False)
    middle.set_xlim(-0.6, len(PAPER_FITS) - 0.4)
    middle.set_ylim(0.0, 18.0)
    middle.set_yticks(range(0, 16, 5))
    corner_label(middle, "Geochemical Earth", size)
    corner_label(lower, "Free radiogenic power", size)
    middle.tick_params(axis="x", labelbottom=False)
    lower.set_xticks(range(len(PAPER_FITS)))
    lower.set_xticklabels([names[site] for site, _ in PAPER_FITS])
    lower.tick_params(axis="x", labelsize=19)
    # The two true lifetimes, once, in the bottom panel, empty above its target: with the power
    # assumed, all four's points could reach up there.
    marks = [Line2D([], [], marker=m, ms=10, mew=1.6, mec=EXCLUDES, mfc=EXCLUDES, ls="none") for m in "os"]
    lower.legend(marks, ["Expected limit", rf"Expected limit$\,/\,{10 ** 0.5:.1f}$"], title="True lifetime",
                 title_fontsize=20, fontsize=20, loc="upper right", bbox_to_anchor=(0.985, 0.97), borderaxespad=0.0,
                 edgecolor="k", framealpha=1.0, fancybox=False)
    # One y-axis label for the two panels, centred on them.
    label = r"Limits that exclude the true lifetime [\%]"
    middle.set_ylabel(label)
    top_, bottom_ = middle.get_position().y1, lower.get_position().y0
    middle.yaxis.label.set_y((0.5 * (top_ + bottom_) - middle.get_position().y0) / middle.get_position().height)
    return fig


# The data, read from the release's data/ files.
from . import data as _data  # noqa: E402


def at_expected_limit(*args, **kwargs):
    """The paper's ``at_expected_limit``, read from the data file of the figure in use."""
    return getattr(_data.current(), 'at_expected_limit')(*args, **kwargs)
