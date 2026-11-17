# The paper's plotting code, its functions that read the analysis's results replaced by ones that read
# the release's data/ files (at the end of this file).
"""Figure 11 of the paper: every lower limit from data, for every channel, experiment and Earth model.

One block per channel and lifetime, down the page: the two chains limit two lifetimes each, so the
eight channels make ten blocks, inverted ordering first.  Each block has a row for each experiment,
then one for each combination.  In a row, a bar runs from the left edge, the excluded side, to the
90% limit from data with the Earth's radiogenic power free, and ends in a filled diamond; an open
marker gives the limit on each of the three Earth models.  A row that sets a limit on some Earth
model and none with the power free has markers and no bar.  The top axis gives the decay length at
2 MeV, c E tau/m.

The notebook draws the figure with :func:`figure`; :func:`limit` reads each limit from
``data/limits.json``.
"""

from __future__ import annotations

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle


from .analysis import dl
from . import style as figure_style

#: The blocks: each channel's case, the lifetime limited, the channel and the lifetime as
#: labelled, and for a chain the chain its two blocks share.
BLOCKS = (
    ("io_n1_n3", "log10_tau1_over_m1", r"$\nu_1\to\nu_3$", r"$\tau_1/m_1$", None),
    ("io_n2_n3", "log10_tau2_over_m2", r"$\nu_2\to\nu_3$", r"$\tau_2/m_2$", None),
    ("io_n2_n1", "log10_tau2_over_m2", r"$\nu_2\to\nu_1$", r"$\tau_2/m_2$", None),
    ("io_n2_n1_n3", "log10_tau2_over_m2", r"$\nu_2\to\nu_1\to\nu_3$", r"$\tau_2/m_2$", "io_chain"),
    ("io_n2_n1_n3", "log10_tau1_over_m1", r"$\nu_2\to\nu_1\to\nu_3$", r"$\tau_1/m_1$", "io_chain"),
    ("no_n2_n1", "log10_tau2_over_m2", r"$\nu_2\to\nu_1$", r"$\tau_2/m_2$", None),
    ("no_n3_n1", "log10_tau3_over_m3", r"$\nu_3\to\nu_1$", r"$\tau_3/m_3$", None),
    ("no_n3_n2", "log10_tau3_over_m3", r"$\nu_3\to\nu_2$", r"$\tau_3/m_3$", None),
    ("no_n3_n2_n1", "log10_tau3_over_m3", r"$\nu_3\to\nu_2\to\nu_1$", r"$\tau_3/m_3$", "no_chain"),
    ("no_n3_n2_n1", "log10_tau2_over_m2", r"$\nu_3\to\nu_2\to\nu_1$", r"$\tau_2/m_2$", "no_chain"),
)
#: How many blocks belong to inverted ordering; the rest are normal ordering.
N_INVERTED = 5
#: The rows of a block, top to bottom: the experiments, each in a colorblind-safe color (Okabe
#: and Ito) and named with the year of the paper whose data are fitted,
#: then the combinations, in gray and black, with taller bars.
ROWS = (("borexino", "Borexino (2020)", "#0072B2", 0.55), ("kamland", "KamLAND (2022)", "#E69F00", 0.55),
        ("snoplus", "SNO+ (2026)", "#009E73", 0.55), ("juno", "JUNO (2026)", "#CC79A7", 0.55),
        ("combined_bx_kl_snoplus", "BX+KL+SNO+", "#6E6E6E", 0.75), ("combined_all", "All four", "#000000", 0.75))
#: The Earth models, each with its marker and how the Earth's power is treated: one of the three
#: models assumed, or the power free, profiled across the range they span (starting from the
#: geochemical Earth).
EARTHS = (("cosmochemical", "Cosmochemical", "v", "assumed"), ("geochemical", "Geochemical", "o", "assumed"),
          ("geodynamical", "Geodynamical", "^", "assumed"), ("geochemical", "Free radiogenic power", "D", "profiled"))
#: The axis range, s/eV (it stops at 1e-7, where no limit reaches, and starts at 1e-12, so that the
#: limits on tau_3/m_3 are no stubs at the edge), and
#: where a row whose data set no limit on any Earth model says so.
LEFT, RIGHT = 1e-12, 1e-7
NO_LIMIT_X = 1.4e-12
#: The decay lengths marked on the top axis, with a dashed line each: the Earth's radius and diameter.
EARTH_MARKS = ((r"$R_\oplus$", 1.0), (r"$2R_\oplus$", 2.0))
#: How close, in decades, a decade's label may come to an Earth mark's before it is left out.
MARK_CLEARANCE = 0.15
#: The legends' distance from the panel's lower right corner, and between one another, in axes fractions.
LEGEND_MARGIN, LEGEND_GAP = 0.012, 0.01
#: The decay length at 2 MeV per unit of tau/m: c times E, in km per (s/eV).
KM_PER_S_PER_EV = 2.998e5 * 2e6
EARTH_RADIUS_KM = 6371.0
#: The bars' opacity.
DARK = 0.55
#: A chain's brace: where its tips sit and how deep it is, in fractions of the axes' width,
#: left of the axis and its lifetime labels.
BRACE_X, BRACE_DEPTH = -0.075, 0.018
#: Spacing, in rows: between the experiments and the combinations of a block, between
#: blocks, and the header above each ordering.
GAP_COMBINED, GAP_BLOCK, GAP_HEADER = 0.5, 1.0, 2.4


def brace(ax, top: float, bottom: float, label: str) -> None:
    """A curly brace joining a chain's two blocks, with the chain to its left.

    Drawn left of the axis, outside it, in the axes' width for x and in rows
    for y; its arms are two pairs of logistic curves, which give a brace's
    rounded tips and its point.

    Args:
        ax: The axes.
        top: The first block's top, in rows.
        bottom: The second block's bottom, in rows.
        label: The chain.
    """
    from matplotlib.transforms import blended_transform_factory
    trans = blended_transform_factory(ax.transAxes, ax.transData)
    y = np.linspace(top, bottom, 401)
    half = y[:201]
    sharp = 40.0 / (bottom - top)
    arm = 1 / (1 + np.exp(-sharp * (half - half[0]))) + 1 / (1 + np.exp(-sharp * (half - half[-1])))
    shape = np.concatenate([arm, arm[-2::-1]])
    shape = (shape - shape.min()) / (shape.max() - shape.min())
    ax.plot(BRACE_X - BRACE_DEPTH * shape, y, transform=trans, color="black", lw=1.3, clip_on=False)
    ax.text(BRACE_X - BRACE_DEPTH - 0.022, 0.5 * (top + bottom), label, transform=trans, ha="center",
            va="center", rotation=90, fontsize=17)


def table() -> dict:
    """Every stored limit the figure shows.

    Returns:
        ``{(block, row): {Earth model's name: limit}}``, as :func:`limit` returns them.
    """
    return {(k, j): {name: limit(site, case, param, bse, treatment) for bse, name, _, treatment in EARTHS}
            for k, (case, param, *_) in enumerate(BLOCKS) for j, (site, *_) in enumerate(ROWS)}


def limited(found: dict) -> bool:
    """Whether a row sets a limit on some Earth model, or with the power free.

    Args:
        found: A row of :func:`table`.

    Returns:
        True if any of its fits sets a 90% limit.
    """
    return any(found.values())


def drawn(found: dict) -> bool:
    """Whether a row's limit reaches the plotted range, so that the row is drawn.

    A limit below :data:`LEFT` is named in the block's note instead.

    Args:
        found: A row of :func:`table`.

    Returns:
        True if its strongest limit lies at or above :data:`LEFT`.
    """
    return limited(found) and max(v for v in found.values() if v) >= LEFT


def layout(limits: dict):
    """Where each block, row, note and header sits, in rows from the top.

    A block draws the rows that set a limit, experiments then combinations, and then one line
    naming those that set none, one such line per block.

    Args:
        limits: From :func:`table`.

    Returns:
        ``(rows, blocks, headers, notes, bottom)``: the y of every drawn row as
        ``rows[(k, j)]``, each block's ``(top, bottom, center)``, each ordering's header
        ``(y, text)``, each block's note ``notes[k] = (y, [rows without a limit], [rows whose
        limit lies below the plotted range])``, and the y
        of the last row's edge.
    """
    rows, blocks, headers, notes = {}, [], [], {}
    y = 0.0
    for k in range(len(BLOCKS)):
        if k in (0, N_INVERTED):
            headers.append((y + 0.75, "Inverted ordering" if k == 0 else "Normal ordering"))
            y += GAP_HEADER
        top = y
        shown = [j for j in range(len(ROWS)) if drawn(limits[(k, j)])]
        for i, j in enumerate(shown):
            if i and shown[i - 1] < 4 <= j:
                y += GAP_COMBINED
            rows[(k, j)] = y + 0.5
            y += 1.0
        missing = [j for j in range(len(ROWS)) if not limited(limits[(k, j)])]
        below = [j for j in range(len(ROWS)) if limited(limits[(k, j)]) and not drawn(limits[(k, j)])]
        if missing or below:
            # A block of the note alone keeps two rows, room for its two-line label.
            pad = 0.5 if not shown else 0.0
            notes[k] = (y + pad + 0.5, missing, below)
            y += 1.0 + 2 * pad
        blocks.append((top, y, 0.5 * (top + y)))
        y += GAP_BLOCK
    return rows, blocks, headers, notes, y - GAP_BLOCK


def draw(ax):
    """Draw every block and its rows, the headers and the decay-length axis.

    The lines naming the rows without a limit and the legends are placed afterwards, by
    :func:`finish`, once the figure's layout is fixed: both are aligned with things whose size
    is known only when drawn.

    Args:
        ax: The axes to draw on.

    Returns:
        ``(top, headers, notes)``: the decay-length axis, the headers' text artists, and each
        line to place as ``(y, text)``.
    """
    limits = table()
    rows, blocks, headers, notes, bottom = layout(limits)
    for k, (top, end, _) in enumerate(blocks):
        if k % 2 == 0:
            ax.axhspan(top - 0.25, end + 0.25, color="0.96", lw=0, zorder=0)
    header_texts = [ax.text(1.5 * LEFT, y, text + r" $\blacktriangledown$", ha="left", va="center", fontsize=20,
                            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="black", lw=1.0))
                    for y, text in headers]
    for _, multiple in EARTH_MARKS:
        ax.axvline(multiple * EARTH_RADIUS_KM / KM_PER_S_PER_EV, color="0.45", lw=1.0, ls="--", zorder=1)
    names = [name.split(" (")[0] for _, name, _, _ in ROWS]
    def note(missing, below):
        """The line under a block: who sets no limit, and whose limit lies below the range."""
        if len(missing) == len(ROWS):
            return r"No 90\% C.L.~limit from any experiment or combination"
        parts = ([r"No 90\% C.L.~limit: " + ", ".join(names[j] for j in missing)] if missing else []) + \
                ([rf"Limit below $10^{{{int(round(np.log10(LEFT)))}}}$~s~eV$^{{-1}}$: " + ", ".join(names[j] for j in below)]
                 if below else [])
        return "; ".join(parts)
    lines = [(y, note(missing, below)) for y, missing, below in notes.values()]
    for (k, j), y in rows.items():
        _, _, color, height = ROWS[j]
        found = limits[(k, j)]
        # The bar reads the limit with the power free; the three Earth models are markers alone.
        free = next(found[name] for _, name, _, treatment in EARTHS if treatment == "profiled")
        if free:
            ax.add_patch(Rectangle((LEFT, y - height / 2), free - LEFT, height, color=color, alpha=DARK, lw=0,
                                   zorder=2))
        for _, name, marker, treatment in EARTHS:
            v = found[name]
            if v:
                ax.plot(v, y, marker=marker, ms=9, mfc=color if treatment == "profiled" else "white", mec=color,
                        mew=1.6, ls="none", zorder=5 if treatment == "profiled" else 4)
    ax.set_xscale("log")
    ax.set_xlim(LEFT, RIGHT)
    ax.set_ylim(bottom + 0.4, -0.9)
    ax.set_yticks([c for _, _, c in blocks])
    ax.set_yticklabels([life if chain else f"{channel}\n{life}" for _, _, channel, life, chain in BLOCKS],
                       fontsize=17, multialignment="center")
    for chain in dict.fromkeys(c for *_, c in BLOCKS if c):
        members = [k for k, block in enumerate(BLOCKS) if block[-1] == chain]
        brace(ax, blocks[members[0]][0], blocks[members[-1]][1], BLOCKS[members[0]][2])
    ax.tick_params(axis="y", which="both", left=False, right=False)
    ax.set_xlabel(r"Neutrino lifetime, $\tau_i/m_i$ [s~eV$^{-1}$]")
    top = ax.secondary_xaxis("top", functions=(lambda x: x * KM_PER_S_PER_EV, lambda d: d / KM_PER_S_PER_EV))
    # The decades, with the Earth's radius and diameter among them; a decade's label too close to
    # one of theirs is left out, its tick kept.
    marks = {label: multiple * EARTH_RADIUS_KM for label, multiple in EARTH_MARKS}
    decades = [10.0 ** n for n in range(-2, 9) if LEFT * KM_PER_S_PER_EV <= 10.0 ** n <= RIGHT * KM_PER_S_PER_EV]
    ticks = decades + list(marks.values())
    labels = ["" if any(abs(np.log10(d / m)) < MARK_CLEARANCE for m in marks.values())
              else rf"$10^{{{int(round(np.log10(d)))}}}$" for d in decades] + list(marks)
    top.set_xticks(ticks)
    top.set_xticklabels(labels)
    # As large as the lower axis's label.
    top.set_xlabel("Decay length at 2 MeV [km]", fontsize=plt.rcParams["axes.labelsize"], labelpad=18)
    top.tick_params(axis="x", which="major", pad=4)
    # The rc draws every axis's ticks on the top edge too; there the lifetime's
    # decades would sit 0.22 of a decade from the decay length's, so the top edge
    # carries the decay length's ticks alone.
    ax.tick_params(axis="x", which="both", top=False)
    return top, header_texts, lines


def finish(ax, header_texts, lines) -> None:
    """Place what is aligned with drawn sizes: the lines without a limit, and the legends.

    Each line naming the rows without a limit starts at the left border of the headers'
    rounded boxes.  The two legends keep their own sizes and stand in the empty lower right
    of the panel, below the normal-ordering header and right of its short bars, stacked and
    right-aligned.

    Args:
        ax: The axes, laid out.
        header_texts: The headers, whose boxes the lines start at.
        lines: Each line as ``(y, text)``.
    """
    from matplotlib.legend import Legend
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    left_px = header_texts[0].get_bbox_patch().get_window_extent(renderer).x0
    x = ax.transData.inverted().transform((left_px, 0.0))[0]
    for y, text in lines:
        ax.text(x, y, text, ha="left", va="center", fontsize=15, color="k", zorder=4)
    style = dict(frameon=True, edgecolor="black", framealpha=1.0, fancybox=False, fontsize=14, title_fontsize=15,
                 borderaxespad=0.0, borderpad=0.6)
    specs = (([Patch(color=c, alpha=DARK, lw=0) for _, _, c, _ in ROWS], [name for _, name, _, _ in ROWS],
              dict(title=r"Lifetime disfavored (90\% C.L.) by", ncol=2)),
             ([Line2D([], [], marker=m, ms=9, mfc="0.3" if treatment == "profiled" else "white", mec="0.3", mew=1.6,
                      ls="none") for _, _, m, treatment in EARTHS], [name for _, name, _, _ in EARTHS],
              dict(title="Earth composition model", ncol=2)))
    to_axes = ax.transAxes.inverted()

    def place(spec, anchor, **where):
        handles, labels, extra = spec
        legend = Legend(ax, handles, labels, bbox_to_anchor=anchor, bbox_transform=ax.transAxes, **where, **style,
                        **extra)
        legend.set_zorder(10)
        ax.add_artist(legend)
        return legend

    trial = [place(spec, (1.0, 2.0), loc="lower right") for spec in specs]
    fig.canvas.draw()
    sizes = []
    for legend in trial:
        box = legend.get_window_extent(renderer)
        (x0, y0), (x1, y1) = to_axes.transform([(box.x0, box.y0), (box.x1, box.y1)])
        sizes.append((x1 - x0, y1 - y0))
        legend.remove()
    # All as wide as the widest, stacked from the panel's lower right
    # corner up, the first legend on top.
    width = max(w for w, _ in sizes)
    margin, gap = LEGEND_MARGIN, LEGEND_GAP
    y = margin
    for spec, (_, height) in reversed(list(zip(specs, sizes))):
        place(spec, (1.0 - margin - width, y, width, 0.0), loc="lower left", mode="expand")
        y += height + gap


def figure():
    """The overview, in the repository's style.

    Returns:
        The matplotlib figure.
    """
    figure_style.use()
    fig, ax = plt.subplots(figsize=(14.5, 17.0))
    _top, header_texts, lines = draw(ax)
    fig.tight_layout()
    finish(ax, header_texts, lines)
    return fig


# The data, read from the release's data/ files.
from . import data as _data  # noqa: E402


def limit(site, case, param, bse, treatment="assumed"):
    """One 90% limit from data, from data/limits.json: 0.0 where the fit sets none, None where it was not run."""
    earth = "free_radiogenic_power" if treatment == "profiled" else bse
    try:
        value = _data.load("limits")["lower_limit_s_per_eV"][case][param.replace("log10_", "")][site][earth]["90"]["observed"]
    except KeyError:
        return None
    return 0.0 if value is None else value
