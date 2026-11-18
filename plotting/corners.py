# The paper's plotting code, its functions that read the analysis's results replaced by ones that read
# the release's data/ files (at the end of this file).
"""Figure E.3 of the paper: how the lifetime trades off against the parameters that move most with it.

A row per fit: the test statistic against tau_1/m_1, then the lifetime against each of the four
parameters that move most with it.  A pair's 68% and 95% C.L. regions come from its profiled
likelihood on a grid of 13 x 13 values, every other parameter fitted at each point, drawn through a
spline that passes through every grid point; the cross marks the data's best fit.

The notebook draws the figure with :func:`figure_rows`.
"""

from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FuncFormatter, MaxNLocator
from scipy.interpolate import CubicSpline, RectBivariateSpline
from scipy.stats import chi2


from . import style as figure_style
from . import names as pn
from . import coverage as plot_coverage
from . import impacts as plot_impacts
from .analysis import ss

#: The records drawn: the study's, with the mixing the paper is drawn from (``dl.PAPER_MIXING``).
MIXING = ss.dl.PAPER_MIXING
ss.use_mixing(MIXING)
#: Panels in a row of :func:`figure_rows`: the lifetime's profile and its four pairs.
CORNER_COLUMNS = 1 + ss.CORNER_PARAMETERS
#: The NLL above its minimum at 68% and 95% C.L., for two parameters and for one.
LEVELS_2D = tuple(0.5 * chi2.ppf(cl, 2) for cl in (0.68, 0.95))
LEVELS_1D = tuple(0.5 * chi2.ppf(cl, 1) for cl in (0.68, 0.95))
#: Points per axis of the spline the regions are drawn through.
FINE = 241
#: The regions' opacity, 68% then 95%, and the diagonal's lines at the same two levels.
ALPHAS = (0.7, 0.28)
STYLES_1D = ("--", "-.")
#: Each fit's color, as the overview and the coverage figure draw it.
COLORS = plot_coverage.COLORS
#: Each fit's name, as the impact figures title it.
TITLES = plot_impacts.TITLES


def records(site: str) -> tuple[list[str], dict[str, dict]]:
    """A fit's corner records.

    Args:
        site: One of ``ss.CORNER_FITS``.

    Returns:
        ``(names, panels)``: the plot's parameters, the lifetime first, and each
        panel's record by its subject, ``a`` or ``a|b``.

    Raises:
        FileNotFoundError: If a panel's record is not on file.
    """
    names = ss.corner_parameters(site, "real")
    subjects = names + [f"{a}|{b}" for i, a in enumerate(names) for b in names[i + 1:]]
    panels = {}
    for subject in subjects:
        record = ss.load("corner", site, "real", subject)
        if record is None:
            raise FileNotFoundError(f"{site}: the corner panel {subject} is not on file")
        panels[subject] = record
    return names, panels


def scale(name: str) -> tuple[float, str]:
    """How an axis is printed: a factor to multiply its values by, and a unit for its label.

    Args:
        name: A parameter's name in the fit.

    Returns:
        The factor and the unit's text, empty for a quantity of order one.
    """
    base = name.partition("@")[0]
    if base in ("osc_dm2_21", "dm2_21_ev2"):
        return 1e5, r" [$10^{-5}$ eV$^2$]"
    if base in ("osc_dm2_32", "dm2_31_ev2"):
        return 1e3, r" [$10^{-3}$ eV$^2$]"
    return 1.0, ""


def ticks(axis: list[float], name: str):
    """An axis's major ticks, none at either end, where the next panel's labels sit.

    The lifetime's at every half or whole unit of log10 tau/m; the others at
    about four round values.

    Args:
        axis: The grid along the axis.
        name: The parameter.

    Returns:
        A matplotlib locator.
    """
    factor = scale(name)[0]
    lo, hi = axis[0] * factor, axis[-1] * factor
    margin = 0.06 * (hi - lo)
    if name == ss.PARAM:
        step = 0.5 if hi - lo < 3.0 else 1.0
        values = np.arange(np.ceil(lo / step) * step, hi, step)
    else:
        values = MaxNLocator(4).tick_values(lo, hi)
    return FixedLocator([v for v in values if lo + margin <= v <= hi - margin])


def axis_label(name: str, combined: bool) -> str:
    """A parameter's axis label: its plain name and its unit.

    Args:
        name: A parameter's name in the fit.
        combined: Whether the fit is a combination.

    Returns:
        The label.
    """
    return pn.label(name, combined) + scale(name)[1]


def pair(ax, record: dict, x_name: str, y_name: str, color: str, best: dict, edge: bool = False) -> None:
    """One pair's regions at 68% and 95% C.L., through a spline of its grid.

    Args:
        ax: The axes.
        record: The panel's record, its grid over ``x_name`` then ``y_name``.
        x_name: The parameter along the horizontal axis, the record's first.
        y_name: The parameter along the vertical axis, the record's second.
        color: The fit's color.
        best: The data's best fit, every parameter by name.
        edge: Whether a best fit beyond the end of the horizontal axis, a lifetime
            without decay, is marked on the axis's end, at its value of the other.
    """
    fx, fy = scale(x_name)[0], scale(y_name)[0]
    x = np.array(record["axes"][x_name]) * fx
    y = np.array(record["axes"][y_name]) * fy
    z = np.array(record["nll"], dtype=float) - record["global_nll"]
    spline = RectBivariateSpline(x, y, z, kx=3, ky=3, s=0)
    gx, gy = np.linspace(x[0], x[-1], FINE), np.linspace(y[0], y[-1], FINE)
    fine = spline(gx, gy).T
    # The 68% region darker, the 95% one lighter: two fills, the inner over the outer.
    ax.contourf(gx, gy, fine, levels=[-np.inf, LEVELS_2D[1]], colors=[color], alpha=ALPHAS[1])
    ax.contourf(gx, gy, fine, levels=[-np.inf, LEVELS_2D[0]], colors=[color], alpha=ALPHAS[0])
    ax.contour(gx, gy, fine, levels=list(LEVELS_2D), colors=[color], linewidths=1.0)
    bx, by = best.get(x_name), best.get(y_name)
    if bx is not None and by is not None and edge and bx * fx > x[-1]:
        # No decay: drawn on the right edge, unclipped, where the lifetime axis stops acting.
        ax.plot(x[-1], by * fy, "x", color="k", ms=CROSS, mew=2.0, clip_on=False, zorder=10)
    elif bx is not None and by is not None and x[0] <= bx * fx <= x[-1] and y[0] <= by * fy <= y[-1]:
        ax.plot(bx * fx, by * fy, "x", color="k", ms=CROSS if edge else 7, mew=2.0 if edge else 1.6)


#: The rows' font sizes: tick labels, axis labels, the fits' names, the levels' names, the common label;
#: and the best fit's cross.  The figure is 18 inches wide and printed at the paper's width, about 0.39 of it.
ROW_FONTS = {"ticks": 17, "labels": 17, "names": 20, "levels": 15, "common": 24}
CROSS = 10
#: The regions' names' font size; how far, in axes fraction, they keep from the panel's frame, and from a
#: marker's center (its half-size, about 0.04, and a gap).
REGION_FONT = 14
FRAME, CLEAR = 0.05, 0.08
#: Axis labels longer than this many characters go on two lines, at the comma, colon or parenthesis
#: nearest their middle.
WRAP = 26
#: The lifetime axis the rows of :func:`figure_rows` share, log10 tau/m: from below the lowest start of
#: any fit's grid (KamLAND's, at -9.63) to where the grids end.
ROWS_LIFETIME = (-9.7, ss.LIFETIME_TOP)
#: Each fit's name as its row is tagged; the combination's on two lines, to fit its panel.
ROW_NAMES = {"borexino": "Borexino (2020)", "kamland": "KamLAND (2022)", "snoplus": "SNO+ (2026)",
             "juno": "JUNO (2026)", "combined_all": "Borexino, KamLAND,\nSNO+, and JUNO"}


def wrapped(text: str) -> str:
    """A long label on two lines, broken at the comma, colon or parenthesis nearest its middle.

    Args:
        text: The label.

    Returns:
        The label, on two lines if it is longer than :data:`WRAP` characters and has a place to break;
        a break at a comma or colon keeps it, one before a parenthesis moves the parenthesis down.
    """
    if len(text) <= WRAP:
        return text
    breaks = [i + 1 for i, c in enumerate(text) if c in ",:" and text[i + 1:i + 2] == " "]
    breaks += [i for i, c in enumerate(text) if c == "(" and text[i - 1:i] == " "]
    if not breaks:
        return text
    i = min(breaks, key=lambda b: abs(b - len(text) / 2))
    return text[:i].rstrip() + "\n" + text[i:].lstrip()


def lifetime_profile(ax, record: dict, profile: dict, color: str, labels: tuple[str, ...], side: str = "left") -> None:
    """The lifetime's test statistic, q = 2 Delta NLL, with the levels its one-sided limits are read at.

    Args:
        ax: The axes.
        record: The lifetime's diagonal panel.
        profile: The fit's ``profile`` record, which holds the thresholds.
        color: The fit's color.
        labels: The levels to name on the panel, among ``"68"``, ``"90"`` and ``"95"``.
        side: Where the names go, ``"left"`` or ``"right"``: the side the profile leaves clear.
    """
    x = np.array(record["axes"][ss.PARAM])
    y = 2.0 * (np.array(record["nll"], dtype=float) - record["global_nll"])
    fine = np.linspace(x[0], x[-1], FINE)
    ax.plot(fine, CubicSpline(x, y)(fine), color=color, lw=2.4)
    for level, style in (("68", STYLES_1D[0]), ("90", "-"), ("95", STYLES_1D[1])):
        threshold = 2.0 * profile["thresholds_nll"][level]
        ax.axhline(threshold, color="0.45", lw=1.1 if level != "90" else 1.3, ls=style)
        if level in labels:
            # Just above its line, on the side the profile leaves clear.
            ax.text(0.03 if side == "left" else 0.97, threshold, f"${level}\\%$ C.L.", transform=ax.get_yaxis_transform(),
                    ha=side, va="bottom", fontsize=ROW_FONTS["levels"], color="0.25")
    ax.set_ylim(0.0, 1.15 * float(np.nanmax(y)))
    ax.yaxis.set_major_locator(MaxNLocator(4, integer=True))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.set_ylabel(r"Test statistic, $q$", fontsize=ROW_FONTS["labels"])


def place_label(ax, gx: np.ndarray, gy: np.ndarray, zone: np.ndarray, texts: tuple[str, ...], color: str,
                avoid: tuple = ()) -> None:
    """Write a label where its whole box fits inside a zone of a panel, as far from the zone's edges as it can.

    Each wording is tried in turn, the first that fits wins.  For a wording, every candidate center is tried
    on a grid over the panel; one is kept only if the label's box, with a margin, lies inside the zone and
    clear of the points to avoid.  Of those, the one deepest inside the zone, measured in the panel's own
    proportions, wins.  The panel's ranges must be final.

    Args:
        ax: The axes.
        gx: The horizontal grid the zone is sampled on, in data units.
        gy: The vertical grid, in data units.
        zone: Booleans over ``(gy, gx)``: where the label may go.
        texts: The wordings, most preferred first, e.g. on one line and then on two.
        color: The label's color.
        avoid: Points in axes fraction the box must keep clear of, such as a marker.

    Raises:
        ValueError: If no wording fits anywhere.
    """
    from scipy.interpolate import RegularGridInterpolator
    from scipy.ndimage import distance_transform_edt
    renderer = ax.figure.canvas.get_renderer()
    frame = ax.get_window_extent(renderer)
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    u, v = (gx - x0) / (x1 - x0), (gy - y0) / (y1 - y0)
    inside = RegularGridInterpolator((v, u), zone.astype(float), bounds_error=False, fill_value=0.0)
    # Depth in the panel's proportions: a step across is frame.width wide, a step up frame.height tall.
    depth = RegularGridInterpolator((v, u), distance_transform_edt(
        zone, sampling=(np.mean(np.diff(v)) * frame.height, np.mean(np.diff(u)) * frame.width)),
        bounds_error=False, fill_value=0.0)
    margin_u, margin_v = 0.025, 0.022
    for text in texts:
        label = ax.text(0.5, 0.5, text, transform=ax.transAxes, ha="center", va="center", ma="center",
                        color=color, fontsize=REGION_FONT, linespacing=1.0)
        box = label.get_window_extent(renderer)
        w, h = box.width / frame.width, box.height / frame.height
        best, at = -1.0, None
        for vc in np.linspace(0.04, 0.96, 93):
            for uc in np.linspace(0.04, 0.96, 93):
                us = np.linspace(uc - w / 2 - margin_u, uc + w / 2 + margin_u, 9)
                vs = np.linspace(vc - h / 2 - margin_v, vc + h / 2 + margin_v, 5)
                if inside(np.array([(b, a) for b in vs for a in us])).min() < 0.999:
                    continue
                # Clear of the panel's frame, and of each marker by its own size and a gap.
                if not (FRAME <= uc - w / 2 and uc + w / 2 <= 1 - FRAME and FRAME <= vc - h / 2 and vc + h / 2 <= 1 - FRAME):
                    continue
                if any(abs(uc - au) < w / 2 + CLEAR and abs(vc - av) < h / 2 + CLEAR for au, av in avoid):
                    continue
                score = float(depth([(vc, uc)])[0])
                if score > best:
                    best, at = score, (uc, vc)
        if at is not None:
            label.set_position(at)
            return
        label.remove()
    raise ValueError(f"no place inside its zone for any of {texts!r}")


def label_regions(ax, record: dict, x_name: str, y_name: str, best: dict, texts: tuple[str, ...] = ("one",)) -> bool:
    """Name a pair's two regions inside themselves: the 68% one in white, the 95% one in its band.

    Args:
        ax: The axes, its ranges final.
        record: The panel's record, its grid over ``x_name`` then ``y_name``.
        x_name: The parameter along the horizontal axis.
        y_name: The parameter along the vertical axis.
        best: The data's best fit, whose cross the labels keep clear of.
        texts: The layouts to try, ``"one"`` line and then ``"two"``, in order.

    Returns:
        Whether both names fit; where one does not, neither is drawn.
    """
    fx, fy = scale(x_name)[0], scale(y_name)[0]
    x = np.array(record["axes"][x_name]) * fx
    y = np.array(record["axes"][y_name]) * fy
    z = np.array(record["nll"], dtype=float) - record["global_nll"]
    spline = RectBivariateSpline(x, y, z, kx=3, ky=3, s=0)
    gx, gy = np.linspace(x[0], x[-1], FINE), np.linspace(y[0], y[-1], FINE)
    fine = spline(gx, gy).T
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    bx = min(best[x_name] * fx, x[-1])
    avoid = (((bx - x0) / (x1 - x0), (best[y_name] * fy - y0) / (y1 - y0)),)
    drawn = len(ax.texts)
    try:
        for level, zone, color in (("68", fine <= LEVELS_2D[0], "white"),
                                   ("95", (fine > LEVELS_2D[0]) & (fine <= LEVELS_2D[1]), "0.15")):
            wordings = {"one": f"${level}\\%$ C.L.", "two": f"${level}\\%$\nC.L."}
            place_label(ax, gx, gy, zone, tuple(wordings[t] for t in texts), color, avoid)
    except ValueError:
        for text in ax.texts[drawn:]:
            text.remove()
        return False
    return True


def figure_rows(sites: tuple = None):
    """The left column of each fit's corner plot, as one row per fit.

    Each row is the fit's test statistic for the lifetime, then the lifetime
    against each of the fit's four parameters that move most with it, every
    panel on the same lifetime axis, so that the fits read down a column.  The
    test statistic carries the levels its one-sided limits are read at (68%,
    90% and 95% C.L.); a pair's regions are those of two parameters.  The
    cross is the best fit; where the data prefer no decay, it sits on the
    right edge.

    Args:
        sites: The fits, one per row; ``ss.CORNER_FITS`` if None.

    Returns:
        The matplotlib figure.
    """
    figure_style.use()
    sites = tuple(sites or ss.CORNER_FITS)
    n = len(sites)
    fig, axes = plt.subplots(n, CORNER_COLUMNS, figsize=(18.0, 3.6 * n), sharex=True,
                             gridspec_kw={"wspace": 0.68, "hspace": 0.12})
    for k, (row, site) in enumerate(zip(axes, sites)):
        names, panels = records(site)
        combined = site in ss.dl.COMBINATIONS
        color = COLORS[site]
        best = ss.best_fit(site, "real")
        profile = ss.load("profile", site, "real")
        # The top row names the levels, at the left, where its profile has not begun; the rows below carry
        # the same lines unnamed.
        lifetime_profile(row[0], panels[ss.PARAM], profile, color, ("68", "90", "95") if k == 0 else ())
        # The fit's name to the left of its row, outside the panels, in a rounded box.
        row[0].annotate(ROW_NAMES[site], xy=(-0.55, 0.5), xycoords="axes fraction", rotation=90, ha="center",
                        va="center", fontsize=ROW_FONTS["names"], ma="center",
                        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="k", lw=0.8))
        for ax, name in zip(row[1:], names[1:]):
            pair(ax, panels[f"{ss.PARAM}|{name}"], ss.PARAM, name, color, best, edge=True)
            y = np.array(panels[name]["axes"][name]) * scale(name)[0]
            ax.set_ylim(y[0], y[-1])
            ax.yaxis.set_major_locator(ticks(panels[name]["axes"][name], name))
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))
            ax.set_ylabel(wrapped(axis_label(name, combined)), fontsize=ROW_FONTS["labels"])
        for ax in row:
            ax.tick_params(labelsize=ROW_FONTS["ticks"])
    for ax in axes[-1]:
        ax.set_xlim(*ROWS_LIFETIME)
        ax.xaxis.set_major_locator(FixedLocator([-9.0, -8.0, -7.0, -6.0]))
        ax.xaxis.set_minor_locator(FixedLocator([d + np.log10(k) for d in range(-10, -6) for k in range(2, 10)]))
        # Every decade ticked; the last, at the axis's end, unlabeled, so that the labels do not crowd.
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _pos: "" if v > -6.5 else f"$10^{{{int(round(v))}}}$"))
    # The regions named inside themselves, once every range is final: in the first panel, reading row by
    # row, where both names fit on one line with room around them; on two lines if none has the room.
    done = False
    for texts in (("one",), ("two",)):
        for k, site in enumerate(sites):
            names, panels = records(site)
            for j, name in enumerate(names[1:], start=1):
                if label_regions(axes[k, j], panels[f"{ss.PARAM}|{name}"], ss.PARAM, name,
                                 ss.best_fit(site, "real"), texts):
                    done = True
                    break
            if done:
                break
        if done:
            break
    # One label for every column, centered on the panels rather than the figure (the row names widen it on
    # the left), set just below the tick labels.
    left, right = axes[-1, 0].get_position().x0, axes[-1, -1].get_position().x1
    below = axes[-1, 0].get_position().y0 - 0.5 / fig.get_figheight()
    fig.supxlabel(r"Neutrino lifetime ($\nu_1 \to \nu_3$, IO), $\tau_1/m_1$ [s~eV$^{-1}$]",
                  fontsize=ROW_FONTS["common"], x=0.5 * (left + right), y=below, va="top")
    return fig


# The data, read from the release's data/ files.
from . import data as _data  # noqa: E402
