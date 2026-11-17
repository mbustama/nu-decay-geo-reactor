"""The style the paper's figures are drawn in: ``matplotlibrc`` in this folder.

Text is set with LaTeX in Palatino, so drawing the figures needs a LaTeX installation with the
psnfss fonts and the amssymb package.
"""

from __future__ import annotations

import os

import matplotlib

#: The matplotlibrc of the paper's figures.
RC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "matplotlibrc")


def use() -> None:
    """Draw every later figure in the paper's style."""
    matplotlib.rc_file(RC)
    # The AMS symbols, for the few labels that need one (a black triangle in the overview's headers).
    matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amssymb}"
