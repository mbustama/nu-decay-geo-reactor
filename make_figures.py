"""Redraw every figure of the paper: run notebooks/figures.ipynb in place, its outputs stored.

    python make_figures.py

The figures go to figures/ as PDF files.
"""

import os

import nbformat
from nbclient import NotebookClient

HERE = os.path.dirname(os.path.abspath(__file__))
NOTEBOOK = os.path.join(HERE, "notebooks", "figures.ipynb")


def main() -> int:
    nb = nbformat.read(NOTEBOOK, as_version=4)
    NotebookClient(nb, timeout=1200, resources={"metadata": {"path": os.path.dirname(NOTEBOOK)}}).execute()
    nbformat.write(nb, NOTEBOOK)
    print("redrew every figure into", os.path.join(HERE, "figures"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
