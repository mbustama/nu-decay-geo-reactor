"""Read the release's data files.

``data/limits.json`` holds every limit of the paper; :func:`limit` reads one. Every other file
holds the numbers of one figure: the results of the functions that the paper's own code calls to
draw it, each stored with the arguments it was called with. :func:`figure_data` returns an object
whose methods have the names of those functions and return the stored results, so that the
notebook draws each figure with the paper's own plotting code.

How the files store Python objects that JSON has no type for:

* an array is a list (``null`` for a value that is not a number);
* a dictionary whose keys are not all strings, such as one keyed by lifetimes, is
  ``{"__keys__": [...], "__values__": [...]}``, the keys and the values in the same order.
"""

from __future__ import annotations

import json
import os
from functools import lru_cache

import numpy as np

#: The data files, next to this package.
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")


@lru_cache(maxsize=None)
def _read(name: str) -> str:
    path = os.path.join(DATA, name + ".json")
    if not os.path.exists(path) and os.path.exists(path + ".gz"):      # the largest files are compressed
        import gzip
        with gzip.open(path + ".gz", "rt") as fh:
            return fh.read()
    with open(path) as fh:
        return fh.read()


def load(name: str) -> dict:
    """One data file, as stored.

    Args:
        name: The file's name without ``.json``, e.g. ``"limits"``.

    Returns:
        Its contents.
    """
    return json.loads(_read(name))


def _numbers(items: list) -> bool:
    """Whether a list holds numbers only (and ``null``, a value that is not a number)."""
    return bool(items) and all(isinstance(v, (int, float)) and not isinstance(v, bool) or v is None for v in items)


def decode(x):
    """A stored object as the paper's code returned it: lists of numbers as arrays, keyed dictionaries and tuples restored."""
    if isinstance(x, dict):
        if set(x) == {"__keys__", "__values__"}:
            return {_key(decode(k)): decode(v) for k, v in zip(x["__keys__"], x["__values__"])}
        if set(x) == {"__tuple__"}:
            return tuple(v.item() if isinstance(v, np.generic) else v for v in (decode(w) for w in x["__tuple__"]))
        return {k: decode(v) for k, v in x.items()}
    if isinstance(x, list):
        if _numbers(x):
            return np.array([np.nan if v is None else v for v in x], dtype=float)
        if x and all(isinstance(v, list) and _numbers(v) for v in x) and len({len(v) for v in x}) == 1:
            return np.array([[np.nan if w is None else w for w in v] for v in x], dtype=float)
        return [decode(v) for v in x]
    return x


def _key(k):
    """A dictionary key as it was: a list or an array back to a tuple, which can key a dictionary."""
    if isinstance(k, np.ndarray):
        return tuple(k.tolist())
    return tuple(_key(v) for v in k) if isinstance(k, (list, tuple)) else k


def encode(x):
    """The inverse of :func:`decode`, for arguments: what a stored call's arguments look like."""
    if isinstance(x, np.ndarray):
        return encode(x.tolist())
    if isinstance(x, np.generic):
        return x.item()
    if isinstance(x, dict):
        if all(isinstance(k, str) for k in x):
            return {k: encode(v) for k, v in x.items()}
        return {"__keys__": [encode(k) for k in x], "__values__": [encode(v) for v in x.values()]}
    if isinstance(x, tuple):
        return {"__tuple__": [encode(v) for v in x]}
    if isinstance(x, list):
        return [encode(v) for v in x]
    return x


#: Significant digits of the stored numbers: arguments are compared at this precision.
DIGITS = 6


def _rounded(x):
    if isinstance(x, float):
        return None if x != x or x in (float("inf"), float("-inf")) else float(f"{x:.{DIGITS}g}")
    if isinstance(x, dict):
        return {k: _rounded(v) for k, v in x.items()}
    if isinstance(x, list):
        return [_rounded(v) for v in x]
    return x


def call_key(function: str, args, kwargs) -> str:
    """The key a stored call is found by: the function and its arguments, as JSON."""
    return json.dumps(_rounded([function, encode(list(args)), encode(dict(kwargs))]), sort_keys=True)


class FigureData:
    """The stored calls of one figure, with the paper's function names as methods.

    ``FigureData("neutrino_spectrum").neutrino_spectrum("borexino", (-8.0, -9.0))`` returns what
    the paper's ``neutrino_spectrum("borexino", (-8.0, -9.0))`` returned when the figure was drawn.
    """

    def __init__(self, name: str):
        stored = load(name)
        self.name = name
        self.content = stored.get("content")
        self._calls = {json.dumps(_rounded([c["function"], c["args"], c["kwargs"]]), sort_keys=True): c["result"]
                       for c in stored["calls"]}
        self._constants = stored.get("constants", {})

    def __getattr__(self, function: str):
        if function.startswith("_"):
            raise AttributeError(function)
        if function in self._constants:
            return decode(self._constants[function])

        def stored(*args, **kwargs):
            key = call_key(function, args, kwargs)
            if key not in self._calls:
                raise KeyError(f"data/{self.name}.json holds no call {function}{args}{kwargs or ''}")
            return decode(self._calls[key])

        return stored


def figure_data(name: str) -> FigureData:
    """The stored calls of figure ``name``; see :class:`FigureData`."""
    return FigureData(name)


_CURRENT = [None]


def use(name: str) -> FigureData:
    """Make figure ``name``'s stored calls the ones the plotting modules read, and return them."""
    _CURRENT[0] = FigureData(name)
    return _CURRENT[0]


def current() -> FigureData:
    """The stored calls set by :func:`use`."""
    if _CURRENT[0] is None:
        raise RuntimeError("no figure's data in use: call plotting.data.use(name) first")
    return _CURRENT[0]


def limit(channel: str, lifetime: str, fit: str, earth: str, cl: str = "90", kind: str = "observed"):
    """One lower limit on tau_i/m_i from ``data/limits.json``.

    Args:
        channel: E.g. ``"io_n1_n3"``.
        lifetime: ``"tau1_over_m1"``, ``"tau2_over_m2"`` or ``"tau3_over_m3"``.
        fit: An experiment or a combination, e.g. ``"combined_all"``.
        earth: ``"cosmochemical"``, ``"geochemical"``, ``"geodynamical"`` or ``"free_radiogenic_power"``.
        cl: ``"90"`` or ``"95"``.
        kind: ``"observed"`` (from data) or ``"expected"`` (Asimov data set).

    Returns:
        The limit in s/eV, or None where that fit sets none.
    """
    return load("limits")["lower_limit_s_per_eV"][channel][lifetime][fit][earth][cl][kind]
