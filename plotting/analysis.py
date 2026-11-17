"""Stand-ins for the names of the paper's analysis package that the plotting modules use.

Constants and the functions ``band`` and ``coverage`` are copied from the paper's analysis; the data
functions return what the paper's returned when the figure was drawn (``plotting/data.py``).
"""

from types import SimpleNamespace

import numpy as np

from . import data as _data


def _stored(name):
    def call(*args, **kwargs):
        return getattr(_data.current(), name)(*args, **kwargs)
    call.__name__ = name
    return call


def band(records: list, cl: str = "90") -> dict:
    """The spread of simulated limits at one level, and why some have none.

    Quantiles are order statistics (``inverted_cdf``), so each is one of the
    simulated limits; a simulated experiment that sets no limit, its profile
    not reaching the threshold even at complete decay, counts as the weakest.

    Args:
        records: The simulated experiments' records.
        cl: The level, as the records key it.

    Returns:
        ``n``, ``median``, ``68`` and ``95`` (each a pair; None where the
        quantile falls on a simulated experiment without a limit), ``no_limit``,
        ``no_decay_excluded`` (no decay itself excluded at this level) and
        ``flagged`` (a fit that failed its convergence check).
    """
    values = np.array([r["limits_s_per_ev"][cl] or 0.0 for r in records])
    q = np.quantile(values, [0.025, 0.16, 0.5, 0.84, 0.975], method="inverted_cdf")
    shown = [None if v == 0.0 else float(v) for v in q]
    return {"n": len(records), "median": shown[2], "68": (shown[1], shown[3]), "95": (shown[0], shown[4]),
            "no_limit": int(np.sum(values == 0.0)),
            "no_decay_excluded": sum(bool(r.get("no_decay_excluded", {}).get(cl)) for r in records),
            "flagged": sum(not (r["global_fit_ok"] and r["profile_at_minimum"]) for r in records)}


def coverage(records: list, truth: float, cl: str = "90") -> tuple[int, int]:
    """How many simulated experiments' limits at one level exclude the true lifetime, and of how many.

    A lower limit excludes the truth when it lies above it.  A simulated
    experiment that sets no limit excludes no lifetime, the truth included.

    Args:
        records: The records of simulated experiments made at the truth.
        truth: log10 of the true tau/m, in s/eV.
        cl: The level, as the records key it.

    Returns:
        ``(excluded, n)``.
    """
    excluded = sum(1 for r in records if r["limits_s_per_ev"][cl] and np.log10(r["limits_s_per_ev"][cl]) > truth)
    return excluded, len(records)


BSE_PROFILE_RANGE_TW = (9.4, 39.6)
dl = SimpleNamespace(PAPER_MIXING='varied', COMBINATIONS=('combined_bx_kl_snoplus', 'combined_all'),
                     COVERAGE_EXPECTED={'68': 0.32, '90': 0.1, '95': 0.05}, band=band, coverage=coverage,
                     current=lambda record: True)
ss = SimpleNamespace(dl=dl, PARAM='log10_tau1_over_m1', LIFETIME_TOP=-6.0, CORNER_FITS=('borexino', 'kamland', 'snoplus', 'juno', 'combined_all'),
                     CORNER_PARAMETERS=4, use_mixing=lambda mixing: None,
                     use_treatment=lambda treatment: None, corner_parameters=_stored("corner_parameters"),
                     load=_stored("load"), best_fit=_stored("best_fit"))
rct = SimpleNamespace(FITS=(('borexino', 'Borexino'), ('kamland', 'KamLAND'), ('snoplus', 'SNO+'), ('juno', 'JUNO'), ('combined_bx_kl_snoplus', 'BX+KL+SNO+'), ('combined_all', 'All four')), CASE='io_n1_n3', PARAM='log10_tau1_over_m1', MIXING='varied',
                      BSE='geochemical', points=_stored("points"))
