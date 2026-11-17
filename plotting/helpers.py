"""The drawing helpers of the paper's figures, and constants that the figure cells use.

Copied from the paper's own code.
"""

import numpy as np
from scipy.stats import chi2

R_EARTH_KM = 6371.0
SCHEMATIC_SITE = (42.5, 13.6, 1.4)
CRUST_R_INNER_KM = 6336.0
LOCAL_BOXES = {'borexino': (40.0, 44.0, 10.0, 16.0), 'kamland': (34.0, 36.0, 134.0, 140.0, 36.0, 38.0, 136.0, 142.0), 'snoplus': (43.0, 47.0, -84.0, -78.0), 'juno': (20.0, 24.0, 109.0, 115.0)}
LOCAL_CRUST_MAX_KM = 350.0
MANTLE_R_INNER_KM = 3480.0


def poisson_interval(counts: np.ndarray, level: float = 0.6827) -> tuple[np.ndarray, np.ndarray]:
    """The central Poisson interval on each count (Garwood), as distances below and above it.

    Args:
        counts: Observed counts.
        level: Probability content.

    Returns:
        ``(below, above)``, each as long as ``counts``.
    """
    from scipy.stats import chi2
    n = np.asarray(counts, dtype=float)
    tail = 0.5 * (1.0 - level)
    low = np.where(n > 0, 0.5 * chi2.ppf(tail, 2.0 * np.maximum(n, 1e-12)), 0.0)
    high = 0.5 * chi2.ppf(1.0 - tail, 2.0 * n + 2.0)
    return n - low, high - n


def ternary_xy(points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Flavor fractions (f_e, f_mu, f_tau) to the plane of the ternary figure.

    The orientation is python-ternary's: f_e = 1 at the bottom right, f_mu = 1 at the top,
    f_tau = 1 at the bottom left, with sides of unit length.

    Args:
        points: Array of shape (n, 3), or (3,) for one point.

    Returns:
        The x and y coordinates.
    """
    p = np.atleast_2d(points)
    return p[:, 0] + 0.5 * p[:, 1], 0.5 * np.sqrt(3.0) * p[:, 1]


def globe_rings(rings: list, latitude_deg: float, longitude_deg: float) -> list:
    """Land outlines on a globe that shows the site at the top of its limb.

    An orthographic projection centred 90 degrees south of the site, on its meridian, so that
    the site sits at the top of the disk.  Points on the far side are pushed onto the limb, so
    that a ring that crosses it is closed along it.  Rings with fewer than three points on the
    near side are dropped.

    Args:
        rings: As :func:`land_outlines` returns them.
        latitude_deg: The site's latitude.
        longitude_deg: The site's longitude.

    Returns:
        One ``(n, 2)`` array of ``(x, y)`` per ring kept, in units of the Earth's radius.
    """
    lat0, lon0 = np.deg2rad(latitude_deg - 90.0), np.deg2rad(longitude_deg)
    out = []
    for ring in rings:
        a = np.deg2rad(np.asarray(ring, dtype=float))
        lon, lat = a[:, 0], a[:, 1]
        cosc = np.sin(lat0) * np.sin(lat) + np.cos(lat0) * np.cos(lat) * np.cos(lon - lon0)
        if (cosc > 0).sum() < 3:
            continue
        x = np.cos(lat) * np.sin(lon - lon0)
        y = np.cos(lat0) * np.sin(lat) - np.sin(lat0) * np.cos(lat) * np.cos(lon - lon0)
        far = cosc <= 0
        if far.any():
            n = np.hypot(x, y)
            n[n == 0] = 1.0
            x, y = np.where(far, x / n, x), np.where(far, y / n, y)
        out.append(np.column_stack([x, y]))
    return out


def production_point(depth_km: float, distance_km: float, detector_depth_km: float = SCHEMATIC_SITE[2]) -> tuple:
    """Where a production point sits in the plane of the cut, with the detector at the top.

    Args:
        depth_km: The production point's depth.
        distance_km: Its straight-line distance from the detector, the chord's length.
        detector_depth_km: The detector's depth.

    Returns:
        ``(x, y)`` in units of the Earth's radius, the Earth's centre at the origin and the
        detector at ``(0, 1 - detector_depth / R)``.
    """
    r_det, r_src = R_EARTH_KM - detector_depth_km, R_EARTH_KM - depth_km
    angle = np.arccos(np.clip((r_det ** 2 + r_src ** 2 - distance_km ** 2) / (2.0 * r_det * r_src), -1.0, 1.0))
    return r_src * np.sin(angle) / R_EARTH_KM, r_src * np.cos(angle) / R_EARTH_KM

