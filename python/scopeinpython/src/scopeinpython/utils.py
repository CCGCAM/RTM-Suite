"""Small numerical utilities shared across the RTMo pipeline.

Direct ports of ``SCOPEinR/R/Sint.R``, ``e2phot.R``, ``ephoton.R`` and
``Planck.R``.
"""
from __future__ import annotations

import numpy as np

from ._data import constants as _constants

__all__ = ["sint", "get_ephoton", "get_e2phot", "get_planck", "satvap"]


def sint(y: np.ndarray, x: np.ndarray) -> float:
    """Trapezoidal integration of ``y`` over ``x``.

    Direct port of ``SCOPEinR::Sint``.
    """
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    step = np.diff(x)
    mean = 0.5 * (y[:-1] + y[1:])
    return float(np.sum(mean * step))


def get_ephoton(lambda_m: np.ndarray, const: dict | None = None) -> np.ndarray:
    """Energy content (J) of one photon at wavelength(s) ``lambda_m`` (m).

    Direct port of ``SCOPEinR::get.ephoton``.
    """
    const = const or _constants()
    h = const["h"]
    c = const["c"]
    return h * c / np.asarray(lambda_m, dtype=float)


def get_e2phot(lambda_m: np.ndarray, E: np.ndarray, const: dict | None = None) -> np.ndarray:
    """Number of moles of photons corresponding to E Joules of energy at
    wavelength(s) lambda (m).

    Direct port of ``SCOPEinR::get.e2phot``.
    """
    const = const or _constants()
    e = get_ephoton(lambda_m, const)
    photons = np.asarray(E, dtype=float) / e
    return photons / const["A"]


def get_planck(wl_nm: np.ndarray, Tb: np.ndarray, em: np.ndarray | None = None) -> np.ndarray:
    """Blackbody (or greybody) spectral radiance, Planck's law.

    Direct port of ``SCOPEinR::get.Planck``.

    Parameters
    ----------
    wl_nm : array_like
        Wavelength(s), nm.
    Tb : array_like
        Temperature(s), K.
    em : array_like, optional
        Emissivity (default 1).
    """
    wl_nm = np.asarray(wl_nm, dtype=float)
    Tb = np.asarray(Tb, dtype=float)
    c1 = 1.191066e-22
    c2 = 14388.33
    if em is None:
        em = 1.0
    Lb = em * c1 * (wl_nm * 1e-9) ** (-5) / (np.exp(c2 / (wl_nm * 1e-3 * Tb)) - 1)
    return Lb


def satvap(Temp: np.ndarray) -> np.ndarray:
    """Saturated vapour pressure at temperature ``Temp`` (deg C), in hPa/mbar.

    Direct port of ``SCOPEinR::satvap``.
    """
    Temp = np.asarray(Temp, dtype=float)
    b = 237.3
    return 6.107 * 10 ** (7.5 * Temp / (b + Temp))
