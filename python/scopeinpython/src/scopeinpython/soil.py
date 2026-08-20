"""Brightness-Shape-Moisture (BSM) soil reflectance model.

Direct, function-by-function port of ``SCOPEinR/R/BSM.R``.

References
----------
Original Matlab/R version: Christiaan van der Tol.
Verhoef & Bach (2007); Jacquemoud, Baret & Hanocq (1992).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import poisson
from toolsrtm import calctav

from ._data import BSMOptipar, bsm_optipar

__all__ = ["SoilParams", "WettingParams", "get_bsm", "soilwat"]


@dataclass
class SoilParams:
    """Soil brightness/shape parameters for BSM (``soilpar`` in the R code)."""

    BSMBrightness: float  # overall soil brightness
    BSMlat: float  # spectral shape "latitude", typical range 20-40 deg
    BSMlon: float  # spectral shape "longitude", typical range 45-65 deg


@dataclass
class WettingParams:
    """Empirical wetting parameters for BSM (``emp`` in the R code)."""

    SMp: float  # soil moisture volume percentage, range 5-55
    SMC: float = 25.0  # soil moisture capacity, recommended 0.25*100=25
    film: float = 0.015  # effective optical thickness of a single water film, recommended 0.015


def soilwat(rdry: np.ndarray, nw: np.ndarray, kw: np.ndarray, SMp: float, SMC: float, film: float) -> np.ndarray:
    """Wet soil reflectance from dry soil reflectance, following a Poisson
    process for the number of water film layers (Verhoef 2012; Yang 2020).

    Direct port of ``SCOPEinR::soilwat``.

    Parameters
    ----------
    rdry : array_like, shape (nwl,)
        Dry soil reflectance spectrum.
    nw : array_like, shape (nwl,)
        Refraction index of water spectrum.
    kw : array_like, shape (nwl,)
        Absorption coefficient of water spectrum.
    SMp : float
        Soil moisture volume percentage (5-55).
    SMC : float
        Soil moisture capacity (recommended 25, i.e. 0.25 in [0,1] units
        used consistently with the R default of 25).
    film : float
        Effective optical thickness of a single water film (recommended 0.015).

    Returns
    -------
    numpy.ndarray, shape (nwl,)
        Wet soil reflectance spectrum (equals ``rdry`` if ``mu <= 0``).
    """
    rdry = np.asarray(rdry, dtype=float)
    nw = np.asarray(nw, dtype=float)
    kw = np.asarray(kw, dtype=float)

    k = np.arange(0, 7)
    nk = len(k)

    mu = (SMp - 5) / SMC

    if mu <= 0:
        return rdry.copy()

    # Lekner & Dorf (1988) modified soil background reflectance for soil
    # refraction index = 2.0; uses the tav-function of PROSPECT (calctav)
    rbac = 1 - (1 - rdry) * (rdry * calctav(90, 2.0 / nw) / calctav(90, 2.0) + 1 - rdry)

    # total reflectance at bottom of water film surface: rho21, water to air, diffuse
    p = 1 - calctav(90, nw) / nw**2

    # reflectance of water film top surface, 40 deg incidence, as in PROSPECT: rho12, air to water, direct
    Rw = 1 - calctav(40, nw)

    # probability of k water film layers (Poisson), P(0)=dry area, P(1)=single film, ...
    # NB: mirrors R's `stats::dpois(round(mu), k)` exactly -- note the R call evaluates
    # the Poisson pmf AT x=round(mu) as k (0:6) is varied as the RATE parameter, i.e.
    # scipy's poisson.pmf(round(mu), mu=k), not poisson.pmf(k, mu=round(mu)).
    fmul = poisson.pmf(round(mu), k)

    # two-way transmittance, exp(-2*kw*k*film)
    tw = np.exp(-2 * np.outer(kw, film * k))  # (nwl, nk)

    Rwet_k = Rw[:, None] + (1 - Rw[:, None]) * (1 - p[:, None]) * tw * rbac[:, None] / (
        1 - p[:, None] * tw * rbac[:, None]
    )
    rwet = rdry * fmul[0] + Rwet_k[:, 1:nk] @ fmul[1:nk]

    return rwet


def get_bsm(soilpar: SoilParams, emp: WettingParams, optipar: BSMOptipar | None = None) -> np.ndarray:
    """Brightness-Shape-Moisture soil reflectance spectrum (400-2400 nm).

    Direct port of ``SCOPEinR::getBSM``.

    Parameters
    ----------
    soilpar : SoilParams
        BSMBrightness, BSMlat, BSMlon.
    emp : WettingParams
        SMp, SMC, film.
    optipar : BSMOptipar, optional
        GSV/Kw/nw spectral basis. Defaults to the bundled
        ``optipar2017.ProspectD`` subset (:func:`scopeinpython._data.bsm_optipar`).

    Returns
    -------
    numpy.ndarray, shape (2001,)
        Wet soil reflectance spectrum, 400-2400 nm (1 nm step).
    """
    if optipar is None:
        optipar = bsm_optipar()

    GSV = optipar.GSV
    kw = optipar.Kw
    nw = optipar.nw

    B = soilpar.BSMBrightness
    lat = soilpar.BSMlat
    lon = soilpar.BSMlon

    f1 = B * np.sin((np.pi / 180) * lat)
    f2 = B * np.cos((np.pi / 180) * lat) * np.sin((np.pi / 180) * lon)
    f3 = B * np.cos((np.pi / 180) * lat) * np.cos((np.pi / 180) * lon)

    rdry = f1 * GSV[:, 0] + f2 * GSV[:, 1] + f3 * GSV[:, 2]

    rwet = soilwat(rdry, nw, kw, emp.SMp, emp.SMC, emp.film)
    return rwet
