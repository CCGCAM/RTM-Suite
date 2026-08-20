"""Scalar-per-timestep thermal/aerodynamic building blocks used by the
SCOPE energy-balance loop (:mod:`scopeinpython.ebal`, not yet ported).

Direct ports of ``SCOPEinR/R/Monin_ObuKhov.R`` (``get.Monin.Obukhov``),
``SCOPEinR/R/resistances.R`` (``get.resistances`` + its ``psim``/``psih``/
``phstar`` stability-correction helpers), and ``SCOPEinR/R/heatfluxes.R``
(``get.heatfluxes``). Unlike the spectral/canopy RTM modules elsewhere in
this package, these operate on scalars (one timestep/pixel at a time),
matching the R source exactly.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._data import constants as _constants

__all__ = [
    "stefan_boltzmann",
    "monin_obukhov",
    "psim",
    "psih",
    "phstar",
    "ResistanceParams",
    "Resistances",
    "get_resistances",
    "HeatFluxes",
    "get_heatfluxes",
]


def stefan_boltzmann(T_C: np.ndarray) -> np.ndarray:
    """Blackbody radiant exitance (W/m2) at temperature(s) ``T_C``
    (degrees Celsius), via the Stefan-Boltzmann law. Direct port of
    ``SCOPEinR::get.Stefan_Boltzmann``."""
    const = _constants()
    C2K, sigmaSB = const["C2K"], const["sigmaSB"]
    return sigmaSB * (np.asarray(T_C, dtype=float) + C2K) ** 4


def monin_obukhov(ustar: float, Ta: float, H: float) -> float:
    """Monin-Obukhov length (m), a stability parameter used to correct
    aerodynamic resistances above the canopy for non-neutral conditions.

    Direct port of ``SCOPEinR::get.Monin.Obukhov``. Returns ``-1e6``
    (near-neutral stability) where the raw computation is undefined (e.g.
    ``H == 0``), matching R's ``L[is.na(L)] <- -1e6``.
    """
    const = _constants()
    rhoa, cp, kappa, g = const["rhoa"], const["cp"], const["kappa"], const["g"]
    with np.errstate(divide="ignore", invalid="ignore"):
        L = -rhoa * cp * ustar**3 * (Ta + 273.15) / (kappa * g * H)
    L = np.asarray(L, dtype=float)
    L = np.where(np.isnan(L), -1e6, L)
    return float(L) if L.ndim == 0 else L


def psim(z: float, L: float, unstable: bool, stable: bool, x: float) -> float:
    """Stability correction function for momentum transfer (Paulson 1970).
    Direct port of ``SCOPEinR::get.psim``. ``0`` under neutral conditions."""
    if unstable:
        return 2 * np.log((1 + x) / 2) + np.log((1 + x**2) / 2) - 2 * np.arctan(x) + np.pi / 2
    if stable:
        return -5 * z / L
    return 0.0


def psih(z: float, L: float, unstable: bool, stable: bool, x: float) -> float:
    """Stability correction function for heat transfer (Paulson 1970).
    Direct port of ``SCOPEinR::get.psih``."""
    if unstable:
        return 2 * np.log((1 + x**2) / 2)
    if stable:
        return -5 * z / L
    return 0.0


def phstar(z: float, zR: float, d: float, L: float, stable: bool, unstable: bool, x: float) -> float:
    """Stability correction function for the roughness sublayer (Paulson
    1970). Direct port of ``SCOPEinR::get.phstar``."""
    if unstable:
        return (z - d) / (zR - d) * (x**2 - 1) / (x**2 + 1)
    if stable:
        return -5 * z / L
    return 0.0


@dataclass
class ResistanceParams:
    """Inputs to :func:`get_resistances` (``data.soil``/``data.canopy``/
    ``data.meteo`` subsets in R)."""

    rbs: float  # soil boundary-layer resistance, s/m
    Cd: float  # leaf drag coefficient
    LAI: float
    rwc: float  # within-canopy aerodynamic resistance, s/m
    z0m: float  # roughness length for momentum, m (`zo` in R)
    d: float  # zero-plane displacement height, m
    hc: float  # vegetation height, m
    leafwidth: float
    z: float  # measurement height, m
    u: float  # wind speed at z, m/s
    L: float  # Monin-Obukhov length, m


@dataclass
class Resistances:
    ustar: float
    uz0: float
    Kh: float
    rai: float
    rar: float
    rac: float
    rws: float
    raa: float
    rawc: float
    raws: float


def get_resistances(p: ResistanceParams) -> Resistances:
    """Aerodynamic/boundary-layer resistances between soil, canopy and
    reference height (Wallace & Verhoef 2000 two-layer scheme, with a
    Monin-Obukhov stability correction). Direct port of
    ``SCOPEinR::get.resistances``.
    """
    const = _constants()
    kappa = const["kappa"]

    Cd, LAI, rwc = p.Cd, p.LAI, p.rwc
    z0m, d, h, w = p.z0m, p.d, p.hc, p.leafwidth
    z = p.z
    u = max(0.3, p.u)
    L = p.L
    rbs = p.rbs

    zr = 2.5 * h
    n = Cd * LAI / (2 * kappa**2)

    unstable = bool(L < 0 and L > -500)
    stable = bool(L > 0 and L < 500)
    x = (1 - 16 * z / L) ** 0.25 if unstable else float("nan")

    pm_z = psim(z - d, L, unstable, stable, x)
    ph_z = psih(z - d, L, unstable, stable, x)
    pm_h = psim(h - d, L, unstable, stable, x)

    ph_zr = psih(zr - d, L, unstable, stable, x) if z >= zr else ph_z

    phs_zr = phstar(zr, zr, d, L, stable, unstable, x)
    phs_h = phstar(h, zr, d, L, stable, unstable, x)

    ustar = max(0.001, kappa * u / (np.log((z - d) / z0m) - pm_z))
    Kh = kappa * ustar * (zr - d)

    if unstable:
        Kh_out = Kh * (1 - 16 * (h - d) / L) ** 0.5
    elif stable:
        Kh_out = Kh * (1 + 5 * (h - d) / L) ** -1
    else:
        Kh_out = Kh

    uh1 = ustar / kappa * (np.log((h - d) / z0m) - pm_h)
    uh1 = 0.0 if np.isnan(uh1) else uh1
    uh = max(uh1, 0.01)
    uz0 = uh * np.exp(n * ((z0m + d) / h - 1))

    rai = (1 / (kappa * ustar) * (np.log((z - d) / (zr - d)) - ph_z + ph_zr)) if z > zr else 0.0
    rar = 1 / (kappa * ustar) * ((zr - h) / (zr - d)) - phs_zr + phs_h
    # NB: R uses the RAW (stability-uncorrected) `Kh` here, not the
    # stability-corrected value returned as the `Kh` output field below --
    # confirmed by direct reading of SCOPEinR::get.resistances (the local
    # `Kh` variable is never reassigned after the corrected value is
    # written into `resist_out[['Kh']]`). Reproduced exactly, not "fixed",
    # since there's no independent way to tell if this is intentional.
    rac = (h * np.sinh(n) / (n * Kh)
           * (np.log((np.exp(n) - 1) / (np.exp(n) + 1))
              - np.log((np.exp(n * (z0m + d) / h) - 1) / (np.exp(n * (z0m + d) / h) + 1))))
    rws = (h * np.sinh(n) / (n * Kh)
           * (np.log((np.exp(n * (z0m + d) / h) - 1) / (np.exp(n * (z0m + d) / h) + 1))
              - np.log((np.exp(n * 0.01 / h) - 1) / (np.exp(n * 0.01 / h) + 1))))

    raa = rai + rar + rac
    rawc = rwc
    raws = rws + rbs

    return Resistances(
        ustar=ustar, uz0=uz0, Kh=Kh_out, rai=rai, rar=rar, rac=rac, rws=rws,
        raa=raa, rawc=rawc, raws=raws,
    )


@dataclass
class HeatFluxes:
    lambda_: np.ndarray  # latent heat of vaporization, J/kg
    s: np.ndarray  # slope of the saturated vapour pressure curve, hPa/degC
    lE: np.ndarray  # latent heat flux, W/m2
    H: np.ndarray  # sensible heat flux, W/m2
    ec: np.ndarray  # vapour pressure at the leaf surface, hPa
    Cc: np.ndarray  # CO2 concentration at the leaf surface, umol/m3


def get_heatfluxes(
    ra: np.ndarray, rs: np.ndarray, Tc: np.ndarray, ea: np.ndarray, Ta: np.ndarray,
    e_to_q: float, Ca: np.ndarray, Ci: np.ndarray,
) -> HeatFluxes:
    """Latent and sensible heat flux of a leaf (or soil surface, called
    with soil-specific ``ra``/``rs``). Direct port of
    ``SCOPEinR::get.heatfluxes``.
    """
    const = _constants()
    rhoa, cp = const["rhoa"], const["cp"]

    Tc = np.asarray(Tc, dtype=float)
    lambda_ = (2.501 - 0.002361 * Tc) * 1e6

    ei = 6.107 * 10 ** (7.5 * Tc / (237.3 + Tc))
    s = ei * 2.3026 * 7.5 * 237.3 / (237.3 + Tc) ** 2

    qi = ei * e_to_q
    qa = np.asarray(ea, dtype=float) * e_to_q

    lE = rhoa / (np.asarray(ra, dtype=float) + np.asarray(rs, dtype=float)) * lambda_ * (qi - qa)
    H = (rhoa * cp) / np.asarray(ra, dtype=float) * (Tc - np.asarray(Ta, dtype=float))
    ec = ea + (ei - ea) * ra / (ra + rs)
    Cc = Ca - (Ca - Ci) * ra / (ra + rs)

    return HeatFluxes(lambda_=lambda_, s=s, lE=lE, H=H, ec=ec, Cc=Cc)
