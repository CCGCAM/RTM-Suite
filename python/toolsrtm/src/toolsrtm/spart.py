"""SPART: top-of-canopy (TOC) and top-of-atmosphere (TOA) reflectance/radiance.

Direct port of ``ToolsRTM::SPART`` -- BSM soil reflectance
(``getBSM.toolsRTM``/``soilwat``) + fourSAIL canopy BRDF + ``Compute_BRF``
(:func:`spart_toc`), plus the SMAC atmospheric TOA path (:func:`spart_toa`,
:mod:`toolsrtm.smac`) that turns TOC reflectance into top-of-atmosphere
reflectance/radiance for a specific sensor.

**Scope of the TOA path**: only Sentinel-2A (MSI)'s sensor coefficients are
bundled so far, out of the 9 the R package ships (see
:mod:`toolsrtm.smac`'s module docstring for why that's a data-export
exercise, not a code-scope limitation, and how to add the other 8).
"""
from __future__ import annotations

import functools
from dataclasses import dataclass
from importlib import resources
from typing import Literal

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.stats import poisson

from .canopy import foursail
from .inform import compute_brf
from .leaf import calctav
from .smac import SmacSensor, get_smac, sentinel2a_msi, spectral_convolution

__all__ = ["spart_toc", "SpartToaResult", "spart_toa"]


@dataclass
class _SpartBSMOptipar:
    wl: np.ndarray
    GSV1: np.ndarray
    GSV2: np.ndarray
    GSV3: np.ndarray
    Kw: np.ndarray
    nw: np.ndarray


@functools.lru_cache(maxsize=None)
def _spart_bsm_optipar() -> _SpartBSMOptipar:
    with resources.files("toolsrtm.data").joinpath("optipar_spart_bsm.csv").open("r", encoding="utf-8") as f:
        f.readline()
        d = np.loadtxt(f, delimiter=",")
    return _SpartBSMOptipar(wl=d[:, 0], GSV1=d[:, 1], GSV2=d[:, 2], GSV3=d[:, 3], Kw=d[:, 4], nw=d[:, 5])


def _soilwat(rdry: np.ndarray, nw: np.ndarray, kw: np.ndarray, SMp: float, SMC: float, film: float) -> np.ndarray:
    """Direct port of ``ToolsRTM::soilwat`` (same algorithm as
    ``scopeinpython.soilwat``, duplicated here since ``toolsrtm`` can't depend on
    ``scopeinpython`` -- see module docstring)."""
    rdry = np.asarray(rdry, dtype=float)
    nw = np.asarray(nw, dtype=float)
    kw = np.asarray(kw, dtype=float)

    k = np.arange(0, 7)
    nk = len(k)
    mu = (SMp - 5) / SMC
    if mu <= 0:
        return rdry.copy()

    rbac = 1 - (1 - rdry) * (rdry * calctav(90, 2.0 / nw) / calctav(90, 2.0) + 1 - rdry)
    p = 1 - calctav(90, nw) / nw**2
    Rw = 1 - calctav(40, nw)

    # mirrors R's stats::dpois(round(mu), k): x=round(mu) fixed, rate=k varied
    fmul = poisson.pmf(round(mu), k)

    tw = np.exp(-2 * np.outer(kw, film * k))
    Rwet_k = Rw[:, None] + (1 - Rw[:, None]) * (1 - p[:, None]) * tw * rbac[:, None] / (
        1 - p[:, None] * tw * rbac[:, None]
    )
    return rdry * fmul[0] + Rwet_k[:, 1:nk] @ fmul[1:nk]


def _get_bsm_toolsrtm(BSMBrightness: float, BSMlat: float, BSMlon: float, SMp: float, SMC: float, film: float) -> np.ndarray:
    """Direct port of ``ToolsRTM::getBSM.toolsRTM``, using ToolsRTM's own
    bundled ``optipar`` GSV/Kw/nw spectra (400-2400 nm, 2001 pts) --
    a separate dataset from SCOPEinR's own BSM optipar table."""
    op = _spart_bsm_optipar()
    rd = np.pi / 180.0
    f1 = BSMBrightness * np.sin(rd * BSMlat)
    f2 = BSMBrightness * np.cos(rd * BSMlat) * np.sin(rd * BSMlon)
    f3 = BSMBrightness * np.cos(rd * BSMlat) * np.cos(rd * BSMlon)
    rdry = f1 * op.GSV1 + f2 * op.GSV2 + f3 * op.GSV3
    return _soilwat(rdry, op.nw, op.Kw, SMp, SMC, film)


def spart_toc(
    inputLUT: dict,
    leaf_model: Literal["PROSPECT-D", "PROSPECT-PRO"] = "PROSPECT-PRO",
    rsoil: np.ndarray | None = None,
    BSMBrightness: float = 0.5,
    BSMlat: float = 25.0,
    BSMlon: float = 45.0,
    SMp: float = 15.0,
    SMC: float = 25.0,
    film: float = 0.015,
) -> np.ndarray:
    """Top-of-canopy BRDF reflectance, 400-2400 nm -- the TOC-only part of
    ``ToolsRTM::SPART`` (see module docstring for what's not ported).

    Equivalent to R's ``SPART(...)$rfl.toc.brdf$rfl.toc`` (pre-sensor-
    convolution values; the R function also spline-interpolates this onto a
    sensor's band centers as part of ``output$rfl.toc.BRDF`` -- that
    resampling step isn't reproduced here, only the underlying 1 nm spectrum
    both are derived from).

    Parameters
    ----------
    inputLUT : dict
        fourSAIL + leaf-model keys, as in :func:`toolsrtm.canopy.foursail`.
    leaf_model : {'PROSPECT-D', 'PROSPECT-PRO'}
    rsoil : array_like, shape (2001,), optional
        Soil reflectance, 400-2400 nm. If omitted, computed from the BSM
        soil model using ``BSMBrightness``/``BSMlat``/``BSMlon``/``SMp``/
        ``SMC``/``film`` (R's ``SPART()`` defaults: brightness=0.5, lat=25,
        lon=45, SMp=15, SMC=25, film=0.015).
    """
    if rsoil is None:
        rsoil = _get_bsm_toolsrtm(BSMBrightness, BSMlat, BSMlon, SMp, SMC, film)
    else:
        rsoil = np.asarray(rsoil, dtype=float)
        if rsoil.shape[0] != 2001:
            raise ValueError(f"rsoil must have length 2001 (400-2400nm, 1nm step), got {rsoil.shape[0]}")

    sail = foursail(inputLUT, rsoil, leaf_model=leaf_model, spectrum_all=False)
    return compute_brf(sail.rdot, sail.rsot, inputLUT["tts"], short_waves=True)


@functools.lru_cache(maxsize=None)
def _extraterrestrial_irradiance() -> np.ndarray:
    """SCOPE/ToolsRTM's bundled extraterrestrial (top-of-atmosphere)
    irradiance spectrum, 400-2400nm, W m-2 nm-1 (``ToolsRTM::Extraterrestrial_irradiance``).
    Default ``df.irradiance`` in R's ``SPART()`` when the caller omits it."""
    with resources.files("toolsrtm.data").joinpath("extraterrestrial_irradiance.csv").open("r", encoding="utf-8") as f:
        f.readline()
        d = np.loadtxt(f, delimiter=",")
    return d[:, 1]


def _interp_to_bands(wl_src: np.ndarray, y: np.ndarray, wl_bands: np.ndarray) -> np.ndarray:
    """Matches R's ``signal::interp1(wl_src, y, wl_bands, 'spline', 1E-4)``
    -- ``not-a-knot`` cubic spline, same documented approximation used
    elsewhere in this port (e.g. ``scopeinpython.rtmf``'s upsampling) vs
    R's ``fmm``-method spline."""
    cs = CubicSpline(wl_src, y, bc_type="not-a-knot", extrapolate=False)
    return cs(wl_bands)


@dataclass
class SpartToaResult:
    """Sensor-band-resolved output of :func:`spart_toa` (R's ``SPART()$output``)."""

    wl_smac: np.ndarray  # sensor band centers, (nbands,)
    rfl_toa: np.ndarray  # TOA reflectance, per band
    rad_toa: np.ndarray  # TOA radiance, per band (W m-2 sr-1 nm-1)
    rfl_toc: np.ndarray  # TOC reflectance (SMAC-combined direct+diffuse), per band
    rfl_toc_brdf: np.ndarray  # TOC BRDF reflectance (spart_toc, interpolated to bands)


def spart_toa(
    inputLUT: dict,
    sensor: SmacSensor | None = None,
    leaf_model: Literal["PROSPECT-D", "PROSPECT-PRO"] = "PROSPECT-PRO",
    rsoil: np.ndarray | None = None,
    BSMBrightness: float = 0.5,
    BSMlat: float = 25.0,
    BSMlon: float = 45.0,
    SMp: float = 15.0,
    SMC: float = 25.0,
    film: float = 0.015,
    irradiance: np.ndarray | None = None,
) -> SpartToaResult:
    """Top-of-atmosphere reflectance/radiance for a specific sensor -- the
    full ``ToolsRTM::SPART()`` pipeline (TOC BRDF + SMAC atmospheric
    correction), sensor-band-resolved. Direct port of ``ToolsRTM::SPART``
    (against the fixed R source; see ``ToolsRTM/R/spart.R``'s own comment
    on the ``optipar2021.Pro.CX``-not-in-ToolsRTM bug already fixed there).

    Parameters
    ----------
    inputLUT : dict
        Everything :func:`spart_toc` needs (fourSAIL + leaf-model keys),
        plus the atmosphere keys ``Pa`` (hPa), ``aot550``, ``uo3`` (atm-cm),
        ``uh2o`` (g/cm2).
    sensor : SmacSensor, optional
        Defaults to :func:`toolsrtm.smac.sentinel2a_msi` (the only sensor
        bundled so far -- see module docstring).
    leaf_model, rsoil, BSMBrightness, BSMlat, BSMlon, SMp, SMC, film
        As in :func:`spart_toc`.
    irradiance : array_like, shape (2001,), optional
        Extraterrestrial irradiance, 400-2400nm, W m-2 nm-1. Defaults to
        the bundled ``ToolsRTM::Extraterrestrial_irradiance``.
    """
    if sensor is None:
        sensor = sentinel2a_msi()

    if rsoil is None:
        rsoil = _get_bsm_toolsrtm(BSMBrightness, BSMlat, BSMlon, SMp, SMC, film)
    else:
        rsoil = np.asarray(rsoil, dtype=float)
        if rsoil.shape[0] != 2001:
            raise ValueError(f"rsoil must have length 2001 (400-2400nm, 1nm step), got {rsoil.shape[0]}")

    tts = inputLUT["tts"]
    sail = foursail(inputLUT, rsoil, leaf_model=leaf_model, spectrum_all=False)
    rfl_canopy_brdf = compute_brf(sail.rdot, sail.rsot, tts, short_waves=True)

    wlO = np.arange(400, 2401, 1, dtype=float)  # 400-2400nm, 1nm step -- matches rsoil/sail's grid
    wl_smac = sensor.wl_smac

    rv_so = _interp_to_bands(wlO, sail.rsot, wl_smac)
    rv_do = _interp_to_bands(wlO, sail.rdot, wl_smac)
    rv_dd = _interp_to_bands(wlO, sail.rddt, wl_smac)
    rv_sd = _interp_to_bands(wlO, sail.rsdt, wl_smac)
    rfl_toc_brdf_sensor = _interp_to_bands(wlO, rfl_canopy_brdf, wl_smac)

    atm = get_smac(
        sensor, tts=tts, tto=inputLUT["tto"], psi=inputLUT["psi"],
        Pa=inputLUT["Pa"], taup550=inputLUT["aot550"], uo3=inputLUT["uo3"], uh2o=inputLUT["uh2o"],
    )

    R_TOC = (atm.Ta_ss * rv_so + atm.Ta_sd * rv_do) / (atm.Ta_ss + atm.Ta_sd)

    rtoa1 = (atm.Ta_sd * rv_do + atm.Ta_ss * rv_sd * atm.Ra_dd * rv_do) * atm.Ta_oo / (1 - rv_dd * atm.Ra_dd)
    rtoa2 = (atm.Ta_ss * rv_sd + atm.Ta_sd * rv_dd) * atm.Ta_do / (1 - rv_dd * atm.Ra_dd)
    rtoa0 = atm.Ra_so + atm.Ta_ss * rv_so * atm.Ta_oo
    R_TOA = atm.Tg * (rtoa0 + rtoa1 + rtoa2)

    if irradiance is None:
        irradiance = _extraterrestrial_irradiance()
    Ea_bands = spectral_convolution(wlO, irradiance, sensor)
    L_TOA = Ea_bands * R_TOA

    return SpartToaResult(
        wl_smac=wl_smac, rfl_toa=R_TOA, rad_toa=L_TOA, rfl_toc=R_TOC, rfl_toc_brdf=rfl_toc_brdf_sensor,
    )
