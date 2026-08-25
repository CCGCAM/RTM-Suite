"""Loaders for the bundled physical-constants table and BSM optical
parameter table (ported from SCOPEinR/data/*.rda: ``constants``,
``optipar2017.ProspectD`` (GSV/Kw/nw columns only, used by BSM)).

Exported once from the original R package to CSV (see
python/scratch/scratch_export.R in the repo root); bundled as package data so this
port has no runtime dependency on R.
"""
from __future__ import annotations

import functools
from dataclasses import dataclass
from importlib import resources

import numpy as np


@functools.lru_cache(maxsize=None)
def constants() -> dict:
    """Physical constants table (SCOPEinR::constants), as a dict of scalars
    keyed by the R ``constant`` column (e.g. ``constants()['sigmaSB']``)."""
    out = {}
    with resources.files("scopeinpython.data").joinpath("constants.csv").open("r", encoding="utf-8") as f:
        f.readline()
        for line in f:
            parts = [p.strip().strip('"') for p in line.strip().split(",")]
            if len(parts) < 2:
                continue
            out[parts[0]] = float(parts[1])
    return out


@dataclass
class BSMOptipar:
    """GSV / water-absorption / water-refraction-index spectra used by the
    BSM soil model (subset of SCOPEinR::optipar2017.ProspectD), 400-2400 nm."""

    wl: np.ndarray
    GSV: np.ndarray  # (nwl, 3)
    Kw: np.ndarray
    nw: np.ndarray


@functools.lru_cache(maxsize=None)
def bsm_optipar() -> BSMOptipar:
    with resources.files("scopeinpython.data").joinpath("optipar_bsm.csv").open("r", encoding="utf-8") as f:
        f.readline()
        d = np.loadtxt(f, delimiter=",")
    return BSMOptipar(wl=d[:, 0], GSV=d[:, 1:4], Kw=d[:, 4], nw=d[:, 5])


@functools.lru_cache(maxsize=None)
def soil_scope_spectra() -> np.ndarray:
    """The 3 reference dry-soil reflectance spectra bundled with SCOPE
    itself (``SCOPEinR/inst/input/soil_spectra/soil_scope.txt``), used when
    ``options.soilspectrum == 0`` (soil-from-file, not BSM). Shape
    ``(2001, 3)``, 400-2400 nm, columns selected by the LUT's ``spectrum``
    field (1-indexed, matching R)."""
    with resources.files("scopeinpython.data").joinpath("soil_scope.csv").open("r", encoding="utf-8") as f:
        f.readline()
        d = np.loadtxt(f, delimiter=",")
    return d[:, 1:4]


@dataclass
class DefaultIrradiance:
    wave: np.ndarray
    Esun_: np.ndarray
    Esky_: np.ndarray


@functools.lru_cache(maxsize=None)
def default_irradiance() -> DefaultIrradiance:
    """SCOPE's bundled default top-of-atmosphere irradiance
    (``SCOPEinR::Esun_[[1]]``/``Esky_[[1]]``), used by ``get_scope`` when
    the caller doesn't supply their own ``Esun_``/``Esky_`` (matching R's
    ``options.irradiance == 0`` path, i.e. everything except the
    measurements-file and MODTRAN-atmosphere-file branches -- see
    :mod:`scopeinpython.scope`). Same wavelength grid as
    ``SpectralConfig.wlS`` (2162 points)."""
    with resources.files("scopeinpython.data").joinpath("default_irradiance.csv").open("r", encoding="utf-8") as f:
        f.readline()
        d = np.loadtxt(f, delimiter=",")
    return DefaultIrradiance(wave=d[:, 0], Esun_=d[:, 1], Esky_=d[:, 2])
