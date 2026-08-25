"""SCOPE spectral region definitions.

Direct port of ``SCOPEinR::get.spectra.SCOPE`` (SCOPEinR/R/define_bands.R).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = ["SpectralConfig", "get_spectra_scope"]


@dataclass
class SpectralConfig:
    reg1: np.ndarray  # 400:2400 step 1 (optical)
    reg2: np.ndarray  # 2500:15000 step 100
    reg3: np.ndarray  # 16000:50000 step 1000
    wlS: np.ndarray  # concatenation of reg1, reg2, reg3 -- all SCOPE wavelengths
    wlP: np.ndarray  # PROSPECT data range == reg1
    wlO: np.ndarray  # optical part == reg1
    wlT: np.ndarray  # thermal part == reg2 + reg3
    wlPAR: np.ndarray  # PAR range, 400-700 nm subset of wlS
    IwlP: np.ndarray  # 0-based indices of reg1 within wlS
    IwlT: np.ndarray  # 0-based indices of the thermal part within wlS
    wlIrrad: np.ndarray  # same grid as wlS (irradiance file wavelengths)
    wlE: np.ndarray = field(default_factory=lambda: np.arange(400, 751, 1))
    wlF: np.ndarray = field(default_factory=lambda: np.arange(640, 851, 1))


def get_spectra_scope() -> SpectralConfig:
    """Return the SCOPE spectral region definitions (wavelength grids and
    band indices) used throughout the RTMo pipeline.

    Direct port of ``SCOPEinR::get.spectra.SCOPE(getSpectral = TRUE)``.
    """
    reg1 = np.arange(400, 2401, 1, dtype=float)
    reg2 = np.arange(2500, 15001, 100, dtype=float)
    reg3 = np.arange(16000, 50001, 1000, dtype=float)
    wlS = np.concatenate([reg1, reg2, reg3])

    wlP = reg1
    wlO = reg1
    wlT = np.concatenate([reg2, reg3])
    wlPAR = wlS[(wlS >= 400) & (wlS <= 700)]

    IwlP = np.arange(len(reg1))  # 0-based, matches R's 1:length(reg1) minus 1
    IwlT = np.arange(len(reg1), len(reg1) + len(reg2) + len(reg3))

    wlIrrad = wlS.copy()

    return SpectralConfig(
        reg1=reg1, reg2=reg2, reg3=reg3, wlS=wlS, wlP=wlP, wlO=wlO, wlT=wlT,
        wlPAR=wlPAR, IwlP=IwlP, IwlT=IwlT, wlIrrad=wlIrrad,
    )
