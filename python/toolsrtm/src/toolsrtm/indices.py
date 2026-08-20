"""Spectral vegetation indices, direct port of ``ToolsRTM::getIndices``.

The R function takes a dataframe with named reflectance columns
(``R.400``, ``R.670``, ...) and a regex pattern to find them; the Python
port instead takes ``wavelengths``/``reflectance`` arrays directly (the
numerically load-bearing part -- linear interpolation onto an integer nm
grid, then per-index algebra on named wavelengths -- is unchanged, only the
dataframe-column-parsing convenience layer is dropped, since it doesn't
change any numbers).

A handful of formulas below are reproduced **as literally written in the R
source**, bugs included, rather than "fixed" -- see the inline notes at
``CIgreen`` and ``TCARI/OSAVI.1510``, and the ``GnyLi``/``CI2`` overwrite
sequences -- so results match R exactly for anyone cross-checking against
the original package.
"""
from __future__ import annotations

from typing import Literal

import numpy as np
from scipy.interpolate import interp1d

__all__ = ["get_indices"]


def _lookup(r_grid: np.ndarray, grid_index: dict, wl) -> np.ndarray:
    """r[wl] lookup against the interpolated grid; NaN if wl isn't an
    integer key on the grid (mirrors R silently returning NA for a
    nonexistent named-vector element, e.g. ``r['B3']``)."""
    j = grid_index.get(wl)
    if j is None:
        return np.full(r_grid.shape[0], np.nan)
    return r_grid[:, j]


def get_indices(
    wavelengths: np.ndarray,
    reflectance: np.ndarray,
    spectral_domain: Literal["VNIR", "SWIR", "VNIR-SWIR"] = "VNIR",
    factor: float | None = None,
) -> dict:
    """Common spectral vegetation indices from reflectance spectra.

    Direct port of ``ToolsRTM::getIndices`` (the reflectance-column
    parsing via ``pattern.rfl`` is replaced by passing ``wavelengths``/
    ``reflectance`` directly -- see module docstring).

    Parameters
    ----------
    wavelengths : array_like, shape (nwl,)
        Wavelengths (nm) of the ``reflectance`` columns.
    reflectance : array_like, shape (nwl,) or (nsamples, nwl)
        Reflectance spectra (0-1), one or more rows.
    spectral_domain : {'VNIR', 'SWIR', 'VNIR-SWIR'}
        'VNIR': indices using 400-850 nm bands. 'SWIR': indices using
        800-2550 nm bands. 'VNIR-SWIR': both.
    factor : float, optional
        Scaling factor applied to reflectance (e.g. 1/10000 for digital
        numbers). Default 1 (no scaling).

    Returns
    -------
    dict
        ``{index_name: numpy.ndarray of shape (nsamples,)}``.
    """
    wavelengths = np.asarray(wavelengths, dtype=float)
    reflectance = np.atleast_2d(np.asarray(reflectance, dtype=float))
    if factor is None:
        factor = 1.0
    values = reflectance * factor

    if spectral_domain not in ("VNIR", "SWIR", "VNIR-SWIR"):
        raise ValueError("spectral_domain must be 'VNIR', 'SWIR', or 'VNIR-SWIR'")

    # R: min_ <- round(min(wavelengths)) + 10; max_ <- round(max(wavelengths)) + 10
    min_ = round(float(wavelengths.min())) + 10
    max_ = round(float(wavelengths.max())) + 10
    grid = np.arange(min_, max_ + 1)  # R's min_:max_ is inclusive of max_
    grid_index = {int(w): j for j, w in enumerate(grid)}

    n = values.shape[0]
    r_grid = np.empty((n, len(grid)))
    for i in range(n):
        f = interp1d(wavelengths, values[i], kind="linear", fill_value="extrapolate")
        r_grid[i] = f(grid)

    def r(wl):
        return _lookup(r_grid, grid_index, wl)

    idx: dict[str, np.ndarray] = {}

    if spectral_domain == "VNIR":
        # ---- Structural indices ----
        idx["NDVI"] = (r(800) - r(670)) / (r(800) + r(670))
        idx["RDVI"] = (r(800) - r(670)) / (r(800) + r(670)) ** 0.5
        idx["SR"] = r(800) / r(670)
        idx["MSR"] = (r(800) / r(670) - 1) / ((r(800) / r(670)) ** 0.5 + 1)
        idx["OSAVI"] = (1 + 0.16) * (r(800) - r(670)) / (r(800) + r(670) + 0.16)
        idx["MSAVI"] = 0.5 * (2 * r(800) + 1 - np.sqrt((2 * r(800) + 1) ** 2 - 8 * (r(800) - r(670))))
        idx["MTVI1"] = 1.2 * (1.2 * (r(800) - r(550)) - 2.5 * (r(670) - r(550)))
        idx["MTVI2"] = (1.5 * (1.2 * (r(800) - r(550)) - 2.5 * (r(670) - r(550)))) / np.sqrt(
            (2 * r(800) + 1) ** 2 - (6 * r(800) - 5 * np.sqrt(r(670))) - 0.5
        )
        idx["MCARI"] = ((r(700) - r(670)) - 0.2 * (r(700) - r(550))) * (r(700) / r(670))
        idx["MCARI1"] = 1.5 * (2.5 * (r(800) - r(670)) - 1.3 * (r(800) - r(550)))
        idx["MCARI2"] = (1.5 * (2.5 * (r(800) - r(670)) - 1.3 * (r(800) - r(550)))) / np.sqrt(
            (2 * r(800) + 1) ** 2 - (6 * r(800) - 5 * np.sqrt(r(670))) - 0.5
        )
        idx["EVI"] = 2.5 * (r(800) - r(670)) / (r(800) + 6 * r(670) - 7.5 * r(400) + 1)
        idx["LIC1"] = (r(800) - r(680)) / (r(800) + r(680))

        # ---- Pigment indices ----
        idx["VOG"] = r(740) / r(720)
        idx["VOG2"] = (r(734) - r(747)) / (r(715) + r(726))
        idx["VOG3"] = (r(734) - r(747)) / (r(715) + r(720))
        idx["GM1"] = r(750) / r(550)
        idx["GM2"] = r(750) / r(700)
        idx["TCARI"] = 3 * ((r(700) - r(670)) - 0.2 * (r(700) - r(550)) * (r(700) / r(670)))
        idx["T.O"] = idx["TCARI"] / idx["OSAVI"]
        idx["CI"] = r(750) / r(710)
        idx["TVI"] = 0.5 * (120 * (r(750) - r(550)) - 200 * (r(670) - r(550)))
        idx["SRPI"] = r(430) / r(680)
        idx["NPQI"] = (r(415) - r(435)) / (r(415) + r(435))
        idx["NPCI"] = (r(680) - r(430)) / (r(680) + r(430))
        idx["CTR1"] = r(695) / r(420)
        idx["CAR"] = r(515) / r(570)
        idx["DCabxc"] = r(672) / (r(550) * (3 * r(708)))
        idx["DNCabxc"] = r(860) / (r(550) * r(708))
        idx["SIPI"] = (r(800) - r(445)) / (r(800) + r(680))
        idx["CRI550"] = (1 / r(510)) - (1 / r(550))
        idx["CRI700"] = (1 / r(510)) - (1 / r(700))
        idx["CRI550m"] = (1 / r(515)) - (1 / r(550))
        idx["CRI700m"] = (1 / r(515)) - (1 / r(700))
        idx["RCRI550"] = (1 / r(510)) - (1 / r(550)) * r(770)
        idx["RCRI700"] = (1 / r(510)) - (1 / r(700)) * r(770)
        idx["PSRI"] = (r(680) - r(500)) / r(750)
        idx["LIC3"] = r(440) / r(740)

        # ---- Red-edge indices ----
        idx["CIre"] = (r(782) / r(705)) - 1
        idx["CIrededge"] = (r(800) / r(705)) - 1
        # R source: r['B3'] -- 'B3' isn't a wavelength key (this is the R
        # source's own bug: NA in R, NaN here), reproduced as-is.
        idx["CIgreen"] = (r(800) / r("B3")) - 1
        idx["Chlred.edge"] = (r(780) / r(705)) ** (-1)
        idx["CVI"] = (r(800) * r(665)) / r(665) ** 2
        idx["IRECI"] = (r(780) - r(665)) / (r(705) / r(740))
        idx["REP"] = 700 + 40 * (((r(665) + r(780)) / 2) - r(705)) / (r(740) - r(705))
        idx["RVI"] = r(800) / r(665)
        idx["RedEg1"] = r(705) / r(665)
        idx["RedEg2"] = (r(705) - r(665)) / (r(705) + r(665))

        # ---- PRI indices ----
        idx["PRI"] = (r(570) - r(530)) / (r(570) + r(530))
        idx["PRI515"] = (r(515) - r(530)) / (r(515) + r(530))
        idx["PRIM1"] = (r(512) - r(531)) / (r(512) + r(531))
        idx["PRIM2"] = (r(600) - r(531)) / (r(600) + r(531))
        idx["PRIM3"] = (r(670) - r(531)) / (r(670) + r(531))
        idx["PRIM4"] = (r(570) - r(531) - r(670)) / (r(570) + r(531) + r(670))
        idx["PRIn"] = idx["PRI"] / (idx["RDVI"] * r(700) / r(670))
        idx["PRI_CI"] = idx["PRI"] * ((r(760) / r(700)) - 1)

        # ---- BGR indices ----
        idx["B"] = r(450) / r(490)
        idx["G"] = r(550) / r(670)
        idx["R"] = r(700) / r(670)
        idx["BGI1"] = r(400) / r(550)
        idx["BGI2"] = r(450) / r(550)
        idx["BF1"] = r(400) / r(410)
        idx["BF2"] = r(400) / r(420)
        idx["BF3"] = r(400) / r(430)
        idx["BF4"] = r(400) / r(440)
        idx["BF5"] = r(400) / r(450)
        idx["BRI1"] = r(400) / r(690)
        idx["BRI2"] = r(450) / r(690)
        idx["RGI"] = r(690) / r(550)
        idx["RARS"] = r(746) / r(513)
        idx["LIC2"] = r(440) / r(690)
        idx["HI"] = (r(534) - r(698)) / (r(534) + r(698)) - 0.5 * r(704)
        idx["CUR"] = (r(675) * r(690)) / r(683) ** 2

        # ---- NIR-VIS indices ----
        idx["PSSRa"] = r(800) / r(680)
        idx["PSSRb"] = r(800) / r(635)
        idx["PSSRc"] = r(800) / r(470)
        idx["PSNDc"] = (r(800) - r(470)) / (r(800) + r(470))

        # ---- Continuum-removal red-edge indices (J.B. Feret CR_SWIR config) ----
        w06, w_800, w05, w07, w_762, w04 = 740, 800, 704, 782, 762, 665
        idx["CR.red.nir.1"] = r(740) / (r(800) + (w06 - w_800) * (r(782) - r(800)) / (w_762 - w_800))
        idx["CR.red.nir.2"] = r(704) / (r(782) + (w05 - w07) * (r(740) - r(782)) / (w06 - w_762))
        idx["CR.red.nir.3"] = r(704) / (r(800) + (w05 - w_800) * (r(740) - r(800)) / (w06 - w_800))
        idx["CR.red.nir.4"] = r(665) / (r(800) + (w04 - w_800) * (r(704) - r(800)) / (w06 - w_800))
        idx["CR.red.nir.5"] = r(740) / (r(800) + (w06 - w_800) * (r(704) - r(800)) / (w05 - w_800))
        idx["CR.red.nir"] = r(782) / (r(800) + (w07 - w_800) * (r(704) - r(800)) / (w05 - w_800))
        idx["CR.red.nir.7"] = r(762) / (r(800) + (w_762 - w_800) * (r(704) - r(800)) / (w05 - w_800))

    if spectral_domain in ("SWIR", "VNIR-SWIR"):
        # ---- SWIR indices ----
        # GnyLi is computed three times at 900/850/950 nm in the R source,
        # each reassignment overwriting 'GnyLi' -- only 'GnyLi.w850' (850 nm)
        # survives under its own name; final 'GnyLi' == 'GnyLi.w950' (950 nm).
        # The 900 nm computation is dead code, reproduced faithfully (dropped).
        idx["GnyLi.w850"] = ((r(850) * r(1050)) - (r(955) * r(1220))) / ((r(850) * r(1050)) + (r(955) * r(1220)))
        idx["GnyLi.w950"] = ((r(950) * r(1050)) - (r(955) * r(1220))) / ((r(950) * r(1050)) + (r(955) * r(1220)))
        idx["GnyLi"] = idx["GnyLi.w950"]

        idx["CI1"] = ((r(736) - r(735)) / 1) * (r(990) / r(720))
        # CI2 is computed at 900 then 950 nm in the R source (900 nm dead code
        # once overwritten); final CI2 == the 950 nm version.
        idx["CI2"] = ((r(736) - r(735)) / 1) * (r(950) / r(720))

        idx["MCARI.1510"] = ((r(700) - r(1510)) - 0.2 * (r(700) - r(550))) * (r(700) / r(1510))
        idx["TCARI.1510"] = 3 * ((r(700) - r(1510)) - 0.2 * (r(700) - r(550)) * (r(700) / r(1510)))
        idx["OSAVI.1510"] = (1 + 0.16) * (r(800) - r(1510)) / (r(800) + r(1510) + 0.16)
        # R source references indices['TCARI 1510']/indices['OSAVI 1510']
        # (space, not the actual '.1510'-suffixed keys defined two lines
        # above) -- neither key exists, so this is NA in R; reproduced as NaN.
        idx["TCARI/OSAVI.1510"] = np.full(n, np.nan)
        idx["NRI.1510"] = (r(1510) - r(660)) / (r(1510) + r(660))
        idx["RSI.990.720"] = r(990) / r(720)
        idx["NDNI"] = (np.log10(1 / r(1510)) - np.log10(1 / r(1680))) / (np.log10(1 / r(1510)) + np.log10(1 / r(1680)))
        idx["S1080"] = (r(1080) - r(660)) / (r(1080) + r(660))
        idx["S1260"] = (r(1260) - r(660)) / (r(1260) + r(660))
        idx["N1645"] = (r(1645) - r(1715)) / (r(1645) + r(1715))
        idx["N870"] = (r(870) - r(1450)) / (r(870) + r(1450))
        idx["N850.1510"] = (r(850) - r(1510)) / (r(850) + r(1510))
        idx["NN1510"] = r(1510) / r(850)

        w11, w8A, w12 = 1614, 865, 2202
        idx["CR.SWIR"] = r(1614) / (r(865) + (w11 - w8A) * (r(2202) - r(865)) / (w12 - w8A))

    # R: drop columns that are all-NA
    idx = {k: v for k, v in idx.items() if not np.all(np.isnan(v))}
    return idx
