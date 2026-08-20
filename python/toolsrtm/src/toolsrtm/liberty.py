"""LIBERTY leaf model (Dawson, Curran & Plummer 1998) -- a Kubelka-Munk-style
model designed for conifer needles, using cell-diameter/intercellular-air-
space/baseline-absorption/leaf-thickness parameters instead of PROSPECT's
N-layer plate model.

Direct port of ``ToolsRTM/R/liberty.R``.
"""
from __future__ import annotations

import functools
from dataclasses import dataclass
from importlib import resources

import numpy as np
from scipy.interpolate import interp1d

__all__ = ["LibertyResult", "liberty"]

_NWL = 420  # ToolsRTM::dataspec.liberty has 421 rows; R's `for (i in seq(1,420))` only uses the first 420


@functools.lru_cache(maxsize=None)
def _dataspec_liberty() -> np.ndarray:
    with resources.files("toolsrtm.data").joinpath("dataspec_liberty.csv").open("r", encoding="utf-8") as f:
        f.readline()
        d = np.loadtxt(f, delimiter=",")
    return d[:_NWL, :]  # (420, 5): k_chl, k_water, ke, k_ligcell, k_proteins


@dataclass
class LibertyResult:
    lambda_: np.ndarray  # 400-2500 nm (1 nm step)
    refl: np.ndarray
    tran: np.ndarray
    RR: np.ndarray  # infinite-medium reflectance (R in the R source), before finite-thickness stacking


def liberty(
    cell_d: float, inter_c: float, baseline_abs: float, leaf_thick: float,
    albino_abs: float, Cab: float, EWT: float, lign_cell: float, Nitrogen: float,
) -> LibertyResult:
    """LIBERTY leaf reflectance/transmittance. Direct port of ``ToolsRTM::liberty``.

    Parameters
    ----------
    cell_d : float
        Average leaf cell diameter (m-6). Typical 40 (20-100).
    inter_c : float
        Intercellular air space fraction. Typical 0.045 (0.01-0.1).
    baseline_abs : float
        Wavelength-independent baseline absorption. Typical 0.0004-0.0006.
    leaf_thick : float
        Leaf thickness (arbitrary units, via Benford layer stacking). Typical 1.6 (1-10).
    albino_abs : float
        Visible-region absorption due to lignin (albino leaf proxy). Typical 2 (0-4).
    Cab : float
        Chlorophyll content, ug/cm2 (converted internally to mg/m2 as R does: ``Cab*10``).
    EWT : float
        Equivalent water thickness, g/cm2 (converted internally to g/m2: ``EWT*10000``).
    lign_cell : float
        Combined lignin+cellulose content, ug/cm2 (converted internally: ``lign_cell*10``).
    Nitrogen : float
        Nitrogen content, g/m2. Typical 1 (0.3-2.0).

    Returns
    -------
    LibertyResult
        ``lambda_``, refl, tran, RR -- all 400-2500 nm, 1 nm step (linearly
        interpolated/extrapolated from the model's native 5 nm grid).
    """
    Cab_ = Cab * 10.0
    EWT_ = EWT * 10000.0
    lign_cell_ = lign_cell * 10.0

    d = _dataspec_liberty()
    k_chl, k_water, ke, k_ligcell, k_proteins = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4]

    transR = np.zeros(_NWL)
    reflR = np.zeros(_NWL)
    wavelength = np.zeros(_NWL)
    RR = np.zeros(_NWL)

    N0 = 1.0
    whole = np.trunc(leaf_thick)
    fraction = leaf_thick - whole

    for i in range(_NWL):
        coeff = cell_d * (
            baseline_abs + k_chl[i] * Cab_ + k_water[i] * EWT_ + ke[i] * albino_abs
            + k_ligcell[i] * lign_cell_ + k_proteins[i] * Nitrogen
        )

        N1 = 1.4891 - 0.0005 * i  # i here is R's (i-1), i.e. 0-indexed

        # me: average Fresnel reflectivity integrated over the full incident
        # hemisphere (0-90 deg), 1-degree steps
        width = np.pi / 180.0
        pp = np.arange(1, 91)
        alpha = pp * width
        beta = np.arcsin(N0 / N1 * np.sin(alpha))
        plus, dif = alpha + beta, alpha - beta
        refl_ang = 0.5 * ((np.sin(dif)) ** 2 / (np.sin(plus)) ** 2 + (np.tan(dif)) ** 2 / (np.tan(plus)) ** 2)
        me = 2 * np.sum(refl_ang * np.sin(alpha) * np.cos(alpha) * width)

        # mi: same integral but only over the total-internal-reflection cone
        critical_deg = np.arcsin(N0 / N1) * 180.0 / np.pi
        jj = np.arange(1, int(np.floor(critical_deg)) + 1)
        alpha_i = jj * width
        beta_i = np.arcsin(N0 / N1 * np.sin(alpha_i))
        plus_i, dif_i = alpha_i + beta_i, alpha_i - beta_i
        refl_i = 0.5 * ((np.sin(dif_i)) ** 2 / (np.sin(plus_i)) ** 2 + (np.tan(dif_i)) ** 2 / (np.tan(plus_i)) ** 2)
        mint = np.sum(refl_i * np.sin(alpha_i) * np.cos(alpha_i) * width)
        critical_rad = critical_deg * np.pi / 180.0
        mi = (1 - np.sin(critical_rad) ** 2) + 2 * mint

        M = 2 * (1 - (coeff + 1) * np.exp(-coeff)) / coeff**2
        Tt = ((1 - mi) * M) / (1 - (mi * M))
        x = inter_c / (1 - (1 - 2 * inter_c) * Tt)
        a = me * Tt + x * Tt - me - Tt - x * me * Tt
        b = 1 + x * me * Tt - 2 * x**2 * me**2 * Tt
        c = 2 * me * x**2 * Tt - x * Tt - 2 * x * me

        R = 0.5
        for _ in range(50):
            R = -(a * R**2 + c) / b

        rb = (2 * x * me) + (x * Tt) - (x * Tt * 2 * x * me)
        tb = np.sqrt(((R - rb) * (1 - (R * rb))) / R)

        top = tb ** (1 + fraction) * ((((1 + tb) ** 2) - rb**2) ** (1 - fraction))
        bot1 = (1 + tb) ** (2 * (1 - fraction)) - rb**2
        bot2 = 1 + (64.0 / 3.0) * fraction * (fraction - 0.5) * (fraction - 1) * 0.001
        tif = top / (bot1 * bot2)
        rif = (1 + rb**2 - tb**2 - np.sqrt((1 + rb**2 - tb**2) ** 2 - 4 * rb**2 * (1 - tif**2))) / (2 * rb)

        if whole >= 2:
            prev_t, prev_r = 1.0, 0.0
            for _ in range(int(whole) - 1):
                cur_t = (prev_t * tb) / (1 - (prev_r * rb))
                cur_r = prev_r + ((prev_t * prev_t) * rb) / (1 - (prev_r * rb))
                prev_t, prev_r = cur_t, cur_r
        else:
            cur_t, cur_r = 1.0, 0.0

        trans = (cur_t * tif) / (1 - (rif * cur_r))
        refl = cur_r + ((cur_t**2 * rif) / (1 - (rif * cur_r)))

        transR[i] = trans
        reflR[i] = refl
        wavelength[i] = 400 + i * 5
        RR[i] = R

    wave_1nm = np.arange(400, 2501)
    refl_f = interp1d(wavelength, reflR, kind="linear", fill_value="extrapolate")
    tran_f = interp1d(wavelength, transR, kind="linear", fill_value="extrapolate")
    RR_f = interp1d(wavelength, RR, kind="linear", fill_value="extrapolate")

    return LibertyResult(
        lambda_=wave_1nm, refl=refl_f(wave_1nm), tran=tran_f(wave_1nm), RR=RR_f(wave_1nm),
    )
