"""RTMt (Stefan-Boltzmann variant): total thermal-infrared outgoing
radiation and net radiation per leaf/soil component, given already-solved
leaf/soil temperatures.

Direct, partial port of ``SCOPEinR::get.RTMt.sb`` (``SCOPEinR/R/RTMt.sb.R``)
-- the spectrally-integrated (Stefan-Boltzmann) thermal RTM, cheaper than
the per-wavelength ``RTMt_planck.R`` (not ported). Leaf/soil temperatures
(``Tcu``/``Tch``/``Tsu``/``Tsh``) are inputs here, not solved for -- that
is :mod:`scopeinpython.ebal`'s job (not yet ported).

**Only the "SCOPE-lite" scalar-per-layer branch is ported** (R's
``is.matrix(Hcsu3) == FALSE`` path, i.e. ``Tcu``/``Tch`` given as plain
length-``nl`` vectors, one temperature per canopy layer -- not the full
``(13, 36, nl)`` per-leaf-angle-class array). This matches every
reference case built during this whole port (``data.opts`` always sets
``lite = 1``) and the full-array branch has its own unresolved question
(see below), so is intentionally left for later.

**Not ported: the ``obsdir`` (observation-direction brightness
temperature) branch.** Not needed for :mod:`scopeinpython.ebal`'s
convergence loop, which only depends on the core net-radiation output
below.

**A structural note, not a bug**: in this formulation, the "direct thermal
emission" term ``Es_`` is initialized to zero and never given a nonzero
seed anywhere in the recursion (unlike the analogous ``Es_``/``Esun_``
in RTMo/RTMf, which are seeded from real solar irradiance) -- so every
``Xsd``/``R_sd``-weighted term in the layer recursion below is always
exactly zero. Kept in the port for structural fidelity with R (and in
case a future ``obsdir``/full port needs it), not simplified away.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .rtmo import RTMoResult
from .thermal import stefan_boltzmann

__all__ = ["RTMtSbResult", "rtmt_sb"]


@dataclass
class RTMtSbResult:
    Emint: np.ndarray  # downward diffuse thermal flux per layer boundary, (nl+1,)
    Eplut: np.ndarray  # upward diffuse thermal flux per layer boundary, (nl+1,)
    Eoutte: float  # TOC outgoing thermal flux, W/m2
    Rnuct: np.ndarray  # net thermal radiation, sunlit leaves, per layer, (nl,)
    Rnhct: np.ndarray  # net thermal radiation, shaded leaves, per layer, (nl,)
    Rnust: float  # net thermal radiation, sunlit soil
    Rnhst: float  # net thermal radiation, shaded soil


def rtmt_sb(
    rtmo: RTMoResult,
    nl: int,
    LAI: float,
    rho_thermal: float,
    tau_thermal: float,
    rs_thermal: float,
    Tcu: np.ndarray,
    Tch: np.ndarray,
    Tsu: float,
    Tsh: float,
) -> RTMtSbResult:
    """Direct port of the scalar/"lite" branch of ``SCOPEinR::get.RTMt.sb``
    (see module docstring for what's not ported).

    Parameters
    ----------
    rtmo : RTMoResult
        From :func:`scopeinpython.rtmo.run_rtmo`, called with the same
        canopy as here. Only the thermal-region (last-column) values of
        ``Xdd``/``rho_dd``/``tau_dd``/``R_dd`` and the (wavelength-
        independent) ``Xss`` are used.
    Tcu, Tch : array_like, shape (nl,)
        Sunlit / shaded leaf temperature per canopy layer, deg C.
    Tsu, Tsh : float
        Sunlit / shaded soil temperature, deg C.
    """
    Ps = rtmo.Ps
    iLAI = LAI / nl

    Xdd = rtmo.Xdd[:, -1]
    Xsd = rtmo.Xsd[:, -1]
    Xss = rtmo.Xss
    R_dd = rtmo.R_dd[:, -1]
    R_sd = rtmo.R_sd[:, -1]
    rho_dd = rtmo.rho_dd[:, -1]
    tau_dd = rtmo.tau_dd[:, -1]

    epsc = 1 - rho_thermal - tau_thermal
    epss = 1 - rs_thermal

    Hcsu = epsc * stefan_boltzmann(np.asarray(Tcu, dtype=float))  # (nl,)
    Hcsh = epsc * stefan_boltzmann(np.asarray(Tch, dtype=float))
    Hssu = epss * stefan_boltzmann(Tsu)
    Hssh = epss * stefan_boltzmann(Tsh)

    Hcsu_ = Hcsu * Ps[:nl]
    Hcsh_ = Hcsh * (1 - Ps[:nl])
    Hc = Hcsu_ + Hcsh_  # hemispherical emittance by leaf layer, (nl,)
    Hs = Hssu * Ps[nl] + Hssh * (1 - Ps[nl])  # hemispherical emittance by soil

    U = np.zeros(nl + 1)
    U[nl] = Hs
    Es_ = np.zeros(nl + 1)  # always stays zero, see module docstring
    Emin = np.zeros(nl + 1)
    Eplu = np.zeros(nl + 1)
    Y = np.zeros(nl)

    for j in range(nl - 1, -1, -1):
        Y[j] = (rho_dd[j] * U[j + 1] + Hc[j] * iLAI) / (1 - rho_dd[j] * R_dd[j + 1])
        U[j] = tau_dd[j] * (R_dd[j + 1] * Y[j] + U[j + 1]) + Hc[j] * iLAI

    for j in range(nl):
        Es_[j + 1] = Xss[j] * Es_[j]
        Emin[j + 1] = Xsd[j] * Es_[j] + Xdd[j] * Emin[j] + Y[j]
        Eplu[j] = R_sd[j] * Es_[j] + R_dd[j] * Emin[j] + U[j]

    Eplu[nl] = R_sd[nl - 1] * Es_[nl - 1] + R_dd[nl - 1] * Emin[nl - 1] + Hs
    Eoutte = float(Eplu[0])

    Rnuc = Emin[:nl] + Eplu[1:nl + 1] - 2 * Hcsu
    Rnhc = Emin[:nl] + Eplu[1:nl + 1] - 2 * Hcsh
    Rnus = float(Emin[nl] - Hssu)
    Rnhs = float(Emin[nl] - Hssh)

    return RTMtSbResult(
        Emint=Emin, Eplut=Eplu, Eoutte=Eoutte,
        Rnuct=Rnuc, Rnhct=Rnhc, Rnust=Rnus, Rnhst=Rnhs,
    )
