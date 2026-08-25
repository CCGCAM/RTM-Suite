"""FLUSPECT-B-Cx leaf model, SCOPE variant: like ``toolsrtm.fluspect_cx``,
but with a caller-configurable ``step`` (nm) for the excitation/emission
wavelength grids used to build the ``Mb``/``Mf`` fluorescence matrices
(default step=5 -> 53x71 matrices, vs. the non-SCOPE version's fixed 1 nm
-> 211x351), and a single combined ``Mb``/``Mf`` pair (using one ``phi``
spectrum) rather than separate PSI/PSII matrices.

Direct port of ``SCOPEinR/R/fluspect_Cx_forSCOPE.R`` (``getFluspect.Cx.SCOPE``),
the function SCOPE's own leaf-optics pipeline (``fluspect_mSCOPE.R``)
actually calls -- **not** the same function as ``ToolsRTM::getFluspect.Cx``
(``toolsrtm.fluspect_cx``), despite very similar code: the two diverge in
the wavelength-grid construction (``step``-parameterized here) and the SIF
scaling constant (``step`` here vs. a fixed ``int=5`` there). Confirmed via
direct inspection this SCOPE variant uses ``match()`` for ``Iwlf`` (not the
``intersect()`` bug in ``ToolsRTM::getFluspect.Cx``), so it does not
reproduce that bug -- it doesn't have it in the first place.

Reuses :func:`toolsrtm.fluspect._prospect_mesophyll` for the shared
PROSPECT-with-interfaces-removed core (identical math to the non-SCOPE
Cx path).
"""
from __future__ import annotations

import functools
from dataclasses import dataclass
from importlib import resources

import numpy as np
from toolsrtm.fluspect import _prospect_mesophyll

__all__ = ["FluspectCxScopeResult", "get_fluspect_cx_scope"]

_WLP = np.arange(400, 2401)


def _r_seq(from_: float, to: float, by: float) -> np.ndarray:
    """Mimic R's ``seq(from, to, by)``: from, from+by, ... up to the last
    value <= to (inclusive if exactly reachable)."""
    n = int(np.floor((to - from_) / by + 1e-9)) + 1
    return from_ + by * np.arange(n)


@dataclass
class _ScopeFluspectOptipar:
    wl: np.ndarray
    nr: np.ndarray
    Kab: np.ndarray
    Kca: np.ndarray
    Ks: np.ndarray
    Kw: np.ndarray
    Kdm: np.ndarray
    KcaV: np.ndarray
    KcaZ: np.ndarray
    Kant: np.ndarray
    Kp: np.ndarray
    Kcbc: np.ndarray
    phi: np.ndarray
    phiI: np.ndarray
    phiII: np.ndarray


@functools.lru_cache(maxsize=None)
def _scope_fluspect_optipar() -> _ScopeFluspectOptipar:
    with resources.files("scopeinpython.data").joinpath("optipar_fluspect_scope.csv").open("r", encoding="utf-8") as f:
        f.readline()
        d = np.loadtxt(f, delimiter=",")
    cols = ("wl", "nr", "Kab", "Kca", "Ks", "Kw", "Kdm", "KcaV", "KcaZ", "Kant", "Kp", "Kcbc", "phi", "phiI", "phiII")
    return _ScopeFluspectOptipar(**{c: d[:, i] for i, c in enumerate(cols)})


@dataclass
class FluspectCxScopeResult:
    lambda_: np.ndarray
    refl: np.ndarray
    tran: np.ndarray
    kChlrel: np.ndarray
    kCarrel: np.ndarray
    Mb: np.ndarray  # (len(wlf), len(wle)); 53x71 at the default step=5
    Mf: np.ndarray


def get_fluspect_cx_scope(
    Cab: float, Car: float, EWT: float, LMA: float, Cs: float, N: float, fqe: float, Cx: float,
    Prot: float, CBC: float, Anth: float, step: float = 5.0,
) -> FluspectCxScopeResult | None:
    """FLUSPECT-B-Cx (SCOPE variant). Direct port of
    ``SCOPEinR::getFluspect.Cx.SCOPE``. Returns ``None`` if ``fqe <= 0``
    (same reasoning as ``toolsrtm.fluspect_cx``: the R source's ``Mb``/``Mf``
    -- and everything else -- are only ever returned inside its
    ``if (fqe_ > 0) {...}`` block).
    """
    op = _scope_fluspect_optipar()
    if LMA > 0 and (Prot > 0 or CBC > 0):
        LMA = 0.0

    if Cx == -999:
        Kca = op.Kca
    else:
        Kca = (1 - Cx) * op.KcaV + Cx * op.KcaZ

    Kall = (Cab * op.Kab + Car * Kca + LMA * op.Kdm + EWT * op.Kw + Cs * op.Ks
            + Anth * op.Kant + Prot * op.Kp + CBC * op.Kcbc) / N

    j = Kall > 0
    kChlrel = np.where(j, Cab * op.Kab / (np.where(j, Kall, 1.0) * N), 0.0)
    kCarrel = np.where(j, Car * Kca / (np.where(j, Kall, 1.0) * N), 0.0)

    core = _prospect_mesophyll(Kall, N, op.nr)
    refl, tran = core["refl"], core["tran"]

    if fqe <= 0:
        return None

    fqe_vec = np.array([fqe / 5, fqe])
    talf, r21, t21 = core["talf"], core["r21"], core["t21"]
    rho, tau, k, s = core["rho"], core["tau"], core["k"], core["s"]
    kChl = kChlrel * k

    wle = _r_seq(400.0, 750.0, step) if step != 1 else _WLP[:351].astype(float)
    wlf = _r_seq(640.0, 850.0, step - 1) if step != 1 else _WLP[240:451].astype(float)
    Iwle = (wle - 400).astype(int)
    Iwlf = (wlf - 400).astype(int)

    k_e, s_e, kChl_e = k[Iwle], s[Iwle], kChl[Iwle]
    k_f, s_f = k[Iwlf], s[Iwlf]

    ndub = 15
    eps = 2.0 ** (-ndub)
    te = 1 - (k_e + s_e) * eps
    tf = 1 - (k_f + s_f) * eps
    re = s_e * eps
    rf = s_f * eps

    sigmoid = 1.0 / (1.0 + np.outer(np.exp(-wlf / 10), np.exp(wle / 10)))  # (nwlf, nwle)
    phi_plot = op.phi[Iwlf]

    Mf = Mb = step * fqe_vec[1] * np.outer(0.5 * phi_plot * eps, kChl_e) * sigmoid

    for _ in range(ndub):
        xe = te / (1 - re * re)
        ten = te * xe
        ren = re * (1 + ten)
        xf = tf / (1 - rf * rf)
        tfn = tf * xf
        rfn = rf * (1 + tfn)

        A11 = xf[:, None] + xe[None, :]
        A12 = (xf[:, None] * xe[None, :]) * (rf[:, None] + re[None, :])
        A21 = 1 + (xf[:, None] * xe[None, :]) * (1 + rf[:, None] * re[None, :])
        A22 = (xf * rf)[:, None] + (xe * re)[None, :]

        Mf, Mb = Mf * A11 + Mb * A12, Mb * A21 + Mf * A22
        te, re, tf, rf = ten, ren, tfn, rfn

    Rb = rho + tau**2 * r21 / (1 - rho * r21)
    Xe = talf[Iwle] / (1 - r21[Iwle] * Rb[Iwle])
    Xf = t21[Iwlf] / (1 - r21[Iwlf] * Rb[Iwlf])
    Ye = tau[Iwle] * r21[Iwle] / (1 - rho[Iwle] * r21[Iwle])
    Yf = tau[Iwlf] * r21[Iwlf] / (1 - rho[Iwlf] * r21[Iwlf])

    A = Xe[None, :] * (1 + Ye[None, :] * Yf[:, None]) * Xf[:, None]
    B = Xe[None, :] * (Ye[None, :] + Yf[:, None]) * Xf[:, None]

    Mb_n = A * Mb + B * Mf
    Mf_n = A * Mf + B * Mb

    return FluspectCxScopeResult(
        lambda_=op.wl.copy(), refl=refl, tran=tran, kChlrel=kChlrel, kCarrel=kCarrel,
        Mb=Mb_n, Mf=Mf_n,
    )
