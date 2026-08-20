"""FLUSPECT-B / FLUSPECT-B-Cx leaf model: PROSPECT reflectance/transmittance
plus chlorophyll-fluorescence excitation-emission matrices, via Verhoef's
doubling method applied to the leaf mesophyll layer alone (leaf-air
interfaces removed before doubling, re-added after).

Direct port of ``ToolsRTM/R/fluspect_B.R`` (``getFluspect.B``) and
``ToolsRTM/R/fluspect_Cx.R`` (``getFluspect.Cx``). These are two separate,
duplicated R functions (not a shared core with a mode flag), so they're
ported as two separate Python functions here too, since their
fluorescence-matrix output shapes genuinely differ (``MbI``/``MbII``/
``MfI``/``MfII`` split by photosystem for the B version, vs a single
``Mb``/``Mf`` pair for the Cx version).

Reuses :func:`toolsrtm.leaf._one_layer` / :func:`toolsrtm.leaf._stokes_n_layers`
for the PROSPECT part (identical math, ``alpha`` fixed at 59 degrees here
instead of a caller-supplied angle, matching ``ToolsRTM::calctav(59, nr)``
in both R functions).
"""
from __future__ import annotations

import functools
from dataclasses import dataclass
from importlib import resources
from typing import Literal

import numpy as np

from .leaf import _one_layer, _stokes_n_layers, calctav

__all__ = ["FluspectBResult", "FluspectCxResult", "fluspect_b", "fluspect_cx"]

# SCOPE's define_bands(): wlP = 400:2400 (PROSPECT range), wlE = 400:750
# (excitation), wlF = 640:850 (fluorescence emission). All 1 nm step, all
# subsets of wlP starting at the same 400 nm origin, so "index of wlE/wlF
# within wlP" is a fixed static slice, computed once here.
_WLP = np.arange(400, 2401)
_WLE = np.arange(400, 751)
_WLF = np.arange(640, 851)
_IWLE = slice(0, len(_WLE))  # position of wlE's start/end within wlP (both start at 400)
_IWLF_CORRECT = slice(_WLF[0] - _WLP[0], _WLF[-1] - _WLP[0] + 1)  # R's `which(wlp>=640 & wlp<=850)`


@dataclass
class _FluspectOptipar:
    wl: np.ndarray
    nr: np.ndarray
    Kab: np.ndarray
    Kca: np.ndarray
    Ks: np.ndarray
    Kw: np.ndarray
    Kdm: np.ndarray
    phiI: np.ndarray
    phiII: np.ndarray
    KcaV: np.ndarray
    KcaZ: np.ndarray
    Kant: np.ndarray
    Kp: np.ndarray
    Kcbc: np.ndarray


@functools.lru_cache(maxsize=None)
def _fluspect_optipar() -> _FluspectOptipar:
    with resources.files("toolsrtm.data").joinpath("optipar_fluspect.csv").open("r", encoding="utf-8") as f:
        f.readline()
        d = np.loadtxt(f, delimiter=",")
    cols = ("wl", "nr", "Kab", "Kca", "Ks", "Kw", "Kdm", "phiI", "phiII", "KcaV", "KcaZ", "Kant", "Kp", "Kcbc")
    return _FluspectOptipar(**{c: d[:, i] for i, c in enumerate(cols)})


def _kca_from_cx(op: _FluspectOptipar, Cx: float) -> np.ndarray:
    if Cx == -999:
        return op.Kca
    return (1 - Cx) * op.KcaV + Cx * op.KcaZ


def _prospect_mesophyll(Kall, N, nr):
    """Shared PROSPECT-with-interfaces-removed core (both fluspect_b/_cx):
    single-layer + N-layer Stokes combination, then strip the leaf-air
    interfaces to get the bare mesophyll rho/tau used by the doubling
    routine, plus the Kubelka-Munk k/s coefficients. Returns a dict of all
    the intermediate spectra later functions need (talf, r21, t21, rho,
    tau, k, s, ...).
    """
    Ra, Ta, r, t = _one_layer(Kall, nr, alpha=59)
    Rsub, Tsub = _stokes_n_layers(r, t, N)

    denom = 1 - Rsub * r
    tran = Ta * Tsub / denom
    refl = Ra + (Ta * Rsub * t) / denom

    talf = calctav(59, nr)
    ralf = 1 - talf
    t21 = calctav(90, nr) / (nr**2)
    r21 = 1 - t21

    # Remove the top (leaf-air) interface, isolating the mesophyll layer
    Rb = (refl - ralf) / (talf * t21 + (refl - ralf) * r21)
    Z = tran * (1 - Rb * r21) / (talf * t21)
    rho = (Rb - r21 * Z**2) / (1 - (r21 * Z) ** 2)
    tau = (1 - Rb * r21) / (1 - (r21 * Z) ** 2) * Z
    tt = tau
    rr = np.maximum(rho, 0.0)

    I_rt = (rr + tt) < 1
    D = np.zeros_like(rr)
    D[I_rt] = np.sqrt(
        (1 + rr[I_rt] + tt[I_rt]) * (1 + rr[I_rt] - tt[I_rt]) * (1 - rr[I_rt] + tt[I_rt]) * (1 - rr[I_rt] - tt[I_rt])
    )
    a = np.ones_like(rr)
    b = np.ones_like(rr)
    a[I_rt] = (1 + rr[I_rt] ** 2 - tt[I_rt] ** 2 + D[I_rt]) / (2 * rr[I_rt])
    b[I_rt] = (1 - rr[I_rt] ** 2 + tt[I_rt] ** 2 + D[I_rt]) / (2 * tt[I_rt])

    s = rr / tt
    I_a = (a > 1) & np.isfinite(a) & (a != np.inf)
    I_na = np.isnan(s) | np.isnan(a) | np.isnan(b)
    mask = I_a & ~I_na
    s = np.where(mask, 2 * a / (a**2 - 1) * np.log(b), s)

    k = np.log(b)
    mask_k = I_a & ~I_na
    k = np.where(mask_k, (a - 1) / (a + 1) * np.log(b), k)

    return dict(refl=refl, tran=tran, talf=talf, r21=r21, t21=t21, rho=rho, tau=tau, k=k, s=s)


@dataclass
class FluspectBResult:
    lambda_: np.ndarray  # 400-2400 nm
    refl: np.ndarray
    tran: np.ndarray
    kChlrel: np.ndarray
    MbI: np.ndarray  # (211, 351): backward-scattering fluorescence matrix, PSI
    MbII: np.ndarray  # PSII
    MfI: np.ndarray  # forward-scattering, PSI
    MfII: np.ndarray  # PSII


def fluspect_b(
    Cab: float, Car: float, EWT: float, LMA: float, Cs: float, N: float, fqe: float, Cx: float,
    Prot: float | None = None, CBC: float | None = None, Anth: float | None = None,
) -> FluspectBResult:
    """FLUSPECT-B leaf model. Direct port of ``ToolsRTM::getFluspect.B``.

    Parameters
    ----------
    Cab, Car : float
        Chlorophyll a+b, carotenoid content (ug/cm2).
    EWT, LMA : float
        Equivalent water thickness, leaf mass per area (g/cm2 both).
    Cs : float
        Senescent/brown pigment content.
    N : float
        Leaf structure parameter.
    fqe : float
        Fluorescence quantum efficiency. If ``fqe <= 0``, ``MbI``/``MbII``/
        ``MfI``/``MfII`` come back all-zero (no fluorescence) -- ``refl``/
        ``tran`` are still computed and returned normally either way, since
        they never depend on fqe at all.
    Cx : float
        Violaxanthin-zeaxanthin transition state, 0-1 (or -999 to use a
        fixed carotenoid absorption spectrum instead of the Cx-interpolated
        one).
    Prot, CBC, Anth : float, optional
        If all three are given, matches R's ``"Prot" %in% colnames(inputsLeaf)``
        auto-detection: the PROSPECT-PRO ``Kall`` formula (with Anth/Prot/CBC
        terms) is used instead of the plain PROSPECT-D one. Note this
        auto-detection *overrides* any notion of a ``version`` argument in
        R too -- there is no way to force PROSPECT-D once these are supplied,
        matching R exactly.

    Returns
    -------
    FluspectBResult

    Notes
    -----
    R's ``ToolsRTM::getFluspect.B`` used to have ``return(LRT)`` *inside*
    its ``if (fqe_ > 0) {...}`` block -- since refl/tran are computed
    unconditionally above that block but the function's only return
    statement lived inside it, ``fqe_ == 0`` made the whole function
    silently return R's ``NULL`` instead of a proper result (an implicit
    fall-through, not an error), which propagated into a zero-length
    rdot/rsot several calls downstream instead of failing loudly. This port
    originally mirrored that by returning ``None`` and required callers to
    check for it (:func:`toolsrtm.canopy._leaf_optics` used to raise
    ``ValueError`` on ``None`` rather than silently propagating a NULL like
    R did). Now that the R source computes and returns refl/tran/zero-Mb/Mf
    unconditionally (see ``ToolsRTM/R/fluspect_B.R``'s own comment), this
    port matches that instead of raising.
    """
    op = _fluspect_optipar()
    is_pro = Prot is not None and CBC is not None and Anth is not None
    if is_pro and LMA > 0 and (Prot > 0 or CBC > 0):
        LMA = 0.0

    Kca = _kca_from_cx(op, Cx)
    if is_pro:
        Kall = (Cab * op.Kab + Car * Kca + LMA * op.Kdm + EWT * op.Kw + Cs * op.Ks
                + Anth * op.Kant + Prot * op.Kp + CBC * op.Kcbc) / N
    else:
        Kall = (Cab * op.Kab + Car * Kca + LMA * op.Kdm + EWT * op.Kw + Cs * op.Ks) / N

    j = Kall > 0
    kChlrel = np.where(j, Cab * op.Kab / (np.where(j, Kall, 1.0) * N), 0.0)

    core = _prospect_mesophyll(Kall, N, op.nr)
    refl, tran = core["refl"], core["tran"]

    if fqe <= 0:
        zeros = np.zeros((len(_WLF), len(_WLE)))
        return FluspectBResult(
            lambda_=op.wl.copy(), refl=refl, tran=tran, kChlrel=kChlrel,
            MbI=zeros, MbII=zeros.copy(), MfI=zeros.copy(), MfII=zeros.copy(),
        )

    fqe_vec = np.array([fqe / 5, fqe])
    talf, r21, t21 = core["talf"], core["r21"], core["t21"]
    rho, tau, k, s = core["rho"], core["tau"], core["k"], core["s"]
    kChl = kChlrel * k

    wle = op.wl[_IWLE]
    wlf = op.wl[_IWLF_CORRECT]
    k_e, s_e, kChl_e = k[_IWLE], s[_IWLE], kChl[_IWLE]
    k_f, s_f = k[_IWLF_CORRECT], s[_IWLF_CORRECT]

    ndub = 15
    eps = 2.0 ** (-ndub)
    te = 1 - (k_e + s_e) * eps
    tf = 1 - (k_f + s_f) * eps
    re = s_e * eps
    rf = s_f * eps

    sigmoid = 1.0 / (1.0 + np.outer(np.exp(-wlf / 10), np.exp(wle / 10)))  # (nwlf, nwle)

    MfI = MbI = fqe_vec[0] * np.outer(0.5 * op.phiI[_IWLF_CORRECT] * eps, kChl_e) * sigmoid
    MfII = MbII = fqe_vec[1] * np.outer(0.5 * op.phiII[_IWLF_CORRECT] * eps, kChl_e) * sigmoid

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

        MfI, MbI = MfI * A11 + MbI * A12, MbI * A21 + MfI * A22
        MfII, MbII = MfII * A11 + MbII * A12, MbII * A21 + MfII * A22

        te, re, tf, rf = ten, ren, tfn, rfn

    Rb = rho + tau**2 * r21 / (1 - rho * r21)
    Xe = talf[_IWLE] / (1 - r21[_IWLE] * Rb[_IWLE])
    Xf = t21[_IWLF_CORRECT] / (1 - r21[_IWLF_CORRECT] * Rb[_IWLF_CORRECT])
    Ye = tau[_IWLE] * r21[_IWLE] / (1 - rho[_IWLE] * r21[_IWLE])
    Yf = tau[_IWLF_CORRECT] * r21[_IWLF_CORRECT] / (1 - rho[_IWLF_CORRECT] * r21[_IWLF_CORRECT])

    A = Xe[None, :] * (1 + Ye[None, :] * Yf[:, None]) * Xf[:, None]
    B = Xe[None, :] * (Ye[None, :] + Yf[:, None]) * Xf[:, None]

    MbI_n = A * MbI + B * MfI
    MfI_n = A * MfI + B * MbI
    MbII_n = A * MbII + B * MfII
    MfII_n = A * MfII + B * MbII

    return FluspectBResult(
        lambda_=op.wl.copy(), refl=refl, tran=tran, kChlrel=kChlrel,
        MbI=MbI_n, MbII=MbII_n, MfI=MfI_n, MfII=MfII_n,
    )


@dataclass
class FluspectCxResult:
    lambda_: np.ndarray
    refl: np.ndarray
    tran: np.ndarray
    kChlrel: np.ndarray
    kCarrel: np.ndarray
    Mb: np.ndarray  # (211, 351)
    Mf: np.ndarray


def fluspect_cx(
    Cab: float, Car: float, EWT: float, LMA: float, Cs: float, N: float, fqe: float, Cx: float,
    Prot: float, CBC: float, Anth: float,
) -> FluspectCxResult:
    """FLUSPECT-B-Cx leaf model. Direct port of ``ToolsRTM::getFluspect.Cx``.

    Always uses the PROSPECT-PRO ``Kall`` formula (Prot/CBC/Anth required,
    no auto-detection like :func:`fluspect_b`). If ``fqe <= 0``, ``Mb``/
    ``Mf`` come back all-zero (no fluorescence) -- ``refl``/``tran`` are
    still computed and returned normally either way (see :func:`fluspect_b`'s
    own docstring for why this isn't a ``None`` return, like it used to be).

    **History**: R's ``ToolsRTM::getFluspect.Cx`` used to have a real
    indexing bug here (``Iwlf <- intersect(wlp, wlf)`` returned the
    *wavelength values* 640:850 used positionally against ``wlp``
    (400:2400), reading Kubelka-Munk coefficients at 1039-1249 nm instead
    of the intended 640-850 nm fluorescence emission range -- confirmed via
    a standalone R repro). ``ToolsRTM::getFluspect.B`` (:func:`fluspect_b`)
    never had this bug (uses ``which()``), nor does SCOPEinR's own
    ``fluspect_Cx_forSCOPE.R`` (uses ``match()``). Fixed in the R source
    (``Iwlf <- which(wlp >= min(wlf) & wlp <= max(wlf))``) and this port
    updated to match; this function now computes the physically-intended
    640-850 nm range (``_IWLF_CORRECT``), not the historical buggy one.
    """
    op = _fluspect_optipar()
    if LMA > 0 and (Prot > 0 or CBC > 0):
        LMA = 0.0

    Kca = _kca_from_cx(op, Cx)
    Kall = (Cab * op.Kab + Car * Kca + LMA * op.Kdm + EWT * op.Kw + Cs * op.Ks
            + Anth * op.Kant + Prot * op.Kp + CBC * op.Kcbc) / N

    j = Kall > 0
    kChlrel = np.where(j, Cab * op.Kab / (np.where(j, Kall, 1.0) * N), 0.0)
    kCarrel = np.where(j, Car * Kca / (np.where(j, Kall, 1.0) * N), 0.0)

    core = _prospect_mesophyll(Kall, N, op.nr)
    refl, tran = core["refl"], core["tran"]

    if fqe <= 0:
        zeros = np.zeros((len(_WLF), len(_WLE)))
        return FluspectCxResult(
            lambda_=op.wl.copy(), refl=refl, tran=tran, kChlrel=kChlrel, kCarrel=kCarrel,
            Mb=zeros, Mf=zeros.copy(),
        )

    fqe_vec = np.array([fqe / 5, fqe])
    talf, r21, t21 = core["talf"], core["r21"], core["t21"]
    rho, tau, k, s = core["rho"], core["tau"], core["k"], core["s"]
    kChl = kChlrel * k

    wle = op.wl[_IWLE]
    k_e, s_e, kChl_e = k[_IWLE], s[_IWLE], kChl[_IWLE]
    k_f, s_f = k[_IWLF_CORRECT], s[_IWLF_CORRECT]
    wlf_labels = _WLF  # sigmoid/output still keyed by the *intended* 640-850 labels

    ndub = 15
    int_factor = 5
    eps = 2.0 ** (-ndub)
    te = 1 - (k_e + s_e) * eps
    tf = 1 - (k_f + s_f) * eps
    re = s_e * eps
    rf = s_f * eps

    sigmoid = 1.0 / (1.0 + np.outer(np.exp(-wlf_labels / 10), np.exp(wle / 10)))

    Mf = Mb = int_factor * fqe_vec[0] * np.outer(0.5 * op.phiI[_IWLF_CORRECT] * eps, kChl_e) * sigmoid

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
    Xe = talf[_IWLE] / (1 - r21[_IWLE] * Rb[_IWLE])
    Xf = t21[_IWLF_CORRECT] / (1 - r21[_IWLF_CORRECT] * Rb[_IWLF_CORRECT])
    Ye = tau[_IWLE] * r21[_IWLE] / (1 - rho[_IWLE] * r21[_IWLE])
    Yf = tau[_IWLF_CORRECT] * r21[_IWLF_CORRECT] / (1 - rho[_IWLF_CORRECT] * r21[_IWLF_CORRECT])

    A = Xe[None, :] * (1 + Ye[None, :] * Yf[:, None]) * Xf[:, None]
    B = Xe[None, :] * (Ye[None, :] + Yf[:, None]) * Xf[:, None]

    Mb_n = A * Mb + B * Mf
    Mf_n = A * Mf + B * Mb

    return FluspectCxResult(
        lambda_=op.wl.copy(), refl=refl, tran=tran, kChlrel=kChlrel, kCarrel=kCarrel,
        Mb=Mb_n, Mf=Mf_n,
    )
