"""Leaf biochemistry: Farquhar-von Caemmerer-Berry photosynthesis (Collatz
C4 variant) coupled with a Ball-Berry/Leuning stomatal-conductance model and
the van der Tol et al. (2014) fluorescence yield model.

Direct port of ``SCOPEinR/R/biochemical.R`` (``get.biochemical``) and its
helpers in ``SCOPEinR/R/Biochemical_functions.R``. This is the leaf-level
photosynthesis+fluorescence solver called *inside* SCOPE's energy-balance
iteration (``ebal.R``, not ported) to get ``A``/``rcw``/``eta`` at a given
leaf temperature -- it does not itself iterate on temperature, so it can be
called and verified standalone given an assumed leaf micro-environment
(matching how the R function itself works: ``data.meteo$Temp`` is an
input, not something this function solves for).

Only the ``tempcor=1`` (temperature-corrected) C3 path and the
``BallBerry0 != 0`` (iterative Ci) path are ported in full generality here;
the ``BallBerry0 == 0`` closed-form Ci path and the C4/no-temperature-
correction paths are ported too but exercised less by the reference tests
-- see ``python/README.md``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq

from ._data import constants as _constants
from .utils import satvap

__all__ = [
    "LeafBio",
    "MeteoLeaf",
    "BiochemResult",
    "sel_root",
    "get_gs_fun",
    "get_ball_berry",
    "get_temperature_function_c3",
    "get_high_temp_inhibtion_c3",
    "get_fluorescence_model",
    "get_ci_next",
    "get_compute_a",
    "get_biochemical",
]


def sel_root(a, b, c, dsign):
    """Root of least magnitude of ``a*x^2 + b*x + c = 0``. Direct port of
    ``SCOPEinR::sel_root``. ``dsign``: -1/0 picks the smaller root, +1 the
    larger (per quadratic-formula sign convention on the discriminant)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    c = np.asarray(c, dtype=float)
    dsign = np.asarray(dsign, dtype=float).copy()
    dsign = np.where(dsign == 0, -1.0, dsign)
    disc = np.sqrt(b**2 - 4 * a * c)
    x = (-b + dsign * disc) / (2 * a)
    return np.where(a == 0, -c / b, x)


def get_gs_fun(Cs, RH, A, BallBerrySlope, BallBerry0):
    """Ball-Berry stomatal conductance. Direct port of ``SCOPEinR::get.gsFun``."""
    Cs = np.asarray(Cs, dtype=float)
    gs = np.maximum(BallBerry0, BallBerrySlope * A * RH / (Cs + 1e-9) + BallBerry0)
    gs = np.where(np.isnan(Cs) | np.isinf(Cs), np.nan, gs)
    return gs


def get_ball_berry(Cs, RH, A, BallBerrySlope, BallBerry0, minCi, Ci_input=None):
    """Ball-Berry/Leuning Ci and (optionally) gs. Direct port of
    ``SCOPEinR::get.BallBerry``. Returns ``(gs, Ci)`` (``gs`` is ``None``
    when not computable, matching R's ``NULL``)."""
    Cs = np.asarray(Cs, dtype=float)
    if Ci_input is not None:
        Ci = np.asarray(Ci_input, dtype=float)
        gs = get_gs_fun(Cs, RH, A, BallBerrySlope, BallBerry0) if A is not None else None
        return gs, Ci
    if np.all(np.asarray(BallBerry0) == 0) or A is None:
        Ci = np.maximum(minCi * Cs, Cs * (1 - 1.6 / (BallBerrySlope * RH)))
        return None, Ci
    gs = get_gs_fun(Cs, RH, A, BallBerrySlope, BallBerry0)
    Ci = np.maximum(minCi * Cs, Cs - 1.6 * A / gs)
    return gs, Ci


def get_temperature_function_c3(Tref, R, Temp, deltaHa):
    """Arrhenius temperature correction factor. Direct port of
    ``SCOPEinR::get.temperature.functionC3``."""
    return np.exp(deltaHa / (Tref * R) * (1 - Tref / Temp))


def get_high_temp_inhibtion_c3(Tref, R, T, deltaS, deltaHd):
    """High-temperature inhibition factor. Direct port of
    ``SCOPEinR::get.high.temp.inhibtionC3``."""
    num = 1 + np.exp((Tref * deltaS - deltaHd) / (Tref * R))
    den = 1 + np.exp((deltaS * T - deltaHd) / (R * T))
    return num / den


def get_fluorescence_model(ps, x, Kp, Kf, Kd, Knparams):
    """van der Tol et al. (2014) fluorescence-yield model. Direct port of
    ``SCOPEinR::get.Fluorescence.model``. Returns a dict with eta, qE, qQ,
    fs, fo, fm, fo0, fm0, Kn."""
    ps = np.asarray(ps, dtype=float)
    x = np.asarray(x, dtype=float)
    Kno, alpha, beta = Knparams

    x_alpha = np.exp(np.log(x) * alpha)
    Kn = Kno * (1 + beta) * x_alpha / (beta + x_alpha)

    fo0 = Kf / (Kf + Kp + Kd)
    fo = Kf / (Kf + Kp + Kd + Kn)
    fm = Kf / (Kf + Kd + Kn)
    fm0 = Kf / (Kf + Kd)
    fs = fm * (1 - ps)
    eta = fs / fo0
    qQ = 1 - (fs - fo) / (fm - fo)
    qE = 1 - (fm - fo) / (fm0 - fo0)

    return dict(eta=eta, qE=qE, qQ=qQ, fs=fs, fo=fo, fm=fm, fo0=fo0, fm0=fm0, Kn=Kn)


def get_compute_a(Ci, Type, g_m, Vs_C3, MM_consts, Rd, Vcmax, Gamma_star, Je, effcon, atheta, kpepcase):
    """Farquhar (C3) / Collatz (C4) net CO2 assimilation. Direct port of
    ``SCOPEinR::get.computeA``. Returns a dict with A, Ag, Vc, Vs, Ve,
    CO2_per_electron (``fcount`` -- a debug iteration counter via R's
    ``<<-`` -- is not reproduced; it has no effect on the physics)."""
    Ci = np.asarray(Ci, dtype=float)
    if Type == "C3":
        Vs = Vs_C3
        if np.any(np.asarray(g_m) < np.inf):
            Vc = sel_root(1 / g_m, -(MM_consts + Ci + (Rd + Vcmax) / g_m),
                           Vcmax * (Ci - Gamma_star + Rd / g_m), -1)
            Ve = sel_root(1 / g_m, -(Ci + 2 * Gamma_star + (Rd + Je * effcon) / g_m),
                           Je * effcon * (Ci - Gamma_star + Rd / g_m), -1)
            CO2_per_electron = Ve / Je
        else:
            Vc = Vcmax * (Ci - Gamma_star) / (MM_consts + Ci)
            CO2_per_electron = (Ci - Gamma_star) / (Ci + 2 * Gamma_star) * effcon
            Ve = Je * CO2_per_electron
    else:  # C4
        Vc = Vcmax
        Vs = kpepcase * Ci
        CO2_per_electron = effcon
        Ve = Je * CO2_per_electron

    V = sel_root(atheta, -(Vc + Ve), Vc * Ve, np.sign(-Vc))
    Ag = sel_root(0.98, -(V + Vs), V * Vs, -1)
    A = Ag - Rd

    return dict(A=A, Ag=Ag, Vc=Vc, Vs=Vs, Ve=Ve, CO2_per_electron=CO2_per_electron)


def get_ci_next(Ci_in, Cs, RH, minCi, BallBerrySlope, BallBerry0, A_fun, ppm2bar):
    """Ci fixed-point residual (Ball-Berry Ci minus guessed Ci_in), used as
    the objective for the Brent root-finder in :func:`get_biochemical`.
    Direct port of ``SCOPEinR::get.Ci.next``."""
    av = A_fun(Cs)
    A_bar = None if av["A"] is None else av["A"] * ppm2bar
    gs, Ci_out = get_ball_berry(Cs, RH, A_bar, BallBerrySlope, BallBerry0, minCi)
    err = Ci_out - Ci_in
    return err, (gs, Ci_out)


@dataclass
class LeafBio:
    """Leaf biochemical parameters (``data.leafbio`` in R)."""

    Type: str  # 'C3' or 'C4'
    stressfactor: float
    Vcmax25: float
    BallBerry0: float
    BallBerrySlope: float
    Rdparam: float
    Kn0: float
    Knalpha: float
    Knbeta: float
    g_m: float | None = None  # mol m-2 s-1 bar-1; None -> Inf (no mesophyll-conductance effect)
    TDP: dict = field(default_factory=dict)  # temperature-dependence params, see get_biochemical


@dataclass
class MeteoLeaf:
    """Leaf micro-environment (``data.meteo`` in R)."""

    Q: float  # absorbed PAR, umol photons m-2 s-1
    Cs: float  # CO2 at the leaf boundary layer, ppm
    Temp: float  # leaf temperature, deg C or K
    eb: float  # vapour pressure in the leaf boundary layer, hPa
    Oa: float  # O2 concentration, mmol/mol
    p: float  # air pressure, hPa


@dataclass
class BiochemResult:
    A: np.ndarray
    Ci: np.ndarray
    Cc: np.ndarray | None
    rcw: np.ndarray
    gs: np.ndarray
    RH: np.ndarray
    Vcmax: np.ndarray
    Rd: np.ndarray
    Ja: np.ndarray
    ps: np.ndarray
    ps_rel: np.ndarray
    Kd: np.ndarray
    Kn: np.ndarray
    NPQ: np.ndarray
    Kf: float
    Kp0: float
    Kp: np.ndarray
    eta: np.ndarray
    qE: np.ndarray
    fs: np.ndarray
    SIF: np.ndarray
    fo0: np.ndarray
    fm0: np.ndarray
    fo: np.ndarray
    fm: np.ndarray
    qQ: np.ndarray
    Phi_N: np.ndarray


def get_biochemical(leafbio: LeafBio, meteo: MeteoLeaf, temp_correction: bool, fV: float = 1.0) -> BiochemResult:
    """Leaf-level photosynthesis (Farquhar/Collatz) + fluorescence yield
    (van der Tol et al. 2014). Direct port of ``SCOPEinR::get.biochemical``.

    Parameters
    ----------
    leafbio : LeafBio
    meteo : MeteoLeaf
    temp_correction : bool
        Whether to apply temperature correction to Vcmax/Rd/Kc/Ko/Gamma_star
        (matches R's ``data.opts`` row-7 ``tempcor`` flag). If True,
        ``leafbio.TDP`` must contain the relevant temperature-dependence
        parameters (C3: ``delHaV``/``delSV``/``delHdV``/``delHaR``/``delSR``/
        ``delHdR``/``delHaKc``/``delHaKo``/``delHaT``; C4: ``Q10``/``s1``-``s6``).
    fV : float, default 1.0
        Scaling factor on ``Vcmax25`` (e.g. a canopy N/Vcmax profile factor).

    Returns
    -------
    BiochemResult
    """
    const = _constants()
    rhoa, Mair, R = const["rhoa"], const["Mair"], const["R"]

    Q = np.asarray(meteo.Q, dtype=float)
    Cs = np.asarray(meteo.Cs, dtype=float)
    Temp = np.asarray(meteo.Temp, dtype=float)
    T_k = np.where(Temp < 200, Temp + 273.15, Temp)
    eb = np.asarray(meteo.eb, dtype=float)
    O = np.asarray(meteo.Oa, dtype=float)
    p = np.asarray(meteo.p, dtype=float)

    Type = leafbio.Type
    stressfactor = leafbio.stressfactor
    Vcmax25 = fV * leafbio.Vcmax25
    BallBerry0 = leafbio.BallBerry0
    BallBerrySlope = leafbio.BallBerrySlope
    RdPerVcmax25 = leafbio.Rdparam
    effcon = 1 / 5 if Type == "C3" else 1 / 6

    Tref = 25 + 273.15
    Kc25 = 405e-6
    Ko25 = 279e-3
    spfy25 = 2444

    ppm2bar = 1e-6 * (p * 1e-3)
    Cs_bar = Cs * ppm2bar
    O_bar = (O * 1e-3) * (p * 1e-3) * (1.0 if Type == "C3" else 0.0)
    Gamma_star25 = 0.5 * O_bar / spfy25
    Rd25 = RdPerVcmax25 * Vcmax25

    atheta = 0.8
    g_m = np.inf if leafbio.g_m is None else leafbio.g_m * 1e6

    Knparams = (leafbio.Kn0, leafbio.Knalpha, leafbio.Knbeta)
    Kf = 0.05
    Kd = np.maximum(0.8738, 0.0301 * (T_k - 273.15) + 0.0773)
    Kp = 4.0

    fl = dict(Vcmax=1.0, Rd=1.0, TPU=1.0, Kc=1.0, Ko=1.0, Gamma_star=1.0)
    Ke = 1.0

    if temp_correction:
        tdp = leafbio.TDP
        if Type == "C4":
            Q10, s1, s2, s3, s4 = tdp["Q10"], tdp["s1"], tdp["s2"], tdp["s3"], tdp["s4"]
            s5, s6 = tdp["s5"], tdp["s6"]
            fHTv = 1 + np.exp(s1 * (T_k - s2))
            fLTv = 1 + np.exp(s3 * (s4 - T_k))
            Vcmax = (Vcmax25 * Q10 ** (0.1 * (T_k - Tref))) / (fHTv * fLTv)
            fHTv = 1 + np.exp(s5 * (T_k - s6))
            Rd = (Rd25 * Q10 ** (0.1 * (T_k - Tref))) / fHTv
            Ke25 = 20000 * Vcmax25
            Ke = Ke25 * Q10 ** (0.1 * (T_k - Tref))
        elif Type == "C3":
            fTv = get_temperature_function_c3(Tref, R, T_k, tdp["delHaV"])
            fHTv = get_high_temp_inhibtion_c3(Tref, R, T_k, tdp["delSV"], tdp["delHdV"])
            fl["Vcmax"] = fTv * fHTv

            fTv = get_temperature_function_c3(Tref, R, T_k, tdp["delHaR"])
            fHTv = get_high_temp_inhibtion_c3(Tref, R, T_k, tdp["delSR"], tdp["delHdR"])
            fl["Rd"] = fTv * fHTv

            fl["Kc"] = get_temperature_function_c3(Tref, R, T_k, tdp["delHaKc"])
            fl["Ko"] = get_temperature_function_c3(Tref, R, T_k, tdp["delHaKo"])
            fl["Gamma_star"] = get_temperature_function_c3(Tref, R, T_k, tdp["delHaT"])
            Ke = 1.0

    if Type == "C3":
        Vcmax = Vcmax25 * fl["Vcmax"] * stressfactor
        Rd = Rd25 * fl["Rd"] * stressfactor
        Kc = Kc25 * fl["Kc"]
        Ko = Ko25 * fl["Ko"]
    Gamma_star = Gamma_star25 * fl["Gamma_star"]

    po0 = Kp / (Kf + Kd + Kp)
    Je = 0.5 * po0 * Q

    if Type == "C3":
        MM_consts = Kc * (1 + O_bar / Ko)
        Vs_C3 = Vcmax / 2
        minCi = 0.3
    else:
        MM_consts = 0.0
        Vs_C3 = 0.0
        minCi = 0.1

    RH = np.minimum(1.0, eb / satvap(T_k - 273.15))

    def compute_a_fun(x):
        return get_compute_a(x, Type, g_m, Vs_C3, MM_consts, Rd, Vcmax, Gamma_star, Je, effcon, atheta, Ke)

    Cs_bar_arr = np.atleast_1d(Cs_bar)
    RH_arr = np.broadcast_to(np.atleast_1d(RH), Cs_bar_arr.shape)

    if np.all(np.asarray(BallBerry0) == 0):
        _, Ci = get_ball_berry(Cs_bar_arr, RH_arr, None, BallBerrySlope, BallBerry0, minCi)
    else:
        Ci = np.empty_like(Cs_bar_arr)
        for i in range(Cs_bar_arr.size):
            lower, upper = Cs_bar_arr.flat[i] - 0.001, Cs_bar_arr.flat[i] + 0.001

            def obj(x, i=i):
                err, _ = get_ci_next(x, Cs_bar_arr.flat[i], RH_arr.flat[i], minCi,
                                      BallBerrySlope, BallBerry0, compute_a_fun, ppm2bar)
                return err
            Ci.flat[i] = brentq(obj, lower, upper, xtol=1e-7, maxiter=1000)
        Ci = Ci.reshape(Cs_bar_arr.shape)
        if np.isscalar(Cs_bar) or np.ndim(Cs_bar) == 0:
            Ci = Ci.reshape(())

    params = compute_a_fun(Ci)
    A = params["A"]
    Ag = params["Ag"]
    CO2_per_electron = params["CO2_per_electron"]

    gs = np.maximum(0.0, 1.6 * A * ppm2bar / (Cs_bar_arr.reshape(np.shape(A)) - Ci))
    Ja = Ag / CO2_per_electron
    rcw = (rhoa / (Mair * 1e-3)) / gs

    ps = po0 * Ja / Je
    ps = np.where(np.isnan(ps), po0, ps)
    ps_rel = np.maximum(0.0, 1 - ps / po0)

    fluo = get_fluorescence_model(ps, ps_rel, Kp, Kf, Kd, Knparams)
    eta, qE, qQ = fluo["eta"], fluo["qE"], fluo["qQ"]
    fs, fo, fm, fm0, fo0 = fluo["fs"], fluo["fo"], fluo["fm"], fluo["fm0"], fluo["fo0"]
    Kn = fluo["Kn"]

    Kpa = ps / fs * Kf

    Cc = None
    if g_m is not None:
        Cc = (Ci - A / g_m) / ppm2bar
    Ci_ppm = Ci / ppm2bar

    return BiochemResult(
        A=A, Ci=Ci_ppm, Cc=Cc, rcw=rcw, gs=gs, RH=RH, Vcmax=Vcmax, Rd=Rd, Ja=Ja,
        ps=ps, ps_rel=ps_rel, Kd=Kd, Kn=Kn, NPQ=Kn / (Kf + Kd), Kf=Kf, Kp0=Kp, Kp=Kpa,
        eta=eta, qE=qE, fs=fs, SIF=fs * Q, fo0=fo0, fm0=fm0, fo=fo, fm=fm, qQ=qQ,
        Phi_N=Kn / (Kn + Kp + Kf + Kd),
    )
