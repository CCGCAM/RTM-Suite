"""RTMo: optical top-of-canopy bidirectional reflectance pipeline.

Direct, function-by-function port of the *optical BRDF* portion of
``SCOPEinR/R/RTMo.R`` and its helpers in ``SCOPEinR/R/RTMo_functions.R``
(``get.volscatt.scope``, ``get.Pso``, ``get.reflectances``,
``get.fluxprofile``). Ported outputs: the four-stream TOC reflectance
factors (rdd, rsd, rdo, rso), the apparent TOC reflectance (refl), the
TOC radiance in viewing direction (``Lo_``), the outgoing top-of-canopy flux
(``Eout_``/Eouto/Eoutt/Lot), and the gap probabilities (Ps, Po, Pso, k, K).

Also includes :func:`net_radiation_lite`, a **partial** port of RTMo.R's
section 4 (PAR / net-radiation absorption breakdown) -- specifically just
``Rnuc``/``Rnhc``/``Rnus``/``Rnhs``/``Pnu_Cab``/``Pnh_Cab``, the six
quantities :mod:`scopeinpython.ebal` actually consumes; the "lite"
(scalar-per-layer, not full ``(13,36,nl)`` per-leaf-angle) branch only,
matching every reference case in this port.

**A real bug found and fixed in R here**: the "lite" branch's direct-beam
term (``Rndir``/``Pndir_Cab``/etc in ``RTMo.R``) used to compute
``fs_ * Asun[j]`` where ``j`` was leftover from an unrelated, already-
completed ``for`` loop above it (always equal to ``nl``, the *last*
canopy layer) -- so the direct-beam contribution to sunlit-leaf net
radiation came out **exactly constant across all canopy layers**
(confirmed via direct repro: ``data.rad$Rnuc - data.rad$Rnhc`` had zero
layer-to-layer variation), which is not physically sensible -- direct-
beam absorption should decay with canopy depth, the same as the diffuse
term a few lines below already correctly does. Fixed in
``SCOPEinR/R/RTMo.R`` to use the full per-layer vectors (``Asun`` etc,
not ``Asun[j]``) -- a scalar ``fs_`` broadcasts against them correctly
without indexing. A second, downstream instance of the same class of
mistake (``Rndir[1] + Rndif`` picking only the first layer, harmless
before the fix since ``Rndir`` was a scalar, wrong after) was fixed in
the same edit.

NOT ported (out of scope for this port, see python/README.md):
  - Everything else in RTMo.R's section 4 (``Rnuc_Car``/``Pnuc_Car``/
    ``Rnuc_PAR``/``Rnhc_PAR``, top-of-canopy incident PAR ``P``/``EPAR``,
    the full ``(13,36,nl)`` per-leaf-angle branch) -- not consumed by
    :mod:`scopeinpython.ebal`, so out of scope for now;
  - the MODTRAN-atmospheric-file branch of ``get.calcTOCirr`` (only the
    "precomputed ``Esun_``/``Esky_``" branch, i.e. the default SCOPE example
    irradiance, is ported);
  - mSCOPE per-layer leaf property variation (a single leaf-optics
    spectrum is broadcast to all ``nl`` canopy layers, as in the R
    example script here, i.e. mly$nly == 1).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad

from ._data import constants as _constants
from .spectral import SpectralConfig
from .utils import get_e2phot, sint

__all__ = [
    "get_volscatt_scope",
    "get_pso",
    "get_reflectances",
    "get_fluxprofile",
    "CanopyStructure",
    "RTMoResult",
    "run_rtmo",
    "NetRadiationLite",
    "net_radiation_lite",
]


# ---------------------------------------------------------------------------
# Geometry / gap-probability helpers
# ---------------------------------------------------------------------------


def get_volscatt_scope(tts: float, tto: float, psi: float, ttli: np.ndarray):
    """Volume scattering phase functions and interception coefficients,
    vectorised over the leaf inclination classes ``ttli``.

    Direct port of ``SCOPEinR::get.volscatt.scope`` (distinct from
    ``toolsrtm.canopy.volscatt``, the scalar-``ttl`` variant used by
    fourSAIL; SCOPE's RTMo uses this vectorised formulation instead).

    Returns
    -------
    dict with keys chi_s, chi_o, frho, ftau (arrays, same length as ttli).
    """
    deg2rad = np.pi / 180.0
    ttli = np.asarray(ttli, dtype=float)
    nli = len(ttli)

    psi_rad = psi * deg2rad * np.ones(nli)
    cos_psi = np.cos(psi * deg2rad)

    cos_ttli = np.cos(ttli * deg2rad)
    sin_ttli = np.sin(ttli * deg2rad)

    cos_tts = np.cos(tts * deg2rad)
    sin_tts = np.sin(tts * deg2rad)
    cos_tto = np.cos(tto * deg2rad)
    sin_tto = np.sin(tto * deg2rad)

    Cs = cos_ttli * cos_tts
    Ss = sin_ttli * sin_tts
    Co = cos_ttli * cos_tto
    So = sin_ttli * sin_tto

    As = np.maximum(Ss, Cs)
    Ao = np.maximum(So, Co)

    bts = np.arccos(-Cs / As)
    bto = np.arccos(-Co / Ao)

    chi_o = 2 / np.pi * ((bto - np.pi / 2) * Co + np.sin(bto) * So)
    chi_s = 2 / np.pi * ((bts - np.pi / 2) * Cs + np.sin(bts) * Ss)

    delta1 = np.abs(bts - bto)
    delta2 = np.pi - np.abs(bts + bto - np.pi)

    Tot = psi_rad + delta1 + delta2

    bt1 = np.minimum(psi_rad, delta1)
    bt3 = np.maximum(psi_rad, delta2)
    bt2 = Tot - bt1 - bt3

    T1 = 2 * Cs * Co + Ss * So * cos_psi
    T2 = np.sin(bt2) * (2 * As * Ao + Ss * So * np.cos(bt1) * np.cos(bt3))

    Jmin = bt2 * T1 - T2
    Jplus = (np.pi - bt2) * T1 + T2

    frho = Jplus / (2 * np.pi**2)
    ftau = -Jmin / (2 * np.pi**2)

    frho = np.maximum(0.0, frho)
    ftau = np.maximum(0.0, ftau)

    return {"chi_s": chi_s, "chi_o": chi_o, "frho": frho, "ftau": ftau}


def get_pso(K: float, k: float, LAI: float, q: float, dso: float, xl: float) -> float:
    """Bi-directional gap probability at normalized canopy depth ``xl``.

    Direct port of ``SCOPEinR::get.Pso``.
    """
    if dso != 0:
        alf = (dso / q) * 2 / (k + K)
        return np.exp((K + k) * LAI * xl + np.sqrt(K * k) * LAI / alf * (1 - np.exp(xl * alf)))
    else:
        return np.exp((K + k) * LAI * xl - np.sqrt(K * k) * LAI * xl)


def get_reflectances(tau_ss, tau_sd, tau_dd, rho_dd, rho_sd, rsoil, nl: int, nwl: int):
    """Propagate thin-layer reflectance/transmittance down through ``nl``
    canopy layers to the soil and back, producing the directional-
    hemispherical (R_sd) and hemispherical-hemispherical (R_dd) reflectance
    at the top of each layer (and the soil, layer index nl).

    Direct port of ``SCOPEinR::get.reflectances``. All spectral inputs are
    (nl, nwl) arrays except ``rsoil`` (nwl,); ``tau_ss`` may be (nl, nwl)
    or a scalar/((nl,) broadcastable) as in the R caller (constant per
    layer for a homogeneous canopy).

    Returns
    -------
    dict with R_sd, R_dd (nl+1, nwl), Xss (nl,), Xsd, Xdd (nl, nwl).
    """
    tau_ss = np.broadcast_to(np.asarray(tau_ss, dtype=float), (nl, nwl))
    tau_sd = np.asarray(tau_sd, dtype=float)
    tau_dd = np.asarray(tau_dd, dtype=float)
    rho_dd = np.asarray(rho_dd, dtype=float)
    rho_sd = np.asarray(rho_sd, dtype=float)
    rsoil = np.asarray(rsoil, dtype=float)

    R_sd = np.zeros((nl + 1, nwl))
    R_dd = np.zeros((nl + 1, nwl))
    Xsd = np.zeros((nl, nwl))
    Xdd = np.zeros((nl, nwl))
    Xss = np.zeros(nl)

    R_sd[nl, :] = rsoil
    R_dd[nl, :] = rsoil

    for j in range(nl - 1, -1, -1):
        Xss[j] = tau_ss[j, 0]
        dnorm = 1 - rho_dd[j, :] * R_dd[j + 1, :]
        Xsd[j, :] = (tau_sd[j, :] + tau_ss[j, :] * R_sd[j + 1, :] * rho_dd[j, :]) / dnorm
        Xdd[j, :] = tau_dd[j, :] / dnorm
        R_sd[j, :] = rho_sd[j, :] + tau_dd[j, :] * (R_sd[j + 1, :] * Xss[j] + R_dd[j + 1, :] * Xsd[j, :])
        R_dd[j, :] = rho_dd[j, :] + tau_dd[j, :] * R_dd[j + 1, :] * Xdd[j, :]

    return {"R_sd": R_sd, "R_dd": R_dd, "Xss": Xss, "Xsd": Xsd, "Xdd": Xdd}


def get_fluxprofile(Esun_, Esky_, rsoil, Xss, Xsd, Xdd, R_sd, R_dd, nl: int, nwl: int, rs_thermal: float = 0.06):
    """Propagate top-of-canopy direct/diffuse irradiance down through the
    canopy (and back up) to the vertical flux profile.

    Direct port of ``SCOPEinR::get.fluxprofile`` (the ``nwl==2162``/no
    spectral-padding path; the R function's dim==2001 branch is dead code
    in the standard pipeline used here, since leaf/soil optics are already
    built at full ``nwl`` width before this is called -- see
    :func:`run_rtmo`).

    Returns
    -------
    dict with ``Es_``, ``Emin_``, ``Eplu_`` ((nl+1, nwl) arrays).
    """
    Es_ = np.full((nl + 1, nwl), rs_thermal)
    Emin_ = np.full((nl + 1, nwl), rs_thermal)
    Eplu_ = np.full((nl + 1, nwl), rs_thermal)

    Es_[0, :] = Esun_
    Emin_[0, :] = Esky_

    for j in range(nl):
        Es_[j + 1, :] = Xss[j] * Es_[j, :]
        Emin_[j + 1, :] = Xsd[j, :] * Es_[j, :] + Xdd[j, :] * Emin_[j, :]
        Eplu_[j, :] = R_sd[j, :] * Es_[j, :] + R_dd[j, :] * Emin_[j, :]

    Eplu_[nl, :] = rsoil * (Es_[nl, :] + Emin_[nl, :])

    return {"Es_": Es_, "Emin_": Emin_, "Eplu_": Eplu_}


# ---------------------------------------------------------------------------
# Canopy structure + main pipeline
# ---------------------------------------------------------------------------


@dataclass
class CanopyStructure:
    """Canopy structural inputs for RTMo (subset of ``data.canopy`` in R
    needed by the optical BRDF pipeline)."""

    LAI: float
    lidf: np.ndarray  # length 13, leaf inclination distribution (see toolsrtm.canopy.dladgen/campbell)
    hot: float  # hotspot parameter = leafwidth / hc
    nlayers: int | None = None  # if None, computed as max(2, ceil(10*LAI))
    litab: np.ndarray = None  # 13 standard leaf angle classes (deg); default SCOPE grid if None
    lazitab: np.ndarray = None  # 36 leaf azimuth classes (deg); default SCOPE grid if None
    xl: np.ndarray = None  # nlayers+1 normalized depths [0, -1/nl, ..., -1]; computed if None

    def __post_init__(self):
        if self.nlayers is None:
            self.nlayers = max(2, int(np.ceil(10 * self.LAI)))
        if self.litab is None:
            self.litab = np.concatenate([np.arange(5, 76, 10), np.arange(81, 90, 2)]).astype(float)
        if self.lazitab is None:
            self.lazitab = np.arange(5, 356, 10).astype(float)
        if self.xl is None:
            nl = self.nlayers
            x = np.linspace(-1.0 / nl, -1.0, nl)
            self.xl = np.concatenate([[0.0], x])


@dataclass
class RTMoResult:
    rdd: np.ndarray  # TOC bi-hemispherical reflectance
    rsd: np.ndarray  # TOC directional-hemispherical reflectance (solar incidence)
    rdo: np.ndarray  # TOC hemispherical-directional reflectance (viewing dir)
    rso: np.ndarray  # TOC bi-directional reflectance factor
    refl: np.ndarray  # TOC apparent reflectance
    Lo_: np.ndarray  # TOC radiance in viewing direction (mW m-2 um-1 sr-1)
    Eout_: np.ndarray  # TOC outgoing (upward) flux spectrum (mW m-2 um-1)
    Eouto: float  # TOC outgoing flux, optical range, spectrally integrated (W m-2)
    Eoutt: float  # TOC outgoing flux, thermal range, spectrally integrated (W m-2)
    Lot: float  # TOC radiance, thermal range, spectrally integrated
    Esun_: np.ndarray
    Esky_: np.ndarray
    k: float  # extinction coefficient, solar direction
    K: float  # extinction coefficient, viewing direction
    Ps: np.ndarray  # gap probability, solar direction (length nl+1)
    Po: np.ndarray  # gap probability, viewing direction (length nl+1)
    Pso: np.ndarray  # bi-directional gap probability (length nl+1)
    # Internal per-layer radiative quantities, exposed for RTMf (fluorescence
    # RTM, fluspect_mscope.py's sibling module) -- not part of R's getRTMo
    # *return* value either, but computed as intermediates inside it and
    # needed unchanged by get.RTMf, so kept here rather than recomputed.
    Emin_: np.ndarray  # downward diffuse flux per layer boundary, (nl+1, nwl)
    Eplu_: np.ndarray  # upward diffuse flux per layer boundary, (nl+1, nwl)
    Emins_: np.ndarray  # Emin_'s sun-only component (Esky_=0), (nl+1, nwl)
    Emind_: np.ndarray  # Emin_'s sky-only component (Esun_=0), (nl+1, nwl)
    Eplus_: np.ndarray  # Eplu_'s sun-only component, (nl+1, nwl)
    Eplud_: np.ndarray  # Eplu_'s sky-only component, (nl+1, nwl)
    rho_dd: np.ndarray  # thin-layer diffuse-diffuse reflectance, (nl, nwl)
    R_dd: np.ndarray  # cumulative diffuse-diffuse reflectance, (nl+1, nwl)
    Xdd: np.ndarray  # diffuse-diffuse transmittance factor, (nl, nwl)
    tau_dd: np.ndarray  # thin-layer diffuse-diffuse transmittance, (nl, nwl)
    vb: np.ndarray  # directional backscatter coefficient, (nl, nwl)
    vf: np.ndarray  # directional forward-scatter coefficient, (nl, nwl)
    Xsd: np.ndarray  # direct-to-diffuse transmittance factor, (nl, nwl)
    Xss: np.ndarray  # direct-to-direct transmittance factor, (nl,)
    R_sd: np.ndarray  # cumulative directional-hemispherical reflectance, (nl+1, nwl)


def run_rtmo(
    spectral: SpectralConfig,
    leaf_refl: np.ndarray,
    leaf_tran: np.ndarray,
    rho_thermal: float,
    tau_thermal: float,
    rsoil: np.ndarray,
    canopy: CanopyStructure,
    tts: float,
    tto: float,
    psi: float,
    Esun_: np.ndarray,
    Esky_: np.ndarray,
) -> RTMoResult:
    """Top-of-canopy optical BRDF: leaf optics + soil + canopy structure +
    geometry -> rdd/rsd/rdo/rso/refl (plus ``Lo_``, ``Eout_`` and gap probabilities).

    Direct port of the optical-BRDF portion of ``SCOPEinR:::getRTMo``
    (sections 0-3.3 and the outgoing-radiance block of section 5; the
    thermal energy balance and PAR/net-radiation breakdown, sections 4 and
    the vertical-profile block of section 5, are NOT ported -- see module
    docstring).

    Parameters
    ----------
    spectral : SpectralConfig
        From :func:`scopeinpython.spectral.get_spectra_scope`.
    leaf_refl, leaf_tran : array_like, shape (nl, 2001) or (2001,)
        Leaf hemispherical reflectance/transmittance, 400-2400 nm (e.g.
        from ``toolsrtm.prospect_d``/``prospect_pro``, sliced to
        ``[:2001]``). A 1-D array is broadcast to every canopy layer
        (equivalent to mSCOPE with a single layer, ``mly$nly == 1``).
    rho_thermal, tau_thermal : float
        Leaf reflectance/transmittance in the thermal region (SCOPE
        default 0.01), extended across the 161 thermal bands.
    rsoil : array_like, shape (2001,)
        Soil reflectance, 400-2400 nm (e.g. from
        :func:`scopeinpython.soil.get_bsm`).
    canopy : CanopyStructure
    tts, tto, psi : float
        Solar zenith, viewing zenith, relative azimuth (degrees).
    ``Esun_``, ``Esky_`` : array_like, shape (nwl,)
        Top-of-canopy direct solar / diffuse sky irradiance, on the
        ``spectral.wlS`` grid. This port only supports the "precomputed
        irradiance" mode of ``get.calcTOCirr`` (i.e. ``atmo`` already
        containing ``Esun_``/``Esky_``, as in the default SCOPEinR example
        data ``SCOPEinR::Esun_``/``SCOPEinR::Esky_``); the MODTRAN-atmospheric-
        file branch is not ported.

    Returns
    -------
    RTMoResult
    """
    wl = spectral.wlS
    nwl = len(wl)
    n_thermal = len(spectral.wlT)  # 161

    nl = canopy.nlayers
    litab = canopy.litab
    lazitab = canopy.lazitab
    LAI = canopy.LAI
    lidf = canopy.lidf
    xl = canopy.xl
    dx = 1.0 / nl
    hot = canopy.hot

    # ---- leaf optics: broadcast to nl layers and extend into the thermal region ----
    leaf_refl = np.asarray(leaf_refl, dtype=float)
    leaf_tran = np.asarray(leaf_tran, dtype=float)
    if leaf_refl.ndim == 1:
        leaf_refl = np.broadcast_to(leaf_refl, (nl, leaf_refl.shape[0])).copy()
    if leaf_tran.ndim == 1:
        leaf_tran = np.broadcast_to(leaf_tran, (nl, leaf_tran.shape[0])).copy()

    rho = np.concatenate([leaf_refl, np.full((nl, n_thermal), rho_thermal)], axis=1)
    tau = np.concatenate([leaf_tran, np.full((nl, n_thermal), tau_thermal)], axis=1)

    rsoil_full = np.concatenate([np.asarray(rsoil, dtype=float), np.full(n_thermal, rho_thermal)])

    # ---- 1.1 geometric quantities ----
    deg2rad = np.pi / 180.0
    cos_tts = np.cos(tts * deg2rad)
    tan_tts = np.tan(tts * deg2rad)
    cos_tto = np.cos(tto * deg2rad)
    tan_tto = np.tan(tto * deg2rad)
    psi_ = abs(psi - 360 * round(psi / 360))
    dso = np.sqrt(tan_tts**2 + tan_tto**2 - 2 * tan_tts * tan_tto * np.cos(psi_ * deg2rad))

    # ---- 1.2 extinction & scattering geometric factors ----
    vg = get_volscatt_scope(tts, tto, psi_, litab)
    chi_s, chi_o, frho, ftau = vg["chi_s"], vg["chi_o"], vg["frho"], vg["ftau"]

    cos_ttli = np.cos(litab * deg2rad)

    ksli = chi_s / cos_tts
    koli = chi_o / cos_tto
    sobli = frho * np.pi / (cos_tts * cos_tto)
    sofli = ftau * np.pi / (cos_tts * cos_tto)
    bfli = cos_ttli**2

    k = float(np.sum(ksli * lidf))
    K = float(np.sum(koli * lidf))
    bf = float(np.sum(bfli * lidf))
    sob = float(np.sum(sobli * lidf))
    sof = float(np.sum(sofli * lidf))

    # ---- 1.3 diffuse/direct scattering weights ----
    sdb = 0.5 * (k + bf)
    sdf = 0.5 * (k - bf)
    ddb = 0.5 * (1 + bf)
    ddf = 0.5 * (1 - bf)
    dob = 0.5 * (K + bf)
    dof = 0.5 * (K - bf)

    # ---- 2.1 thin-layer reflectance/transmittance factors ----
    sigb = ddb * rho + ddf * tau
    sigf = ddf * rho + ddb * tau
    sb = sdb * rho + sdf * tau
    sf = sdf * rho + sdb * tau
    vb = dob * rho + dof * tau
    vf = dof * rho + dob * tau
    w = sob * rho + sof * tau
    a = 1 - sigf

    # ---- 3. flux calculation ----
    iLAI = LAI / nl
    tau_ss = 1 - k * iLAI  # scalar, identical for all layers/wavelengths (thin-layer approx)
    tau_dd = 1 - a * iLAI
    tau_sd = sf * iLAI
    rho_sd = sb * iLAI
    rho_dd = sigb * iLAI

    refl_layers = get_reflectances(tau_ss, tau_sd, tau_dd, rho_dd, rho_sd, rsoil_full, nl, nwl)
    R_sd, R_dd = refl_layers["R_sd"], refl_layers["R_dd"]
    Xss, Xsd, Xdd = refl_layers["Xss"], refl_layers["Xsd"], refl_layers["Xdd"]

    rdd = R_dd[0, :]
    rsd = R_sd[0, :]

    Esun_ = np.asarray(Esun_, dtype=float)
    Esky_ = np.asarray(Esky_, dtype=float)

    Eflux1 = get_fluxprofile(Esun_, 0 * Esky_, rsoil_full, Xss, Xsd, Xdd, R_sd, R_dd, nl, nwl, rs_thermal=0.06)
    Emins_, Eplus_ = Eflux1["Emin_"], Eflux1["Eplu_"]

    Eflux2 = get_fluxprofile(0 * Esun_, Esky_, rsoil_full, Xss, Xsd, Xdd, R_sd, R_dd, nl, nwl, rs_thermal=0.06)
    Emind_, Eplud_ = Eflux2["Emin_"], Eflux2["Eplu_"]

    Emin_ = Emins_ + Emind_
    Eplu_ = Eplus_ + Eplud_

    # ---- 1.5 gap probabilities Ps, Po, Pso ----
    Ps = np.exp(k * xl * LAI)
    Po = np.exp(K * xl * LAI)
    Ps[:nl] = Ps[:nl] * (1 - np.exp(-k * LAI * dx)) / (k * LAI * dx)
    Po[:nl] = Po[:nl] * (1 - np.exp(-K * LAI * dx)) / (K * LAI * dx)

    q = hot
    Pso = np.zeros(len(xl))
    for j in range(len(xl)):
        val, _ = quad(lambda y: get_pso(K, k, LAI, q, dso, y), xl[j] - dx, xl[j])
        Pso[j] = val / dx

    both = np.minimum(Po, Ps)
    Pso = np.where(Pso > Po, both, Pso)
    Pso = np.where(Pso > Ps, both, Pso)

    # ---- 3.3 outgoing fluxes (viewing direction) ----
    vb_rs_thermal, vf_rs_thermal, w_rs_thermal = vb, vf, w
    rsoil_thermal = rsoil_full

    vb_Po_Emin_ = vb_rs_thermal * Po[:nl, None] * Emind_[:nl, :]
    vf_Po_Eplud_ = vf_rs_thermal * Po[:nl, None] * Eplud_[:nl, :]
    piLocd_ = (vb_Po_Emin_.sum(axis=0) + vf_Po_Eplud_.sum(axis=0)) * iLAI

    Emin_Po = Emind_[nl - 1, :] * Po[nl - 1]  # NB: R uses index `nl` (1-based) == Python nl-1, see get.fluxprofile note
    piLosd_ = rsoil_thermal * Emin_Po

    vb_Po_Emins_ = vb_rs_thermal * Po[:nl, None] * Emins_[:nl, :]
    vf_Po_Eplus_ = vf_rs_thermal * Po[:nl, None] * Eplus_[:nl, :]
    w_Po_Esun = (w_rs_thermal * Pso[:nl, None]).sum(axis=0) * Esun_
    piLocu_ = (vb_Po_Emins_.sum(axis=0) + vf_Po_Eplus_.sum(axis=0) + w_Po_Esun) * iLAI

    Emins_Po_Esun_ = Emins_[nl - 1, :] * Po[nl - 1] + Esun_ * Pso[nl - 1]
    piLosu_ = rsoil_thermal * Emins_Po_Esun_

    piLod_ = piLocd_ + piLosd_
    piLou_ = piLocu_ + piLosu_
    piLoc_ = piLocu_ + piLocd_
    piLos_ = piLosu_ + piLosd_
    piLo_ = piLoc_ + piLos_
    Lo_ = piLo_ / np.pi

    # Esun_/Esky_ are genuinely near-zero at a handful of water-vapor/O2
    # absorption wavelengths -- rso/rdo/refl are expected to come back
    # NaN/Inf there (refl is masked back to a finite fallback via rso just
    # below); suppress the resulting divide warning rather than let NumPy
    # print it (with this file's own absolute path) on every call.
    with np.errstate(divide="ignore", invalid="ignore"):
        rso = piLou_ / Esun_
        rdo = piLod_ / Esky_
        refl = piLo_ / (Esky_ + Esun_)
    refl = np.where(Esky_ < 1e-4, rso, refl)
    idx = Esky_ < 2e-4 * np.max(Esky_)
    refl = np.where(idx, rso, refl)

    # ---- section 5: TOC outgoing spectrally-integrated quantities ----
    Eout_ = Eplu_[0, :]
    Eouto = 0.001 * sint(Eout_[spectral.IwlP], spectral.wlP)
    Eoutt = 0.001 * sint(Eout_[spectral.IwlT], spectral.wlT)
    Lot = 0.001 * sint(Lo_[spectral.IwlT], spectral.wlT)

    return RTMoResult(
        rdd=rdd, rsd=rsd, rdo=rdo, rso=rso, refl=refl, Lo_=Lo_,
        Eout_=Eout_, Eouto=Eouto, Eoutt=Eoutt, Lot=Lot,
        Esun_=Esun_, Esky_=Esky_,
        k=k, K=K, Ps=Ps, Po=Po, Pso=Pso,
        Emin_=Emin_, Eplu_=Eplu_, rho_dd=rho_dd, R_dd=R_dd, Xdd=Xdd, vb=vb, vf=vf, tau_dd=tau_dd,
        Emins_=Emins_, Emind_=Emind_, Eplus_=Eplus_, Eplud_=Eplud_,
        Xsd=Xsd, Xss=Xss, R_sd=R_sd,
    )


# ---------------------------------------------------------------------------
# Section 4 (partial): PAR / net-radiation absorption breakdown
# ---------------------------------------------------------------------------


@dataclass
class NetRadiationLite:
    Rnuc: np.ndarray  # net radiation, sunlit leaves, per layer (W/m2), (nl,)
    Rnhc: np.ndarray  # net radiation, shaded leaves, per layer (W/m2), (nl,)
    Rnus: float  # net radiation, sunlit soil (W/m2)
    Rnhs: float  # net radiation, shaded soil (W/m2)
    Pnu_Cab: np.ndarray  # net PAR absorbed by Cab, sunlit leaves, per layer (umol m-2 s-1), (nl,)
    Pnh_Cab: np.ndarray  # net PAR absorbed by Cab, shaded leaves, per layer (umol m-2 s-1), (nl,)


def net_radiation_lite(
    spectral: SpectralConfig,
    rtmo: RTMoResult,
    canopy: CanopyStructure,
    tts: float,
    lazitab: np.ndarray,
    leaf_refl: np.ndarray,
    leaf_tran: np.ndarray,
    rho_thermal: float,
    tau_thermal: float,
    rsoil: np.ndarray,
    kChlrel: np.ndarray,
) -> NetRadiationLite:
    """Direct, partial port of RTMo.R's section 4 ("lite" branch only,
    against the fixed R source -- see module docstring). Computes just
    the 6 quantities :mod:`scopeinpython.ebal` needs.

    Parameters
    ----------
    rtmo : RTMoResult
        From :func:`run_rtmo`, called with the same ``canopy``/``tts``/
        ``leaf_refl``/``leaf_tran``/``rho_thermal``/``tau_thermal``/
        ``rsoil`` as here (``tto``/``psi``/``Esky_`` aren't needed here).
    lazitab : array_like, shape (36,)
        Leaf azimuth classes, degrees (same grid as ``canopy.lidf``'s
        13 inclination classes pair with).
    kChlrel : array_like, shape (nl, 2001) or (2001,)
        Relative contribution of chlorophyll to leaf absorption, 400-2400nm
        (from a Fluspect leaf model's ``kChlrel`` output; a 1-D array
        broadcasts to every layer; pass zeros for a plain PROSPECT leaf
        model, matching R's ``data.leafopt$kChlrel`` for that case).
    """
    const = _constants()
    nl = canopy.nlayers
    litab = canopy.litab
    lidf = canopy.lidf
    lazitab = np.asarray(lazitab, dtype=float)

    wl = spectral.wlS
    n_thermal = len(spectral.wlT)
    Ipar = np.arange(len(spectral.wlPAR))
    wlP = spectral.wlP
    IwlP = spectral.IwlP

    leaf_refl = np.asarray(leaf_refl, dtype=float)
    leaf_tran = np.asarray(leaf_tran, dtype=float)
    if leaf_refl.ndim == 1:
        leaf_refl = np.broadcast_to(leaf_refl, (nl, leaf_refl.shape[0])).copy()
    if leaf_tran.ndim == 1:
        leaf_tran = np.broadcast_to(leaf_tran, (nl, leaf_tran.shape[0])).copy()
    rho = np.concatenate([leaf_refl, np.full((nl, n_thermal), rho_thermal)], axis=1)
    tau = np.concatenate([leaf_tran, np.full((nl, n_thermal), tau_thermal)], axis=1)
    epsc = 1 - rho - tau  # (nl, nwl)

    rsoil_full = np.concatenate([np.asarray(rsoil, dtype=float), np.full(n_thermal, rho_thermal)])
    epss = 1 - rsoil_full  # (nwl,)

    kChlrel = np.asarray(kChlrel, dtype=float)
    if kChlrel.ndim == 1:
        kChlrel = np.broadcast_to(kChlrel, (nl, kChlrel.shape[0])).copy()

    Esun_ = rtmo.Esun_
    Emin_ = rtmo.Emin_
    Eplu_ = rtmo.Eplu_

    Asun = np.empty(nl)
    Pnsun_Cab = np.empty(nl)
    for j in range(nl):
        Asun[j] = 0.001 * sint(Esun_ * epsc[j, :], wl)
        Pnsun_Cab[j] = 0.001 * sint(
            get_e2phot(wlP * 1e-9, kChlrel[j, IwlP] * Esun_[IwlP] * epsc[j, IwlP], const), wlP,
        )

    # 4.3 direct-beam contribution (fixed: per-layer, not a stray scalar)
    deg2rad = np.pi / 180.0
    cos_tts, sin_tts = np.cos(tts * deg2rad), np.sin(tts * deg2rad)
    cos_ttli = np.cos(litab * deg2rad)
    sin_ttli = np.sin(litab * deg2rad)
    cos_phils = np.cos(lazitab * deg2rad)
    cds = cos_ttli[:, None] * cos_tts + sin_ttli[:, None] * (sin_tts * cos_phils)[None, :]
    absfs = np.abs(cds / cos_tts)
    fs_ = float(np.sum(lidf * absfs.mean(axis=1)))

    Rndir = fs_ * Asun
    Pndir_Cab = fs_ * Pnsun_Cab

    # 4.4 diffuse contribution
    Rndif = np.empty(nl)
    Pndif_Cab = np.empty(nl)
    for j in range(nl):
        E_j = 0.5 * (Emin_[j, :] + Emin_[j + 1, :] + Eplu_[j, :] + Eplu_[j + 1, :])
        Rndif_row = E_j * epsc[j, :]
        Rndif[j] = 0.001 * sint(Rndif_row, wl)
        Pndif_Cab_row = 0.001 * get_e2phot(wlP * 1e-9, kChlrel[j, IwlP] * Rndif_row[IwlP], const)
        Pndif_Cab[j] = sint(Pndif_Cab_row, wlP)

    Rnhc = Rndif
    Pnhc_Cab = Pndif_Cab
    Rnuc = Rndir + Rndif
    Pnuc_Cab = Pndir_Cab + Pndif_Cab

    Rndirsoil = 0.001 * sint(Esun_ * epss, wl)
    Rndifsoil = 0.001 * sint(Emin_[nl, :] * epss, wl)
    Rnus = float(Rndifsoil + Rndirsoil)
    Rnhs = float(Rndifsoil)

    return NetRadiationLite(
        Rnuc=Rnuc, Rnhc=Rnhc, Rnus=Rnus, Rnhs=Rnhs,
        Pnu_Cab=1e6 * Pnuc_Cab, Pnh_Cab=1e6 * Pnhc_Cab,
    )
