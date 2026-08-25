"""fourSAIL canopy radiative transfer model and its supporting leaf-angle
distribution / scattering-geometry helpers.

Direct, function-by-function port of:
  - ToolsRTM/R/volscatt.R
  - ToolsRTM/R/campbell.R
  - ToolsRTM/R/dladgen.R
  - ToolsRTM/R/dcum.R
  - ToolsRTM/R/Jfunc1.R, Jfunc2.R, Jfunc4.R
  - ToolsRTM/R/NonConservativeScatering.R
  - ToolsRTM/R/foursail.R  (single-layer canopy; all 5 leaf models ToolsRTM
    itself supports -- PROSPECT-D/-PRO, Liberty, Fluspect-B/-B-Cx -- are
    wired in via ``_leaf_optics``, matching R's own dispatch)

References
----------
Verhoef W & Bach H, 2007. Coupled soil-leaf-canopy and atmosphere radiative
transfer modeling ... Remote Sensing of Environment, 109:166-182.
Verhoef, Jia, Xiao & Su, 2007. Unified optical-thermal four-stream
radiative transfer theory for homogeneous vegetation canopies. IEEE TGRS
45:1808-1822.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from .fluspect import fluspect_b, fluspect_cx
from .leaf import prospect_d, prospect_pro
from .liberty import liberty

__all__ = [
    "volscatt",
    "campbell",
    "dladgen",
    "dcum",
    "jfunc1",
    "jfunc2",
    "jfunc4",
    "non_conservative_scattering",
    "conservative_scattering",
    "scattering",
    "LeafAngleDistribution",
    "FourSAILResult",
    "foursail_core",
    "foursail",
    "FourSAIL2Result",
    "foursail2_core",
    "foursail2",
]


# ---------------------------------------------------------------------------
# Leaf angle distribution
# ---------------------------------------------------------------------------


@dataclass
class LeafAngleDistribution:
    lidf: np.ndarray
    litab: np.ndarray


def dcum(a: float, b: float, t: float) -> float:
    """Cumulative leaf inclination distribution value at angle ``t`` (deg).

    Direct port of ``ToolsRTM::dcum``.
    """
    rd = np.pi / 180.0
    if a >= 1:
        f = 1 - np.cos(rd * t)
    else:
        eps = 1e-8
        delx = 1.0
        x = 2 * rd * t
        p = x
        y = x
        while delx >= eps:
            y = a * np.sin(x) + 0.5 * b * np.sin(2.0 * x)
            dx = 0.5 * (y - x + p)
            x = x + dx
            delx = abs(dx)
        f = (2.0 * y + p) / np.pi
    return f


def dladgen(a: float, b: float) -> LeafAngleDistribution:
    """Bimodal (Verhoef) leaf angle distribution function, from parameters
    a (average leaf slope) and b (bimodality). Direct port of
    ``ToolsRTM::dladgen`` (a.k.a. ``SCOPEinR::leafangles``).

    Constraint: ``abs(a) + abs(b) < 1``.
    """
    litab = np.array([5, 15, 25, 35, 45, 55, 65, 75, 81, 83, 85, 87, 89], dtype=float)
    freq = np.zeros(13)
    for i1 in range(1, 9):
        t = i1 * 10
        freq[i1 - 1] = dcum(a, b, t)
    for i2 in range(9, 13):
        t = 80 + (i2 - 8) * 2
        freq[i2 - 1] = dcum(a, b, t)
    freq[12] = 1
    for i in range(12, 0, -1):
        freq[i] = freq[i] - freq[i - 1]
    return LeafAngleDistribution(lidf=freq, litab=litab)


def campbell(ala: float) -> LeafAngleDistribution:
    """Ellipsoidal leaf angle distribution (Campbell, 1986), parametrised
    by the average leaf inclination angle ``ala`` (degrees).

    Direct port of ``ToolsRTM::campbell``.
    """
    tx1 = np.array([10, 20, 30, 40, 50, 60, 70, 80, 82, 84, 86, 88, 90], dtype=float)
    tx2 = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 82, 84, 86, 88], dtype=float)
    litab = (tx2 + tx1) / 2
    n = len(litab)
    tl1 = tx1 * (np.pi / 180)
    tl2 = tx2 * (np.pi / 180)
    excent = np.exp(-1.6184e-5 * ala**3 + 2.1145e-3 * ala**2 - 1.2390e-1 * ala + 3.2491)

    freq = np.zeros(n)
    for i in range(n):
        x1 = excent / np.sqrt(1 + excent**2 * np.tan(tl1[i]) ** 2)
        x2 = excent / np.sqrt(1 + excent**2 * np.tan(tl2[i]) ** 2)
        if excent == 1:
            freq[i] = abs(np.cos(tl1[i]) - np.cos(tl2[i]))
        else:
            alpha = excent / np.sqrt(abs(1 - excent**2))
            alpha2 = alpha**2
            x12 = x1**2
            x22 = x2**2
            if excent > 1:
                alpx1 = np.sqrt(alpha2 + x12)
                alpx2 = np.sqrt(alpha2 + x22)
                dum = x1 * alpx1 + alpha2 * np.log(x1 + alpx1)
                freq[i] = abs(dum - (x2 * alpx2 + alpha2 * np.log(x2 + alpx2)))
            else:
                almx1 = np.sqrt(alpha2 - x12)
                almx2 = np.sqrt(alpha2 - x22)
                dum = x1 * almx1 + alpha2 * np.arcsin(x1 / alpha)
                freq[i] = abs(dum - (x2 * almx2 + alpha2 * np.arcsin(x2 / alpha)))
    freq0 = freq / freq.sum()
    return LeafAngleDistribution(lidf=freq0, litab=litab)


# ---------------------------------------------------------------------------
# Volume scattering geometry
# ---------------------------------------------------------------------------


def volscatt(tts: float, tto: float, psi: float, ttl: float):
    """Volume scattering functions and interception coefficients for given
    solar zenith, viewing zenith, azimuth and (scalar) leaf inclination
    angle. Direct port of ``ToolsRTM::volscatt`` (scalar-``ttl`` variant,
    called once per leaf-angle class inside :func:`foursail_core`).

    Returns
    -------
    tuple(chi_s, chi_o, frho, ftau)
    """
    rd = np.pi / 180.0
    costs = np.cos(rd * tts)
    costo = np.cos(rd * tto)
    sints = np.sin(rd * tts)
    sinto = np.sin(rd * tto)
    cospsi = np.cos(rd * psi)

    costl = np.cos(rd * ttl)
    sintl = np.sin(rd * ttl)

    cs = costl * costs
    co = costl * costo
    ss = sintl * sints
    so = sintl * sinto

    cosbts = 5.0
    if abs(ss) > 1e-6:
        cosbts = -cs / ss

    cosbto = 5.0
    if abs(so) > 1e-6:
        cosbto = -co / so

    if abs(cosbts) < 1:
        bts = np.arccos(cosbts)
        ds = ss
    else:
        bts = np.pi
        ds = cs

    chi_s = 2 / np.pi * ((bts - np.pi * 0.5) * cs + np.sin(bts) * ss)

    if abs(cosbto) < 1:
        bto = np.arccos(cosbto)
        doo = so
    elif tto < 90:
        bto = np.pi
        doo = co
    else:
        bto = 0.0
        doo = -co
    chi_o = 2 / np.pi * ((bto - np.pi * 0.5) * co + np.sin(bto) * so)

    btran1 = abs(bts - bto)
    btran2 = np.pi - abs(bts + bto - np.pi)

    psir = rd * psi
    if psir <= btran1:
        bt1 = psir
        bt2 = btran1
        bt3 = btran2
    else:
        bt1 = btran1
        if psir <= btran2:
            bt2 = psir
            bt3 = btran2
        else:
            bt2 = btran2
            bt3 = psir

    t1 = 2 * cs * co + ss * so * cospsi
    t2 = 0.0
    if bt2 > 0:
        t2 = np.sin(bt2) * (2 * ds * doo + ss * so * np.cos(bt1) * np.cos(bt3))

    denom = 2 * np.pi * np.pi
    frho = ((np.pi - bt2) * t1 + t2) / denom
    ftau = (-bt2 * t1 + t2) / denom

    frho = max(frho, 0.0)
    ftau = max(ftau, 0.0)

    return chi_s, chi_o, frho, ftau


# ---------------------------------------------------------------------------
# J functions (avoid singularities in SAIL solution)
# ---------------------------------------------------------------------------


def jfunc1(k: float, l: np.ndarray, t: float) -> np.ndarray:
    """J1 function with avoidance of singularity problem. Port of
    ``ToolsRTM::Jfunc1`` (``k``, ``t`` scalars; ``l`` array over wavelength)."""
    l = np.asarray(l, dtype=float)
    del_ = (k - l) * t
    out = np.zeros_like(l)
    mask = np.abs(del_) > 1e-3
    out[mask] = (np.exp(-l[mask] * t) - np.exp(-k * t)) / (k - l[mask])
    out[~mask] = 0.5 * t * (np.exp(-k * t) + np.exp(-l[~mask] * t)) * (1 - del_[~mask] * del_[~mask] / 12)
    return out


def jfunc2(k: float, l, t: float):
    """J2 function. Port of ``ToolsRTM::Jfunc2``. ``k``/``l`` may be scalars
    or arrays (broadcastable); ``t`` a scalar."""
    return (1 - np.exp(-(k + l) * t)) / (k + l)


def jfunc4(m: np.ndarray, t: float) -> np.ndarray:
    """J4 function for treating (near) conservative scattering. Port of
    ``ToolsRTM::Jfunc4``."""
    m = np.asarray(m, dtype=float)
    del_ = m * t
    out = np.zeros_like(m)
    mask = del_ > 1e-3
    out[mask] = (1 - np.exp(-del_[mask])) / (m[mask] * (1 + np.exp(-del_[mask])))
    out[~mask] = 0.5 * t * (1.0 - del_[~mask] * del_[~mask] / 12.0)
    return out


def conservative_scattering(m, lai, att, sigb, ks, ko, sf, sb, vf, vb, tss, too):
    """Near/complete conservative scattering solution (m close to 0), used
    by foursail2/INFORM for wavelengths where non_conservative_scattering's
    general exponential-decay formulation is numerically unstable.

    Direct port of ``ToolsRTM::ConservativeScattering``. Same argument/return
    shape as :func:`non_conservative_scattering`.
    """
    j4 = jfunc4(m, lai)
    amsig = att - sigb
    apsig = att + sigb
    rtp = (1 - amsig * j4) / (1 + amsig * j4)
    rtm = (-1 + apsig * j4) / (1 + apsig * j4)
    rdd = 0.5 * (rtp + rtm)
    tdd = 0.5 * (rtp - rtm)

    dns = ks * ks - m * m
    dno = ko * ko - m * m
    cks = (sb * (ks - att) - sf * sigb) / dns
    cko = (vb * (ko - att) - vf * sigb) / dno
    dks = (-sf * (ks + att) - sb * sigb) / dns
    dko = (-vf * (ko + att) - vb * sigb) / dno
    ho = (sf * cko + sb * dko) / (ko + ks)

    rsd = cks * (1 - tss * tdd) - dks * rdd
    rdo = cko * (1 - too * tdd) - dko * rdd
    tsd = dks * (tss - tdd) - cks * tss * rdd
    tdo = dko * (too - tdd) - cko * too * rdd
    rsod = ho * (1 - tss * too) - cko * tsd * too - dko * rsd

    return {
        "tdd": tdd, "rdd": rdd, "tsd": tsd, "rsd": rsd,
        "tdo": tdo, "rdo": rdo, "rsod": rsod,
    }


def scattering(m, lai, att, sigb, ks, ko, sf, sb, vf, vb, tss, too):
    """Dispatch per-wavelength between :func:`conservative_scattering`
    (``m <= 0.01``) and :func:`non_conservative_scattering` (``m > 0.01``),
    matching ``foursail2.R``/``inform.R``'s ``f_ConS``/``f_Non_ConS`` index
    split -- computed for every wavelength with both formulations, then
    selected with ``np.where`` (simpler and equally correct vs. subsetting
    arrays and reassembling, since both branches are vectorised already).
    plain ``foursail()`` never needs this: it always uses the non-conservative
    formulation unconditionally, matching ``ToolsRTM::foursail.R`` itself.
    """
    m = np.asarray(m, dtype=float)
    nc = non_conservative_scattering(m, lai, att, sigb, ks, ko, sf, sb, vf, vb, tss, too)
    cs = conservative_scattering(m, lai, att, sigb, ks, ko, sf, sb, vf, vb, tss, too)
    is_cons = m <= 0.01
    return {k: np.where(is_cons, cs[k], nc[k]) for k in nc}


def non_conservative_scattering(m, lai, att, sigb, ks, ko, sf, sb, vf, vb, tss, too):
    """Non-conservative scattering solution of the SAIL 4-stream equations.

    Direct port of ``ToolsRTM::NonConservativeScattering``. All spectral
    quantities (``m``, ``att``, ``sigb``, ``sf``, ``sb``, ``vf``, ``vb``)
    are arrays over wavelength; ``lai``, ``ks``, ``ko``, ``tss``, ``too``
    are scalars.

    Returns
    -------
    dict with keys tdd, rdd, tsd, rsd, tdo, rdo, rsod (arrays over wavelength).
    """
    e1 = np.exp(-m * lai)
    e2 = e1 * e1
    rinf = (att - m) / sigb
    rinf2 = rinf * rinf
    re = rinf * e1
    denom = 1 - rinf2 * e2

    J1ks = jfunc1(ks, m, lai)
    J2ks = jfunc2(ks, m, lai)
    J1ko = jfunc1(ko, m, lai)
    J2ko = jfunc2(ko, m, lai)

    Ps = (sf + sb * rinf) * J1ks
    Qs = (sf * rinf + sb) * J2ks
    Pv = (vf + vb * rinf) * J1ko
    Qv = (vf * rinf + vb) * J2ko

    tdd = (1 - rinf2) * e1 / denom
    rdd = rinf * (1 - e2) / denom
    tsd = (Ps - re * Qs) / denom
    rsd = (Qs - re * Ps) / denom
    tdo = (Pv - re * Qv) / denom
    rdo = (Qv - re * Pv) / denom

    z = jfunc2(ks, ko, lai)

    g1 = (z - J1ks * too) / (ko + m)
    g2 = (z - J1ko * tss) / (ks + m)

    Tv1 = (vf * rinf + vb) * g1
    Tv2 = (vf + vb * rinf) * g2

    T1 = Tv1 * (sf + sb * rinf)
    T2 = Tv2 * (sf * rinf + sb)
    T3 = (rdo * Qs + tdo * Ps) * rinf

    rsod = (T1 + T2 - T3) / (1 - rinf2)

    return {
        "tdd": tdd,
        "rdd": rdd,
        "tsd": tsd,
        "rsd": rsd,
        "tdo": tdo,
        "rdo": rdo,
        "rsod": rsod,
    }


# ---------------------------------------------------------------------------
# fourSAIL
# ---------------------------------------------------------------------------


@dataclass
class FourSAILResult:
    rdot: np.ndarray  # hemispherical-directional reflectance factor (viewing dir)
    rsot: np.ndarray  # bi-directional reflectance factor
    rddt: np.ndarray  # bi-hemispherical reflectance factor
    rsdt: np.ndarray  # directional-hemispherical reflectance factor (solar)


def _foursail_scattering_core(
    rho: np.ndarray,
    tau: np.ndarray,
    rsoil: np.ndarray,
    lidf: np.ndarray,
    litab: np.ndarray,
    lai: float,
    hotspot: float,
    tts: float,
    tto: float,
    psi: float,
) -> FourSAILResult:
    """SAIL geometry + non-conservative-scattering + hotspot + soil
    interaction, given a precomputed leaf-angle distribution
    (``lidf``/``litab``). This is sections 1.2 onward of ``ToolsRTM::foursail``
    (everything after LIDF acquisition), factored out so
    :mod:`toolsrtm.inform`'s internal helpers -- which source their LIDF
    differently (a hardcoded partial table for ``foursail.inform``/
    ``foursail_t_s``/``foursail_t_o``, see that module) -- can reuse the
    exact same physics as :func:`foursail_core`.
    """
    rho = np.asarray(rho, dtype=float)
    tau = np.asarray(tau, dtype=float)
    rsoil = np.asarray(rsoil, dtype=float)

    rd = np.pi / 180.0
    cts = np.cos(rd * tts)
    cto = np.cos(rd * tto)
    ctscto = cts * cto
    ttans = np.tan(rd * tts)
    ttano = np.tan(rd * tto)
    cospsi = np.cos(rd * psi)
    dso = np.sqrt(ttans * ttans + ttano * ttano - 2 * ttans * ttano * cospsi)

    ks = ko = bf = sob = sof = 0.0
    na = len(litab)
    for i in range(na):
        ttl = litab[i]
        ctl = np.cos(rd * ttl)
        chi_s, chi_o, frho, ftau = volscatt(tts, tto, psi, ttl)

        ksli = chi_s / cts
        koli = chi_o / cto
        sobli = frho * np.pi / ctscto
        sofli = ftau * np.pi / ctscto
        bfli = ctl * ctl

        ks += ksli * lidf[i]
        ko += koli * lidf[i]
        bf += bfli * lidf[i]
        sob += sobli * lidf[i]
        sof += sofli * lidf[i]

    sdb = 0.5 * (ks + bf)
    sdf = 0.5 * (ks - bf)
    ddb = 0.5 * (1.0 + bf)
    ddf = 0.5 * (1.0 - bf)
    dob = 0.5 * (ko + bf)
    dof = 0.5 * (ko - bf)

    sigb = ddb * rho + ddf * tau
    sigf = ddf * rho + ddb * tau
    att = 1 - sigf
    m2 = (att + sigb) * (att - sigb)
    m2 = np.where(m2 <= 0, 0.0, m2)
    m = np.sqrt(m2)

    sb = sdb * rho + sdf * tau
    sf = sdf * rho + sdb * tau
    vb = dob * rho + dof * tau
    vf = dof * rho + dob * tau
    w = sob * rho + sof * tau

    if lai < 0:
        nwl = rho.shape[0]
        rddt = rsoil.copy()
        rsdt = rsoil.copy()
        rdot = rsoil.copy()
        rsot = rsoil.copy()
        return FourSAILResult(rdot=rdot, rsot=rsot, rddt=rddt, rsdt=rsdt)

    tss = np.exp(-ks * lai)
    too = np.exp(-ko * lai)

    scat = non_conservative_scattering(m, lai, att, sigb, ks, ko, sf, sb, vf, vb, tss, too)
    tdd, rdd = scat["tdd"], scat["rdd"]
    tsd, rsd = scat["tsd"], scat["rsd"]
    tdo, rdo = scat["tdo"], scat["rdo"]
    rsod = scat["rsod"]

    alf = 1e6
    if hotspot > 0:
        alf = (dso / hotspot) * 2 / (ks + ko)
    if alf > 200:
        alf = 200.0

    if alf == 0:
        tsstoo = tss
        sumint = (1 - tss) / (ks * lai)
    else:
        fhot = lai * np.sqrt(ko * ks)
        x1 = 0.0
        y1 = 0.0
        f1 = 1.0
        fint = (1.0 - np.exp(-alf)) * 0.05
        sumint = 0.0
        for i in range(1, 21):
            if i < 20:
                x2 = -np.log(1 - i * fint) / alf
            else:
                x2 = 1.0
            y2 = -(ko + ks) * lai * x2 + fhot * (1 - np.exp(-alf * x2)) / alf
            f2 = np.exp(y2)
            sumint += (f2 - f1) * (x2 - x1) / (y2 - y1)
            x1, y1, f1 = x2, y2, f2
        tsstoo = f1

    rsos = w * lai * sumint
    rso = rsos + rsod

    dn = 1 - rsoil * rdd
    rddt = rdd + tdd * rsoil * tdd / dn
    rsdt = rsd + (tsd + tss) * rsoil * tdd / dn
    rdot = rdo + tdd * rsoil * (tdo + too) / dn
    rsodt = rsod + ((tss + tsd) * tdo + (tsd + tss * rsoil * rdd) * too) * rsoil / dn
    rsost = rsos + tsstoo * rsoil
    rsot = rsost + rsodt

    return FourSAILResult(rdot=rdot, rsot=rsot, rddt=rddt, rsdt=rsdt)


def foursail_core(
    rho: np.ndarray,
    tau: np.ndarray,
    rsoil: np.ndarray,
    LIDFa: float,
    LIDFb: float,
    TypeLidf: Literal[1, 2],
    lai: float,
    hotspot: float,
    tts: float,
    tto: float,
    psi: float,
) -> FourSAILResult:
    """fourSAIL canopy bidirectional reflectance, given precomputed leaf
    reflectance/transmittance spectra ``rho``/``tau`` and soil reflectance
    ``rsoil`` (all same length, over wavelength).

    This is the model-agnostic core of ``ToolsRTM::foursail`` (the part
    after leaf-model dispatch): sections 1.2 onward of the R function.

    Parameters
    ----------
    rho, tau : array_like
        Leaf hemispherical reflectance / transmittance spectra.
    rsoil : array_like
        Soil reflectance spectrum, same length as rho/tau.
    LIDFa, LIDFb : float
        Leaf inclination distribution parameters (see :func:`dladgen`,
        :func:`campbell`).
    TypeLidf : {1, 2}
        1: use :func:`dladgen` (LIDFa, LIDFb shape parameters);
        2: use :func:`campbell` (LIDFa = average leaf angle).
    lai : float
        Leaf area index. If negative, canopy is bare soil (LAI=0 case).
    hotspot : float
        Hot-spot size parameter.
    tts, tto, psi : float
        Solar zenith, viewing zenith, relative azimuth (degrees).

    Returns
    -------
    FourSAILResult
        rdot, rsot, rddt, rsdt spectra.
    """
    if TypeLidf == 1:
        ld = dladgen(LIDFa, LIDFb)
    elif TypeLidf == 2:
        ld = campbell(LIDFa)
    else:
        raise ValueError("TypeLidf must be 1 (dladgen) or 2 (campbell)")
    return _foursail_scattering_core(rho, tau, rsoil, ld.lidf, ld.litab, lai, hotspot, tts, tto, psi)


_LeafModel = Literal["PROSPECT-PRO", "PROSPECT-D", "Liberty", "Fluspect-B", "Fluspect-B-Cx"]


def _leaf_optics(inputLUT: dict, leaf_model: _LeafModel) -> tuple[np.ndarray, np.ndarray, bool]:
    """Dispatch to any of the 5 ported leaf models, matching the pattern
    shared by ``ToolsRTM::foursail``, ``get.foursail2.leafopt`` (used by
    ``foursail2``), and ``inform``. Returns ``(refl, tran, force_2001)`` --
    ``force_2001`` is True for the two Fluspect leaf models, whose
    ``optipar`` table only spans 400-2400 nm (2001 pts), so callers must
    truncate ``rsoil``/output to that range regardless of their own
    ``spectrum_all`` flag, matching R's unconditional
    ``rsoil <- rsoil[1:2001]`` in those two branches of ``foursail.R``.

    ``Fluspect-B``'s ``Prot``/``CBC``/``Anth`` are passed through only if
    all three keys are present in ``inputLUT`` -- matching
    ``getFluspect.B``'s own column-presence auto-detection of the
    PROSPECT-PRO-style Kall formula (see :func:`toolsrtm.fluspect.fluspect_b`).
    """
    if leaf_model == "PROSPECT-PRO":
        lrt = prospect_pro(
            inputLUT["N"], inputLUT["Cab"], inputLUT["Car"], inputLUT["Anth"], inputLUT["Cbrown"],
            inputLUT["EWT"], inputLUT["LMA"], inputLUT["alpha"], inputLUT["Prot"], inputLUT["CBC"],
        )
        return lrt.refl, lrt.tran, False
    elif leaf_model == "PROSPECT-D":
        lrt = prospect_d(
            inputLUT["N"], inputLUT["Cab"], inputLUT["Car"], inputLUT["Anth"], inputLUT["Cbrown"],
            inputLUT["EWT"], inputLUT["LMA"], inputLUT["alpha"],
        )
        return lrt.refl, lrt.tran, False
    elif leaf_model == "Liberty":
        lrt = liberty(
            inputLUT["cell.d"], inputLUT["inter.c"], inputLUT["baseline.abs"], inputLUT["leaf.thick"],
            inputLUT["albino.abs"], inputLUT["Cab"], inputLUT["EWT"], inputLUT["lign.cell"], inputLUT["Nitrogen"],
        )
        return lrt.refl, lrt.tran, False
    elif leaf_model == "Fluspect-B":
        pro_kwargs = {}
        if all(k in inputLUT for k in ("Prot", "CBC", "Anth")):
            pro_kwargs = dict(Prot=inputLUT["Prot"], CBC=inputLUT["CBC"], Anth=inputLUT["Anth"])
        lrt = fluspect_b(inputLUT["Cab"], inputLUT["Car"], inputLUT["EWT"], inputLUT["LMA"],
                          inputLUT["Cs"], inputLUT["N"], inputLUT["fqe"], inputLUT["Cx"], **pro_kwargs)
        return lrt.refl, lrt.tran, True
    elif leaf_model == "Fluspect-B-Cx":
        lrt = fluspect_cx(inputLUT["Cab"], inputLUT["Car"], inputLUT["EWT"], inputLUT["LMA"],
                           inputLUT["Cs"], inputLUT["N"], inputLUT["fqe"], inputLUT["Cx"],
                           inputLUT["Prot"], inputLUT["CBC"], inputLUT["Anth"])
        return lrt.refl, lrt.tran, True
    else:
        raise ValueError(
            f"leaf_model={leaf_model!r} not supported; use one of 'PROSPECT-PRO', 'PROSPECT-D', "
            "'Liberty', 'Fluspect-B', 'Fluspect-B-Cx'"
        )


def foursail(
    inputLUT: dict,
    rsoil: np.ndarray,
    leaf_model: _LeafModel = "PROSPECT-PRO",
    spectrum_all: bool = True,
) -> FourSAILResult:
    """fourSAIL simulation for a single set of input parameters, dispatching
    to a leaf model exactly as ``ToolsRTM::foursail`` does -- all 5 leaf
    models ToolsRTM itself supports are ported and wired in here.

    Parameters
    ----------
    inputLUT : dict
        Scalar input parameters. Must contain the fourSAIL geometry/canopy
        keys (``LIDFa``, ``LIDFb``, ``TypeLidf``, ``LAI``, ``hspot``,
        ``tts``, ``tto``, ``psi``) plus the leaf-model parameters: for
        PROSPECT-D/-PRO, ``N``, ``Cab``, ``Car``, ``Anth``, ``Cbrown``,
        ``EWT``, ``LMA``, ``alpha`` (PRO also ``Prot``, ``CBC``); for
        Liberty, ``cell.d``, ``inter.c``, ``baseline.abs``, ``leaf.thick``,
        ``albino.abs``, ``Cab``, ``EWT``, ``lign.cell``, ``Nitrogen``; for
        Fluspect-B/-B-Cx, ``Cab``, ``Car``, ``EWT``, ``LMA``, ``Cs``, ``N``,
        ``fqe``, ``Cx`` (Cx also ``Prot``, ``CBC``, ``Anth``). See
        :func:`toolsrtm.liberty.liberty`/:func:`toolsrtm.fluspect.fluspect_b`/
        :func:`toolsrtm.fluspect.fluspect_cx` for parameter meanings.
    rsoil : array_like
        Soil reflectance spectrum, 400-2500 nm (2101 values) if
        ``spectrum_all`` else 400-2400 nm (2001 values). Always truncated
        to 2001 for the Fluspect leaf models regardless of ``spectrum_all``
        (see :func:`_leaf_optics`).
    leaf_model : {'PROSPECT-PRO', 'PROSPECT-D', 'Liberty', 'Fluspect-B', 'Fluspect-B-Cx'}
    spectrum_all : bool
        True: full 400-2500 nm PROSPECT range (2101 pts). False: truncate
        to 400-2400 nm (2001 pts), as used inside SCOPE-style pipelines.

    Returns
    -------
    FourSAILResult
    """
    rho, tau, force_2001 = _leaf_optics(inputLUT, leaf_model)
    rsoil = np.asarray(rsoil, dtype=float)
    if not spectrum_all or force_2001:
        rho = rho[:2001]
        tau = tau[:2001]
        rsoil = rsoil[:2001]

    return foursail_core(
        rho, tau, rsoil,
        LIDFa=inputLUT["LIDFa"], LIDFb=inputLUT["LIDFb"], TypeLidf=inputLUT["TypeLidf"],
        lai=inputLUT["LAI"], hotspot=inputLUT["hspot"],
        tts=inputLUT["tts"], tto=inputLUT["tto"], psi=inputLUT["psi"],
    )


# ---------------------------------------------------------------------------
# foursail2 (two-layer green/brown canopy)
# ---------------------------------------------------------------------------


@dataclass
class FourSAIL2Result:
    rdot: np.ndarray
    rsot: np.ndarray
    rddt: np.ndarray
    rsdt: np.ndarray
    alfast: np.ndarray  # canopy absorptance, direct solar incident flux
    alfadt: np.ndarray  # canopy absorptance, hemispherical diffuse incident flux


def foursail2_core(
    rho_green: np.ndarray, tau_green: np.ndarray,
    rho_brown: np.ndarray, tau_brown: np.ndarray,
    rsoil: np.ndarray,
    LIDFa: float, LIDFb: float, TypeLidf: Literal[1, 2],
    lai: float, hotspot: float, tts: float, tto: float, psi: float,
    fraction_brown: float, diss: float, Cv: float, Zeta: float,
) -> FourSAIL2Result:
    """Two-layer (green over brown/senescent) canopy bidirectional
    reflectance, given precomputed leaf optics for each layer.

    Direct port of ``ToolsRTM::foursail2``'s model-agnostic core (the part
    after leaf-model dispatch) -- non-Lambertian soil is not supported here
    either, matching the R version's own documented limitation.
    """
    rho_green = np.asarray(rho_green, dtype=float); tau_green = np.asarray(tau_green, dtype=float)
    rho_brown = np.asarray(rho_brown, dtype=float); tau_brown = np.asarray(tau_brown, dtype=float)
    rsoil = np.asarray(rsoil, dtype=float)
    rddsoil = rdosoil = rsdsoil = rsosoil = rsoil

    rd = np.pi / 180.0
    if lai <= 0:
        nwl = rsoil.shape[0]
        z = rsoil.copy()
        zeros = np.zeros(nwl)
        return FourSAIL2Result(rdot=z, rsot=z, rddt=z, rsdt=z, alfast=zeros, alfadt=zeros)

    if TypeLidf == 1:
        ld = dladgen(LIDFa, LIDFb)
    elif TypeLidf == 2:
        ld = campbell(LIDFa)
    else:
        raise ValueError("TypeLidf must be 1 (dladgen) or 2 (campbell)")
    lidf, litab = ld.lidf, ld.litab

    cts = np.cos(rd * tts); cto = np.cos(rd * tto); ctscto = cts * cto
    ttans = np.tan(rd * tts); ttano = np.tan(rd * tto); cospsi = np.cos(rd * psi)
    dso = np.sqrt(ttans * ttans + ttano * ttano - 2 * ttans * ttano * cospsi)

    # Crown/vegetation clumping (FLIM-style), same as INFORM
    Cs = Co = 1.0
    if Cv <= 1.0:
        Cs = 1.0 - (1.0 - Cv) ** (1.0 / cts)
        Co = 1.0 - (1.0 - Cv) ** (1.0 / cto)
    Overlap = 0.0
    if Zeta > 0.0:
        Overlap = min(Cs * (1.0 - Co), Co * (1.0 - Cs)) * np.exp(-dso / Zeta)
    Fcd = Cs * Co + Overlap
    Fcs = (1.0 - Cs) * Co - Overlap
    Fod = Cs * (1.0 - Co) - Overlap
    Fos = (1.0 - Cs) * (1.0 - Co) + Overlap
    Fcdc = 1.0 - (1.0 - Fcd) ** (0.5 / cts + 0.5 / cto)

    # Green/brown leaf-optics mixing within each layer
    leafgreen_r, leafgreen_t = rho_green, tau_green
    leafbrown_r, leafbrown_t = rho_brown, tau_brown
    fb = fraction_brown
    if fraction_brown == 0.0:
        fb = 0.5
        leafbrown_r, leafbrown_t = leafgreen_r, leafgreen_t
    if fraction_brown == 1.0:
        fb = 0.5
        leafgreen_r, leafgreen_t = leafbrown_r, leafbrown_t
    s = (1.0 - diss) * fb * (1.0 - fb)
    rho1 = ((1 - fb - s) * leafgreen_r + s * leafbrown_r) / (1 - fb)
    tau1 = ((1 - fb - s) * leafgreen_t + s * leafbrown_t) / (1 - fb)
    rho2 = (s * leafgreen_r + (fb - s) * leafbrown_r) / fb
    tau2 = (s * leafgreen_t + (fb - s) * leafbrown_t) / fb

    ks = ko = bf = sob = sof = 0.0
    for i in range(len(litab)):
        ttl = litab[i]
        ctl = np.cos(rd * ttl)
        chi_s, chi_o, frho, ftau = volscatt(tts, tto, psi, ttl)
        ks += (chi_s / cts) * lidf[i]
        ko += (chi_o / cto) * lidf[i]
        bf += (ctl * ctl) * lidf[i]
        sob += (frho * np.pi / ctscto) * lidf[i]
        sof += (ftau * np.pi / ctscto) * lidf[i]

    sdb = 0.5 * (ks + bf); sdf = 0.5 * (ks - bf)
    dob = 0.5 * (ko + bf); dof = 0.5 * (ko - bf)
    ddb = 0.5 * (1.0 + bf); ddf = 0.5 * (1.0 - bf)

    lai1 = (1 - fb) * lai
    lai2 = fb * lai

    tss_full = np.exp(-ks * lai)
    ck = np.exp(-ks * lai1)
    alf = 1e6
    if hotspot > 0.0:
        alf = (dso / hotspot) * 2.0 / (ks + ko)
    alf = min(alf, 200.0)

    if alf == 0.0:
        tsstoo = tss_full
        s1 = (1 - ck) / (ks * lai)
        s2 = (ck - tss_full) / (ks * lai)
    else:
        fhot = lai * np.sqrt(ko * ks)
        ca = np.exp(alf * (fb - 1.0))
        x1 = y1 = 0.0; f1 = 1.0
        fint = (1.0 - ca) * 0.05
        s1 = 0.0
        for istep in range(1, 21):
            x2 = -np.log(1.0 - istep * fint) / alf if istep < 20 else 1.0 - fb
            y2 = -(ko + ks) * lai * x2 + fhot * (1.0 - np.exp(-alf * x2)) / alf
            f2 = np.exp(y2)
            s1 += (f2 - f1) * (x2 - x1) / (y2 - y1)
            x1, y1, f1 = x2, y2, f2
        fint = (ca - np.exp(-alf)) * 0.05
        s2 = 0.0
        for istep in range(1, 21):
            x2 = -np.log(ca - istep * fint) / alf if istep < 20 else 1.0
            y2 = -(ko + ks) * lai * x2 + fhot * (1.0 - np.exp(-alf * x2)) / alf
            f2 = np.exp(y2)
            s2 += (f2 - f1) * (x2 - x1) / (y2 - y1)
            x1, y1, f1 = x2, y2, f2
        tsstoo = f1

    # --- Bottom layer (index 2: brown-heavier) ---
    tss = np.exp(-ks * lai2); too = np.exp(-ko * lai2)
    sb = sdb * rho2 + sdf * tau2; sf = sdf * rho2 + sdb * tau2
    vb = dob * rho2 + dof * tau2; vf = dof * rho2 + dob * tau2
    w2 = sob * rho2 + sof * tau2
    sigb = ddb * rho2 + ddf * tau2; sigf = ddf * rho2 + ddb * tau2
    att = 1.0 - sigf
    m2 = np.where((att + sigb) * (att - sigb) < 0, 0.0, (att + sigb) * (att - sigb))
    m = np.sqrt(m2)
    scat = scattering(m, lai2, att, sigb, ks, ko, sf, sb, vf, vb, tss, too)
    tdd, rdd, tsd, rsd, tdo, rdo, rsod = (scat[k] for k in ("tdd", "rdd", "tsd", "rsd", "tdo", "rdo", "rsod"))

    rddb, rsdb, rdob, rsodb = rdd, rsd, rdo, rsod
    tddb, tsdb, tdob, toob, tssb = tdd, tsd, tdo, too, tss

    # --- Top layer (index 1: green-heavier) ---
    tss = np.exp(-ks * lai1); too = np.exp(-ko * lai1)
    sb = sdb * rho1 + sdf * tau1; sf = sdf * rho1 + sdb * tau1
    vb = dob * rho1 + dof * tau1; vf = dof * rho1 + dob * tau1
    w1 = sob * rho1 + sof * tau1
    sigb = ddb * rho1 + ddf * tau1; sigf = ddf * rho1 + ddb * tau1
    att = 1.0 - sigf
    m2 = np.where((att + sigb) * (att - sigb) < 0, 0.0, (att + sigb) * (att - sigb))
    m = np.sqrt(m2)
    scat = scattering(m, lai1, att, sigb, ks, ko, sf, sb, vf, vb, tss, too)
    tdd, rdd, tsd, rsd, tdo, rdo, rsod = (scat[k] for k in ("tdd", "rdd", "tsd", "rsd", "tdo", "rdo", "rsod"))

    # Combine layers (adding method)
    rn = 1.0 - rdd * rddb
    tup = (tss * rsdb + tsd * rddb) / rn
    tdn = (tsd + tss * rsdb * rdd) / rn
    rsdt = rsd + tup * tdd
    rdot = rdo + tdd * (rddb * tdo + rdob * too) / rn
    rsodt = rsod + (tss * rsodb + tdn * rdob) * too + tup * tdo
    rsost = (w1 * s1 + w2 * s2) * lai
    rsot = rsost + rsodt

    rddt_t = rdd + tdd * rddb * tdd / rn
    rddt_b = rddb + tddb * rdd * tddb / rn
    tsst = tss * tssb; toot = too * toob
    tsdt = tss * tsdb + tdn * tddb
    tdot = tdob * too + tddb * (tdo + rdd * rdob * too) / rn
    tddt = tdd * tddb / rn

    # Apply clumping
    rddcb = Cv * rddt_b; rddct = Cv * rddt_t
    tddc = 1 - Cv + Cv * tddt
    rsdc = Cs * rsdt; tsdc = Cs * tsdt
    rdoc = Co * rdot; tdoc = Co * tdot
    tssc = 1 - Cs + Cs * tsst; tooc = 1 - Co + Co * toot

    rsoc = Fcdc * rsot
    tssooc = Fcd * tsstoo + Fcs * toot + Fod * tsst + Fos
    alfas = 1.0 - tssc - tsdc - rsdc
    alfad = 1.0 - tddc - rddct

    rn = 1 - rddcb * rddsoil
    tup = (tssc * rsdsoil + tsdc * rddsoil) / rn
    tdn = (tsdc + tssc * rsdsoil * rddcb) / rn
    rddt = rddct + tddc * rddsoil * tddc / rn
    rsdt = rsdc + tup * tddc
    rdot = rdoc + tddc * (rddsoil * tdoc + rdosoil * tooc) / rn
    rsot = rsoc + tssooc * rsosoil + tdn * rdosoil * tooc + tup * tdoc

    alfast = alfas + tup * alfad
    alfadt = alfad * (1.0 + tddc * rddsoil / rn)

    return FourSAIL2Result(rdot=rdot, rsot=rsot, rddt=rddt, rsdt=rsdt, alfast=alfast, alfadt=alfadt)


# Illustrative default senescent (brown) leaf spectrum -- matches
# ToolsRTM::foursail2's own LUT_GB row-2 default (used when the caller
# doesn't supply a real field-observed brown-leaf spectrum), extended with
# the Fluspect/Liberty columns that same default LUT_GB row also carries
# (foursail2.R lines ~105-112: "covering the columns needed by any of the
# 5 leaf models").
_BROWN_LEAF_DEFAULT_D = dict(N=2.0, Cab=5.0, Car=5.0, Anth=0.0, Cbrown=0.0, EWT=0.005, LMA=0.008, alpha=40.0)
_BROWN_LEAF_DEFAULT_PRO = dict(**_BROWN_LEAF_DEFAULT_D, Prot=0.0, CBC=0.0)
_BROWN_LEAF_DEFAULT_EXTRA = dict(
    Cs=0.1, fqe=0.01, Cx=0.0,
    **{"cell.d": 40.0, "inter.c": 0.045, "baseline.abs": 0.0012, "leaf.thick": 1.6,
       "albino.abs": 2.0, "lign.cell": 6.0, "Nitrogen": 0.5},
)
_BROWN_LEAF_DEFAULT_LIBERTY = {**_BROWN_LEAF_DEFAULT_D, **_BROWN_LEAF_DEFAULT_EXTRA}
_BROWN_LEAF_DEFAULT_FLUSPECT = {**_BROWN_LEAF_DEFAULT_PRO, **_BROWN_LEAF_DEFAULT_EXTRA}


def foursail2(
    inputLUT: dict,
    rsoil: np.ndarray,
    leaf_model: _LeafModel = "PROSPECT-PRO",
    spectrum_all: bool = True,
) -> FourSAIL2Result:
    """foursail2 (two-layer green/brown canopy) for a single set of input
    parameters, dispatching to a leaf model exactly as
    ``ToolsRTM::foursail2``'s ``get.foursail2.leafopt`` does -- all 5 leaf
    models are wired in (see :func:`foursail` for the full parameter list
    per leaf model). Green-vegetation leaf optics come from ``inputLUT``
    itself (matching ``ToolsRTM::foursail2``'s default ``FieldObserv=NULL``
    behaviour); brown-vegetation leaf optics use the illustrative default
    senescent spectrum (``_BROWN_LEAF_DEFAULT_*``) the R function also falls
    back to when no field-observed ``LUT_GB`` is supplied.

    Extra ``inputLUT`` keys beyond plain :func:`foursail`: ``fraction_brown``
    (0-1), ``diss`` (layer dissociation factor), ``Cv`` (vertical crown
    cover fraction), ``Zeta`` (tree shape factor, crown diameter / height).
    """
    brown_defaults = {
        "PROSPECT-PRO": _BROWN_LEAF_DEFAULT_PRO,
        "PROSPECT-D": _BROWN_LEAF_DEFAULT_D,
        "Liberty": _BROWN_LEAF_DEFAULT_LIBERTY,
        "Fluspect-B": _BROWN_LEAF_DEFAULT_FLUSPECT,
        "Fluspect-B-Cx": _BROWN_LEAF_DEFAULT_FLUSPECT,
    }.get(leaf_model)
    if brown_defaults is None:
        raise ValueError(
            f"leaf_model={leaf_model!r} not supported; use one of 'PROSPECT-PRO', 'PROSPECT-D', "
            "'Liberty', 'Fluspect-B', 'Fluspect-B-Cx'"
        )

    rho_g, tau_g, force_2001 = _leaf_optics(inputLUT, leaf_model)
    rho_b, tau_b, _ = _leaf_optics(brown_defaults, leaf_model)
    rsoil = np.asarray(rsoil, dtype=float)
    if not spectrum_all or force_2001:
        rho_g, tau_g, rho_b, tau_b = rho_g[:2001], tau_g[:2001], rho_b[:2001], tau_b[:2001]
        rsoil = rsoil[:2001]

    return foursail2_core(
        rho_g, tau_g, rho_b, tau_b, rsoil,
        LIDFa=inputLUT["LIDFa"], LIDFb=inputLUT["LIDFb"], TypeLidf=inputLUT["TypeLidf"],
        lai=inputLUT["LAI"], hotspot=inputLUT["hspot"],
        tts=inputLUT["tts"], tto=inputLUT["tto"], psi=inputLUT["psi"],
        fraction_brown=inputLUT["fraction_brown"], diss=inputLUT["diss"],
        Cv=inputLUT["Cv"], Zeta=inputLUT["Zeta"],
    )
