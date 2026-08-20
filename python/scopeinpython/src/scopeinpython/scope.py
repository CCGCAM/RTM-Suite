"""End-to-end SCOPE simulation wrapper: leaf optics -> soil -> optical BRDF
-> thermal energy balance -> [fluorescence] -> [zeaxanthin], driven by a
single LUT input row.

Direct, partial port of ``SCOPEinR::get.SCOPE`` (``SCOPEinR/R/get.SCOPE.R``),
"SCOPE-lite" only (matches every other module in this port). Ties together
:func:`~scopeinpython.fluspect_mscope.fluspect_mscope`,
:func:`~scopeinpython.soil.get_bsm` (or the bundled reference soil spectra,
``SCOPEinR``'s own default), :func:`~scopeinpython.rtmo.run_rtmo`,
:func:`~scopeinpython.ebal.ebal`, :func:`~scopeinpython.rtmf.rtmf` and
:func:`~scopeinpython.rtmz.rtmz`.

**Not ported** (same scope as documented for their own modules, or newly
scoped out here):

- ``get.SCOPE.parallel`` -- R's parallel-backend variant (``foreach``/
  ``doParallel``) of the same per-row loop below. Not a separate function
  here: parallelize :func:`get_scope` yourself (``multiprocessing``,
  ``joblib``, a plain loop, ...) -- R's parallel backend has no 1:1 Python
  equivalent worth porting.
- ``options.calc_directional`` (full BRDF over many angles, ``get.brdf``),
  ``options.calc_spectrum_planck`` (``RTMt_planck.R``, the per-wavelength
  thermal RTM), ``options.mSCOPE`` with more than one profile layer
  (``run_rtmo`` only ever sees one leaf-optics spectrum per canopy layer --
  see its own docstring), ``options.simulation`` time-series mode,
  ``options.LIDF`` from an angle file (only LIDFa/LIDFb-derived LIDF is
  supported), and ``options.irradiance`` measurement-file /
  MODTRAN-atmosphere-file modes (only the bundled default spectrum, or a
  caller-supplied ``Esun_``/``Esky_``) -- none of these are exposed as
  parameters here at all, matching how e.g. :func:`scopeinpython.rtmt_sb.rtmt_sb`
  simply omits its own unported ``obsdir`` branch rather than accepting
  and silently ignoring the argument.
- The canopy-level "derived data products" section of ``get.SCOPE.R``
  beyond what's listed on :class:`ScopeResult` -- ``Pnsun_Car``/
  ``Rnsun_Cab``/``Rnsun_PAR``/``LST``/etc all need net-radiation or
  radiance breakdowns (``Rnuc_Car``, the directional-brightness-temperature
  ``Lote``, ...) that :func:`scopeinpython.rtmo.net_radiation_lite` and
  :func:`scopeinpython.rtmt_sb.rtmt_sb` don't compute yet -- see their own
  module docstrings for exactly what's missing and why.
- ``options.soil_heat_method`` 0/1 and ``options.calc_rss_rbs == 1``
  (recomputing ``rss``/``rbs`` from ``SMC``/``LAI`` via ``calc_rssrbs`` --
  not ported) -- :func:`get_scope` always uses the LUT's own ``rss``/
  ``rbs`` columns directly (R's ``calc.rss_rbs == 0`` default) and the
  simple ``G = 0.35*Rn`` soil-heat method (R's ``soil_heat_method == 2``
  default), matching :func:`scopeinpython.ebal.ebal`'s own scope.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
from toolsrtm import dladgen

from ._data import constants as _constants
from ._data import default_irradiance, soil_scope_spectra
from .biochemical import LeafBio
from .ebal import EbalCanopyParams, EbalMeteo, EbalResult, EbalSoilParams, aggregator_ebal, ebal
from .fluspect import _scope_fluspect_optipar
from .fluspect_mscope import MultiLayerLeafBio, fluspect_mscope
from .rtmf import RTMfResult, rtmf
from .rtmo import CanopyStructure, RTMoResult, net_radiation_lite, run_rtmo
from .rtmz import RTMzResult, rtmz
from .soil import SoilParams, WettingParams, get_bsm
from .spectral import SpectralConfig, get_spectra_scope
from .utils import get_ephoton, sint

__all__ = ["TDP_DEFAULT", "get_zo_and_d", "ScopeOptions", "ScopeResult", "get_scope"]


TDP_DEFAULT: dict = {
    "delHaV": 65330.0, "delSV": 485.0, "delHdV": 149250.0,
    "delHaJ": 43540.0, "delSJ": 495.0, "delHdJ": 152040.0,
    "delHaP": 53100.0, "delSP": 490.0, "delHdP": 150650.0,
    "delHaR": 46390.0, "delSR": 490.0, "delHdR": 150650.0,
    "delHaKc": 79430.0, "delHaKo": 36380.0, "delHaT": 37830.0,
    "Q10": 2.0, "s1": 0.3, "s2": 313.15, "s3": 0.2, "s4": 288.15, "s5": 1.3, "s6": 328.15,
}
"""Fixed temperature-response coefficients for
:func:`~scopeinpython.biochemical.get_biochemical` (``leafbio.TDP``).
Direct port of ``SCOPEinR::define_temp_response_biochem`` -- always these
exact constants, nothing here is derived from caller input."""


def get_zo_and_d(
    CR: float, CSSOIL: float, CD1: float, Psicor: float, LAI: float, hc: float, kappa: float,
) -> tuple[float, float]:
    """Roughness length for momentum (``zom``) and zero-plane displacement
    height (``d``), from vegetation height and LAI (Verhoef, McNaughton &
    Jacobs 1997). Direct port of ``SCOPEinR::get.zo_and_d`` (against the
    fixed R source: an undefined-variable bug in its degenerate-canopy
    branch, ``zo_and_d$d <- d`` with no ``d`` ever assigned on that branch
    -- errors, or silently picks up a stale ``d`` left over from a previous
    call in the same R session -- is fixed in ``SCOPEinR/R/zo_and_d.R`` to
    the evidently-intended ``d <- 0``).

    Returns
    -------
    (zom, d) : tuple[float, float]
    """
    sq = np.sqrt(CD1 * LAI / 2)
    G1 = max(3.3, (CSSOIL + CR * LAI / 2) ** (-0.5))
    if LAI > 1e-7 and hc > 1e-7:
        d = hc * (1 - (1 - np.exp(-sq)) / sq)
    else:
        d = 0.0
    zom = (hc - d) * np.exp(-kappa * G1 + Psicor)
    return float(zom), float(d)


@dataclass
class ScopeOptions:
    """Subset of R's ``options.SCOPE`` this port actually implements (see
    module docstring for what isn't exposed at all)."""

    calc_fluor: bool = True  # options.calc_fluorescence
    calc_xanthophyllabs: bool = True  # options.calc_xanthophyllabs
    apply_t_corr: bool = True  # options.applTcorr
    use_monin_obukhov: bool = True  # options.MoninObukhov
    use_bsm_soil: bool = False  # options.soilspectrum: False = bundled reference soil file (R's own shipped default), True = BSM
    k_maxit: int = 100
    maxEBer: float = 1.0


@dataclass
class ScopeResult:
    rtmo: RTMoResult
    ebal: EbalResult
    rtmf: RTMfResult | None  # None unless options.calc_fluor
    rtmz: RTMzResult | None  # None unless options.calc_xanthophyllabs
    nlayers: int
    LAIsunlit: float
    LAIshaded: float
    Pnsun_Cab: float  # net PAR absorbed by Cab, sunlit canopy total, umol m-2 s-1
    Pnsha_Cab: float  # ... shaded canopy total
    Pntot_Cab: float  # Pnsun_Cab + Pnsha_Cab
    Ja: float  # canopy-total electron transport rate, umol m-2 s-1
    PNPQ: float  # canopy-total non-photochemical-quenching energy (photon-flux form), umol m-2 s-1
    fqe: float | None  # canopy-level reabsorption-corrected apparent fluorescence quantum efficiency; None unless calc_fluor


def get_scope(
    lut: Mapping,
    options: ScopeOptions | None = None,
    spectral: SpectralConfig | None = None,
    Esun_: np.ndarray | None = None,
    Esky_: np.ndarray | None = None,
    rsoil: np.ndarray | None = None,
) -> ScopeResult:
    """Run one full SCOPE simulation for a single LUT row. Direct port of
    the per-row body of ``SCOPEinR::get.SCOPE`` (see module docstring for
    exact scope).

    Parameters
    ----------
    lut : Mapping
        One row of SCOPE's ``LUT_input.csv`` layout (a ``dict`` or
        ``pandas.Series`` with at least the columns used below -- see
        ``SCOPEinR/inst/input/LUT_input.csv`` for the full reference set).
    options : ScopeOptions
    spectral : SpectralConfig, optional
        Defaults to :func:`scopeinpython.spectral.get_spectra_scope`.
    Esun_, Esky_ : array_like, shape (2162,), optional
        Top-of-atmosphere direct/diffuse irradiance, on ``spectral.wlS``.
        Defaults to SCOPE's own bundled example spectrum (R's
        ``options.irradiance == 0`` path).
    rsoil : array_like, shape (2001,), optional
        Soil reflectance, 400-2400nm. If omitted, computed from
        ``options.use_bsm_soil`` and the LUT's soil columns.
    """
    const = _constants()
    if options is None:
        options = ScopeOptions()
    if spectral is None:
        spectral = get_spectra_scope()
    if Esun_ is None or Esky_ is None:
        irr = default_irradiance()
        if Esun_ is None:
            Esun_ = irr.Esun_
        if Esky_ is None:
            Esky_ = irr.Esky_

    def f(key: str) -> float:
        return float(lut[key])

    # ---- leaf optics + biochemistry ----
    N, Cab, Car, EWT, LMA = f("N"), f("Cab"), f("Car"), f("EWT"), f("LMA")
    Prot, CBC, Cx, Cs = f("Prot"), f("CBC"), f("Cx"), f("Cs")
    rho_thermal, tau_thermal, fqe_leaf = f("rho_thermal"), f("tau_thermal"), f("fqe")
    Anth = f("Anth")

    leafbio = LeafBio(
        Type=str(lut["Type"]), stressfactor=f("stressfactor"), Vcmax25=f("Vcmax25"),
        BallBerry0=f("BallBerry0"), BallBerrySlope=f("BallBerrySlope"), Rdparam=f("Rdparam"),
        Kn0=f("Kn0"), Knalpha=f("Knalpha"), Knbeta=f("Knbeta"),
        TDP=dict(TDP_DEFAULT) if options.apply_t_corr else {},
    )
    leaf_emis = 1 - rho_thermal - tau_thermal

    # ---- canopy structure ----
    LAI, hc = f("LAI"), f("hc")
    LIDFa, LIDFb = f("LIDFa"), f("LIDFb")
    leafwidth = f("leafwidth")
    hot = leafwidth / hc
    lidf = dladgen(LIDFa, LIDFb).lidf
    Rin = f("Rin")
    nlayers = max(2, int(np.ceil(10 * LAI)) + (60 if (Rin < 200 and options.use_monin_obukhov) else 0))
    canopy = CanopyStructure(LAI=LAI, lidf=lidf, hot=hot, nlayers=nlayers)
    nl = nlayers

    zom, d = get_zo_and_d(f("CR"), f("CSSOIL"), f("CD1"), f("Psicor"), LAI, hc, const["kappa"])
    ebal_canopy = EbalCanopyParams(Cd=f("Cd"), rwc=f("rwc"), z0m=zom, d=d, hc=hc, leafwidth=leafwidth, kV=f("kV"))

    # ---- meteo ----
    meteo = EbalMeteo(Ta=f("Ta"), ea=f("ea"), Ca=f("Ca"), p=f("p"), u=f("u"), z=f("z"))

    # ---- soil ----
    if rsoil is None:
        if options.use_bsm_soil:
            rsoil = get_bsm(
                SoilParams(BSMBrightness=f("BSMBrightness"), BSMlat=f("BSMlat"), BSMlon=f("BSMlon")),
                WettingParams(SMp=15.0, SMC=25.0, film=0.015),
            )
        else:
            spectrum_idx = int(f("spectrum")) - 1  # R's `spectrum + 1` 1-indexed column pick, 0-indexed here
            rsoil = soil_scope_spectra()[:, spectrum_idx]
    rsoil = np.asarray(rsoil, dtype=float)
    rs_thermal = f("rs_thermal")
    soil = EbalSoilParams(rbs=f("rbs"), rss=f("rss"), rs_thermal=rs_thermal)

    # ---- leaf optics (Fluspect-Cx, single mSCOPE profile layer) ----
    mly = MultiLayerLeafBio(
        nly=1, pLAI=np.array([LAI]), pCab=np.array([Cab]), pEWT=np.array([EWT]),
        pCar=np.array([Car]), pLMA=np.array([LMA]), pCs=np.array([Cs]), pN=np.array([N]),
    )
    leafopt = fluspect_mscope(mly, spectral, nl, Cx=0.0, fqe=fqe_leaf, Prot=Prot, CBC=CBC, Anth=Anth, step=5.0)

    tts, tto, psi = f("tts"), f("tto"), f("psi")

    rtmo = run_rtmo(
        spectral=spectral, leaf_refl=leafopt.refl, leaf_tran=leafopt.tran,
        rho_thermal=rho_thermal, tau_thermal=tau_thermal, rsoil=rsoil, canopy=canopy,
        tts=tts, tto=tto, psi=psi, Esun_=Esun_, Esky_=Esky_,
    )

    ebal_res = ebal(
        spectral=spectral, rtmo=rtmo, canopy=canopy, ebal_canopy=ebal_canopy, meteo=meteo, soil=soil,
        rho_thermal=rho_thermal, tau_thermal=tau_thermal, leaf_refl=leafopt.refl, leaf_tran=leafopt.tran,
        rsoil=rsoil, kChlrel=leafopt.kChlrel, leafbio=leafbio, leaf_emis=leaf_emis, tts=tts,
        lazitab=canopy.lazitab, use_monin_obukhov=options.use_monin_obukhov,
        k_maxit=options.k_maxit, maxEBer=options.maxEBer,
    )

    Ps = rtmo.Ps[:nl]
    Ph = 1 - Ps

    LAIsunlit = LAI * float(np.mean(Ps))
    LAIshaded = LAI - LAIsunlit

    # ebal() computes net_radiation_lite internally but doesn't expose it on
    # EbalResult -- recompute (cheap, no iteration) for the canopy-level
    # PAR-absorption aggregates below.
    net_rad = net_radiation_lite(
        spectral=spectral, rtmo=rtmo, canopy=canopy, tts=tts, lazitab=canopy.lazitab,
        leaf_refl=leafopt.refl, leaf_tran=leafopt.tran, rho_thermal=rho_thermal, tau_thermal=tau_thermal,
        rsoil=rsoil, kChlrel=leafopt.kChlrel,
    )

    Pnsun_Cab = LAI * float(np.sum(Ps * net_rad.Pnu_Cab) / nl)
    Pnsha_Cab = LAI * float(np.sum(Ph * net_rad.Pnh_Cab) / nl)
    Pntot_Cab = Pnsun_Cab + Pnsha_Cab

    Ja = aggregator_ebal(LAI, ebal_res.bcu.Ja, ebal_res.bch.Ja, Ps, nl)
    PNPQ = aggregator_ebal(
        LAI, net_rad.Pnu_Cab * ebal_res.bcu.Phi_N, net_rad.Pnh_Cab * ebal_res.bch.Phi_N, Ps, nl,
    )

    rtmf_res = None
    fqe_canopy = None
    if options.calc_fluor:
        etau_full = np.broadcast_to(ebal_res.bcu.eta[:, None, None], (nl, 13, 36)).copy()
        rtmf_res = rtmf(
            spectral=spectral, rtmo=rtmo, canopy=canopy, tts=tts, tto=tto, psi=psi, rsoil=rsoil,
            Mb=leafopt.Mb, Mf=leafopt.Mf, etau=etau_full, etah=ebal_res.bch.eta,
        )
        # R also computes `PoutFrc <- leafbio$fqe * aPAR_Cab_eta` and a
        # `sigmaF` diagnostic here, but neither feeds `data.canopy$fqe`
        # below (only `EoutFrc`/`Pntot_Cab` does) -- not reproduced.
        op = _scope_fluspect_optipar()
        iw_wlF = (spectral.wlF - 400).astype(int)  # wlF is always a subset of wlP's 400-2400 grid
        ep = const["A"] * get_ephoton(spectral.wlF * 1e-9, const)
        EoutFrc_ = 1e-3 * ep * op.phi[iw_wlF]
        EoutFrc = 1e-3 * sint(EoutFrc_, spectral.wlF)
        fqe_canopy = float(EoutFrc / Pntot_Cab) if Pntot_Cab != 0 else float("nan")

    rtmz_res = None
    if options.calc_xanthophyllabs:
        leafoptZ = fluspect_mscope(mly, spectral, nl, Cx=1.0, fqe=fqe_leaf, Prot=Prot, CBC=CBC, Anth=Anth, step=5.0)
        Knu_full = np.broadcast_to(ebal_res.bcu.Kn[:, None, None], (nl, 13, 36)).copy()
        rtmz_res = rtmz(
            spectral=spectral, rtmo=rtmo, canopy=canopy, tts=tts, tto=tto, psi=psi, rsoil=rsoil,
            refl=leafopt.refl, tran=leafopt.tran, reflZ=leafoptZ.refl, tranZ=leafoptZ.tran,
            Knu=Knu_full, Knh=ebal_res.bch.Kn,
        )

    return ScopeResult(
        rtmo=rtmo, ebal=ebal_res, rtmf=rtmf_res, rtmz=rtmz_res, nlayers=nl,
        LAIsunlit=LAIsunlit, LAIshaded=LAIshaded, Pnsun_Cab=Pnsun_Cab, Pnsha_Cab=Pnsha_Cab,
        Pntot_Cab=Pntot_Cab, Ja=Ja, PNPQ=PNPQ, fqe=fqe_canopy,
    )
