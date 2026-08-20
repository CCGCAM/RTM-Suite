"""SCOPE energy-balance closure loop: iterates leaf (sunlit/shaded) and
soil (sunlit/shaded) temperature until sensible+latent heat flux matches
net radiation, coupling radiative transfer (:func:`~scopeinpython.rtmo.net_radiation_lite`,
:func:`~scopeinpython.rtmt_sb.rtmt_sb`), aerodynamics
(:func:`~scopeinpython.thermal.get_resistances`,
:func:`~scopeinpython.thermal.monin_obukhov`), photosynthesis/fluorescence
(:func:`~scopeinpython.biochemical.get_biochemical`) and heat-flux
partitioning (:func:`~scopeinpython.thermal.get_heatfluxes`).

Direct, **"SCOPE-lite"-only** port of ``SCOPEinR::get.ebal``
(``SCOPEinR/R/ebal.R``). Matches every reference case built during this
whole port:

- Only the empirical/default fluorescence-model branch (``get.biochemical``,
  not the alternative ``get.biochemical.MD12`` -- not ported, only used
  when ``options.Fluorescence_model$Value == 1``).
- Only the simple ground-heat-flux method (R's ``SoilHeatMethod == 2``,
  ``G = 0.35 * Rn`` -- always the case for a single-timestep/non-timeseries
  run, i.e. ``options.simulation$Value != 1``); the two time-series-history
  soil-inertia methods (0/1, needing a rolling ``Tsold`` state across
  timesteps) are not ported.
- Only ``meanleaf.v2``'s ``'layers'`` aggregation mode (:func:`aggregator_ebal`)
  -- the "lite" pipeline's ``Tcu``/``Rnuc``/etc are plain length-``nl``
  vectors throughout (one value per canopy layer, no per-leaf-angle-class
  detail), so the full ``(13,36,nl)`` ``'angles'``/``'angles_and_layers'``
  modes never apply.

See :mod:`scopeinpython.rtmo` (``net_radiation_lite``) and
:mod:`scopeinpython.rtmt_sb` for two bugs found and fixed in R's own
dependencies of this loop before this port was written.
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass

import numpy as np

from ._data import constants as _constants
from .biochemical import BiochemResult, LeafBio, MeteoLeaf, get_biochemical
from .rtmo import CanopyStructure, NetRadiationLite, RTMoResult, net_radiation_lite
from .rtmt_sb import rtmt_sb
from .spectral import SpectralConfig
from .thermal import ResistanceParams, get_heatfluxes, get_resistances, monin_obukhov

__all__ = ["EbalMeteo", "EbalCanopyParams", "EbalSoilParams", "EbalResult", "aggregator_ebal", "ebal"]


def _biochemical_per_layer(
    leafbio: LeafBio, Temp: np.ndarray, eb: np.ndarray, Cs: np.ndarray, Q: np.ndarray,
    Oa: float, p: float, fV: np.ndarray, temp_correction: bool,
) -> BiochemResult:
    """Call :func:`get_biochemical` once per canopy layer (scalar
    micro-environment each time) and stack the results into per-layer
    arrays. ``get_biochemical`` itself is only proven correct for scalar
    inputs (its internal Brent Ci-solver indexes ``Cs``/``RH`` per-element
    but not the leaf-level photosynthesis parameters derived from
    ``Vcmax``/``Rd``/etc, which stay full-length arrays if given as such
    -- discovered while building this function, not chased further since
    per-layer scalar calls are already exactly how the rest of this port
    exercises ``get_biochemical``, and matches ``fV`` genuinely varying by
    layer as well, same as R's ``fV`` vertical Vcmax profile)."""
    nl = len(Temp)
    per_layer = [
        get_biochemical(
            leafbio, MeteoLeaf(Q=float(Q[i]), Cs=float(Cs[i]), Temp=float(Temp[i]), eb=float(eb[i]), Oa=Oa, p=p),
            temp_correction=temp_correction, fV=float(fV[i]),
        )
        for i in range(nl)
    ]
    kwargs = {}
    for f in dataclasses.fields(BiochemResult):
        values = [getattr(r, f.name) for r in per_layer]
        if values[0] is None:
            kwargs[f.name] = None
        else:
            kwargs[f.name] = np.array([float(np.atleast_1d(v).flatten()[0]) for v in values])
    return BiochemResult(**kwargs)


def aggregator_ebal(LAI: float, sunlit_flux: np.ndarray, shaded_flux: np.ndarray, Fs: np.ndarray, nl: int) -> float:
    """LAI-scaled canopy-integrated total flux, combining sunlit (weighted
    by ``Fs``) and shaded (weighted by ``1 - Fs``) leaf-scale contributions.
    Direct port of ``SCOPEinR::get.aggregator.ebal`` + ``meanleaf.v2``,
    ``'layers'`` mode only (see module docstring)."""
    Fs = np.asarray(Fs, dtype=float)
    sunlit_mean = np.sum(Fs * np.asarray(sunlit_flux, dtype=float)) / nl
    shaded_mean = np.sum((1 - Fs) * np.asarray(shaded_flux, dtype=float)) / nl
    return float(LAI * (sunlit_mean + shaded_mean))


@dataclass
class EbalMeteo:
    Ta: float  # air temperature, deg C
    ea: float  # air vapour pressure, hPa
    Ca: float  # ambient CO2, umol/mol
    p: float  # air pressure, hPa
    u: float  # wind speed at z, m/s
    z: float  # measurement height, m


@dataclass
class EbalCanopyParams:
    """Aerodynamic/structural canopy properties beyond
    :class:`~scopeinpython.rtmo.CanopyStructure` needed by the resistance
    scheme and the Vcmax vertical profile."""

    Cd: float  # leaf drag coefficient
    rwc: float  # within-canopy aerodynamic resistance, s/m
    z0m: float  # roughness length for momentum, m
    d: float  # zero-plane displacement height, m
    hc: float  # vegetation height, m
    leafwidth: float
    kV: float  # Vcmax25 vertical decay exponent


@dataclass
class EbalSoilParams:
    rbs: float  # soil boundary-layer resistance, s/m
    rss: float  # soil surface resistance for vapour transport, s/m
    rs_thermal: float  # soil thermal-region reflectance


@dataclass
class EbalResult:
    counter: int
    Tcu: np.ndarray
    Tch: np.ndarray
    Tsu: float
    Tsh: float
    bcu: BiochemResult
    bch: BiochemResult
    canopyemis: float
    Rnctot: float
    lEctot: float
    Hctot: float
    Actot: float
    Tcave: float
    Rnstot: float
    lEstot: float
    Hstot: float
    Gtot: float
    Tsave: float
    Rntot: float
    lEtot: float
    Htot: float
    maxEBercu: float
    maxEBerch: float
    maxEBers: float


def ebal(
    spectral: SpectralConfig,
    rtmo: RTMoResult,
    canopy: CanopyStructure,
    ebal_canopy: EbalCanopyParams,
    meteo: EbalMeteo,
    soil: EbalSoilParams,
    rho_thermal: float,
    tau_thermal: float,
    leaf_refl: np.ndarray,
    leaf_tran: np.ndarray,
    rsoil: np.ndarray,
    kChlrel: np.ndarray,
    leafbio: LeafBio,
    leaf_emis: float,
    tts: float,
    lazitab: np.ndarray,
    use_monin_obukhov: bool = True,
    k_maxit: int = 100,
    maxEBer: float = 1.0,
) -> EbalResult:
    """Direct port of ``SCOPEinR::get.ebal`` ("SCOPE-lite" only, see
    module docstring).

    Parameters
    ----------
    rtmo : RTMoResult
        From :func:`scopeinpython.rtmo.run_rtmo`, same canopy/geometry as
        elsewhere here.
    leaf_refl, leaf_tran, rsoil : array_like
        Same leaf/soil optics passed to ``run_rtmo`` (400-2400nm, 2001
        points; ``leaf_refl``/``leaf_tran`` may be ``(nl, 2001)`` or
        ``(2001,)``).
    kChlrel : array_like, shape (nl, 2001) or (2001,)
        See :func:`scopeinpython.rtmo.net_radiation_lite`.
    leafbio : LeafBio
        ``TDP`` must be populated for temperature-corrected biochemistry
        (matches how every biochemistry reference case in this port is
        built).
    leaf_emis : float
        Leaf thermal-IR emissivity (``data.leafbio$emis`` in R).
    k_maxit : int
        Maximum number of energy-balance iterations (``k.maxit`` in R).
    maxEBer : float
        Convergence threshold, maximum acceptable energy-balance error
        for any component (W/m2).
    """
    const = _constants()
    rhoa, cp, sigmaSB = const["rhoa"], const["cp"], const["sigmaSB"]
    MH2O, Mair = const["MH2O"], const["Mair"]

    nl = canopy.nlayers
    LAI = canopy.LAI
    Ps = rtmo.Ps

    Ta, ea, Ca, p = meteo.Ta, meteo.ea, meteo.Ca, meteo.p

    # Static optical-RTM contribution to net radiation (computed once;
    # RTMt.sb's thermal contribution is added to this fresh every
    # iteration below, matching R's `data.rad$Rnuc + data.rad$Rnuct` etc).
    net_rad: NetRadiationLite = net_radiation_lite(
        spectral=spectral, rtmo=rtmo, canopy=canopy, tts=tts, lazitab=lazitab,
        leaf_refl=leaf_refl, leaf_tran=leaf_tran, rho_thermal=rho_thermal,
        tau_thermal=tau_thermal, rsoil=rsoil, kChlrel=kChlrel,
    )

    e_to_q = MH2O / Mair / p
    Fc = Ps[:nl]  # sunlit leaf-area fraction per layer, used to weight canopy aggregation
    Fs = np.array([1 - Ps[nl], Ps[nl]])  # [shaded, sunlit] weight for soil (2,)
    Oa = 209.0  # O2 concentration, mmol/mol (SCOPE default)

    fV = np.exp(ebal_canopy.kV * canopy.xl[:nl])  # vertical Vcmax25 profile

    ech = np.full(nl, ea)
    Cch = np.full(nl, Ca)
    ecu = np.full(nl, ea)
    Ccu = np.full(nl, Ca)

    Ts = np.array([Ta + 3.0, Ta + 3.0])  # [shaded, sunlit] soil temperature
    Tch = np.full(nl, Ta + 0.1)
    Tcu = np.full(nl, Ta + 0.3)
    L = -1e6

    counter = 0
    Wc = 1.0
    CONT = True
    maxEBercu = maxEBerch = maxEBers = np.nan

    resistance_base = dict(
        rbs=soil.rbs, Cd=ebal_canopy.Cd, LAI=LAI, rwc=ebal_canopy.rwc,
        z0m=ebal_canopy.z0m, d=ebal_canopy.d, hc=ebal_canopy.hc, leafwidth=ebal_canopy.leafwidth,
        z=meteo.z, u=meteo.u,
    )

    lEch = Hch = lEcu = Hcu = lEs = Hs = G = None
    data_bch = data_bcu = None

    while CONT:
        rtmt = rtmt_sb(
            rtmo=rtmo, nl=nl, LAI=LAI, rho_thermal=rho_thermal, tau_thermal=tau_thermal,
            rs_thermal=soil.rs_thermal, Tcu=Tcu, Tch=Tch, Tsu=Ts[1], Tsh=Ts[0],
        )

        Rnhc = net_rad.Rnhc + rtmt.Rnhct
        Rnuc = net_rad.Rnuc + rtmt.Rnuct
        Rnhs = net_rad.Rnhs + rtmt.Rnhst
        Rnus = net_rad.Rnus + rtmt.Rnust
        Rns = np.array([Rnhs, Rnus])

        data_bch = _biochemical_per_layer(leafbio, Tch, ech, Cch, net_rad.Pnh_Cab, Oa, p, fV, temp_correction=True)
        data_bcu = _biochemical_per_layer(leafbio, Tcu, ecu, Ccu, net_rad.Pnu_Cab, Oa, p, fV, temp_correction=True)

        resist_out = get_resistances(ResistanceParams(**resistance_base, L=L))
        raa, rawc, raws = resist_out.raa, resist_out.rawc, resist_out.raws
        rac = (LAI + 1) * (raa + rawc)
        ras = (LAI + 1) * (raa + raws)

        out_ch = get_heatfluxes(ra=rac, rs=data_bch.rcw, Tc=Tch, ea=ea, Ta=Ta, e_to_q=e_to_q, Ca=Ca, Ci=data_bch.Ci)
        lEch, Hch, ech, Cch = out_ch.lE, out_ch.H, out_ch.ec, out_ch.Cc
        lambdah, sh = out_ch.lambda_, out_ch.s

        out_cu = get_heatfluxes(ra=rac, rs=data_bcu.rcw, Tc=Tcu, ea=ea, Ta=Ta, e_to_q=e_to_q, Ca=Ca, Ci=data_bcu.Ci)
        lEcu, Hcu, ecu, Ccu = out_cu.lE, out_cu.H, out_cu.ec, out_cu.Cc
        lambdau, su = out_cu.lambda_, out_cu.s

        out_s = get_heatfluxes(ra=ras, rs=soil.rss, Tc=Ts, ea=ea, Ta=Ta, e_to_q=e_to_q, Ca=Ca, Ci=Ca)
        lEs, Hs = out_s.lE, out_s.H
        lambdas, ss = out_s.lambda_, out_s.s

        Hstot = float(np.sum(Fs * Hs))
        Hctot = aggregator_ebal(LAI, Hcu, Hch, Ps[:nl], nl)
        Htot = Hstot + Hctot

        if use_monin_obukhov:
            L = float(monin_obukhov(resist_out.ustar, Ta, Htot))

        G = 0.35 * Rns
        dG = 4 * (1 - soil.rs_thermal) * sigmaSB * (Ts + 273.15) ** 3 * 0.35

        EBerch = Rnhc - lEch - Hch
        EBercu = Rnuc - lEcu - Hcu
        EBers = Rns - lEs - Hs - G

        counter += 1
        maxEBercu = float(np.max(np.abs(EBercu)))
        maxEBerch = float(np.max(np.abs(EBerch)))
        maxEBers = float(np.max(np.abs(EBers)))

        CONT = (maxEBercu > maxEBer or maxEBerch > maxEBer or maxEBers > maxEBer) and (counter < k_maxit + 1)
        if not CONT:
            break

        if counter == 10:
            Wc = 0.8
        if counter == 20:
            Wc = 0.6

        Tch = Tch + Wc * EBerch / (
            (rhoa * cp) / rac
            + rhoa * lambdah * e_to_q * sh / (rac + data_bch.rcw)
            + 4 * leaf_emis * sigmaSB * (Tch + 273.15) ** 3
        )
        Tcu = Tcu + Wc * EBercu / (
            (rhoa * cp) / rac
            + rhoa * lambdau * e_to_q * su / (rac + data_bcu.rcw)
            + 4 * leaf_emis * sigmaSB * (Tcu + 273.15) ** 3
        )
        Ts = Ts + Wc * EBers / (
            rhoa * cp / ras
            + rhoa * lambdas * e_to_q * ss / (ras + soil.rss)
            + 4 * (1 - soil.rs_thermal) * sigmaSB * (Ts + 273.15) ** 3
            + dG
        )

        Tch = np.where(np.abs(Tch) > 100, Ta, Tch)
        Tcu = np.where(np.abs(Tcu) > 100, Ta, Tcu)

    # emissivity (real leaf/soil vs. a hypothetical black surface at the same temperatures)
    rtmt_final = rtmt_sb(
        rtmo=rtmo, nl=nl, LAI=LAI, rho_thermal=rho_thermal, tau_thermal=tau_thermal,
        rs_thermal=soil.rs_thermal, Tcu=Tcu, Tch=Tch, Tsu=Ts[1], Tsh=Ts[0],
    )
    rtmt_black = rtmt_sb(
        rtmo=rtmo, nl=nl, LAI=LAI, rho_thermal=0.0, tau_thermal=0.0,
        rs_thermal=0.0, Tcu=Tcu, Tch=Tch, Tsu=Ts[1], Tsh=Ts[0],
    )
    canopyemis = float(rtmt_final.Eoutte / rtmt_black.Eoutte)

    Rnctot = aggregator_ebal(LAI, Rnuc, Rnhc, Fc, nl)
    lEctot = aggregator_ebal(LAI, lEcu, lEch, Fc, nl)
    Hctot_final = aggregator_ebal(LAI, Hcu, Hch, Fc, nl)
    Actot = aggregator_ebal(LAI, data_bcu.A, data_bch.A, Fc, nl)
    Tcave = aggregator_ebal(1.0, Tcu, Tch, Fc, nl)
    Rnstot = float(np.sum(Fs * Rns))
    lEstot = float(np.sum(Fs * lEs))
    Hstot_final = float(np.sum(Fs * Hs))
    Gtot = float(np.sum(Fs * G))
    Tsave = float(np.sum(Fs * Ts))
    Rntot = Rnctot + Rnstot
    lEtot = lEctot + lEstot
    Htot_final = Hctot_final + Hstot_final

    return EbalResult(
        counter=counter, Tcu=Tcu, Tch=Tch, Tsu=float(Ts[1]), Tsh=float(Ts[0]),
        bcu=data_bcu, bch=data_bch, canopyemis=canopyemis,
        Rnctot=Rnctot, lEctot=lEctot, Hctot=Hctot_final, Actot=Actot, Tcave=Tcave,
        Rnstot=Rnstot, lEstot=lEstot, Hstot=Hstot_final, Gtot=Gtot, Tsave=Tsave,
        Rntot=Rntot, lEtot=lEtot, Htot=Htot_final,
        maxEBercu=maxEBercu, maxEBerch=maxEBerch, maxEBers=maxEBers,
    )
