"""RTMf: SCOPE canopy fluorescence radiative transfer model.

Direct port of ``SCOPEinR::get.RTMf`` (``SCOPEinR/R/RTMf.R``), computing
the TOC fluorescence radiance in the observation direction and the TOC
hemispherical upward fluorescence flux, given per-canopy-layer
fluorescence excitation-emission matrices (``Mb``/``Mf``, from
:func:`scopeinpython.fluspect_mscope.fluspect_mscope`) and per-layer
fluorescence quantum efficiencies (``etau`` for sunlit leaves, ``etah``
for shaded leaves -- from :func:`scopeinpython.biochemical.get_biochemical`,
called once per leaf micro-environment; NOT computed here, matching how
the R function itself takes them as plain arguments rather than an
iterative energy-balance dependency, so this is composable without the
SCOPE thermal energy-balance loop, ``ebal.R``, being ported).

**Ported against a fixed R source.** ``SCOPEinR::get.RTMf`` had two real
bugs, found and fixed during this port (see ``python/docs/verification.rst``
for the full history):

- A column-recycling issue in ~18 expressions computing
  ``wfEs``/``sfEs``/``sbEs``/``vfEplu_*``/``vbEmin_*``/``sigfEmin_*``/
  ``sigbEmin_*``/``sigfEplu_*``/``sigbEplu_*``/``Femmin``/``Femplu``:
  R's ``vector * matrix`` only broadcasts correctly per matrix-column
  when ``nrow(matrix) == length(vector)``, which didn't hold here (the
  matrices have ``nf`` (53) rows, the vectors have ``nl`` (canopy-layer
  count) elements). Found in two passes -- the first covered the 16
  expressions building ``wfEs`` through ``sigbEplu_u/h``; a second,
  closer look (triggered by this Python port's own reference values
  disagreeing with an R run at the time, which is how the *first* pass's
  fix was confirmed complete but this remaining pair was found to still
  need it) turned up two more instances two lines later, in
  ``Femmin``/``Femplu`` themselves (``Qs * Fsmin`` etc). All fixed in R
  with ``sweep()``. NumPy needs no equivalent fix anywhere:
  ``matrix * vector_of_length_ncol`` already broadcasts correctly by
  construction (NumPy aligns trailing axes; R recycles a flat vector), so
  the natural direct translation below is correct without any
  bug-workaround helper.
- An ``absfs_nl <- c(absfsfo)`` copy-paste mixup (should have been
  ``c(absfs)``), fixed in R to ``absfs_nl <- c(absfs)``.

Both fixes are reflected directly in the derivation below (i.e. this
module was written against the corrected R source, not translated from
the buggy version and then patched).

**A real, deliberate numerical approximation, not an exact match**: R's
``signal::interp1(..., method='spline')`` calls ``stats::splinefun()``,
whose default method (``"fmm"``, Forsythe-Malcolm-Moler) estimates
end-point derivatives via the unique cubic through the first/last 4
points -- a different construction from any boundary condition
:class:`scipy.interpolate.CubicSpline` offers directly. This port uses
``bc_type='not-a-knot'``, which is close (confirmed via a standalone
R-vs-Python comparison on synthetic data: the discrepancy is localized to
within a few nm of each end of the native 640-850 nm, 4 nm-step
fluorescence grid, decaying to ~0 a few points in; ~1.5e-3 absolute worst
case against a spectrum spanning a range of ~2) but not bit-identical.
Every other quantity in this module (everything computed *before* the
spline step -- ``piLo1``..``piLo4``, the native-grid ``LoF_``/``Fhem_``,
etc.) matches R to the same floating-point precision as the rest of this
port; only the final upsampling from the 53-point native fluorescence
grid to the 211-point ``spectral.wlF`` display grid carries this
localized approximation. R's own ``extrap=NA`` default (unused here
since ``signal::interp1`` is called without an explicit ``extrap``
argument) means R itself leaves the last 2 of ``spectral.wlF``'s 211
points (849, 850 nm, beyond the native grid's 848 nm endpoint) as ``NaN``
-- reproduced here via ``CubicSpline(..., extrapolate=False)`` for those
two points specifically (see :func:`_interp_wlf`).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline

from ._data import constants as _constants
from .fluspect import _r_seq
from .rtmo import CanopyStructure, RTMoResult
from .spectral import SpectralConfig
from .utils import get_e2phot, get_ephoton, sint

__all__ = ["RTMfResult", "rtmf"]

_WLF_NATIVE = _r_seq(640.0, 850.0, 4.0)  # 53 pts -- matches fluspect_mscope(step=5)'s Mb/Mf wlf grid
_WLE_NATIVE = _r_seq(400.0, 750.0, 5.0)  # 71 pts -- matches its wle grid
_IWLFI = (_WLE_NATIVE - 400).astype(int)  # 0-based positions of wlE within wlS (both start at 400nm)
_IWLFO = (_WLF_NATIVE - 400).astype(int)  # 0-based positions of wlF within wlS


def _interp_wlf(y_native: np.ndarray, wlf_full: np.ndarray) -> np.ndarray:
    """Upsample a spectrum on ``_WLF_NATIVE`` (640-850nm, 4nm step) onto the
    full ``wlf_full`` grid (640-850nm, 1nm step), matching R's
    ``signal::interp1(_WLF_NATIVE, y_native, xi=wlf_full, method='spline')``
    -- including its default ``extrap=NA`` behaviour (points of ``wlf_full``
    outside ``[_WLF_NATIVE.min(), _WLF_NATIVE.max()]`` become NaN, not
    extrapolated). See module docstring for the ``not-a-knot`` vs R's
    ``fmm`` spline approximation."""
    cs = CubicSpline(_WLF_NATIVE, y_native, bc_type="not-a-knot", extrapolate=False)
    return cs(wlf_full)


@dataclass
class RTMfResult:
    LoF_: np.ndarray  # TOC fluorescence radiance in viewing direction, on spectral.wlF grid (mW m-2 um-1 sr-1)
    EoutF_: np.ndarray  # TOC hemispherical fluorescence flux spectrum, on spectral.wlF grid
    LoF_sunlit: np.ndarray
    LoF_shaded: np.ndarray
    LoF_scattered: np.ndarray
    LoF_soil: np.ndarray
    EoutF: float  # spectrally-integrated hemispherical flux (W m-2)
    LoutF: float  # spectrally-integrated radiance
    Femliave_: np.ndarray
    F685: float  # peak fluorescence value in the 685nm band, or NaN if the peak is at the band edge
    wl685: float
    F740: float
    wl740: float
    F684: float  # fluorescence value at exactly 684nm
    F761: float  # fluorescence value at exactly 761nm


def rtmf(
    spectral: SpectralConfig,
    rtmo: RTMoResult,
    canopy: CanopyStructure,
    tts: float,
    tto: float,
    psi: float,
    rsoil: np.ndarray,
    Mb: np.ndarray,
    Mf: np.ndarray,
    etau: np.ndarray,
    etah: np.ndarray,
) -> RTMfResult:
    """Direct port of ``SCOPEinR::get.RTMf`` (against the fixed R source,
    see module docstring).

    Parameters
    ----------
    rtmo : RTMoResult
        From :func:`scopeinpython.rtmo.run_rtmo`, called with the same
        ``canopy``/``tts``/``tto``/``psi``/``rsoil`` as here.
    rsoil : array_like, shape (2001,)
        Same 400-2400nm soil reflectance array passed to ``run_rtmo``.
    Mb, Mf : array_like, shape (53, 71, nl)
        From :func:`scopeinpython.fluspect_mscope.fluspect_mscope` called
        with ``step=5`` (its default) and this same ``nl``.
    etau : array_like, shape (nl, 13, 36)
        Sunlit-leaf fluorescence quantum efficiency (``eta``, from
        :func:`scopeinpython.biochemical.get_biochemical`) per canopy
        layer x leaf-inclination class x leaf-azimuth class.
    etah : array_like, shape (nl,)
        Shaded-leaf fluorescence quantum efficiency per canopy layer.
    """
    const = _constants()
    nf = len(_IWLFO)  # 53
    nl = canopy.nlayers
    LAI = canopy.LAI
    litab = canopy.litab
    lazitab = canopy.lazitab
    lidf = canopy.lidf
    nlazi = len(lazitab)  # 36
    nlinc = len(litab)  # 13
    nlori = nlinc * nlazi  # 468

    Ps, Po, Pso = rtmo.Ps, rtmo.Po, rtmo.Pso
    Qs = Ps[:nl]

    Esunf_ = rtmo.Esun_[_IWLFI]  # (71,)
    Eminf_ = rtmo.Emin_[:nl, _IWLFI].T  # (71, nl) -- Emin_ at the TOP of each of the nl layers
    Epluf_ = rtmo.Eplu_[:nl, _IWLFI].T  # (71, nl)

    iLAI = LAI / nl

    Xdd = rtmo.Xdd[:, _IWLFO]  # (nl, 53)
    rho_dd = rtmo.rho_dd[:, _IWLFO]  # (nl, 53)
    R_dd = rtmo.R_dd[:, _IWLFO]  # (nl+1, 53)
    tau_dd = rtmo.tau_dd[:, _IWLFO]  # (nl, 53)
    vb = rtmo.vb[:, _IWLFO]  # (nl, 53)
    vf = rtmo.vf[:, _IWLFO]  # (nl, 53)

    deg2rad = np.pi / 180.0
    rs = np.asarray(rsoil, dtype=float)[_IWLFO]  # (53,)
    cos_tto, sin_tto = np.cos(tto * deg2rad), np.sin(tto * deg2rad)
    cos_tts, sin_tts = np.cos(tts * deg2rad), np.sin(tts * deg2rad)
    cos_ttli = np.cos(litab * deg2rad)  # (13,)
    sin_ttli = np.sin(litab * deg2rad)
    cos_phils = np.cos(lazitab * deg2rad)  # (36,)
    cos_philo = np.cos((lazitab - psi) * deg2rad)

    cds = cos_ttli[:, None] * cos_tts + sin_ttli[:, None] * (sin_tts * cos_phils)[None, :]  # (13, 36)
    cdo = cos_ttli[:, None] * cos_tto + sin_ttli[:, None] * (sin_tto * cos_philo)[None, :]

    fs = cds / cos_tts
    absfs = np.abs(fs)
    fo = cdo / cos_tto
    absfo = np.abs(fo)
    fsfo = fs * fo
    absfsfo = np.abs(fsfo)
    foctl = fo * cos_ttli[:, None]
    fsctl = fs * cos_ttli[:, None]
    ctl2 = np.broadcast_to(cos_ttli[:, None] ** 2, (nlinc, nlazi))

    # flatten (nlinc, nlazi) -> (nlori,) column-major (inclination fastest),
    # matching R's `matrix(x, nlori, 1, byrow=F)`
    absfs_nl = absfs.flatten(order="F")
    absfo_nl = absfo.flatten(order="F")
    fsfo_nl = fsfo.flatten(order="F")
    absfsfo_nl = absfsfo.flatten(order="F")
    foctl_nl = foctl.flatten(order="F")
    fsctl_nl = fsctl.flatten(order="F")
    ctl2_nl = ctl2.flatten(order="F")

    Mplu = 0.5 * (Mb + Mf)  # (53, 71, nl)
    Mmin = 0.5 * (Mb - Mf)

    ep = const["A"] * get_ephoton(_WLF_NATIVE * 1e-9, const)  # (53,) energy per mole of photons at each emission wl

    MpluEmin = np.empty((nf, nl))
    MpluEplu = np.empty((nf, nl))
    MminEmin = np.empty((nf, nl))
    MminEplu = np.empty((nf, nl))
    MpluEsun = np.empty((nf, nl))
    MminEsun = np.empty((nf, nl))
    for j in range(nl):
        e2phot_min = get_e2phot(_WLE_NATIVE * 1e-9, Eminf_[:, j], const)
        e2phot_plu = get_e2phot(_WLE_NATIVE * 1e-9, Epluf_[:, j], const)
        e2phot_sun = get_e2phot(_WLE_NATIVE * 1e-9, Esunf_, const)
        MpluEmin[:, j] = ep * (Mplu[:, :, j] @ e2phot_min)
        MpluEplu[:, j] = ep * (Mplu[:, :, j] @ e2phot_plu)
        MminEmin[:, j] = ep * (Mmin[:, :, j] @ e2phot_min)
        MminEplu[:, j] = ep * (Mmin[:, :, j] @ e2phot_plu)
        MpluEsun[:, j] = ep * (Mplu[:, :, j] @ e2phot_sun)
        MminEsun[:, j] = ep * (Mmin[:, :, j] @ e2phot_sun)

    laz = 1.0 / 36
    etau_perm = np.transpose(np.asarray(etau, dtype=float), (1, 2, 0))  # (13, 36, nl)
    etau_reshape = etau_perm.reshape(nlori, nl, order="F")  # (468, nl)
    lidf_laz = np.tile(np.asarray(lidf, dtype=float) * laz, nlazi)  # (468,)

    etau_lidf = etau_reshape * lidf_laz[:, None]  # (468, nl)
    etah_lidf = np.outer(lidf_laz, np.asarray(etah, dtype=float))  # (468, nl)

    # Each `(etau_lidf * w_nl[:, None]).sum(axis=0)` is a length-nl per-layer
    # scalar; multiplying it against an (nf, nl) matrix broadcasts correctly
    # per column (NumPy aligns trailing axes) -- the direct equivalent of
    # the sweep()-fixed R expressions, see module docstring.
    wfEs = MpluEsun * (etau_lidf * absfsfo_nl[:, None]).sum(axis=0) \
        + MminEsun * (etau_lidf * fsfo_nl[:, None]).sum(axis=0)

    sfEs = MpluEsun * (etau_lidf * absfs_nl[:, None]).sum(axis=0) \
        - MminEsun * (etau_lidf * fsctl_nl[:, None]).sum(axis=0)
    sbEs = MpluEsun * (etau_lidf * absfs_nl[:, None]).sum(axis=0) \
        + MminEsun * (etau_lidf * fsctl_nl[:, None]).sum(axis=0)

    vfEplu_h = MpluEplu * (etah_lidf * absfo_nl[:, None]).sum(axis=0) \
        - MminEplu * (etah_lidf * foctl_nl[:, None]).sum(axis=0)
    vfEplu_u = MpluEplu * (etau_lidf * absfo_nl[:, None]).sum(axis=0) \
        - MminEplu * (etau_lidf * foctl_nl[:, None]).sum(axis=0)

    vbEmin_h = MpluEmin * (etah_lidf * absfo_nl[:, None]).sum(axis=0) \
        + MminEmin * (etah_lidf * foctl_nl[:, None]).sum(axis=0)
    vbEmin_u = MpluEmin * (etau_lidf * absfo_nl[:, None]).sum(axis=0) \
        + MminEmin * (etau_lidf * foctl_nl[:, None]).sum(axis=0)

    sigfEmin_h = MpluEmin * etah_lidf.sum(axis=0) - MminEmin * (etah_lidf * ctl2_nl[:, None]).sum(axis=0)
    sigfEmin_u = MpluEmin * etau_lidf.sum(axis=0) - MminEmin * (etau_lidf * ctl2_nl[:, None]).sum(axis=0)
    sigbEmin_h = MpluEmin * etah_lidf.sum(axis=0) + MminEmin * (etah_lidf * ctl2_nl[:, None]).sum(axis=0)
    sigbEmin_u = MpluEmin * etau_lidf.sum(axis=0) + MminEmin * (etau_lidf * ctl2_nl[:, None]).sum(axis=0)

    sigfEplu_h = MpluEplu * etah_lidf.sum(axis=0) - MminEplu * (etah_lidf * ctl2_nl[:, None]).sum(axis=0)
    sigfEplu_u = MpluEplu * etau_lidf.sum(axis=0) - MminEplu * (etau_lidf * ctl2_nl[:, None]).sum(axis=0)
    sigbEplu_h = MpluEplu * etah_lidf.sum(axis=0) + MminEplu * (etah_lidf * ctl2_nl[:, None]).sum(axis=0)
    sigbEplu_u = MpluEplu * etau_lidf.sum(axis=0) + MminEplu * (etau_lidf * ctl2_nl[:, None]).sum(axis=0)

    ################################################################
    #  Emitted fluorescence
    ################################################################
    piLs = wfEs + vfEplu_u + vbEmin_u  # sunlit, per layer
    piLd = vbEmin_h + vfEplu_h  # shaded, per layer

    Fsmin = sfEs + sigfEmin_u + sigbEplu_u
    Fsplu = sbEs + sigbEmin_u + sigfEplu_u
    Fdmin = sigfEmin_h + sigbEplu_h
    Fdplu = sigbEmin_h + sigfEplu_h

    Femmin = iLAI * (Qs * Fsmin) + iLAI * ((1 - Qs) * Fdmin)
    Femplu = iLAI * (Qs * Fsplu) + iLAI * ((1 - Qs) * Fdplu)

    Y = np.zeros((nl, nf))
    U = np.zeros((nl + 1, nf))
    Fmin_ = np.zeros((nl + 1, nf))
    Fplu_ = np.zeros((nl + 1, nf))

    # from bottom to top
    for j in range(nl - 1, -1, -1):
        Y[j, :] = (rho_dd[j, :] * U[j + 1, :] + Femmin[:, j]) / (1 - rho_dd[j, :] * R_dd[j + 1, :])
        U[j, :] = tau_dd[j, :] * (R_dd[j + 1, :] * Y[j, :] + U[j + 1, :]) + Femplu[:, j]

    # from top to bottom
    for j in range(nl):
        Fmin_[j + 1, :] = Xdd[j, :] * Fmin_[j, :] + Y[j, :]
        Fplu_[j, :] = R_dd[j, :] * Fmin_[j, :] + U[j, :]

    piLo1 = iLAI * (piLs @ Pso[:nl])
    piLo2 = iLAI * (piLd @ (Po[:nl] - Pso[:nl]))
    piLo3 = iLAI * ((vb * Fmin_[:nl, :] + vf * Fplu_[:nl, :]).T @ Po[:nl])
    piLo4 = rs * Fmin_[nl, :] * Po[nl]

    piLtot = piLo1 + piLo2 + piLo3 + piLo4

    LoF_native = piLtot / np.pi
    Fhem_native = Fplu_[0, :]

    ###############################################################################################
    # Output
    ###############################################################################################
    wlf_full = spectral.wlF
    LoF_ = _interp_wlf(LoF_native, wlf_full)
    EoutF_ = _interp_wlf(Fhem_native, wlf_full)
    LoF_sunlit = _interp_wlf(piLo1 / np.pi, wlf_full)
    LoF_shaded = _interp_wlf(piLo2 / np.pi, wlf_full)
    LoF_scattered = _interp_wlf(piLo3 / np.pi, wlf_full)
    LoF_soil = _interp_wlf(piLo4 / np.pi, wlf_full)

    EoutF = 0.001 * sint(Fhem_native, _WLF_NATIVE)
    LoutF = 0.001 * sint(LoF_native, _WLF_NATIVE)

    Femliave_ = _interp_wlf((Femmin + Femplu).sum(axis=1), wlf_full)

    F685 = np.nanmax(LoF_[:55])
    iwl685 = int(np.nanargmax(LoF_[:55]))
    wl685 = wlf_full[iwl685]
    if iwl685 == 54:  # R's `iwl685 == 55` (1-based) -> 0-based index 54
        F685 = np.nan
        wl685 = np.nan

    F740 = np.nanmax(LoF_[69:])
    wl740 = wlf_full[int(np.nanargmax(LoF_[69:])) + 69]
    F684 = LoF_[685 - int(wlf_full[0])]
    F761 = LoF_[762 - int(wlf_full[0])]

    return RTMfResult(
        LoF_=LoF_, EoutF_=EoutF_, LoF_sunlit=LoF_sunlit, LoF_shaded=LoF_shaded,
        LoF_scattered=LoF_scattered, LoF_soil=LoF_soil, EoutF=EoutF, LoutF=LoutF,
        Femliave_=Femliave_, F685=F685, wl685=wl685, F740=F740, wl740=wl740,
        F684=F684, F761=F761,
    )
