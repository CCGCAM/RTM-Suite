"""RTMz: SCOPE canopy zeaxanthin ("Z") radiative transfer model.

Direct port of ``SCOPEinR::get.RTMz`` (``SCOPEinR/R/RTMz.R``), computing
the small modification of TOC outgoing radiance due to the
violaxanthin-to-zeaxanthin conversion in leaves (a photoprotection
signal, 500-600 nm). Structurally very similar to
:func:`scopeinpython.rtmf.rtmf` (same geometric-factor setup, same
layer-recursion pattern) but simpler: RTMz works directly in the 500-600
nm wavelength domain (no excitation-emission matrix like ``Mb``/``Mf`` --
the zeaxanthin signal is a per-wavelength reflectance/transmittance
*difference*, ``reflZ - refl`` / ``tranZ - tran``, not a fluorescence
redistribution), and produces a *correction* to be added onto an existing
:class:`~scopeinpython.rtmo.RTMoResult`'s ``Lo_``/``rso``/``rdo``/``rfl``/
``Eout_`` at the 500-600 nm band, rather than a self-contained new
spectrum on its own display grid (RTMf's ``LoF_`` etc).

**Ported against a fixed R source.** ``SCOPEinR::get.RTMz`` had four real
bugs, found and fixed during this port (see ``python/docs/verification.rst``
for the full history):

- The exact same column-recycling issue as ``get.RTMf`` (``vector *
  matrix`` only broadcasts correctly per column when
  ``nrow(matrix) == length(vector)``), across the same family of
  expressions plus its own ``Femmin``/``Femplu``. Fixed with ``sweep()``.
- Its real output (the layer recursion, ``piLo1``..``piLo4``, ``LoF_``)
  was computed entirely inside an *unreachable* ``if (get.plots ==
  TRUE)`` block placed **after** the loop that produced its inputs --
  since that loop's `k`/`piLs`/etc. variables aren't `k`-indexed there,
  the block only ever used the LAST iteration's values, and only ran at
  all when ``get.plots = TRUE``. With the normal ``get.plots = FALSE``
  call, ``get.RTMz()`` silently returned its input **completely
  unmodified**. Fixed by moving the real computation inside the loop,
  unconditionally.
- A ``Po[1:nl+1]`` bug (R parses this as ``(1:nl)+1``, a length-``nl``
  vector, not the intended scalar ``Po[nl+1]``) -- fixed to match the
  equivalent, correctly-written line in ``get.RTMf``.
- A hardcoded ``dim = c(30, 13, 36)`` (should be ``c(nl, 13, 36)``) in the
  ``data.Knu``-is-a-plain-vector branch of its ``etau`` construction --
  silently wrong for any canopy with ``nl != 30``. Fixed to use ``nl``.
- A MATLAB-to-R mistranslation: ``sum(LoF_, 2)`` (meant to sum across the
  sun/sky (``k``) dimension, MATLAB's ``sum(A,2)`` idiom) -- R's
  ``sum()`` has no dimension argument, so this collapsed the whole
  ``(nwlZ, 2)`` matrix to a single scalar (total + 2) and recycled it
  across every wavelength, instead of the intended per-wavelength
  sun+sky total. Fixed to ``rowSums(LoF_)``.

Like :func:`scopeinpython.rtmf.rtmf`, none of these fixes need an
equivalent "bug-workaround" on the Python side: NumPy's default
broadcasting (aligns trailing axes) already does what R needed
``sweep()`` for.

**A likely 7th bug, flagged but NOT fixed here** (uncertain, unlike the
above which were provable from R's documented language semantics):
``vfEplu_u`` uses ``etau_lidf * fsctl_nl`` for its ``Mmin`` term, while
the structurally-parallel ``vfEplu_h`` (and *both* of ``get.RTMf``'s
``vfEplu_h``/``vfEplu_u``) consistently use ``foctl_nl`` there. This
looks like a copy-paste from the ``sfEs``/``sbEs`` block just above (the
one place ``fsctl_nl`` genuinely belongs) that wasn't updated to
``foctl_nl``. Ported here exactly as currently written in R
(``fsctl_nl``, matching the below code) since there's no independent way
to *prove* this one the way the recycling bug was proved from R's own
documented `*` semantics -- worth a second look before relying on RTMz
for real science.

**A design note on ``etau``'s R-side reshape, NOT independently
confirmed as a bug and NOT changed in R**: unlike ``get.RTMf`` (which
permutes its per-layer ``etau`` array to ``(13, 36, nl)`` before
flattening to ``(468, nl)``, ensuring each column is one whole layer's
468 orientation values), ``get.RTMz``'s plain-vector ``data.Knu`` branch
skips that permute step, reshaping the "natural" ``(nl, 13, 36)``
[layer, inclination, azimuth] array directly into ``(468, nl)`` --
which, by the mechanics of R's column-major array storage, does *not*
recover "column = one layer's 468 values" unless the original
``data.Knu`` vector was already constructed in a matching order by its
caller (not verified here; would require reading ``get.SCOPE.R``'s call
site, out of scope for this port). Since a uniform ``etau`` (this port's
only available verification case, absent the not-yet-ported thermal
energy-balance loop that would produce a genuinely varying one) can't
distinguish a correct reshape from a scrambled-but-uniform one, this
is flagged here rather than silently assumed either way. The Python
``etau`` parameter below takes the unambiguous ``(nl, 13, 36)`` shape
directly (same design choice as :func:`scopeinpython.rtmf.rtmf`) and
applies the *same* correctly-permuting reshape ``get.RTMf`` uses, since
that is unambiguously the physically-intended construction regardless of
this R-side question.

**A confirmed, separate issue in R's OTHER ``data.Knu`` branch** (hit
when ``data.Knu`` has a `dim` attribute, e.g. is itself already an
array rather than a plain vector -- discovered while building this
port's own R reference case, which originally passed an array there by
mistake and got visibly wrong numbers as a result): that branch never
calls ``Kn2Cx()`` on ``data.Knu`` at all, unlike the plain-vector branch
-- so the empirical NPQ-to-zeaxanthin conversion is silently skipped
whenever a caller passes an already-shaped array. This port's R
reference case (``python/scratch/scratch_rtmz_export.R``) was fixed to pass a
genuine plain vector (hitting the correctly-behaving, `Kn2Cx`-applying
branch) rather than chase down whether the array branch is fixable too.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .fluspect import _r_seq
from .rtmo import CanopyStructure, RTMoResult
from .spectral import SpectralConfig

__all__ = ["RTMzResult", "rtmz"]

_WLZ_NATIVE = _r_seq(500.0, 600.0, 1.0)  # 101 pts, matches R's `dummy <- c(wlZ[1]:wlZ[2])`
_IWLFI = (_WLZ_NATIVE - 400).astype(int)  # 0-based positions within wlS (and within a (nl,2001) leaf-optics array)


def _kn2cx(Kn: np.ndarray) -> np.ndarray:
    """Empirical Kn (NPQ) -> zeaxanthin-related Cx conversion (Vilfan et
    al. 2018, 2019). Direct port of R's inline ``Kn2Cx``."""
    return 0.3187 * np.asarray(Kn, dtype=float)


@dataclass
class RTMzResult:
    """Corrections to be ADDED to the corresponding 500-600nm slice of an
    existing :class:`~scopeinpython.rtmo.RTMoResult` (indices
    ``_IWLFI`` = ``wlZ_native - 400`` into any ``wlS``-grid array, or
    equivalently ``spectral.wlS[_IWLFI] == _WLZ_NATIVE``)."""

    Lo_delta: np.ndarray  # add to rtmo.Lo_ (well, to data.rad['Lo_'], not currently on RTMoResult) at wlZ positions -- see rtmz() docstring
    rso_delta: np.ndarray  # add to rtmo.rso at wlZ positions
    rdo_delta: np.ndarray  # add to rtmo.rdo at wlZ positions
    Eout_delta: np.ndarray  # add to rtmo.Eout_ at wlZ positions
    LoF_: np.ndarray  # (101, 2) raw sun/sky outgoing zeaxanthin radiance, pi*Lo units (piLtot/pi per k)


def rtmz(
    spectral: SpectralConfig,
    rtmo: RTMoResult,
    canopy: CanopyStructure,
    tts: float,
    tto: float,
    psi: float,
    rsoil: np.ndarray,
    refl: np.ndarray,
    tran: np.ndarray,
    reflZ: np.ndarray,
    tranZ: np.ndarray,
    Knu: np.ndarray,
    Knh: np.ndarray,
) -> RTMzResult:
    """Direct port of ``SCOPEinR::get.RTMz`` (against the fixed R source,
    see module docstring).

    Parameters
    ----------
    rtmo : RTMoResult
        From :func:`scopeinpython.rtmo.run_rtmo`, called with the same
        ``canopy``/``tts``/``tto``/``psi``/``rsoil`` as here.
    rsoil : array_like, shape (2001,)
        Same 400-2400nm soil reflectance array passed to ``run_rtmo``.
    refl, tran : array_like, shape (nl, 2001)
        Baseline (``Cx=0``) leaf reflectance/transmittance, from
        :func:`scopeinpython.fluspect_mscope.fluspect_mscope`.
    reflZ, tranZ : array_like, shape (nl, 2001)
        Full-zeaxanthin (``Cx=1``) leaf reflectance/transmittance, from a
        second :func:`~scopeinpython.fluspect_mscope.fluspect_mscope`
        call with the same ``mly``/``nl``/``step`` but ``Cx=1``.
    Knu : array_like, shape (nl, 13, 36)
        Sunlit-leaf NPQ (``Kn``, from
        :func:`scopeinpython.biochemical.get_biochemical`) per canopy
        layer x leaf-inclination class x leaf-azimuth class.
    Knh : array_like, shape (nl,)
        Shaded-leaf NPQ per canopy layer.
    """
    nl = canopy.nlayers
    LAI = canopy.LAI
    litab = canopy.litab
    lazitab = canopy.lazitab
    lidf = canopy.lidf
    nlazi = len(lazitab)
    nlinc = len(litab)
    nlori = nlinc * nlazi
    nwlZ = len(_IWLFI)  # 101

    Ps, Po, Pso = rtmo.Ps, rtmo.Po, rtmo.Pso
    Qs = Ps[:nl]

    RZ = (np.asarray(reflZ, dtype=float)[:, _IWLFI] - np.asarray(refl, dtype=float)[:, _IWLFI]).T  # (101, nl)
    TZ = (np.asarray(tranZ, dtype=float)[:, _IWLFI] - np.asarray(tran, dtype=float)[:, _IWLFI]).T

    iLAI = LAI / nl

    Xdd = rtmo.Xdd[:, _IWLFI]
    rho_dd = rtmo.rho_dd[:, _IWLFI]
    R_dd = rtmo.R_dd[:, _IWLFI]
    tau_dd = rtmo.tau_dd[:, _IWLFI]
    vb = rtmo.vb[:, _IWLFI]
    vf = rtmo.vf[:, _IWLFI]

    deg2rad = np.pi / 180.0
    rs = np.asarray(rsoil, dtype=float)[_IWLFI]
    cos_tto, sin_tto = np.cos(tto * deg2rad), np.sin(tto * deg2rad)
    cos_tts, sin_tts = np.cos(tts * deg2rad), np.sin(tts * deg2rad)
    cos_ttli = np.cos(litab * deg2rad)
    sin_ttli = np.sin(litab * deg2rad)
    cos_phils = np.cos(lazitab * deg2rad)
    cos_philo = np.cos((lazitab - psi) * deg2rad)

    cds = cos_ttli[:, None] * cos_tts + sin_ttli[:, None] * (sin_tts * cos_phils)[None, :]
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

    absfs_nl = absfs.flatten(order="F")
    absfo_nl = absfo.flatten(order="F")
    fsfo_nl = fsfo.flatten(order="F")
    absfsfo_nl = absfsfo.flatten(order="F")
    foctl_nl = foctl.flatten(order="F")
    fsctl_nl = fsctl.flatten(order="F")
    ctl2_nl = ctl2.flatten(order="F")

    Esunf_ = rtmo.Esun_[_IWLFI]  # (101,)
    Eminf_ = np.stack([rtmo.Emins_[:nl, _IWLFI].T, rtmo.Emind_[:nl, _IWLFI].T], axis=-1)  # (101, nl, 2)
    Epluf_ = np.stack([rtmo.Eplus_[:nl, _IWLFI].T, rtmo.Eplud_[:nl, _IWLFI].T], axis=-1)

    laz = 1.0 / 36
    etah = _kn2cx(Knh)  # (nl,)
    etau_full = _kn2cx(np.asarray(Knu, dtype=float))  # (nl, 13, 36)
    etau_perm = np.transpose(etau_full, (1, 2, 0))  # (13, 36, nl)
    etau_reshape = etau_perm.reshape(nlori, nl, order="F")  # (468, nl)
    lidf_laz = np.tile(np.asarray(lidf, dtype=float) * laz, nlazi)  # (468,)

    etau_lidf = etau_reshape * lidf_laz[:, None]
    etah_lidf = np.outer(lidf_laz, etah)

    LoF_ = np.zeros((nwlZ, 2))
    Fmin_ = np.zeros((nl + 1, nwlZ, 2))
    Fplu_ = np.zeros((nl + 1, nwlZ, 2))

    # Esunf_ has length nwlZ == nrow(RZ), so this is a genuine per-row
    # (per-wavelength) scale -- explicit [:, None] needed since NumPy
    # broadcasting aligns TRAILING axes by default (would otherwise try to
    # match against RZ's nl axis, not its nwlZ axis).
    MpluEsun = RZ * Esunf_[:, None]
    MminEsun = TZ * Esunf_[:, None]

    for k in range(2):
        MpluEmin = RZ * Eminf_[:, :, k]
        MpluEplu = RZ * Epluf_[:, :, k]
        MminEmin = TZ * Eminf_[:, :, k]
        MminEplu = TZ * Epluf_[:, :, k]

        wfEs = MpluEsun * (etau_lidf * absfsfo_nl[:, None]).sum(axis=0) \
            + MminEsun * (etau_lidf * fsfo_nl[:, None]).sum(axis=0)

        sfEs = MpluEsun * (etau_lidf * absfs_nl[:, None]).sum(axis=0) \
            + MminEsun * (etau_lidf * fsctl_nl[:, None]).sum(axis=0)
        sbEs = sfEs  # R: sfEs and sbEs use the identical expression here (both "+")

        vfEplu_h = MpluEplu * (etah_lidf * absfo_nl[:, None]).sum(axis=0) \
            + MminEplu * (etah_lidf * foctl_nl[:, None]).sum(axis=0)
        vfEplu_u = MpluEplu * (etau_lidf * absfo_nl[:, None]).sum(axis=0) \
            + MminEplu * (etau_lidf * fsctl_nl[:, None]).sum(axis=0)

        vbEmin_h = MpluEmin * (etah_lidf * absfo_nl[:, None]).sum(axis=0) \
            + MminEmin * (etah_lidf * foctl_nl[:, None]).sum(axis=0)
        vbEmin_u = MpluEmin * (etau_lidf * absfo_nl[:, None]).sum(axis=0) \
            + MminEmin * (etau_lidf * foctl_nl[:, None]).sum(axis=0)

        sigfEmin_h = MpluEmin * etah_lidf.sum(axis=0) + MminEmin * (etah_lidf * ctl2_nl[:, None]).sum(axis=0)
        sigfEmin_u = MpluEmin * etau_lidf.sum(axis=0) + MminEmin * (etau_lidf * ctl2_nl[:, None]).sum(axis=0)
        sigbEmin_h = sigfEmin_h  # R: identical expressions (both "+")
        sigbEmin_u = sigfEmin_u

        sigbEplu_h = MpluEplu * etah_lidf.sum(axis=0) + MminEplu * (etah_lidf * ctl2_nl[:, None]).sum(axis=0)
        sigfEplu_u = MpluEplu * etau_lidf.sum(axis=0) + MminEplu * (etau_lidf * ctl2_nl[:, None]).sum(axis=0)
        sigfEplu_h = sigbEplu_h  # R: identical expressions
        sigbEplu_u = sigfEplu_u

        piLs = wfEs + vfEplu_u + vbEmin_u
        piLd = vbEmin_h + vfEplu_h
        Fsmin = sfEs + sigfEmin_u + sigbEplu_u
        Fsplu = sbEs + sigbEmin_u + sigfEplu_u
        Fdmin = sigfEmin_h + sigbEplu_h
        Fdplu = sigbEmin_h + sigfEplu_h

        Femmin = iLAI * (Qs * Fsmin) + iLAI * ((1 - Qs) * Fdmin)
        Femplu = iLAI * (Qs * Fsplu) + iLAI * ((1 - Qs) * Fdplu)

        Y = np.zeros((nl, nwlZ))
        U = np.zeros((nl + 1, nwlZ))
        for j in range(nl - 1, -1, -1):
            Y[j, :] = (rho_dd[j, :] * U[j + 1, :] + Femmin[:, j]) / (1 - rho_dd[j, :] * R_dd[j + 1, :])
            U[j, :] = tau_dd[j, :] * (R_dd[j + 1, :] * Y[j, :] + U[j + 1, :]) + Femplu[:, j]
        for j in range(nl):
            Fmin_[j + 1, :, k] = Xdd[j, :] * Fmin_[j, :, k] + Y[j, :]
            Fplu_[j, :, k] = R_dd[j, :] * Fmin_[j, :, k] + U[j, :]

        piLo1 = iLAI * (piLs @ Pso[:nl])
        piLo2 = iLAI * (piLd @ (Po[:nl] - Pso[:nl]))
        piLo3 = iLAI * ((vb * Fmin_[:nl, :, k] + vf * Fplu_[:nl, :, k]).T @ Po[:nl])
        piLo4 = rs * Fmin_[nl, :, k] * Po[nl]

        piLtot = piLo1 + piLo2 + piLo3 + piLo4
        LoF_[:, k] = piLtot / np.pi

    Fhem_ = Fplu_[0, :, :].sum(axis=1)  # colSums(Fplu_[1,,,drop=FALSE]) then rowSums -> sum over k at layer 0

    Lo_delta = LoF_.sum(axis=1)  # fixed: R's `sum(LoF_,2)` bug -> `rowSums(LoF_)`
    rso_delta = LoF_[:, 0] / rtmo.Esun_[_IWLFI]
    rdo_delta = LoF_[:, 1] / rtmo.Esky_[_IWLFI]
    Eout_delta = Fhem_

    return RTMzResult(
        Lo_delta=Lo_delta, rso_delta=rso_delta, rdo_delta=rdo_delta,
        Eout_delta=Eout_delta, LoF_=LoF_,
    )
