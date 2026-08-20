"""Multi-layer (mSCOPE) leaf-optics wrapper around
:func:`scopeinpython.fluspect.get_fluspect_cx_scope`.

Direct port of ``SCOPEinR::get.fluspect_mSCOPE`` (``SCOPEinR/R/fluspect_mSCOPE.R``).
Leaf optical properties (reflectance, transmittance, fluorescence
excitation-emission matrices ``Mb``/``Mf``, and pigment contribution
factors) are computed once per distinct leaf-biochemistry profile layer
(``mly.nly`` layers), then replicated across the ``nl`` canopy layers each
profile layer spans (weighted by ``mly.pLAI``).

**R quirk reproduced exactly, not "fixed"**: at each profile-layer
boundary, the canopy sublayer shared between two consecutive profile
layers (``indStar[i+1]`` in R / ``indStar[i]`` here) is assigned twice --
once by the earlier profile layer's replication, once by the later one --
and the later assignment wins. This mirrors R's own ``rho_temp[in1:in2,]
<- ...`` overlapping-range assignment (the R source's own comment notes an
earlier version of this bug was *more* severe, losing entire profile
layers rather than just double-writing one boundary row per layer).

**A real, confirmed crash in the R source, not offered here**: if
``step`` is omitted in R, ``get.fluspect_mSCOPE`` sets its *internal*
``step_to_model`` to 5 but pre-allocates the ``Mb``/``Mf`` output array
using ``spectral$wlE``/``spectral$wlF`` directly (which are always fixed
1 nm grids, 351 x 211, from ``get.spectra.SCOPE``) -- not the 53x71 shape
``getFluspect.Cx.SCOPE(..., step=5)`` actually returns. The resulting
``leafopt$Mb[,,in1:in2] <- array(...)`` assignment then fails with R's
"number of items to replace is not a multiple of replacement length"
(confirmed via a standalone repro), i.e. calling ``get.fluspect_mSCOPE``
without ``step`` always crashes. Since this code path can never succeed,
``step`` is a required (not optional-with-a-buggy-default) parameter here.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .fluspect import _scope_fluspect_optipar, get_fluspect_cx_scope
from .spectral import SpectralConfig

__all__ = ["MultiLayerLeafBio", "FluspectMScopeResult", "fluspect_mscope"]


@dataclass
class MultiLayerLeafBio:
    """Multi-layer (mSCOPE) leaf-biochemistry profile. One value per
    distinct biochemistry layer (``nly`` of them); ``pLAI`` is that
    layer's share of total LAI (need not sum to 1 -- normalised
    internally, matching R)."""

    nly: int
    pLAI: np.ndarray
    pCab: np.ndarray
    pEWT: np.ndarray
    pCar: np.ndarray
    pLMA: np.ndarray
    pCs: np.ndarray
    pN: np.ndarray


@dataclass
class FluspectMScopeResult:
    refl: np.ndarray  # (nl, nwlP)
    tran: np.ndarray  # (nl, nwlP)
    kChlrel: np.ndarray  # (nl, nwlP)
    kCarrel: np.ndarray  # (nl, nwlP)
    Mb: np.ndarray  # (nwlf, nwle, nl)
    Mf: np.ndarray  # (nwlf, nwle, nl)
    phiI: np.ndarray  # (nwlF,)
    phiII: np.ndarray  # (nwlF,)


def fluspect_mscope(
    mly: MultiLayerLeafBio,
    spectral: SpectralConfig,
    nl: int,
    Cx: float,
    fqe: float,
    Prot: float,
    CBC: float,
    Anth: float,
    step: float = 5.0,
) -> FluspectMScopeResult:
    """Direct port of ``SCOPEinR::get.fluspect_mSCOPE``. ``Cx``/``fqe``/
    ``Prot``/``CBC``/``Anth`` are the baseline leaf properties shared by
    every profile layer (only ``Cab``/``EWT``/``Car``/``LMA``/``Cs``/``N``
    vary per layer, taken from ``mly``). The R function also accepts a
    ``leafopt``-unrelated ``optipar`` parameter and an ``soil``/plotting
    path -- both dropped here: ``optipar`` is documented but never
    actually used in the R source body (it always calls
    ``getFluspect.Cx.SCOPE`` with the hardcoded ``optipar2021.Pro.CX``),
    and the ``get.plots`` branch only produces diagnostic plots, never
    modifying the returned ``leafopt``.
    """
    nwlP = len(spectral.wlP)
    refl = np.empty((nl, nwlP))
    tran = np.empty((nl, nwlP))
    kChlrel = np.empty((nl, nwlP))
    kCarrel = np.empty((nl, nwlP))
    Mb = None
    Mf = None

    pLAI = np.asarray(mly.pLAI, dtype=float)
    frac = pLAI / pLAI.sum()
    indStar = np.concatenate(([1.0], np.floor(np.cumsum(frac) * nl)))

    for i in range(1, mly.nly + 1):
        j = i - 1
        res = get_fluspect_cx_scope(
            Cab=mly.pCab[j], Car=mly.pCar[j], EWT=mly.pEWT[j],
            LMA=mly.pLMA[j], Cs=mly.pCs[j], N=mly.pN[j],
            fqe=fqe, Cx=Cx, Prot=Prot, CBC=CBC, Anth=Anth, step=step,
        )

        if Mb is None:
            Mb = np.empty((res.Mb.shape[0], res.Mb.shape[1], nl))
            Mf = np.empty_like(Mb)

        in1 = int(indStar[i - 1])
        in2 = int(indStar[i])
        sl = slice(in1 - 1, in2)

        refl[sl, :] = res.refl
        tran[sl, :] = res.tran
        kChlrel[sl, :] = res.kChlrel
        kCarrel[sl, :] = res.kCarrel
        Mb[:, :, sl] = res.Mb[:, :, None]
        Mf[:, :, sl] = res.Mf[:, :, None]

    op = _scope_fluspect_optipar()
    iw_coincidents = (spectral.wlF - 400).astype(int)  # wlF is always a subset of wlP's 400-2400 grid
    phiI = op.phiI[iw_coincidents]
    phiII = op.phiII[iw_coincidents]

    return FluspectMScopeResult(
        refl=refl, tran=tran, kChlrel=kChlrel, kCarrel=kCarrel,
        Mb=Mb, Mf=Mf, phiI=phiI, phiII=phiII,
    )
