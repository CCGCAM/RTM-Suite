"""INFORM: the INvertible FOrest Reflectance Model (Atzberger 2000;
Schlerf & Atzberger 2006), a FLIM (Rosema et al. 1992) + fourSAIL + PROSPECT
combination for forest-stand bidirectional reflectance, 400-2500 nm.

Direct port of ``ToolsRTM/R/inform.R`` and its internal helpers
(``foursail.inform.R``, ``foursail.inf.R``, ``foursail_t_s.R``,
``foursail_t_o.R``, ``Compute_BRF.R``), restricted to what ``inform()``
itself actually exercises. All 5 leaf models ToolsRTM itself supports
(PROSPECT-D/-PRO, Liberty, Fluspect-B/-B-Cx) are wired in via
:func:`toolsrtm.canopy._leaf_optics`, the same dispatcher :func:`toolsrtm.canopy.foursail`
uses.

Reproduces the R source's own behavior exactly, so results match R
bit-for-bit -- see inline comments at ``_INFORM_QUIRKY_LIDF`` and
``compute_brf``.
"""
from __future__ import annotations

from typing import Literal

import numpy as np

from ._data import data_spec_pdb
from .canopy import _foursail_scattering_core, _leaf_optics, _LeafModel, campbell, foursail_core
from .leaf import prospect_d

__all__ = ["compute_brf", "foursail_inform", "foursail_inf", "foursail_t_s", "foursail_t_o", "inform"]


# ---------------------------------------------------------------------------
# LIDF quirk shared by foursail.inform / foursail_t_s / foursail_t_o
# ---------------------------------------------------------------------------
#
# For TypeLidf == 1, the R source of these three functions (but NOT
# foursail.inf, which calls dladgen() normally) has the proper
# ``dladgen(LIDFa, LIDFb)`` call commented out and replaced with a
# hardcoded 65-value literal -- of which only the first 13 are ever read,
# because the weighted-sum loop runs over ``length(litab) == 13``. This is
# reproduced verbatim (values copied from the R source) rather than fixed,
# so TypeLidf==1 output matches R exactly.
_INFORM_QUIRKY_LIDF = np.array([
    0.015192247821401716, 0.04511513125626976, 0.07366721933874043,
    0.09998095795183792, 0.12325683701156054, 0.14278760558390557,
    0.1579798610972386, 0.16837196070089022, 0.034475081609994906,
    0.034644632694447175, 0.03477199446646273, 0.03485697201080851,
    0.03489949845644191,
])
_tx1 = np.array([10, 20, 30, 40, 50, 60, 70, 80, 82, 84, 86, 88, 90], dtype=float)
_tx2 = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 82, 84, 86, 88], dtype=float)
_INFORM_QUIRKY_LITAB = (_tx1 + _tx2) / 2.0


def _inform_lidf(LIDFa: float, TypeLidf: Literal[1, 2]):
    """LIDF source used by foursail.inform/foursail_t_s/foursail_t_o (the
    quirky hardcoded table for TypeLidf==1, real campbell() for TypeLidf==2
    -- LIDFb is not used by either branch here, matching the R source)."""
    if TypeLidf == 1:
        return _INFORM_QUIRKY_LIDF, _INFORM_QUIRKY_LITAB
    elif TypeLidf == 2:
        ld = campbell(LIDFa)
        return ld.lidf, ld.litab
    raise ValueError("TypeLidf must be 1 or 2")


def compute_brf(rdot, rsot, tts: float, es=None, ed=None, short_waves: bool = False):
    """Bidirectional reflectance factor from direct/diffuse light mixing.

    Direct port of ``ToolsRTM::Compute_BRF``. That R function has two Es/Ed
    sourcing paths gated on ``is.null(data.light)``: a positional
    ``dataSpec_PDB[,11]``/``[,12]`` path (columns actually named
    ``dry_soil``/``wet_soil`` -- unused dead code, since no call site in
    the package ever omits ``data.light``) and a named ``data.light$
    direct_light``/``data.light$diffuse_light`` path, which is the one
    every real call site (``foursail.inform``, ``foursail.inf``) hits by
    always passing ``data.light=ToolsRTM::dataSpec_PDB`` explicitly. ``es``/
    ``ed`` here default to those ``direct_light``/``diffuse_light`` columns
    to match.

    ``es``/``ed`` are truncated to ``len(rdot)`` directly rather than being
    gated on ``short_waves`` -- ``rdot``/``rsot`` are 2101 long for
    PROSPECT-D/-PRO/Liberty but only 2001 long for Fluspect-B/-B-Cx (its
    native 400-2400nm domain); a caller who passes mismatched-length
    ``rdot``/``short_waves`` combinations (e.g. calling this directly on a
    Fluspect-based ``foursail()`` result without also setting
    ``short_waves=True``) would otherwise hit a hard NumPy broadcasting
    error. ``ToolsRTM::Compute_BRF`` matches Es/Ed to ``rdot``'s length the
    same way this port does; ``short_waves`` is kept only for backward API
    compatibility and no longer has any effect.
    """
    rdot = np.asarray(rdot, dtype=float)
    rsot = np.asarray(rsot, dtype=float)
    rd = np.pi / 180.0

    if es is None or ed is None:
        d = data_spec_pdb()
        es = d.direct_light if es is None else es
        ed = d.diffuse_light if ed is None else ed
    es = np.asarray(es, dtype=float)[: len(rdot)]
    ed = np.asarray(ed, dtype=float)[: len(rdot)]

    sin90tts = np.sin((90.0 - tts) * rd)
    skyl = 0.847 - 1.61 * sin90tts + 1.04 * sin90tts * sin90tts
    PARdiro = (1 - skyl) * es
    PARdifo = skyl * ed
    return (rdot * PARdifo + rsot * PARdiro) / (PARdiro + PARdifo)


def foursail_inform(inputLUT: dict, rsoil: np.ndarray) -> np.ndarray:
    """Understorey reflectance (dicotyledoneae defaults), BRF-blended.

    Direct port of ``ToolsRTM::foursail.inform`` with ``typeLAI='understorey'``
    -- the only ``typeLAI`` value ``inform()`` itself ever calls it with (the
    ``'Inf-crownTree'``/default branches are dead code from the caller's
    perspective and are not ported).
    """
    N, Cab, Car, Anth, Cbrown = 2.0, 30.0, 0.0, 0.0, 0.0
    EWT, LMA = 0.05, 0.05
    alpha = inputLUT["alpha"]
    LIDFa_fixed = 45.0
    TypeLidf = inputLUT["TypeLidf"]
    lai = inputLUT["LAIu"]
    hotspot = inputLUT["hspot"]
    tts, tto, psi = inputLUT["tts"], inputLUT["tto"], inputLUT["psi"]

    leaf = prospect_d(N, Cab, Car, Anth, Cbrown, EWT, LMA, alpha)
    rho, tau = leaf.refl, leaf.tran

    lidf, litab = _inform_lidf(LIDFa_fixed, TypeLidf)
    sail = _foursail_scattering_core(rho, tau, rsoil, lidf, litab, lai, hotspot, tts, tto, psi)
    return compute_brf(sail.rdot, sail.rsot, tts)


def foursail_inf(
    inputLUT: dict, rsoil: np.ndarray, rleaf: np.ndarray, tleaf: np.ndarray, short_waves: bool = False,
) -> np.ndarray:
    """Infinite (very dense, LAI=15) crown reflectance, BRF-blended.

    Direct port of ``ToolsRTM::foursail.inf`` as called by ``inform()``.
    ``short_waves`` defaults to ``False`` (the non-Fluspect leaf models'
    call site omits the argument, so R's ``missing()`` guard resolves it
    to ``FALSE`` -- see :func:`compute_brf`); ``inform()`` passes
    ``short_waves=True`` explicitly for the Fluspect-B/-B-Cx leaf models
    (matching R's own ``foursail.inf(..., short.waves=T)`` call there).
    Unlike :func:`foursail_inform`, TypeLidf==1 here uses a real
    ``dladgen`` call in the R source (not the hardcoded quirky table), so
    this reuses :func:`toolsrtm.canopy.foursail_core` directly.
    """
    lai, hotspot = 15.0, 0.04
    LIDFa, LIDFb, TypeLidf = inputLUT["LIDFa"], inputLUT["LIDFb"], inputLUT["TypeLidf"]
    tts, tto, psi = inputLUT["tts"], inputLUT["tto"], inputLUT["psi"]

    sail = foursail_core(rleaf, tleaf, rsoil, LIDFa, LIDFb, TypeLidf, lai, hotspot, tts, tto, psi)
    return compute_brf(sail.rdot, sail.rsot, tts, short_waves=short_waves)


def _crown_transmittance(inputLUT: dict, rsoil: np.ndarray, rleaf: np.ndarray, tleaf: np.ndarray, tts: float):
    """Shared body of foursail_t_s/foursail_t_o (differ only in which angle
    is passed as ``tts``, and in not calling Compute_BRF at all -- they
    blend rdot/rsot directly with the LUT's ``skyl`` fraction)."""
    lai = inputLUT["LAI"]
    hotspot = 0.0
    LIDFa, TypeLidf = inputLUT["LIDFa"], inputLUT["TypeLidf"]
    tto, psi = inputLUT["tto"], inputLUT["psi"]
    skyl = inputLUT["skyl"]

    lidf, litab = _inform_lidf(LIDFa, TypeLidf)
    sail = _foursail_scattering_core(rleaf, tleaf, rsoil, lidf, litab, lai, hotspot, tts, tto, psi)

    PARdifo = skyl
    PARdiro = 1 - skyl
    return (sail.rdot * PARdiro + sail.rsot * PARdifo) / (PARdiro + PARdifo)


def foursail_t_s(inputLUT: dict, rsoil: np.ndarray, rleaf: np.ndarray, tleaf: np.ndarray) -> np.ndarray:
    """Crown transmittance in the sun direction. Port of ``foursail_t_s.R``."""
    return _crown_transmittance(inputLUT, rsoil, rleaf, tleaf, tts=inputLUT["tts"])


def foursail_t_o(inputLUT: dict, rsoil: np.ndarray, rleaf: np.ndarray, tleaf: np.ndarray) -> np.ndarray:
    """Crown transmittance in the observation direction. Port of
    ``foursail_t_o.R`` -- note the R source passes the *viewing* zenith
    angle (``tto``) in as ``tts`` too, so illumination and viewing
    geometry coincide; reproduced as-is."""
    return _crown_transmittance(inputLUT, rsoil, rleaf, tleaf, tts=inputLUT["tto"])


def inform(
    inputLUT: dict,
    rsoil: np.ndarray,
    leaf_model: _LeafModel = "PROSPECT-PRO",
) -> np.ndarray:
    """INFORM forest-stand reflectance. Direct port of ``ToolsRTM::inform``,
    dispatching to any of the 5 leaf models ToolsRTM itself supports (see
    :func:`toolsrtm.canopy._leaf_optics`) -- same as :func:`toolsrtm.canopy.foursail`.

    Parameters
    ----------
    inputLUT : dict
        Must contain the fourSAIL keys (``LIDFa``, ``LIDFb``, ``TypeLidf``,
        ``LAI``, ``hspot``, ``tts``, ``tto``, ``psi``) plus the leaf-model
        parameters (see :func:`toolsrtm.canopy.foursail` for the full list
        per leaf model), plus the INFORM-specific canopy keys ``LAIu``
        (understorey LAI), ``sd`` (stem density, ha-1), ``cd`` (crown
        diameter, m), ``h`` (tree height, m), and ``skyl`` (diffuse-light
        fraction).
    rsoil : array_like
        Soil reflectance spectrum, 400-2500 nm (2101 values).
    leaf_model : {'PROSPECT-D', 'PROSPECT-PRO', 'Liberty', 'Fluspect-B', 'Fluspect-B-Cx'}

    Returns
    -------
    numpy.ndarray
        Forest reflectance spectrum -- 2101 values (400-2500 nm) for
        PROSPECT-D/-PRO/Liberty, 2001 values (400-2400 nm) for the two
        Fluspect leaf models (matching R's ``r_understorey[1:2001]``
        truncation there and ``foursail.inf(..., short.waves=T)``).
    """
    tts, tto, psi = inputLUT["tts"], inputLUT["tto"], inputLUT["psi"]
    sd, cd, h = inputLUT["sd"], inputLUT["cd"], inputLUT["h"]

    refl_soil = np.asarray(rsoil, dtype=float)

    r_understorey = foursail_inform(inputLUT, refl_soil)

    r_leaf, t_leaf, force_2001 = _leaf_optics(inputLUT, leaf_model)
    if force_2001:
        r_understorey = r_understorey[:2001]

    r_sail_inf = foursail_inf(inputLUT, r_understorey, r_leaf, t_leaf, short_waves=force_2001)

    # Ground coverage (FLIM model, Rosema et al. 1992)
    adapt = 1.0
    k = adapt * (np.pi * (cd / 2.0) ** 2) / 10000.0
    tto_ = tto * np.pi / 180.0
    tts_ = tts * np.pi / 180.0
    psi_ = psi * np.pi / 180.0

    co = 1 - np.exp(-k * sd / np.cos(tto_))
    cs = 1 - np.exp(-k * sd / np.cos(tts_))
    g = ((np.tan(tto_) ** 2 + np.tan(tts_)) ** 2 - 2 * np.tan(tto_) * np.tan(tts_) * np.cos(psi_)) ** 0.5
    p = np.exp(-g * h / cd)

    Fcd = co * cs + p * (co * (1 - co) * cs * (1 - cs)) ** 0.5
    Fcs = co * (1 - cs) - p * (co * (1 - co) * cs * (1 - cs)) ** 0.5
    Fod = (1 - co) * cs - p * (co * (1 - co) * cs * (1 - cs)) ** 0.5
    Fos = (1 - co) * (1 - cs) + p * (co * (1 - co) * cs * (1 - cs)) ** 0.5

    t_s = foursail_t_s(inputLUT, r_understorey, r_leaf, t_leaf)
    t_o = foursail_t_o(inputLUT, r_understorey, r_leaf, t_leaf)

    G = Fcd * t_s * t_o + Fcs * t_o + Fod * t_s + Fos
    C = Fcd * (1 - t_s * t_o)

    return (r_sail_inf * C) + (r_understorey * G)
