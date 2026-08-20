"""Leaf optical property models: calctav, PROSPECT-D, PROSPECT-PRO.

Direct, function-by-function port of:
  - ToolsRTM/R/calctav.R
  - ToolsRTM/R/prospect_DB.R   (PROSPECT-D, exposed here as ``prospect_d``)
  - ToolsRTM/R/prospect_PRO.R  (PROSPECT-PRO, exposed here as ``prospect_pro``)

References
----------
Feret J-B, Gitelson AA, Noble SD & Jacquemoud S, 2017. PROSPECT-D: Towards
modeling leaf optical properties through a complete lifecycle. Remote
Sensing of Environment, 193, 204-215.

Feret, J.B., Berger, K., de Boissieu, F., Malenovsky, Z., 2021. PROSPECT-PRO
for estimating content of nitrogen-containing leaf proteins and other
carbon-based constituents. Remote Sens. Environ. 252.

Stern F. (1964); Allen W.A. (1973): transmissivity of a dielectric surface
(calctav).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import exp1

from ._data import data_spec_pdb, data_spec_pro

__all__ = ["calctav", "LeafOptics", "prospect_d", "prospect_pro"]


def calctav(alpha: float, nr: np.ndarray) -> np.ndarray:
    """Transmissivity of a dielectric plane surface, averaged over all
    directions of incidence and polarizations.

    Direct port of ``ToolsRTM::calctav``.

    Parameters
    ----------
    alpha : float
        Maximum incidence angle (degrees) defining the solid angle of
        incident light.
    nr : array_like
        Refractive index spectrum (or scalar).

    Returns
    -------
    numpy.ndarray
        Transmissivity of a dielectric plane surface, same shape as ``nr``.

    References
    ----------
    Stern F. (1964), Transmission of isotropic radiation across an
    interface between two dielectrics, Appl. Opt., 3(1):111-113.
    Allen W.A. (1973), Transmission of isotropic light across a dielectric
    surface in two and three dimensions, J. Opt. Soc. Am., 63(6):664-666.
    """
    nr = np.asarray(nr, dtype=float)
    rd = np.pi / 180.0
    n2 = nr**2
    np_ = n2 + 1
    nm = n2 - 1
    a = (nr + 1) * (nr + 1) / 2
    k = -(n2 - 1) * (n2 - 1) / 4
    sa = np.sin(alpha * rd)
    b2 = (sa**2) - (np_ / 2)

    if alpha == 90:
        b1 = 0 * b2
    else:
        b1 = np.sqrt((b2**2) + k)

    b = b1 - b2
    b3 = b**3
    a3 = a**3

    ts = (k**2 / (6 * b3) + k / b - b / 2) - (k**2 / (6 * a3) + (k / a) - (a / 2))

    tp1 = -2 * n2 * (b - a) / (np_**2)
    tp2 = -2 * n2 * np_ * np.log(b / a) / (nm**2)
    tp3 = n2 * ((1 / b) - (1 / a)) / 2
    tp4 = (
        16
        * n2**2
        * ((n2**2) + 1)
        * np.log(((2 * np_ * b) - (nm**2)) / ((2 * np_ * a) - (nm**2)))
        / ((np_**3) * (nm**2))
    )
    tp5 = 16 * (n2**3) * (1 / ((2 * np_ * b) - (nm**2)) - (1 / (2 * np_ * a - (nm**2)))) / (np_**3)
    tp = tp1 + tp2 + tp3 + tp4 + tp5
    tav = (ts + tp) / (2 * (sa**2))

    return tav


@dataclass
class LeafOptics:
    """Output of a leaf RT model: wavelength, reflectance and transmittance."""

    lambda_: np.ndarray
    refl: np.ndarray
    tran: np.ndarray


def _stokes_n_layers(r: np.ndarray, t: np.ndarray, N: float):
    """Stokes equations: reflectance/transmittance of N (real-valued) layers
    from the single-layer r, t and the pigment transmittance tau.
    Shared by prospect_d/prospect_pro (both R functions duplicate this block).
    """
    D = np.sqrt((1 + r + t) * (1 + r - t) * (1 - r + t) * (1 - r - t))
    rq = r**2
    tq = t**2
    a = (1 + rq - tq + D) / (2 * r)
    b = (1 - rq + tq + D) / (2 * t)

    bNm1 = b ** (N - 1)
    bN2 = bNm1**2
    a2 = a**2
    denom = a2 * bN2 - 1
    Rsub = a * (bN2 - 1) / denom
    Tsub = bNm1 * (a2 - 1) / denom

    # case of zero absorption (r + t >= 1)
    j = r + t >= 1
    Tsub = np.where(j, t / (t + (1 - t) * (N - 1)), Tsub)
    Rsub = np.where(j, 1 - Tsub, Rsub)

    return Rsub, Tsub


def _one_layer(Kall: np.ndarray, nr: np.ndarray, alpha: float):
    """Non-conservative scattering + reflectance/transmittance of a single
    layer (Allen et al. 1969), shared by prospect_d/prospect_pro."""
    j = Kall > 0
    t1 = (1 - Kall) * np.exp(-Kall)
    t2 = (Kall**2) * exp1(np.where(j, Kall, 1.0))  # exp1 undefined at 0; masked by `j` below
    tau = np.ones_like(Kall)
    tau = np.where(j, t1 + t2, tau)

    talf = calctav(alpha, nr)
    ralf = 1 - talf
    t12 = calctav(90, nr)
    r12 = 1 - t12
    t21 = t12 / (nr**2)
    r21 = 1 - t21

    denom = 1 - (r21 * r21 * (tau**2))
    Ta = (talf * tau * t21) / denom
    Ra = ralf + (r21 * tau * Ta)

    t = t12 * tau * t21 / denom
    r = r12 + (r21 * tau * t)
    return Ra, Ta, r, t


def prospect_d(
    N: float,
    Cab: float,
    Car: float,
    Anth: float,
    Cbrown: float,
    EWT: float,
    LMA: float,
    alpha: float = 40.0,
) -> LeafOptics:
    """PROSPECT-D leaf optical properties model (400-2500 nm, 1 nm step).

    Direct port of ``ToolsRTM::prospect_DB`` (also called PROSPECT-Dynamic
    with brown pigments in the R source).

    Parameters
    ----------
    N : float
        Leaf structure parameter.
    Cab : float
        Chlorophyll a+b content (microg/cm2).
    Car : float
        Carotenoid content (microg/cm2).
    Anth : float
        Anthocyanin content (microg/cm2).
    Cbrown : float
        Brown pigment content (arbitrary units).
    EWT : float
        Equivalent water thickness (g/cm2).
    LMA : float
        Leaf mass per area (g/cm2).
    alpha : float, default 40
        Maximum incidence angle for calctav.

    Returns
    -------
    LeafOptics
        ``lambda_`` (400-2500 nm), refl, tran.
    """
    d = data_spec_pdb()
    Kall = (Cab * d.Kab + Car * d.Kcar + Anth * d.Kant + Cbrown * d.KBrown + EWT * d.Kw + LMA * d.Km) / N

    Ra, Ta, r, t = _one_layer(Kall, d.nr, alpha)
    Rsub, Tsub = _stokes_n_layers(r, t, N)

    denom = 1 - Rsub * r
    tran = Ta * Tsub / denom
    refl = Ra + (Ta * Rsub * t) / denom

    return LeafOptics(lambda_=d.lambda_.copy(), refl=refl, tran=tran)


def prospect_pro(
    N: float,
    Cab: float,
    Car: float,
    Anth: float,
    Cbrown: float,
    EWT: float,
    LMA: float,
    alpha: float,
    Prot: float,
    CBC: float,
) -> LeafOptics:
    """PROSPECT-PRO leaf optical properties model (400-2500 nm, 1 nm step).

    Direct port of ``ToolsRTM::prospect_PRO``.

    Parameters
    ----------
    N, Cab, Car, Anth, Cbrown, EWT, LMA, alpha : see :func:`prospect_d`.
        Note: ``Anth`` is expressed in nmol/cm2 for PROSPECT-PRO (vs
        microg/cm2 for PROSPECT-D), consistent with the R implementation.
    Prot : float
        Protein content (g/cm2).
    CBC : float
        Non-protein carbon-based constituent content (g/cm2).

    Notes
    -----
    Mirrors the R behaviour: if ``LMA > 0`` and (``Prot > 0`` or
    ``CBC > 0``), LMA is forced to 0 (LMA = Prot + CBC), with a message.

    Returns
    -------
    LeafOptics
        ``lambda_`` (400-2500 nm), refl, tran.
    """
    d = data_spec_pro()

    if LMA > 0 and (Prot > 0 or CBC > 0):
        LMA = 0.0

    Kall = (
        Cab * d.Kab
        + Car * d.Kcar
        + Anth * d.Kant
        + Cbrown * d.KBrown
        + EWT * d.Kw
        + LMA * d.Km
        + Prot * d.Kprot
        + CBC * d.Knonprot
    ) / N

    Ra, Ta, r, t = _one_layer(Kall, d.nr, alpha)
    Rsub, Tsub = _stokes_n_layers(r, t, N)

    denom = 1 - Rsub * r
    tran = Ta * Tsub / denom
    refl = Ra + (Ta * Rsub * t) / denom

    return LeafOptics(lambda_=d.lambda_.copy(), refl=refl, tran=tran)
