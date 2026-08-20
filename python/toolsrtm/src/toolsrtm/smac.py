"""SMAC (Simplified Method for Atmospheric Correction, Rahman & Dedieu 1994)
atmospheric radiative transfer, and the sensor coefficient/spectral-response
data it needs -- the piece that turns SPART's top-of-canopy (TOC) BRDF into
top-of-atmosphere (TOA) reflectance/radiance for a specific sensor.

Direct port of ``ToolsRTM::get.smac``/``get.coef.SMAC``/``get.spectral.convolution``
(``ToolsRTM/R/get.smac.R``, ``get.coef.SMAC.R``, ``Spectral.convolution.R``).

**Scope**: all 9 sensors the R package ships are bundled -- Landsat 4/5/7/8,
Sentinel-2A/B, Sentinel-3A/B, Terra/Aqua MODIS (see ``ToolsRTM/data/*.rda``
and :data:`SENSORS`/:func:`get_sensor`). The atmospheric-correction
*physics* (:func:`get_smac`, :func:`spectral_convolution`) is fully general
and sensor-agnostic given a :class:`SmacSensor`; each sensor's coefficient/SRF
tables were exported via ``python/scratch/scratch_export_smac_sensors.R``.
Sentinel-2A verified against a real, unmodified ``ToolsRTM::SPART()`` call
(see :func:`toolsrtm.spart.spart_toa`); the other 8 verified against a real
``ToolsRTM::get.smac()`` call each (see ``tests/test_smac.py``).
"""
from __future__ import annotations

import functools
from dataclasses import dataclass
from importlib import resources

import numpy as np

__all__ = [
    "SmacSensor", "SENSORS", "get_sensor",
    "sentinel2a_msi", "sentinel2b_msi", "sentinel3a_olci", "sentinel3b_olci",
    "landsat4_tm", "landsat5_tm", "landsat7_etm", "landsat8_oli", "terra_aqua_modis",
    "SmacAtmosphere", "get_smac", "spectral_convolution",
]

_COEF_NAMES = (
    "ah2o", "nh2o", "ao3", "no3", "ao2", "no2", "po2", "aco2", "nco2", "pco2",
    "ach4", "nch4", "pch4", "ano2", "nno2", "pno2", "aco", "nco", "pco",
    "a0s", "a1s", "a2s", "a3s", "a0T", "a1T", "a2T", "a3T", "taur", "sr",
    "a0taup", "a1taup", "wo", "gc", "a0P", "a1P", "a2P", "a3P", "a4P",
    "Rest1", "Rest2", "Rest3", "Rest4", "Resr1", "Resr2", "Resr3",
    "Resa1", "Resa2", "Resa3", "Resa4",
)


@dataclass
class SmacSensor:
    """One sensor's SMAC coefficients + spectral-response data (subset of
    R's ``sensor`` list, e.g. ``ToolsRTM::Sentinel2A.MSI``, actually used
    by :func:`get_smac`/:func:`spectral_convolution`)."""

    mission: str
    wl_smac: np.ndarray  # band center wavelengths SMAC's coefficients are defined at, (nbands,)
    coef: dict  # one (nbands,) array per name in _COEF_NAMES
    wl_srf: np.ndarray  # (nsrf, nbands), NaN-padded -- SRF sample wavelengths per band
    p_srf: np.ndarray  # (nsrf, nbands), NaN-padded -- SRF weights per band


def _load_sensor(slug: str, mission: str) -> SmacSensor:
    """Shared loader for every bundled sensor's 4 CSV files
    (``smac_bands_<slug>.csv``, ``smac_coef_<slug>.csv``,
    ``smac_srf_wl_<slug>.csv``, ``smac_srf_weight_<slug>.csv``)."""
    with resources.files("toolsrtm.data").joinpath(f"smac_bands_{slug}.csv").open("r", encoding="utf-8") as f:
        f.readline()
        bands = np.atleast_1d(np.loadtxt(f, delimiter=","))
        wl_smac = bands[:, 0] if bands.ndim == 2 else bands

    coef = {}
    with resources.files("toolsrtm.data").joinpath(f"smac_coef_{slug}.csv").open("r", encoding="utf-8") as f:
        f.readline()
        for line in f:
            parts = line.strip().split(",")
            name = parts[0].strip('"')
            coef[name] = np.array([float(v) for v in parts[1:]])

    with resources.files("toolsrtm.data").joinpath(f"smac_srf_wl_{slug}.csv").open("r", encoding="utf-8") as f:
        f.readline()
        wl_srf = np.genfromtxt(f, delimiter=",")
    with resources.files("toolsrtm.data").joinpath(f"smac_srf_weight_{slug}.csv").open("r", encoding="utf-8") as f:
        f.readline()
        p_srf = np.genfromtxt(f, delimiter=",")

    return SmacSensor(mission=mission, wl_smac=wl_smac, coef=coef, wl_srf=wl_srf, p_srf=p_srf)


@functools.lru_cache(maxsize=None)
def sentinel2a_msi() -> SmacSensor:
    """Sentinel-2A MSI: 13 bands. Direct export of ``ToolsRTM::Sentinel2A.MSI``."""
    return _load_sensor("sentinel2a", "Sentinel-2A")


@functools.lru_cache(maxsize=None)
def sentinel2b_msi() -> SmacSensor:
    """Sentinel-2B MSI: 13 bands. Direct export of ``ToolsRTM::Sentinel2B.MSI``."""
    return _load_sensor("sentinel2b", "Sentinel-2B")


@functools.lru_cache(maxsize=None)
def sentinel3a_olci() -> SmacSensor:
    """Sentinel-3A OLCI: 21 bands. Direct export of ``ToolsRTM::Sentinel3A.OLCI``."""
    return _load_sensor("sentinel3a", "Sentinel-3A")


@functools.lru_cache(maxsize=None)
def sentinel3b_olci() -> SmacSensor:
    """Sentinel-3B OLCI: 21 bands. Direct export of ``ToolsRTM::Sentinel3B.OLCI``."""
    return _load_sensor("sentinel3b", "Sentinel-3B")


@functools.lru_cache(maxsize=None)
def landsat4_tm() -> SmacSensor:
    """Landsat 4 TM: 6 bands. Direct export of ``ToolsRTM::LANDSAT4.TM``."""
    return _load_sensor("landsat4", "Landsat-4")


@functools.lru_cache(maxsize=None)
def landsat5_tm() -> SmacSensor:
    """Landsat 5 TM: 6 bands. Direct export of ``ToolsRTM::LANDSAT5.TM``."""
    return _load_sensor("landsat5", "Landsat-5")


@functools.lru_cache(maxsize=None)
def landsat7_etm() -> SmacSensor:
    """Landsat 7 ETM+: 6 bands. Direct export of ``ToolsRTM::LANDSAT7.ETM``."""
    return _load_sensor("landsat7", "Landsat-7")


@functools.lru_cache(maxsize=None)
def landsat8_oli() -> SmacSensor:
    """Landsat 8 OLI: 9 bands. Direct export of ``ToolsRTM::LANDSAT8.OLI``."""
    return _load_sensor("landsat8", "Landsat-8")


@functools.lru_cache(maxsize=None)
def terra_aqua_modis() -> SmacSensor:
    """Terra/Aqua MODIS: 20 bands. Direct export of ``ToolsRTM::TerraAqua.MODIS``."""
    return _load_sensor("modis", "Terra/Aqua")


#: Every bundled sensor, by the same short name used in R's ``ToolsRTM::get.smac(sensor=...)``.
SENSORS = {
    "Sentinel2A.MSI": sentinel2a_msi,
    "Sentinel2B.MSI": sentinel2b_msi,
    "Sentinel3A.OLCI": sentinel3a_olci,
    "Sentinel3B.OLCI": sentinel3b_olci,
    "LANDSAT4.TM": landsat4_tm,
    "LANDSAT5.TM": landsat5_tm,
    "LANDSAT7.ETM": landsat7_etm,
    "LANDSAT8.OLI": landsat8_oli,
    "TerraAqua.MODIS": terra_aqua_modis,
}


def get_sensor(name: str) -> SmacSensor:
    """Look up a bundled sensor by name (see :data:`SENSORS` for the exact
    keys, matching R's own sensor object names)."""
    try:
        return SENSORS[name]()
    except KeyError:
        raise ValueError(f"Unknown sensor {name!r}. Choose from {list(SENSORS)}.") from None


@dataclass
class SmacAtmosphere:
    """Per-band atmospheric optical quantities, on ``sensor.wl_smac``."""

    Ta_ss: np.ndarray  # directional transmittance, direct incidence (downward)
    Ta_sd: np.ndarray  # hemispherical transmittance, direct incidence (downward)
    Ta_oo: np.ndarray  # directional transmittance, viewing direction (upward)
    Ta_do: np.ndarray  # hemispherical transmittance, viewing direction (upward)
    Ta_s: np.ndarray  # directional transmittance, diffuse light (downward)
    Ta_o: np.ndarray  # hemispherical transmittance, diffuse light (upward)
    Tg: np.ndarray  # total gaseous scattering transmission
    Ra_dd: np.ndarray  # hemispherical atmospheric reflectance, diffuse light
    Ra_so: np.ndarray  # directional atmospheric reflectance, direct incidence


def get_smac(
    sensor: SmacSensor,
    tts: float, tto: float, psi: float,
    Pa: float, taup550: float, uo3: float, uh2o: float,
) -> SmacAtmosphere:
    """Atmospheric transmittance/reflectance terms (SMAC, Rahman & Dedieu
    1994), per sensor band. Direct port of ``ToolsRTM::get.smac``.

    Parameters
    ----------
    sensor : SmacSensor
    tts, tto, psi : float
        Sun zenith, view zenith, relative azimuth (degrees).
    Pa : float
        Surface air pressure (hPa). Use ``ToolsRTM::get.Altitude2Pa``'s
        formula yourself first if you only have altitude -- not ported
        here since every call site in this port supplies ``Pa`` directly.
    taup550 : float
        Aerosol optical thickness at 550nm.
    uo3 : float
        Ozone content (atm-cm).
    uh2o : float
        Water vapour content (g/cm2).
    """
    c = sensor.coef

    def g(name: str) -> np.ndarray:
        return c[name]

    cdr = np.pi / 180.0
    crd = 180.0 / np.pi

    us = np.cos(tts * cdr)
    uv = np.cos(tto * cdr)
    Peq = Pa / 1013.25

    m = 1 / us + 1 / uv

    taup = g("a0taup") + g("a1taup") * taup550

    uo2 = Peq ** g("po2")
    uco2 = Peq ** g("pco2")
    uch4 = Peq ** g("pch4")
    uno2 = Peq ** g("pno2")
    uco = Peq ** g("pco")

    to3 = np.exp(g("ao3") * (uo3 * m) ** g("no3"))
    th2o = np.exp(g("ah2o") * (uh2o * m) ** g("nh2o"))
    to2 = np.exp(g("ao2") * (uo2 * m) ** g("no2"))
    tco2 = np.exp(g("aco2") * (uco2 * m) ** g("nco2"))
    tch4 = np.exp(g("ach4") * (uch4 * m) ** g("nch4"))
    tno2 = np.exp(g("ano2") * (uno2 * m) ** g("nno2"))
    tco = np.exp(g("aco") * (uco * m) ** g("nco"))

    tg = th2o * to3 * to2 * tco2 * tch4 * tco * tno2

    s = g("a0s") * Peq + g("a3s") + g("a1s") * taup550 + g("a2s") * taup550**2

    ttetas = g("a0T") + g("a1T") * taup550 / us + (g("a2T") * Peq + g("a3T")) / (1 + us)
    ttetav = g("a0T") + g("a1T") * taup550 / uv + (g("a2T") * Peq + g("a3T")) / (1 + uv)

    cksi = -(us * uv + np.sqrt(1 - us**2) * np.sqrt(1 - uv**2) * np.cos(psi * cdr))
    cksi = max(cksi, -1.0)
    ksiD = crd * np.arccos(cksi)

    ray_phase = 0.7190443 * (1 + cksi**2) + 0.0412742
    ray_ref = (g("taur") * ray_phase) / (4 * us * uv)
    ray_ref = ray_ref * Pa / 1013.25
    taurz = g("taur") * Peq

    aer_phase = g("a0P") + g("a1P") * ksiD + g("a2P") * ksiD**2 + g("a3P") * ksiD**3 + g("a4P") * ksiD**4
    wo, gc = g("wo"), g("gc")
    ak2 = (1 - wo) * (3 - wo * 3 * gc)
    ak = np.sqrt(ak2)

    e = -3 * us**2 * wo / (4 * (1 - ak2 * us**2))
    f = -(1 - wo) * 3 * gc * us**2 * wo / (4 * (1 - ak2 * us**2))
    dp = e / (3 * us) + us * f
    d = e + f
    b = 2 * ak / (3 - wo * 3 * gc)
    delta = np.exp(ak * taup) * (1 + b) ** 2 - np.exp(-ak * taup) * (1 - b) ** 2
    ww = wo / 4
    ss = us / (1 - ak2 * us**2)
    q1 = 2 + 3 * us + (1 - wo) * 3 * gc * us * (1 + 2 * us)
    q2 = 2 - 3 * us - (1 - wo) * 3 * gc * us * (1 - 2 * us)
    q3 = q2 * np.exp(-taup / us)
    c1 = ((ww * ss) / delta) * (q1 * np.exp(ak * taup) * (1 + b) + q3 * (1 - b))
    c2 = -((ww * ss) / delta) * (q1 * np.exp(-ak * taup) * (1 - b) + q3 * (1 + b))
    cp1 = c1 * ak / (3 - wo * 3 * gc)
    cp2 = -c2 * ak / (3 - wo * 3 * gc)
    z = d - wo * 3 * gc * uv * dp + wo * aer_phase / 4
    x = c1 - wo * 3 * gc * uv * cp1
    y = c2 - wo * 3 * gc * uv * cp2
    aa1 = uv / (1 + ak * uv)
    aa2 = uv / (1 - ak * uv)
    aa3 = us * uv / (us + uv)

    aer_ref1 = x * aa1 * (1 - np.exp(-taup / aa1))
    aer_ref2 = y * aa2 * (1 - np.exp(-taup / aa2))
    aer_ref3 = z * aa3 * (1 - np.exp(-taup / aa3))
    aer_ref = (aer_ref1 + aer_ref2 + aer_ref3) / (us * uv)

    Res_ray = g("Resr1") + g("Resr2") * g("taur") * ray_phase / (us * uv) \
        + g("Resr3") * (g("taur") * ray_phase / (us * uv)) ** 2

    Res_aer = (g("Resa1") + g("Resa2") * (taup * m * cksi) + g("Resa3") * (taup * m * cksi) ** 2) \
        + g("Resa4") * (taup * m * cksi) ** 3

    tautot = taup + taurz
    Res_6s = (g("Rest1") + g("Rest2") * (tautot * m * cksi) + g("Rest3") * (tautot * m * cksi) ** 2) \
        + g("Rest4") * (tautot * m * cksi) ** 3

    atm_ref = ray_ref - Res_ray + aer_ref - Res_aer + Res_6s

    tdir_tts = np.exp(-tautot / us)
    tdir_ttv = np.exp(-tautot / uv)
    tdif_tts = ttetas - tdir_tts
    tdif_ttv = ttetav - tdir_ttv

    return SmacAtmosphere(
        Ta_ss=tdir_tts, Ta_sd=tdif_tts, Ta_oo=tdir_ttv, Ta_do=tdif_ttv,
        Ta_s=ttetas, Ta_o=ttetav, Tg=tg, Ra_dd=s, Ra_so=atm_ref,
    )


def spectral_convolution(wave: np.ndarray, values: np.ndarray, sensor: SmacSensor) -> np.ndarray:
    """Weighted-average a high-resolution spectrum onto a sensor's bands
    using its spectral response function (SRF). Direct port of
    ``ToolsRTM::get.spectral.convolution`` (the non-Sentinel-3/MODIS
    per-detector-averaging branch -- Sentinel-2A has one SRF per band
    already, no ``colMeans`` step needed).

    Parameters
    ----------
    wave : array_like, shape (nwl,)
        Integer-nm wavelength grid ``values`` is defined on (must cover
        every SRF sample wavelength for a band to get a non-NaN result).
    values : array_like, shape (nwl,)
    """
    wave = np.asarray(wave)
    values = np.asarray(values, dtype=float)
    nbands = sensor.wl_srf.shape[1]
    out = np.full(nbands, np.nan)
    for b in range(nbands):
        wl_b = sensor.wl_srf[:, b]
        p_b = sensor.p_srf[:, b]
        valid = ~np.isnan(wl_b) & ~np.isnan(p_b)
        if not np.any(valid):
            continue
        idx = np.searchsorted(wave, wl_b[valid])
        in_range = (idx >= 0) & (idx < len(wave)) & (wave[np.clip(idx, 0, len(wave) - 1)] == wl_b[valid])
        w = p_b[valid][in_range]
        v = values[idx[in_range]]
        if w.sum() == 0:
            continue
        out[b] = float(np.sum(w * v) / np.sum(w))
    return out
