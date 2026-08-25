"""Convolve reflectance onto a sensor using a plain per-band SRF table --
i.e. no SMAC atmospheric-correction coefficients involved at all, just a
spectral response function.

Direct port of ``ToolsRTM::get.spectral.convolution.srf()``
(``ToolsRTM/R/get.spectral.convolution.srf.R``). That R function itself
generalizes what used to be two separate, app-only helper functions
(``convolve_prisma()``/``convolve_smac_sensor()``) in the AEO-Course
PROSAIL Shiny app's own app.R into one shared, package-level
implementation -- this module is the Python side of that same
generalization.

``smac.py``'s :func:`toolsrtm.smac.spectral_convolution` needs a sensor with
SMAC atmospheric-correction coefficients bundled (only Sentinel-2A ships
those in this Python port so far). PRISMA has no SMAC coefficients at all;
Sentinel-2A/B additionally ship a second, plain, publisher-original SRF
table alongside their SMAC bundle. :func:`spectral_convolution_srf` covers
all three from one shared implementation.
"""
from __future__ import annotations

import functools
from dataclasses import dataclass
from importlib import resources

import numpy as np

__all__ = [
    "SrfTable",
    "SrfConvolutionResult",
    "srf_prisma",
    "fwhm_prisma",
    "srf_sentinel2a",
    "srf_sentinel2b",
    "spectral_convolution_srf",
    "enmap_characteristics",
    "sensor_characteristics",
    "spectral_convolution_gaussian",
]


@dataclass
class SrfTable:
    """A plain per-band SRF table: one weight column per sensor band, all
    sampled on the same wavelength grid (the ``wl`` field)."""

    wl: np.ndarray  # (nwl,)
    band_names: list[str]
    weights: np.ndarray  # (nwl, nbands)


@dataclass
class SrfConvolutionResult:
    """One row per SRF band: the ``wl`` field is the SRF-weighted mean
    center wavelength, ``fwhm`` the full width at half maximum, ``rfl`` the
    convolved reflectance -- rows sorted by center wavelength."""

    band_names: list[str]
    wl: np.ndarray
    fwhm: np.ndarray
    rfl: np.ndarray


def _load_srf_csv(filename: str, wl_col: int = 0) -> tuple[np.ndarray, list[str], np.ndarray]:
    with resources.files("toolsrtm.data").joinpath(filename).open("r", encoding="utf-8") as f:
        header = f.readline().strip().split(",")
        data = np.genfromtxt(f, delimiter=",")
    band_names = [h.strip('"') for h in header[1:]]
    wl = data[:, wl_col]
    weights = data[:, 1:]
    return wl, band_names, weights


@functools.lru_cache(maxsize=None)
def srf_prisma() -> SrfTable:
    """PRISMA's 234 hyperspectral bands. Direct export of ``ToolsRTM::srf.prisma``."""
    wl, band_names, weights = _load_srf_csv("srf_prisma.csv")
    return SrfTable(wl=wl, band_names=band_names, weights=weights)


@functools.lru_cache(maxsize=None)
def fwhm_prisma() -> tuple[np.ndarray, np.ndarray]:
    """PRISMA's officially calibrated per-band FWHM (nm), positionally
    aligned with :func:`srf_prisma`'s 234 bands -- more precise than the
    half-max-crossing estimate :func:`spectral_convolution_srf` would
    otherwise derive from the SRF weight profile itself. Direct export of
    ``ToolsRTM::fwhm.prisma``.

    Returns
    -------
    wl : ndarray, shape (234,)
        Nominal band center wavelength (nm).
    fwhm : ndarray, shape (234,)
        Full width at half maximum (nm).
    """
    with resources.files("toolsrtm.data").joinpath("fwhm_prisma.csv").open("r", encoding="utf-8") as f:
        f.readline()
        data = np.genfromtxt(f, delimiter=",")
    return data[:, 0], data[:, 1]


@functools.lru_cache(maxsize=None)
def srf_sentinel2a() -> SrfTable:
    """Sentinel-2A MSI's 13 bands, plain publisher-original SRF (distinct
    from :func:`toolsrtm.smac.sentinel2a_msi`'s SMAC-bundled copy of the
    same physical curve -- see this module's docstring). Direct export of
    ``ToolsRTM::srf.sentinel2a``."""
    wl, band_names, weights = _load_srf_csv("srf_sentinel2a.csv")
    return SrfTable(wl=wl, band_names=band_names, weights=weights)


@functools.lru_cache(maxsize=None)
def srf_sentinel2b() -> SrfTable:
    """Sentinel-2B MSI's 13 bands -- a real, slightly different sensor from
    2A (up to ~17nm band-center difference in B12), not a duplicate. Direct
    export of ``ToolsRTM::srf.sentinel2b``."""
    wl, band_names, weights = _load_srf_csv("srf_sentinel2b.csv")
    return SrfTable(wl=wl, band_names=band_names, weights=weights)


def spectral_convolution_srf(
    wave: np.ndarray,
    values: np.ndarray,
    srf: SrfTable,
    fwhm: tuple[np.ndarray, np.ndarray] | None = None,
) -> SrfConvolutionResult:
    """Weighted-average a high-resolution spectrum onto ``srf``'s bands.

    Parameters
    ----------
    wave : array_like, shape (nwl,)
        Integer-nm wavelength grid ``values`` is defined on.
    values : array_like, shape (nwl,)
        Reflectance (or any other high-resolution spectrum) on ``wave``.
    srf : SrfTable
        e.g. :func:`srf_prisma`, :func:`srf_sentinel2a`, :func:`srf_sentinel2b`.
    fwhm : (wl, fwhm) tuple, optional
        Bundled precise FWHM, positionally aligned with ``srf``'s bands
        (e.g. :func:`fwhm_prisma`'s return value) -- when omitted, FWHM is
        estimated directly from ``srf``'s own sampled weight profile
        (coarser; this is the only option for Sentinel-2A/B, which have no
        separately bundled FWHM table).
    """
    wave = np.asarray(wave)
    values = np.asarray(values, dtype=float)
    nbands = srf.weights.shape[1]
    use_bundled_fwhm = fwhm is not None and len(fwhm[1]) == nbands

    wl_out = np.full(nbands, np.nan)
    fwhm_out = np.full(nbands, np.nan)
    rfl_out = np.full(nbands, np.nan)

    for b in range(nbands):
        p_all = srf.weights[:, b]
        valid = ~np.isnan(p_all) & (p_all > 0)
        if not np.any(valid):
            continue
        wl_v = srf.wl[valid]
        p_v = p_all[valid]

        wl_out[b] = np.sum(wl_v * p_v) / np.sum(p_v)

        if use_bundled_fwhm:
            fwhm_out[b] = fwhm[1][b]
        else:
            half = p_v.max() / 2
            above = wl_v[p_v >= half]
            fwhm_out[b] = above.max() - above.min() if len(above) >= 2 else np.nan

        idx = np.searchsorted(wave, np.round(wl_v))
        in_range = (idx >= 0) & (idx < len(wave)) & (wave[np.clip(idx, 0, len(wave) - 1)] == np.round(wl_v))
        w = p_v[in_range]
        v = values[idx[in_range]]
        if w.sum() == 0:
            continue
        rfl_out[b] = float(np.sum(w * v) / np.sum(w))

    order = np.argsort(wl_out)
    names_out = [srf.band_names[i] for i in order]
    return SrfConvolutionResult(band_names=names_out, wl=wl_out[order], fwhm=fwhm_out[order], rfl=rfl_out[order])


# --- Gaussian convolution from nominal band characteristics (no measured SRF) ---
#
# :func:`spectral_convolution_srf` above needs a real, measured, per-nm SRF
# table -- most sensors don't have one published/available at all. Often all
# that's known (or all a student has, e.g. from an instrument's own ENVI
# header, or their own camera's calibration sheet) is a list of band center
# wavelengths and FWHM. Direct port of
# ``ToolsRTM::get.spectral.convolution.gaussian()``
# (``ToolsRTM/R/get.spectral.convolution.gaussian.R``) -- see that
# function's own docstring for the full rationale.

_SENSOR_CHARACTERISTICS_SENSORS = (
    "ALI", "Hyperion", "Landsat4", "Landsat5", "Landsat7", "Landsat8",
    "MODIS", "Quickbird", "RapidEye", "Sentinel2a", "Sentinel2b",
    "WorldView2-4", "WorldView2-8",
)


@functools.lru_cache(maxsize=None)
def enmap_characteristics() -> tuple[np.ndarray, np.ndarray]:
    """EnMAP's 242 hyperspectral channels: nominal (center, fwhm) in nm,
    already given directly (no band-edge derivation needed). Direct export
    of ``ToolsRTM::EnMap.characteristics``."""
    with resources.files("toolsrtm.data").joinpath("enmap_characteristics.csv").open("r", encoding="utf-8") as f:
        f.readline()
        rows = [line.strip().split(",") for line in f if line.strip()]
    center = np.array([float(r[2]) for r in rows])
    fwhm = np.array([float(r[3]) for r in rows])
    return center, fwhm


@functools.lru_cache(maxsize=None)
def _sensor_characteristics_table() -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    with resources.files("toolsrtm.data").joinpath("sensor_characteristics.csv").open("r", encoding="utf-8") as f:
        f.readline()
        rows = [line.strip().split(",") for line in f if line.strip()]
    table: dict[str, list[list[float]]] = {}
    for r in rows:
        name = r[0].strip('"')
        table.setdefault(name, []).append([float(r[2]), float(r[3]), float(r[4])])  # lb, ub, average
    return {
        name: (np.array([v[2] for v in vals]), np.array([v[0] for v in vals]), np.array([v[1] for v in vals]))
        for name, vals in table.items()
    }


def sensor_characteristics(sensor: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Nominal (center, lb, ub) band characteristics for one bundled sensor
    -- published band edges, not a measured SRF. Direct export of
    ``ToolsRTM::sensor.characteristics``.

    Parameters
    ----------
    sensor : one of ``"ALI"``, ``"Hyperion"``, ``"Landsat4"``,
        ``"Landsat5"``, ``"Landsat7"``, ``"Landsat8"``, ``"MODIS"``,
        ``"Quickbird"``, ``"RapidEye"``, ``"Sentinel2a"``, ``"Sentinel2b"``,
        ``"WorldView2-4"``, ``"WorldView2-8"``. (Sentinel-2A/B and PRISMA
        also have a REAL measured SRF table bundled -- prefer
        :func:`srf_sentinel2a`/:func:`srf_sentinel2b`/:func:`srf_prisma`
        with :func:`spectral_convolution_srf` for those, more accurate than
        this Gaussian approximation.)

    Returns
    -------
    center, lb, ub : ndarray
        Band center wavelength and lower/upper edge, nm.
    """
    table = _sensor_characteristics_table()
    if sensor not in table:
        raise ValueError(
            f'Unknown sensor "{sensor}". Bundled options: {", ".join(_SENSOR_CHARACTERISTICS_SENSORS)}. '
            "For any other sensor (including your own camera), pass `centers` (and optionally `fwhm`) "
            "directly to spectral_convolution_gaussian() instead of `sensor`."
        )
    return table[sensor]


def spectral_convolution_gaussian(
    wave: np.ndarray,
    values: np.ndarray,
    sensor: str | None = None,
    centers: np.ndarray | None = None,
    fwhm: np.ndarray | None = None,
) -> SrfConvolutionResult:
    """Convolve onto a sensor using only nominal band characteristics
    (center + FWHM, approximated as a Gaussian response -- optionally
    truncated to a published band-edge range), when no real measured SRF
    table is available at all.

    Three ways to call this:

    1. ``sensor="EnMAP"`` -- :func:`enmap_characteristics`'s 242 channels
       (center + FWHM already given).
    2. ``sensor="MODIS"`` (or any other name in
       :func:`sensor_characteristics`'s docstring) -- these ship published
       band EDGES, not FWHM directly; FWHM is derived as ``ub - lb`` and
       the Gaussian response is additionally hard-truncated to
       ``[lb, ub]``.
    3. Your OWN sensor/camera: pass `centers` yourself, in nm (e.g. copied
       straight out of an ENVI header's ``wavelength = {...}`` block).
       `fwhm` is optional -- if omitted, each band's width is approximated
       from its distance to its neighboring bands (the standard assumption
       for a CONTIGUOUS pushbroom imaging spectrometer, e.g. a Headwall
       camera, where the true per-band SRF calibration isn't available).
       Pass `fwhm` explicitly (e.g. from a camera datasheet, or an ENVI
       header's own ``fwhm = {...}`` block) for a more accurate result.

    Parameters
    ----------
    wave, values : array_like
        High-resolution wavelength grid (nm) and spectrum on it (e.g. a
        simulated reflectance spectrum).
    sensor : str, optional
        A bundled sensor name (see above). If given, `centers`/`fwhm` are
        looked up automatically and any values also passed for them are
        ignored.
    centers : array_like, optional (required if `sensor` is not given)
        Your own sensor's band center wavelengths, nm.
    fwhm : array_like, optional
        Your own sensor's per-band FWHM, nm, same length/order as
        `centers`. Derived from band spacing if omitted.
    """
    lb = ub = None

    if sensor is not None:
        if sensor == "EnMAP":
            centers, fwhm = enmap_characteristics()
        else:
            center, lb, ub = sensor_characteristics(sensor)
            centers = center
            fwhm = ub - lb

    if centers is None:
        raise ValueError("Supply either `sensor` (a bundled sensor name) or your own `centers`.")

    centers = np.asarray(centers, dtype=float)
    order = np.argsort(centers)
    centers = centers[order]
    if fwhm is not None:
        fwhm = np.asarray(fwhm, dtype=float)[order]
    if lb is not None:
        lb = np.asarray(lb, dtype=float)[order]
        ub = np.asarray(ub, dtype=float)[order]

    if fwhm is None:
        d = np.diff(centers)
        fwhm = np.concatenate(([d[0]], (d[:-1] + d[1:]) / 2, [d[-1]]))

    wave = np.asarray(wave)
    values = np.asarray(values, dtype=float)
    nbands = len(centers)
    rfl_out = np.full(nbands, np.nan)

    for i in range(nbands):
        sigma = fwhm[i] / (2 * np.sqrt(2 * np.log(2)))
        w = np.exp(-0.5 * ((wave - centers[i]) / sigma) ** 2)
        if lb is not None:
            w = np.where((wave < lb[i]) | (wave > ub[i]), 0.0, w)
        if w.sum() == 0:
            continue
        rfl_out[i] = float(np.sum(w * values) / np.sum(w))

    band_names = [f"band{i + 1}" for i in range(nbands)]
    return SrfConvolutionResult(band_names=band_names, wl=centers, fwhm=fwhm, rfl=rfl_out)
