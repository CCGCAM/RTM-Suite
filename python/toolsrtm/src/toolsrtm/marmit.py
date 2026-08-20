"""MARMIT-1 and MARMIT-2 soil reflectance models (ported from
ToolsRTM/R/marmit1.R, marmit2.R and the get.marmit.rsoil() wrapper added to
ToolsRTM this session).

Only the Bablet_2016 soil database (17 IDs) is bundled -- same deliberate
scoping as the R side's ``get.marmit.rsoil()`` (see
``ToolsRTM/R/get.marmit.rsoil.R``): only the driest spectrum per soil ID is
needed (MARMIT computes wet reflectance FROM that one dry reference), so
that's all that's bundled here, exported via
``python/scratch/scratch_export_marmit.py``.
"""
from __future__ import annotations

import functools
from dataclasses import dataclass
from importlib import resources
from typing import Literal

import numpy as np
from scipy.special import exp1


@dataclass
class MarmitSoil:
    """Result of :func:`get_marmit_rsoil`."""

    wavelength: np.ndarray
    rsoil_dry: np.ndarray
    rsoil_wet: np.ndarray
    smc: float


def marmit1(n: np.ndarray, alpha: np.ndarray, rd: np.ndarray, L: float, eps: float) -> np.ndarray:
    """Wet-soil reflectance from a dry reference spectrum (MARMIT-1, Bablet et al. 2018).

    Parameters
    ----------
    n : spectral optical index of water (real refractive index).
    alpha : water absorption spectral coefficient, cm^-1.
    rd : reflectance of the dry soil reference.
    L : thickness of the surface water layer, cm.
    eps : fraction of the soil surface that is wet (0-1).

    Returns
    -------
    Wet soil reflectance, same shape as ``n``/``alpha``/``rd``.
    """
    n = np.asarray(n, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    rd = np.asarray(rd, dtype=float)

    r12_diffuse = (
        (3 * n**2 + 2 * n + 1) / (3 * (n + 1) ** 2)
        - 2 * n**3 * (n**2 + 2 * n - 1) / ((n**2 + 1) ** 2 * (n**2 - 1))
        + n**2 * (n**2 + 1) * np.log(n) / (n**2 - 1) ** 2
        - n**2 * (n**2 - 1) ** 2 * np.log(n * (n + 1) / (n - 1)) / (n**2 + 1) ** 3
    )
    t12_diffuse = 1 - r12_diffuse
    r21_diffuse = 1 - (1 - r12_diffuse) / n**2
    t21_diffuse = 1 - r21_diffuse

    if L > 0:
        # exp1(x) is the exponential integral E1(x) = integral_x^inf e^-t/t dt,
        # matching R's numerically-integrated exp1_base() exactly.
        tw_diffuse = (1 - alpha * L) * np.exp(-alpha * L) + (alpha * L) ** 2 * exp1(alpha * L)
    else:
        tw_diffuse = np.ones_like(n)

    rw = (t12_diffuse * t21_diffuse * rd * tw_diffuse**2) / (1 - r21_diffuse * rd * tw_diffuse**2)
    return eps * rw + (1 - eps) * rd


def marmit2(
    n_w: np.ndarray, alpha_w: np.ndarray, n_i: float, k_i: float,
    rd: np.ndarray, L: float, eps: float, d_i: float, wls: np.ndarray,
) -> np.ndarray:
    """Wet-soil reflectance from a dry reference spectrum (MARMIT-2,
    accounts for soil particle size/refractive index -- generally more
    accurate than :func:`marmit1` for coarser soils).

    Direct port of ``ToolsRTM::get.marmit2``. Differs from :func:`marmit1`
    in two ways: the effective medium's refractive index (``n``/``k``,
    hence the water-layer transmittance ``tw_diffuse``) is a dielectric
    mixture of water (``n_w``/``alpha_w``) and soil particles (``n_i``/
    ``k_i``, weighted by the particle volume fraction ``d_i``), not pure
    water optics; and the wet/dry mixing uses a power-law (Hapke-like)
    rule with exponent 1/2.27 instead of :func:`marmit1`'s linear mixing.

    Parameters
    ----------
    n_w : spectral optical index of water (real refractive index).
    alpha_w : water absorption spectral coefficient, cm^-1.
    n_i : real part of the soil particles' refractive index.
    k_i : imaginary part of the soil particles' refractive index.
    rd : reflectance of the dry soil reference.
    L : thickness of the surface water layer, cm.
    eps : fraction of the soil surface that is wet (0-1).
    d_i : particle volume fraction of the water/particle mixture.
    wls : wavelengths, nm (needed for the water absorption -> imaginary
        refractive index conversion, unlike :func:`marmit1`).

    Returns
    -------
    Wet soil reflectance, same shape as ``n_w``/``alpha_w``/``rd``/``wls``.
    """
    n_w = np.asarray(n_w, dtype=float)
    alpha_w = np.asarray(alpha_w, dtype=float)
    rd = np.asarray(rd, dtype=float)
    wls = np.asarray(wls, dtype=float)

    k_w = alpha_w * wls * 1e-7 / (4 * np.pi)  # imaginary part of water's refractive index
    e_w = (n_w + 1j * k_w) ** 2  # complex permittivity of water
    e_i = (n_i + 1j * k_i) ** 2  # complex permittivity of soil particles
    e = d_i * e_i + (1 - d_i) * e_w  # dielectric average
    n = np.sqrt(e).real  # effective refractive index of the mixture
    k = np.sqrt(e).imag
    alpha = 4 * np.pi * k / (wls * 1e-7)  # effective absorption coefficient

    r12_diffuse = (
        (3 * n**2 + 2 * n + 1) / (3 * (n + 1) ** 2)
        - 2 * n**3 * (n**2 + 2 * n - 1) / ((n**2 + 1) ** 2 * (n**2 - 1))
        + n**2 * (n**2 + 1) * np.log(n) / (n**2 - 1) ** 2
        - n**2 * (n**2 - 1) ** 2 * np.log(n * (n + 1) / (n - 1)) / (n**2 + 1) ** 3
    )
    t12_diffuse = 1 - r12_diffuse
    r21_diffuse = 1 - (1 - r12_diffuse) / n**2
    t21_diffuse = 1 - r21_diffuse

    if L > 0:
        tw_diffuse = (1 - alpha * L) * np.exp(-alpha * L) + (alpha * L) ** 2 * exp1(alpha * L)
    else:
        tw_diffuse = np.ones_like(n)

    rw = t12_diffuse * t21_diffuse * rd * tw_diffuse**2 / (1 - r21_diffuse * rd * tw_diffuse**2)
    return (eps * rw ** (1 / 2.27) + (1 - eps) * rd ** (1 / 2.27)) ** 2.27


def sigmoid_soil(phi: float, k: float, a: float, psi: float) -> float:
    """Soil moisture content (gravimetric %) from the wetness parameter phi = L*eps."""
    return k / (1 + a * np.exp(-psi * phi))


@functools.lru_cache(maxsize=None)
def _load_water_optics() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with resources.files("toolsrtm.data.marmit").joinpath("water_optics.csv").open("r") as f:
        f.readline()
        data = np.loadtxt(f, delimiter=",")
    return data[:, 0], data[:, 1], data[:, 2]  # wl, n, alpha


@functools.lru_cache(maxsize=None)
def _load_bablet_index() -> dict[int, dict]:
    with resources.files("toolsrtm.data.marmit").joinpath("bablet_2016_index.csv").open("r") as f:
        f.readline()
        rows = [line.strip().split(",") for line in f if line.strip()]
    return {int(r[0]): {"name": r[1], "K": float(r[2]), "a": float(r[3]), "psi": float(r[4])} for r in rows}


@functools.lru_cache(maxsize=None)
def _load_bablet_spectra() -> dict[int, tuple[np.ndarray, np.ndarray]]:
    with resources.files("toolsrtm.data.marmit").joinpath("bablet_2016_spectra.csv").open("r") as f:
        f.readline()
        rows = [line.strip().split(",") for line in f if line.strip()]
    by_id: dict[int, list[tuple[float, float]]] = {}
    for sid, wl, r in rows:
        by_id.setdefault(int(sid), []).append((float(wl), float(r)))
    return {sid: (np.array([p[0] for p in pts]), np.array([p[1] for p in pts])) for sid, pts in by_id.items()}


def _load_external_index(db_root: str, database: str):
    """Read ``<db_root>/<database>/<database>.csv`` -- any MARMIT database in
    the full RTM-Suite layout (e.g. from the repo's own ``databases/``
    folder, all 8 databases, ~200MB total, not bundled with the package).
    Extra columns beyond ``ID, Refl_file, SMCg, K, a, psi`` are ignored,
    matching ``ToolsRTM::get.marmit.rsoil()``'s R-side reader.
    """
    import csv
    from pathlib import Path

    root = Path(db_root)
    db_dir = root / database
    index_path = db_dir / f"{database}.csv"
    if not index_path.is_file():
        available = sorted(p.name for p in root.iterdir() if p.is_dir()) if root.is_dir() else []
        raise ValueError(
            f"get_marmit_rsoil(): soil database {database!r} not found under {db_root!r}. "
            f"Available: {', '.join(available)}."
        )
    with open(index_path, "r", encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))
    by_id: dict[int, list[dict]] = {}
    for r in rows:
        by_id.setdefault(int(r["ID"]), []).append(r)
    return by_id, db_dir


def _load_external_spectrum(db_dir, refl_file: str) -> tuple[np.ndarray, np.ndarray]:
    from pathlib import Path

    path = Path(db_dir) / "spectra" / refl_file
    data = np.loadtxt(path, delimiter="\t", skiprows=1)
    return data[:, 0], data[:, 1]


def get_marmit_rsoil(
    soil_id: int = 1,
    L: float = 0.05,
    eps: float = 0.3,
    version: Literal["marmit1", "marmit2"] = "marmit1",
    n_i: float = 1.53,
    k_i: float = 0.001,
    d_i: float = 0.0005,
    wl_out: np.ndarray | None = None,
    database: str = "Bablet_2016",
    db_root: str | None = None,
) -> MarmitSoil:
    """Build a canopy-model-ready soil reflectance spectrum from MARMIT-1
    or MARMIT-2.

    Python port of ``ToolsRTM::get.marmit.rsoil()``. Only the Bablet_2016
    soil database (17 IDs) is bundled with the package, to keep install
    size small. The other 7 MARMIT databases (Dupiau 2020, Humper 2015,
    Lesaignoux 2008, Liu 2002, Lobell 2002, Marcq 2012, Philpot 2014 -- see
    https://pss-gitlab.math.univ-paris-diderot.fr/marmit/marmit) ship in the
    RTM-Suite monorepo's own ``databases/`` folder (repo root, ~200MB
    total, not bundled here either). Point at it directly with ``db_root``,
    e.g. ``get_marmit_rsoil(database="Liu_2002", db_root="databases")`` run
    from the repo root -- no copying required. Any other folder with the
    same layout (an index CSV ``<name>/<name>.csv`` with columns ``ID,
    Refl_file, SMCg, K, a, psi`` -- extra columns ignored -- plus
    ``<name>/spectra/<Refl_file>`` tab-separated ``Wvl,R`` files) works the
    same way. See that R function's docstring for the physical background.

    Parameters
    ----------
    soil_id : soil ID within the database's index (the ``ID`` column). For
        Bablet_2016, 1-17.
    L : thickness of the surface water layer, cm.
    eps : fraction of the soil surface that is wet (0-1).
    version : {'marmit1', 'marmit2'}. MARMIT-2 additionally accounts for
        soil particle size/refractive index (``n_i``/``k_i``/``d_i``) and
        is generally more accurate for coarser soils; MARMIT-1 is simpler
        and matches the original 2018 paper.
    n_i, k_i, d_i : MARMIT-2-only soil-particle parameters (real
        refractive index, imaginary refractive index, particle volume
        fraction). Ignored when ``version='marmit1'``. Defaults match the
        MARMIT Shiny app's defaults.
    wl_out : wavelength grid (nm) to resample/pad onto. Defaults to
        ``np.arange(400, 2501)`` (400-2500nm, 1nm step), matching
        :func:`toolsrtm.canopy.foursail`'s default 2101-point grid.
    database : soil database name. Default ``"Bablet_2016"``, the only one
        bundled with the package (ignored -- always Bablet_2016 -- unless
        ``db_root`` is given).
    db_root : directory containing database subfolders (e.g. ``"databases"``
        at the RTM-Suite repo root, which has all 8 MARMIT databases -- see
        above). When ``None`` (default), only the bundled Bablet_2016
        database is available and ``database`` is ignored.
    """
    if version not in ("marmit1", "marmit2"):
        raise ValueError(f"version must be 'marmit1' or 'marmit2', got {version!r}")

    if db_root is None:
        index = _load_bablet_index()
        if soil_id not in index:
            raise ValueError(f"soil_id {soil_id} not in Bablet_2016 (available: {sorted(index)})")
        meta = index[soil_id]

        spectra = _load_bablet_spectra()
        wl_raw, rd_raw = spectra[soil_id]
    else:
        by_id, db_dir = _load_external_index(db_root, database)
        if soil_id not in by_id:
            raise ValueError(f"soil_id {soil_id} not in {database!r} (available: {sorted(by_id)})")
        rows = by_id[soil_id]
        dry_row = min(rows, key=lambda r: float(r["SMCg"]))
        meta = {"K": float(dry_row["K"]), "a": float(dry_row["a"]), "psi": float(dry_row["psi"])}
        wl_raw, rd_raw = _load_external_spectrum(db_dir, dry_row["Refl_file"])

    wl_native = np.arange(max(wl_raw.min(), 400), wl_raw.max() + 1)
    rd = np.interp(wl_native, wl_raw, rd_raw)

    wl_w, n_w_raw, alpha_w_raw = _load_water_optics()
    n_w = np.interp(wl_native, wl_w, n_w_raw)
    alpha_w = np.interp(wl_native, wl_w, alpha_w_raw)

    if version == "marmit2":
        rw = marmit2(n_w, alpha_w, n_i, k_i, rd, L, eps, d_i, wl_native)
    else:
        rw = marmit1(n_w, alpha_w, rd, L, eps)
    phi = L * eps
    smc = sigmoid_soil(phi, meta["K"], meta["a"], meta["psi"])

    if wl_out is None:
        wl_out = np.arange(400, 2501)
    wl_out = np.asarray(wl_out, dtype=float)
    # np.interp clamps out-of-range x to the boundary y value by default,
    # matching R's approx(..., rule = 2) used on the R side.
    rsoil_wet = np.interp(wl_out, wl_native, rw)
    rsoil_dry = np.interp(wl_out, wl_native, rd)

    return MarmitSoil(wavelength=wl_out, rsoil_dry=rsoil_dry, rsoil_wet=rsoil_wet, smc=float(smc))
