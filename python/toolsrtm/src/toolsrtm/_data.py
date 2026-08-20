"""Loaders for the bundled spectral parameter tables (ported from ToolsRTM/data/*.rda).

The original R package ships these as .rda binary objects (``dataSpec_PDB``,
``dataSpec_PRO``). They were exported once to CSV from the original R
package (see ``python/scratch/scratch_export.R`` in the repo root) and are bundled
here as package data so this port has no runtime dependency on R.

The R code in ``prospect_DB``/``prospect_PRO`` accesses these tables
*positionally* (``data[, 1]``, ``data[, 2]``, ...), so we do the same here
rather than relying on column names (which differ slightly between the two
CSV exports).
"""
from __future__ import annotations

import functools
from dataclasses import dataclass
from importlib import resources

import numpy as np


@functools.lru_cache(maxsize=None)
def _load_csv(filename: str) -> np.ndarray:
    with resources.files("toolsrtm.data").joinpath(filename).open("r", encoding="utf-8") as f:
        f.readline()  # header
        data = np.loadtxt(f, delimiter=",")
    return data


@dataclass
class ProspectDBTable:
    """Columns of dataSpec_PDB.rda: 400-2500 nm (1 nm step).

    The last four columns (``direct_light``/``diffuse_light``/``dry_soil``/
    ``wet_soil``) aren't used by PROSPECT-D itself -- they're bundled here
    because ``ToolsRTM::Compute_BRF`` (used by ``inform``) indexes this same
    table positionally for its Es/Ed illumination spectra.
    """

    lambda_: np.ndarray
    nr: np.ndarray
    Kab: np.ndarray
    Kcar: np.ndarray
    Kant: np.ndarray
    KBrown: np.ndarray
    Kw: np.ndarray
    Km: np.ndarray
    direct_light: np.ndarray
    diffuse_light: np.ndarray
    dry_soil: np.ndarray
    wet_soil: np.ndarray


@dataclass
class ProspectPROTable:
    """Columns of dataSpec_PRO.rda: 400-2500 nm (1 nm step)."""

    lambda_: np.ndarray
    nr: np.ndarray
    Kab: np.ndarray
    Kcar: np.ndarray
    Kant: np.ndarray
    KBrown: np.ndarray
    Kw: np.ndarray
    Km: np.ndarray
    Kprot: np.ndarray
    Knonprot: np.ndarray


@functools.lru_cache(maxsize=None)
def data_spec_pdb() -> ProspectDBTable:
    d = _load_csv("dataSpec_PDB.csv")
    return ProspectDBTable(*(d[:, i] for i in range(12)))


@functools.lru_cache(maxsize=None)
def data_spec_pro() -> ProspectPROTable:
    d = _load_csv("dataSpec_PRO.csv")
    return ProspectPROTable(*(d[:, i] for i in range(10)))
