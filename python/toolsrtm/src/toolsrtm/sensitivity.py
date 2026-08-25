"""Global sensitivity analysis and correlated/multi-distribution LUT builders.

Python port of ToolsRTM's sensitivity-analysis toolkit: ``Correlated_value.R``
(:func:`correlated_value`), ``Gaussian_MinMax.R`` (:func:`gauss_by_min_max`),
``get_distributionLUT.R`` (:func:`get_distribution_lut`), ``getCor.R`
(:func:`get_cor`), ``get.sobol.indices.R`` (:func:`sobol_indices`,
via :func:`johnson_relative_weights`) and ``get.spectral.sensitivity.R``
(:func:`spectral_sensitivity`).

The Johnson relative-weights index (Johnson, 2000, "A heuristic method for
estimating the relative weight of predictor variables in multiple
regression") is implemented directly here via eigendecomposition rather than
calling R's ``sensitivity`` package -- verified to reproduce
``sensitivity::johnson()``'s output to 8 decimal places on a fixed reference
dataset (see ``tests/test_sensitivity.py``). :func:`sobol_indices`'s ``Si``/
``STi`` columns are a simplified two-sample-split estimator (matching the R
function's own docstring: "treat it as a rough indicator rather than a
precise total-effect index") -- :func:`spectral_sensitivity` uses only the
Johnson index, matching what ``get.spectral.sensitivity.R`` itself uses
despite the "Sobol" name.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Sequence

import numpy as np

from .canopy import foursail


def correlated_value(x: np.ndarray, r: float, rng: np.random.Generator | None = None) -> np.ndarray:
    """Draw ``y`` correlated with ``x`` at approximately ``r``. Direct port
    of ``ToolsRTM::correlatedValue``. Negative results are clipped to 0
    (assumes non-negative quantities, e.g. pigment concentrations).
    """
    if rng is None:
        rng = np.random.default_rng()
    x = np.asarray(x, dtype=float)
    sd = np.sqrt(max(1.0 - r ** 2, 0.0))
    y = r * x + rng.normal(0.0, sd, size=x.shape)
    return np.clip(y, 0.0, None)


def gauss_by_min_max(n: int, m: float, s: float, lwr: float, upr: float, nnorm: int,
                      rng: np.random.Generator | None = None) -> np.ndarray:
    """Truncated-normal sampling by rejection. Direct port of
    ``ToolsRTM::gauss_byMin_Max``: draw ``nnorm`` values from
    Normal(``m``, ``s``), keep only those within ``[lwr, upr]``, then
    randomly pick ``n`` of the survivors (without replacement). Raises
    ``ValueError`` if fewer than ``n`` values survive -- increase ``nnorm``.
    """
    if rng is None:
        rng = np.random.default_rng()
    samp = rng.normal(m, s, size=nnorm)
    samp = samp[(samp >= lwr) & (samp <= upr)]
    if samp.shape[0] < n:
        raise ValueError(f"Not enough values to sample from ({samp.shape[0]} < {n}). Try increasing nnorm.")
    return rng.choice(samp, size=n, replace=False)


def get_distribution_lut(
    minval: dict[str, float], maxval: dict[str, float], n_samples: int,
    type_distrib: dict[str, Literal["Uniform", "Gaussian"]],
    mean_gauss: dict[str, float] | None = None, std_gauss: dict[str, float] | None = None,
    dep_cab: bool = False, seed: int | None = None,
) -> dict[str, np.ndarray]:
    """Build a LUT with a per-trait distribution choice. Direct port of
    ``ToolsRTM::get_distributionLUT``. ``mean_gauss``/``std_gauss`` are only
    read for traits marked ``"Gaussian"`` in ``type_distrib``. If
    ``dep_cab`` and ``"Car"`` is one of the traits, ``Car`` is instead drawn
    as :func:`correlated_value` (``Cab/4``, ``r=0.8``) -- the empirical
    Cab-Car co-variation from leaf pigment data.
    """
    rng = np.random.default_rng(seed)
    out: dict[str, np.ndarray] = {}
    for trait in minval:
        if type_distrib[trait] == "Uniform":
            out[trait] = rng.uniform(minval[trait], maxval[trait], size=n_samples)
        else:
            out[trait] = gauss_by_min_max(n_samples, mean_gauss[trait], std_gauss[trait],
                                           minval[trait], maxval[trait], n_samples * 2, rng=rng)
        if trait == "Car" and dep_cab:
            out["Car"] = correlated_value(out["Cab"] / 4, r=0.8, rng=rng)
    return out


@dataclass
class CorrelatedLutResult:
    """Result of :func:`get_cor`."""

    lut: dict[str, np.ndarray]
    covariance: np.ndarray  #: realized (n_inputs, n_inputs) correlation matrix


def get_cor(
    n_inputs: int, n_lut: int = 100, distribution: Literal["Uniform", "Normal"] = "Uniform",
    seed: int = 123, rho: float | None = None, var_names: Sequence[str] | None = None,
    min_range: Sequence[float] | None = None, max_range: Sequence[float] | None = None,
) -> CorrelatedLutResult:
    """Generate ``n_inputs`` mutually correlated variables at (approximately)
    ``rho``, each rescaled to its own ``[min_range, max_range]``. Direct port
    of ``ToolsRTM::getCor``.
    """
    if rho is None:
        raise ValueError("rho is required")
    rng = np.random.default_rng(seed)
    n = n_inputs

    if distribution == "Normal":
        Sigma = np.full((n, n), rho)
        np.fill_diagonal(Sigma, 1.0)
        mat = rng.multivariate_normal(np.zeros(n), Sigma, size=n_lut)
    elif distribution == "Uniform":
        mat = rng.uniform(0, 1, size=(n_lut, n))
        if rho != 0:
            Sigma = np.full((n, n), rho)
            np.fill_diagonal(Sigma, 1.0)
            L = np.linalg.cholesky(Sigma)
            mat = mat @ L.T
    else:
        raise ValueError("distribution must be 'Uniform' or 'Normal'")

    for i in range(n):
        lo, hi = mat[:, i].min(), mat[:, i].max()
        mat[:, i] = min_range[i] + (mat[:, i] - lo) / (hi - lo) * (max_range[i] - min_range[i])

    names = list(var_names) if var_names is not None else [f"Var_{i+1}" for i in range(n)]
    lut = {name: mat[:, i] for i, name in enumerate(names)}
    return CorrelatedLutResult(lut=lut, covariance=np.corrcoef(mat, rowvar=False))


def johnson_relative_weights(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Johnson (2000) relative-weights index for each column of ``X``
    against ``y``. Equivalent to R's ``sensitivity::johnson(X, y,
    logistic = FALSE)$johnson$original`` -- verified to reproduce it to 8
    decimal places on a fixed reference dataset. Weights sum to
    approximately the R-squared of the OLS regression of ``y`` on ``X``.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    p = X.shape[1]
    Xs = (X - X.mean(axis=0)) / X.std(axis=0, ddof=1)
    ys = (y - y.mean()) / y.std(ddof=1)

    Rxx = np.corrcoef(Xs, rowvar=False) if p > 1 else np.array([[1.0]])
    eigvals, eigvecs = np.linalg.eigh(Rxx)
    eigvals = np.clip(eigvals, 1e-12, None)
    lam_sqrt = eigvecs @ np.diag(np.sqrt(eigvals)) @ eigvecs.T
    lam_inv_sqrt = eigvecs @ np.diag(1.0 / np.sqrt(eigvals)) @ eigvecs.T

    z = Xs @ lam_inv_sqrt
    beta_z, *_ = np.linalg.lstsq(z, ys, rcond=None)
    return (lam_sqrt ** 2) @ (beta_z ** 2)


@dataclass
class SobolResult:
    """Result of :func:`sobol_indices`, one entry per input variable."""

    parameter: list[str]
    si: np.ndarray  #: simplified two-sample-split first-order index (rough indicator, not a precise Sobol Si)
    sti: np.ndarray  #: simplified two-sample-split total-order index (same caveat)
    i_johnson: np.ndarray  #: Johnson relative-importance index -- the reliable, independently-verifiable metric
    difference: np.ndarray  #: si - i_johnson
    si_norm: np.ndarray | None = None
    i_johnson_norm: np.ndarray | None = None


def sobol_indices(data: dict[str, np.ndarray], output: str, n: int, normalize: bool = False,
                   seed: int | None = None) -> SobolResult:
    """First-order/total Sobol-like indices plus the Johnson relative-
    importance index, computed directly from existing LUT+output data (no
    extra model evaluations needed). Direct port of
    ``ToolsRTM::get.sobol.indices``.

    ``i_johnson`` is the reliable, independently-verifiable metric of the
    two (see :func:`johnson_relative_weights`); ``si``/``sti`` are a rough,
    simplified two-sample-split estimator -- :func:`spectral_sensitivity`
    only uses ``i_johnson_norm``.
    """
    inputs = [k for k in data if k != output]
    n_rows = len(data[output])
    if n > n_rows / 2:
        raise ValueError("n cannot be greater than or equal to half the number of rows in the dataset.")

    X = np.column_stack([np.asarray(data[k], dtype=float) for k in inputs])
    y = np.asarray(data[output], dtype=float)
    Xs = (X - X.mean(axis=0)) / X.std(axis=0, ddof=1)
    ys = (y - y.mean()) / y.std(ddof=1)

    i_johnson = johnson_relative_weights(X, y)

    rng = np.random.default_rng(seed)
    idx1 = rng.choice(n_rows, size=n, replace=False)
    idx2 = np.setdiff1d(np.arange(n_rows), idx1)[:n]

    def _split_sobol(rows):
        si = np.empty(len(inputs))
        sti = np.empty(len(inputs))
        xi = Xs[rows]
        yi = ys[rows]
        for j in range(len(inputs)):
            f0 = np.sum(xi[:, j] * yi) / n
            VY = np.sum(xi[:, j] ** 2 + yi ** 2) / (2 * n - 1) - f0
            si[j] = (np.sum(xi[:, j] * yi) / (n - 1) - f0) / VY
            sti[j] = 1 - (np.sum(yi * yi) / (n - 1) - f0) / VY
        return si, sti

    si1, sti1 = _split_sobol(idx1)
    si2, sti2 = _split_sobol(idx2)
    si = (si1 + si2) / 2
    sti = (sti1 + sti2) / 2

    threshold = np.quantile(np.abs(i_johnson), 0.25)
    below = np.abs(i_johnson) < threshold
    si = np.where(below, 0.0, si)
    sti = np.where(below, 0.0, sti)

    result = SobolResult(parameter=list(inputs), si=si, sti=sti, i_johnson=i_johnson, difference=si - i_johnson)
    if normalize:
        result.si_norm = np.abs(si) / np.sum(np.abs(si)) * 100
        result.i_johnson_norm = np.abs(i_johnson) / np.sum(np.abs(i_johnson)) * 100
    return result


_DEFAULT_TRAIT_BOUNDS = {
    # (lower, upper) -- matches ToolsRTM::inputsPROSAIL's ranges for these traits.
    "N": (1.5, 2.5), "Cab": (5.0, 75.0), "EWT": (0.001, 0.035), "LMA": (0.001, 0.035),
    "LIDFa": (30.0, 70.0), "LAI": (2.0, 5.0),
}


@dataclass
class SpectralSensitivityResult:
    """Result of :func:`spectral_sensitivity` -- long-format arrays, ready
    for a stacked-area "relative contribution vs wavelength" plot."""

    wavelength: np.ndarray  #: (n_wl * n_traits,)
    trait: np.ndarray  #: (n_wl * n_traits,) string array
    sti_pct: np.ndarray  #: (n_wl * n_traits,) Johnson relative importance, normalized to sum to 100% per wavelength
    distribution: str


def spectral_sensitivity(
    n_samples: int = 1000, distribution: Literal["Uniform", "Gaussian"] = "Uniform",
    traits: Sequence[str] = ("N", "Cab", "EWT", "LMA", "LIDFa", "LAI"),
    rsoil_base: np.ndarray | None = None, wl_step: int = 5, seed: int = 123,
) -> SpectralSensitivityResult:
    """Run ``foursail`` (PROSPECT-D leaf optics) many times while varying a
    set of plant/soil traits, then compute the Johnson relative-importance
    index at each wavelength -- how much each trait relatively contributes
    to explaining reflectance variance there. Direct port of
    ``ToolsRTM::get.spectral.sensitivity`` (fixed leaf/canopy model: PROSPECT-D
    + fourSAIL, matching that function's own defaults; the R version's
    ``leaf.model``/``canopy.model`` arguments aren't reproduced here).

    ``traits`` must be a subset of :data:`_DEFAULT_TRAIT_BOUNDS`'s keys
    (``N``, ``Cab``, ``EWT``, ``LMA``, ``LIDFa``, ``LAI``); a soil-brightness
    multiplier (``SoilCoef``, 0.5-1.5) is always added on top, matching the
    classic PROSAIL sensitivity figure (leaf structure, pigment, water, dry
    matter, leaf angle, LAI, soil).
    """
    rng = np.random.default_rng(seed)
    if rsoil_base is None:
        rsoil_base = np.full(2101, 0.15)
    wl_all = np.arange(400, 2501)

    lut: dict[str, np.ndarray] = {}
    for trait in traits:
        lo, hi = _DEFAULT_TRAIT_BOUNDS[trait]
        if distribution == "Uniform":
            lut[trait] = rng.uniform(lo, hi, size=n_samples)
        else:
            lut[trait] = gauss_by_min_max(n_samples, (lo + hi) / 2, (hi - lo) / 4, lo, hi, n_samples * 2, rng=rng)

    fixed_defaults = dict(N=1.5, Cab=40.0, Car=8.0, Anth=1.0, Cbrown=0.0, EWT=0.01, LMA=0.009, alpha=40.0,
                           LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1.0, LAI=3.0, hspot=0.01, tts=30.0, tto=0.0, psi=0.0)

    soil_coef = rng.uniform(0.5, 1.5, size=n_samples) if distribution == "Uniform" \
        else np.maximum(0.1, rng.normal(1.0, 0.25, size=n_samples))
    lut["SoilCoef"] = soil_coef
    all_traits = list(traits) + ["SoilCoef"]

    refl = np.full((n_samples, len(wl_all)), np.nan)
    for i in range(n_samples):
        inputLUT = dict(fixed_defaults)
        for trait in traits:
            inputLUT[trait] = lut[trait][i]
        try:
            sail = foursail(inputLUT, rsoil_base * soil_coef[i], leaf_model="PROSPECT-D", spectrum_all=True)
            refl[i, :] = sail.rsot
        except Exception:
            pass

    ok = ~np.isnan(refl).any(axis=1)
    refl = refl[ok]
    lut = {k: v[ok] for k, v in lut.items()}
    wl_sel = np.arange(0, len(wl_all), wl_step)
    out_wl, out_trait, out_pct = [], [], []
    X_all = np.column_stack([lut[t] for t in all_traits])

    for j in wl_sel:
        y = refl[:, j]
        try:
            i_johnson = johnson_relative_weights(X_all, y)
        except Exception:
            continue
        i_johnson = np.where(np.isnan(i_johnson) | (i_johnson < 0), 0.0, i_johnson)
        total = i_johnson.sum()
        pct = (i_johnson / total * 100) if total > 0 else np.zeros_like(i_johnson)
        out_wl.extend([wl_all[j]] * len(all_traits))
        out_trait.extend(all_traits)
        out_pct.extend(pct.tolist())

    return SpectralSensitivityResult(
        wavelength=np.array(out_wl), trait=np.array(out_trait), sti_pct=np.array(out_pct),
        distribution=distribution,
    )
