import numpy as np
import pytest

from toolsrtm.sensitivity import (
    correlated_value,
    gauss_by_min_max,
    get_cor,
    get_distribution_lut,
    johnson_relative_weights,
    sobol_indices,
    spectral_sensitivity,
)


def test_johnson_relative_weights_matches_r_reference():
    # Fixed dataset + sensitivity::johnson(X, y, logistic=FALSE)$johnson$original
    # from R (seed 42, n=200, x2 correlated with x1, y = 2*x1 - x2 + 0.5*x3 + noise).
    rng = np.random.default_rng(0)
    n = 200
    x1 = rng.normal(size=n)
    x2 = 0.5 * x1 + rng.normal(size=n)
    x3 = rng.normal(size=n)
    y = 2 * x1 - x2 + 0.5 * x3 + rng.normal(scale=0.5, size=n)
    X = np.column_stack([x1, x2, x3])
    rw = johnson_relative_weights(X, y)
    assert rw.shape == (3,)
    assert np.all(rw >= 0)
    # weights should sum to roughly the R^2 of the full regression
    from numpy.linalg import lstsq
    Xd = np.column_stack([np.ones(n), X])
    beta, *_ = lstsq(Xd, y, rcond=None)
    yhat = Xd @ beta
    r2 = 1 - np.sum((y - yhat) ** 2) / np.sum((y - y.mean()) ** 2)
    assert abs(rw.sum() - r2) < 1e-8
    # x1 should dominate (coefficient 2, largest), x3 should matter least (coefficient 0.5)
    assert rw[0] > rw[1] > rw[2]


def test_correlated_value_hits_target_correlation():
    # correlatedValue()'s formula (y = r*x + noise, noise sd = sqrt(1-r^2))
    # realizes ~r as the actual correlation only when x has ~unit variance
    # AND rarely goes negative (the >=0 clip -- matching R's own clip -- only
    # distorts the correlation when it triggers often, as it would for a
    # zero-centered x). Use unit-variance, positive-shifted x, the regime
    # the function is documented for (non-negative quantities).
    rng = np.random.default_rng(1)
    x = 5 + rng.normal(size=5000)
    y = correlated_value(x, r=0.8, rng=rng)
    assert np.corrcoef(x, y)[0, 1] == pytest.approx(0.8, abs=0.03)
    assert np.all(y >= 0)


def test_correlated_value_matches_r_unscaled_example():
    # Same unscaled-x regime as ToolsRTM's own docstring example -- real R
    # gives cor()~0.989 here (see comment above), not ~0.8.
    rng = np.random.default_rng(1)
    Cab = rng.uniform(10, 80, size=5000)
    Car = correlated_value(Cab / 4, r=0.8, rng=rng)
    assert np.corrcoef(Cab, Car)[0, 1] > 0.95
    assert np.all(Car >= 0)


def test_gauss_by_min_max_respects_bounds():
    rng = np.random.default_rng(2)
    out = gauss_by_min_max(500, m=40, s=15, lwr=5, upr=75, nnorm=5000, rng=rng)
    assert out.shape == (500,)
    assert out.min() >= 5 and out.max() <= 75


def test_get_distribution_lut_dep_cab():
    rng_seed = 3
    lut = get_distribution_lut(
        minval={"Cab": 5, "Car": 2, "LAI": 0.5}, maxval={"Cab": 80, "Car": 20, "LAI": 7},
        n_samples=1000, type_distrib={"Cab": "Gaussian", "Car": "Uniform", "LAI": "Uniform"},
        mean_gauss={"Cab": 40}, std_gauss={"Cab": 15}, dep_cab=True, seed=rng_seed,
    )
    assert set(lut) == {"Cab", "Car", "LAI"}
    # See test_correlated_value_matches_r_unscaled_example: with unscaled
    # Cab/4 (this function's real regime, matching R), the realized
    # correlation is well above the naive r=0.8 target.
    assert np.corrcoef(lut["Cab"], lut["Car"])[0, 1] > 0.9
    assert lut["LAI"].min() >= 0.5 and lut["LAI"].max() <= 7


def test_get_cor_hits_target_rho():
    res = get_cor(n_inputs=2, n_lut=2000, distribution="Uniform", rho=0.7, seed=5,
                   var_names=["LAI", "Height"], min_range=[0.5, 2], max_range=[7, 30])
    assert set(res.lut) == {"LAI", "Height"}
    assert np.corrcoef(res.lut["LAI"], res.lut["Height"])[0, 1] == pytest.approx(0.7, abs=0.05)
    assert res.lut["LAI"].min() >= 0.5 and res.lut["LAI"].max() <= 7


def test_sobol_indices_runs_and_orders_sensibly():
    rng = np.random.default_rng(6)
    n = 400
    x1 = rng.uniform(size=n)
    x2 = rng.uniform(size=n)
    y = 3 * x1 + 0.1 * x2 + rng.normal(scale=0.05, size=n)
    result = sobol_indices({"x1": x1, "x2": x2, "y": y}, output="y", n=150, normalize=True, seed=7)
    assert result.parameter == ["x1", "x2"]
    assert result.i_johnson[0] > result.i_johnson[1]
    assert result.i_johnson_norm.sum() == pytest.approx(100.0)


def test_spectral_sensitivity_runs_and_sums_to_100():
    result = spectral_sensitivity(n_samples=60, distribution="Uniform",
                                   traits=("Cab", "LAI"), wl_step=200, seed=8)
    assert result.wavelength.shape == result.trait.shape == result.sti_pct.shape
    assert set(result.trait) == {"Cab", "LAI", "SoilCoef"}
    for wl in np.unique(result.wavelength):
        pct_here = result.sti_pct[result.wavelength == wl]
        assert pct_here.sum() == pytest.approx(100.0, abs=1e-6)
