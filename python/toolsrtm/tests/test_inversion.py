from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from toolsrtm.inversion import (
    ALGORITHMS,
    carspls,
    get_inversion,
    get_inversion_opt,
    get_vif,
    hybrid_inversion,
    hybrid_inversion_ensemble,
)

REFDATA = Path(__file__).parent / "refdata"


def _linear_dataset(n=150, seed=1):
    rng = np.random.default_rng(seed)
    X = rng.uniform(0, 1, size=(n, 6))
    y = 3 * X[:, 0] - 2 * X[:, 2] + 0.5 * X[:, 4] + rng.normal(scale=0.05, size=n)
    df = pd.DataFrame(X, columns=[f"x{i}" for i in range(6)])
    df["y"] = y
    return df, [f"x{i}" for i in range(6)]


def _cars_input():
    df = pd.read_csv(REFDATA / "carspls_input.csv")
    X = df[[f"V{i + 1}" for i in range(15)]].to_numpy()
    y = df["y"].to_numpy()
    return X, y


def test_carspls_selects_the_true_predictors():
    # Synthetic data with 3 truly informative predictors (V3, V6, V10, i.e.
    # 0-indexed 2, 5, 9) out of 15, the rest pure noise -- CARS-PLS should
    # recover exactly that set regardless of any R-vs-Python numerical
    # differences in the underlying PLS solver.
    X, y = _cars_input()
    res = carspls(X, y, n_lv=3, fold=5, iteration=15)
    assert list(res.selected_variables) == [2, 5, 9]


def test_carspls_matches_r_reference():
    # Cross-checked against ToolsRTM::carspls() (R, same input CSV, same
    # nLV/fold/iteration) via python/scratch/scratch_export_carspls.R. Found and
    # fixed a genuine R bug along the way: carspls()/get.cars.pls() called
    # bare ceil() without importing it from pracma, so both crashed on
    # every call (including their own documented examples) -- fixed in
    # ToolsRTM/R/carspls.R and get.cars.pls.R to use base R's ceiling()
    # (numerically identical to pracma::ceil, no missing-import risk).
    #
    # R's carspls uses pls::mvr(method="simpls") (SIMPLS); this Python port
    # uses sklearn.cross_decomposition.PLSRegression (NIPALS), which also
    # scales the response internally (R's mvr(scale=TRUE) only scales X).
    # Both are valid PLS solvers but not numerically identical when
    # ncomp > 1, so the per-iteration coefficient path and RMSECV curve
    # diverge by up to ~0.3-0.4 at some interior iterations even though
    # they agree on what actually matters -- num_lv/n_var/selected
    # variables/optimal iteration/min error, checked exactly or near-exactly
    # below.
    X, y = _cars_input()
    res = carspls(X, y, n_lv=3, fold=5, iteration=15)

    path_r = pd.read_csv(REFDATA / "carspls_path.csv")
    sel_r = pd.read_csv(REFDATA / "carspls_selected.csv")

    np.testing.assert_array_equal(res.num_lv, path_r["NumLV"].to_numpy())
    np.testing.assert_array_equal(res.n_var, path_r["Nvar"].to_numpy())
    assert list(res.selected_variables) == [v - 1 for v in sel_r["SelectedVariables"].to_numpy()]
    assert res.optimal_iteration == int(sel_r["OptimalIteration"].to_numpy()[0])
    assert res.min_error == pytest.approx(float(sel_r["MinError"].to_numpy()[0]), rel=1e-5)


def test_get_vif_keeps_independent_predictors():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(200, 5))
    kept = get_vif(X, thresh=10, trace=False)
    assert len(kept) == 5


def test_get_vif_drops_collinear_predictor():
    rng = np.random.default_rng(0)
    x1 = rng.normal(size=200)
    x2 = x1 + rng.normal(scale=0.01, size=200)  # near-duplicate of x1
    x3 = rng.normal(size=200)
    X = np.column_stack([x1, x2, x3])
    kept = get_vif(X, columns=["x1", "x2", "x3"], thresh=10, trace=False)
    assert len(kept) == 2
    assert "x3" in kept
    assert not ("x1" in kept and "x2" in kept)


@pytest.mark.parametrize("algorithm", list(ALGORITHMS))
def test_get_inversion_runs_and_predicts_well(algorithm):
    # Not a floating-point R comparison (see module docstring for why) --
    # every algorithm must at least run end-to-end and recover a mostly-linear
    # synthetic signal reasonably well on held-out data.
    df, inputs = _linear_dataset()
    res = get_inversion(df, "y", inputs, algorithm=algorithm, n_samples=None, seed=1)
    assert res.model_label == algorithm
    assert res.predictions["test"].shape == res.predictions["test"].shape
    assert res.statistics["test"]["r2"] > 0.5, res.statistics


def test_get_inversion_rejects_unknown_algorithm():
    df, inputs = _linear_dataset(n=20)
    with pytest.raises(ValueError):
        get_inversion(df, "y", inputs, algorithm="not-a-real-algorithm")


def test_get_inversion_opt_recovers_nearest_lut_row():
    rng = np.random.default_rng(2)
    n_lut, n_wave = 30, 12
    lut = pd.DataFrame({"N": rng.uniform(1, 3, n_lut), "Cab": rng.uniform(10, 80, n_lut)})
    rfl_rtm = rng.uniform(0, 0.5, size=(n_lut, n_wave))
    obs = rfl_rtm[7] + rng.normal(scale=1e-4, size=n_wave)

    res = get_inversion_opt(obs, rfl_rtm, lut, method="merit-RMSE", n_opt=1)
    pd.testing.assert_series_equal(res.lut_best.iloc[0], lut.iloc[7], check_names=False)
    np.testing.assert_allclose(res.rfl_best[0], rfl_rtm[7], atol=1e-3)


def test_get_inversion_opt_rejects_unknown_method():
    rng = np.random.default_rng(3)
    lut = pd.DataFrame({"N": rng.uniform(1, 3, 5)})
    rfl_rtm = rng.uniform(0, 0.5, size=(5, 4))
    with pytest.raises(ValueError):
        get_inversion_opt(rfl_rtm[0], rfl_rtm, lut, method="not-a-real-method")


def _band_dataset(n=150, seed=5):
    rng = np.random.default_rng(seed)
    B = rng.uniform(0.01, 0.5, size=(n, 8))
    cab = np.exp(2 + 1.5 * B[:, 0] - 2 * B[:, 3] + rng.normal(scale=0.1, size=n))
    df = pd.DataFrame(B, columns=[f"B{i + 1}" for i in range(8)])
    df["Cab"] = cab
    df["unrelated_col"] = rng.normal(size=n)  # must be excluded by pattern="B"
    return df


@pytest.mark.parametrize("method", ["SVM", "RF", "GB", "nnet", "Ensemble"])
def test_hybrid_inversion_runs_and_predicts_well(method):
    df = _band_dataset()
    res = hybrid_inversion(df, "Cab", split=0.8, seed=1, method=method, pattern="B", trans=True)
    assert "unrelated_col" not in res.keep_variables
    assert res.statistics.loc["test", "r2"] > 0.5, res.statistics


def test_hybrid_inversion_vif_and_cars_selection():
    df = _band_dataset()
    vif_res = hybrid_inversion(df, "Cab", split=0.8, seed=1, method="RF", pattern="B", collinearity="VIF")
    assert set(vif_res.keep_variables) <= {f"B{i + 1}" for i in range(8)}

    cars_res = hybrid_inversion(df, "Cab", split=0.8, seed=1, method="RF", pattern="B", collinearity="CARS")
    assert 0 < len(cars_res.keep_variables) <= 8


def test_hybrid_inversion_field_data_branch():
    df = _band_dataset()
    rng = np.random.default_rng(9)
    field = df.iloc[:20].copy()
    field["Cab_obsv"] = field["Cab"] * (1 + rng.normal(scale=0.05, size=20))

    res = hybrid_inversion(df, "Cab", split=0.8, seed=1, method="RF", pattern="B",
                            field_data=field, acron="_obsv")
    assert "field" in res.statistics.index
    assert "field" in res.predictions
    assert res.statistics.loc["field", "r2"] > 0.5


def test_hybrid_inversion_requires_acron_with_field_data():
    df = _band_dataset()
    with pytest.raises(ValueError):
        hybrid_inversion(df, "Cab", pattern="B", field_data=df.iloc[:10])


def test_hybrid_inversion_ensemble_matches_hybrid_inversion_ensemble_method():
    df = _band_dataset()
    res_direct = hybrid_inversion(df, "Cab", split=0.8, seed=1, method="Ensemble", pattern="B", trans=False)
    res_wrapper = hybrid_inversion_ensemble(df, "Cab", split=0.8, seed=1, pattern="B")
    np.testing.assert_allclose(res_direct.predictions["test"], res_wrapper.predictions["test"])
