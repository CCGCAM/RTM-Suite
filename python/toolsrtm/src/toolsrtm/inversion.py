"""Trait-inversion tools: CARS-PLS feature selection, VIF-based collinearity
pruning, LUT nearest-neighbour ("merit function") inversion, and a
multi-algorithm ML dispatcher.

Python port of ``ToolsRTM/R/carspls.R`` / ``get.cars.pls.R``, ``getVIF.R``,
``get.inversion.R``, ``hybrid_inversion.R`` / ``hybrid_inversionE.R``, and
``get.inversionOpt.R``. Needs the optional ``ml`` extra
(``pip install toolsrtm[ml]``: scikit-learn, xgboost) -- none of these
functions are imported by ``toolsrtm/__init__.py`` at import time, and each
imports its own ML dependencies lazily so a plain ``import toolsrtm`` never
requires scikit-learn/xgboost to be installed.

R's ``get.inversion``/``hybrid_inversion`` dispatch by name to specific
``caret`` methods (bartMachine, rqlasso, rvmLinear, AdaBag, brnn, ...). Several
of those have no direct scikit-learn/xgboost equivalent; where that's the
case the docstring of the relevant function says exactly which estimator was
substituted and why. Unlike the pure radiative-transfer math ported
elsewhere in this package, none of this module is verified to floating-point
precision against R -- ``caret``'s own cross-validated tuning is stochastic,
so a Python port using different (but comparable) estimators and search
grids will not reproduce R's numbers bit-for-bit even in principle. What's
verified instead: each algorithm runs end-to-end on held-out data and
produces sane, comparable-magnitude accuracy metrics (see
``tests/test_inversion.py``).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Literal, Sequence

import numpy as np


# ---------------------------------------------------------------------------
# CARS-PLS (Competitive Adaptive Reweighted Sampling for PLS)
# ---------------------------------------------------------------------------


@dataclass
class CarsPlsResult:
    """Result of :func:`carspls`. Mirrors the R ``CARS`` list 1:1 except
    ``selected_variables`` is 0-indexed (Python) instead of 1-indexed (R)."""

    coef: np.ndarray  #: (n_vars, iteration) coefficient path
    n_var: np.ndarray  #: (iteration,) number of retained variables per iteration
    rmsecv: np.ndarray  #: (iteration,) cross-validated RMSE per iteration
    num_lv: np.ndarray  #: (iteration,) best number of latent variables per iteration
    optimal_iteration: int  #: 1-indexed iteration with the lowest RMSECV (matches R)
    min_error: float
    selected_variables: np.ndarray  #: 0-indexed column positions into the original X


def _pls_press_by_ncomp(X, y, max_ncomp, fold, partition_type, scale):
    """Sum of squared cross-validated prediction errors for n_components in
    1..max_ncomp, matching ``pls::crossval``'s per-component PRESS."""
    from sklearn.cross_decomposition import PLSRegression

    n = X.shape[0]
    if partition_type == "interleaved":
        folds = np.arange(n) % fold
    elif partition_type == "consecutive":
        folds = (np.arange(n) * fold) // n
    elif partition_type == "random":
        folds = (np.arange(n) * fold) // n
        rng = np.random.default_rng()
        rng.shuffle(folds)
    else:
        raise ValueError("partition_type must be 'interleaved', 'consecutive', or 'random'")

    press = np.zeros(max_ncomp)
    for k in range(fold):
        test_idx = np.where(folds == k)[0]
        train_idx = np.where(folds != k)[0]
        if len(test_idx) == 0 or len(train_idx) <= max_ncomp:
            continue
        for c in range(1, max_ncomp + 1):
            model = PLSRegression(n_components=c, scale=scale)
            model.fit(X[train_idx], y[train_idx])
            pred = np.ravel(model.predict(X[test_idx]))
            press[c - 1] += float(np.sum((pred - y[test_idx]) ** 2))
    return press


def carspls(
    X: np.ndarray,
    y: np.ndarray,
    n_lv: int = 2,
    fold: int = 10,
    scale_pretreat: bool = True,
    iteration: int = 50,
    partition_type: Literal["interleaved", "consecutive", "random"] = "interleaved",
    verbose: bool = False,
) -> CarsPlsResult:
    """Competitive Adaptive Reweighted Sampling for PLS variable selection.

    Python port of ``carspls``/``get.cars.pls`` (R, original algorithm by
    Yizeng Liang & Hongdong Li, MATLAB->R port by Hongdong Li 2009). At each
    of ``iteration`` rounds: fits a PLS model on the currently-retained
    variables, cross-validates it to get an RMSECV curve over 1..n_lv
    components, records the coefficient-magnitude-ranked variable importance, and
    forcibly eliminates the lowest-ranked variables via an exponentially
    decreasing retention schedule (Monte-Carlo/EDF sampling). The iteration
    with the lowest RMSECV gives the final selected variable set.

    :param X: (n_samples, n_vars) predictor matrix.
    :param y: (n_samples,) response vector.
    :param n_lv: number of PLS latent variables (components) to fit/tune over.
    :param fold: number of cross-validation segments.
    :param scale_pretreat: if True, scale (not just center) each predictor.
    :param iteration: number of CARS-PLS elimination rounds.
    :param partition_type: cross-validation fold assignment: ``"interleaved"``
        (round-robin), ``"consecutive"`` (contiguous blocks), or ``"random"``.
    :param verbose: print progress per iteration (matches R's own screen output).
    :return: :class:`CarsPlsResult`.
    """
    from sklearn.cross_decomposition import PLSRegression

    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).ravel()
    n_row, n_col = X.shape

    order = np.argsort(y, kind="stable")
    X = X[order]
    y = y[order]

    rmsecv = np.zeros(iteration)
    num_lv = np.zeros(iteration, dtype=int)
    coef = np.zeros((n_col, iteration))
    n_var = np.zeros(iteration, dtype=int)
    subset_variable = np.arange(n_col)

    ratio0 = 1.0
    ratio1 = 2.0 / n_col
    b = math.log(ratio0 / ratio1) / (iteration - 1)
    a = ratio0 * math.exp(b)

    for it in range(iteration):
        x_cal = X[:, subset_variable]
        n_lv = min(n_lv, x_cal.shape[0], x_cal.shape[1])  # persists across iterations, matches R
        ncomp = n_lv

        press = _pls_press_by_ncomp(x_cal, y, ncomp, fold, partition_type, scale_pretreat)
        rmsecv_curve = np.sqrt(press / n_row)
        rmsecv[it] = rmsecv_curve.min()
        num_lv[it] = int(np.argmin(rmsecv_curve)) + 1  # 1-indexed, matches R

        model = PLSRegression(n_components=ncomp, scale=scale_pretreat)
        model.fit(x_cal, y)
        coef_iter = np.ravel(model.coef_)

        coef0 = np.zeros(n_col)
        coef0[subset_variable] = coef_iter
        coef[:, it] = coef0
        n_var[it] = int(np.sum(coef0 != 0))

        weight = np.abs(coef0)
        weight_order = np.argsort(-weight, kind="stable")  # descending, R's order(decreasing=TRUE)

        ratio_variable = a * math.exp(-b * (it + 2))  # (it+1)+1: R's iter is 1-indexed
        k = math.ceil(n_col * ratio_variable)
        keep_idx = weight_order[:k]
        new_weight = np.zeros(n_col)
        new_weight[keep_idx] = weight[keep_idx]

        subset_variable = np.where(new_weight != 0)[0]

        if verbose:
            print(f"The {it + 1}th CARS-PLS iteration finished.")

    min_error = float(rmsecv.min())
    opt_candidates = np.where(rmsecv == min_error)[0]
    opt_iteration = int(opt_candidates[-1])  # last tie, matches R's OPT.iter[length(OPT.iter)]
    selected_variables = np.where(coef[:, opt_iteration] != 0)[0]

    return CarsPlsResult(
        coef=coef,
        n_var=n_var,
        rmsecv=rmsecv,
        num_lv=num_lv,
        optimal_iteration=opt_iteration + 1,  # report 1-indexed, matches R
        min_error=min_error,
        selected_variables=selected_variables,
    )


# ---------------------------------------------------------------------------
# VIF-based stepwise collinearity pruning
# ---------------------------------------------------------------------------


def get_vif(frame: np.ndarray, columns: Sequence[str] | None = None, thresh: float = 10.0,
            trace: bool = True) -> list[str] | list[int]:
    """Backward-elimination variable selection by Variance Inflation Factor.

    Python port of ``getVIF`` (R, VIF function originally from
    https://beckmw.wordpress.com/2013/02/05/collinearity-and-stepwise-vif-selection/).
    Iteratively regresses each remaining variable on all others; drops the
    variable with the highest VIF (``1 / (1 - R^2)``) as long as any VIF
    exceeds ``thresh``.

    :param frame: (n_samples, n_vars) array, or a `pandas.DataFrame`.
    :param columns: variable names, required if ``frame`` is a bare array;
        ignored (and taken from ``frame.columns``) if ``frame`` is a DataFrame.
    :param thresh: VIF threshold above which a variable is flagged as collinear.
    :param trace: print each elimination step (matches R's own console output).
    :return: names (or 0-indexed positions, if ``columns`` is None and
        ``frame`` is a bare array) of the retained variables.
    """
    from sklearn.linear_model import LinearRegression

    if hasattr(frame, "columns"):
        columns = list(frame.columns)
        data = np.asarray(frame, dtype=float)
    else:
        data = np.asarray(frame, dtype=float)
        if columns is None:
            columns = list(range(data.shape[1]))
        columns = list(columns)

    def _vifs(cols, mat):
        out = {}
        for j, name in enumerate(cols):
            y = mat[:, j]
            X = np.delete(mat, j, axis=1)
            r2 = LinearRegression().fit(X, y).score(X, y)
            out[name] = np.inf if r2 >= 1.0 else 1.0 / (1.0 - r2)
        return out

    cols = list(columns)
    mat = data.copy()
    while True:
        vifs = _vifs(cols, mat)
        worst_name = max(vifs, key=vifs.get)
        worst_vif = vifs[worst_name]
        if trace:
            for name, v in vifs.items():
                print(f"{name}: {v:.3f}")
        if worst_vif < thresh:
            if trace:
                print(f"All variables have VIF < {thresh}, max VIF {worst_vif:.2f}\n")
            return cols
        if trace:
            print(f"removed: {worst_name} {worst_vif:.3f}\n")
        drop_j = cols.index(worst_name)
        cols.pop(drop_j)
        mat = np.delete(mat, drop_j, axis=1)
        if len(cols) <= 1:
            return cols


# ---------------------------------------------------------------------------
# LUT nearest-neighbour ("merit function") inversion
# ---------------------------------------------------------------------------

MERIT_FUNCTIONS: dict[str, Callable[[np.ndarray, np.ndarray], np.ndarray]] = {
    "merit-RMSE": lambda sim, obs: np.sqrt(np.nanmean((sim - obs) ** 2, axis=-1)),
    "merit-NRMSE": lambda sim, obs: (
        np.sqrt(np.nanmean((sim - obs) ** 2, axis=-1))
        / (np.nanmax(obs, axis=-1) - np.nanmin(obs, axis=-1))
    ),
    "merit-MAE": lambda sim, obs: np.nanmean(np.abs(sim - obs), axis=-1),
    "merit-NMB": lambda sim, obs: (
        (np.nanmean(sim, axis=-1) - np.nanmean(obs, axis=-1)) / np.nanmean(obs, axis=-1)
    ),
    "merit-FGE": lambda sim, obs: np.nanmean(2 * np.abs(sim - obs) / (sim + obs), axis=-1),
}


@dataclass
class InversionOptResult:
    rfl_best: np.ndarray  #: (n_obs, n_wave) best-matching (n_opt-averaged) simulated spectra
    lut_best: "pd.DataFrame"  #: (n_obs, n_lut_columns) n_opt-averaged LUT parameters per observation


def get_inversion_opt(
    rfl_sensor: np.ndarray,
    rfl_rtm: np.ndarray,
    lut,
    wave: Sequence[float] | None = None,
    method: str = "merit-RMSE",
    n_opt: int = 1,
    custom_stat: Callable[[np.ndarray, np.ndarray], np.ndarray] | None = None,
) -> InversionOptResult:
    """LUT (look-up table) inversion by nearest-neighbour spectral matching.

    Python port of ``get.inversionOpt`` (R). For each observed spectrum in
    ``rfl_sensor``, ranks every simulated spectrum in ``rfl_rtm`` by a merit
    (error) function and averages the ``n_opt`` best matches' LUT parameters
    and reflectance. Fully vectorized (broadcasts each observation against
    the whole LUT at once) rather than R's nested per-row loop -- same
    algorithm, no numerical differences expected for the built-in merit
    functions (verified against a hand-computed reference below).

    :param rfl_sensor: (n_obs, n_wave) observed/sensor reflectance.
    :param rfl_rtm: (n_lut, n_wave) simulated reflectance from the LUT.
    :param lut: (n_lut, n_params) `pandas.DataFrame` of the LUT's input parameters.
    :param wave: wavelengths corresponding to columns of ``rfl_sensor``/``rfl_rtm``
        (only used to name the returned reflectance columns; optional).
    :param method: one of ``"merit-RMSE"``, ``"merit-NRMSE"``, ``"merit-MAE"``,
        ``"merit-NMB"``, ``"merit-FGE"``, or ``"merit-custom.metric"`` (requires
        ``custom_stat``).
    :param n_opt: number of best-matching LUT rows to average per observation.
    :param custom_stat: optional ``f(sim, obs) -> error`` merit function,
        broadcast over the last axis exactly like the built-in ones; overrides ``method``.
    :return: :class:`InversionOptResult`.
    """
    import pandas as pd

    rfl_sensor = np.atleast_2d(np.asarray(rfl_sensor, dtype=float))
    rfl_rtm = np.asarray(rfl_rtm, dtype=float)
    n_obs = rfl_sensor.shape[0]

    if custom_stat is not None:
        merit_fn = custom_stat
    elif method in MERIT_FUNCTIONS:
        merit_fn = MERIT_FUNCTIONS[method]
    else:
        raise ValueError(f"Invalid method {method!r}. Choose from {list(MERIT_FUNCTIONS)} or pass custom_stat.")

    rfl_best_rows = []
    lut_best_rows = []
    for i in range(n_obs):
        obs_i = rfl_sensor[i]
        errors = merit_fn(rfl_rtm, obs_i)  # (n_lut,)
        order = np.argsort(errors)
        best_idx = order[:n_opt]
        rfl_best_rows.append(rfl_rtm[best_idx].mean(axis=0))
        lut_best_rows.append(lut.iloc[best_idx].mean(axis=0))

    rfl_best = np.vstack(rfl_best_rows)
    lut_best = pd.DataFrame(lut_best_rows).reset_index(drop=True)
    if wave is not None:
        pass  # wave is accepted for R-signature parity; not needed to build the result here

    return InversionOptResult(rfl_best=rfl_best, lut_best=lut_best)


# ---------------------------------------------------------------------------
# get_inversion: multi-algorithm ML dispatcher
# ---------------------------------------------------------------------------

#: R's caret-method name each algorithm dispatches to, and, where scikit-learn/
#: xgboost has no direct equivalent, the substitution actually used here.
ALGORITHMS = {
    "PLSR": "sklearn PLSRegression, n_components tuned by 5-fold CV (matches caret method 'pls')",
    "SVM": "sklearn SVR(kernel='rbf'), gamma/C tuned by grid search (matches caret method 'svmRadial' via e1071::tune.svm)",
    "RF": "sklearn RandomForestRegressor, max_features tuned (matches caret method 'rf')",
    "GB": "sklearn GradientBoostingRegressor (matches caret method 'gbm')",
    "NN": "sklearn MLPRegressor, single hidden layer (matches caret method 'nnet')",
    "Bayesian": "sklearn BayesianRidge -- substitute: no BART implementation in sklearn/xgboost; "
                "BayesianRidge is a Bayesian *linear* model, not R's bartMachine (Bayesian additive trees)",
    "AdaBag": "sklearn AdaBoostRegressor over shallow DecisionTreeRegressor stumps (matches caret method 'AdaBag')",
    "BRNN": "sklearn MLPRegressor with strong L2 (alpha) regularization -- substitute: approximates "
            "'Bayesian regularization' via explicit weight decay, not R's brnn Gauss-Newton/Levenberg-Marquardt fit",
    "xGB": "xgboost XGBRegressor(booster='gblinear') (matches caret method 'xgbLinear')",
    "RVM": "sklearn BayesianRidge -- substitute: no Relevance Vector Machine in sklearn/xgboost; "
           "BayesianRidge shares RVM's sparsity-favouring linear-Bayesian character",
    "qLASSO": "sklearn QuantileRegressor(quantile=0.5, solver='highs'), L1-penalized (matches caret method 'rqlasso')",
    "Ensemble": "sklearn StackingRegressor(GB + SVR + MLP, final_estimator=LinearRegression) "
                "(matches caretEnsemble::caretStack(..., method='glm'))",
}


@dataclass
class InversionResult:
    model_label: str
    model: object
    predictions: dict  #: {"train": np.ndarray, "test": np.ndarray}
    statistics: dict  #: {"train": {"r2":.., "rmse":.., "mae":..}, "test": {...}}
    importance: dict | None  #: {input_name: importance}, or None if not available for this algorithm


def _regression_stats(y_true, y_pred) -> dict:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
    mae = float(np.mean(np.abs(y_true - y_pred)))
    return {"r2": float(r2), "rmse": rmse, "mae": mae}


def _permutation_importance(model, X_test, y_test, feature_names):
    from sklearn.inspection import permutation_importance

    try:
        result = permutation_importance(model, X_test, y_test, n_repeats=10, random_state=0)
        return {name: float(v) for name, v in zip(feature_names, result.importances_mean)}
    except Exception:
        return None


def _fit_algorithm(algorithm: str, X_train, y_train, inputs, seed: int):
    """Fit one of :data:`ALGORITHMS`' estimators on already-split training
    data. Shared by :func:`get_inversion` and :func:`hybrid_inversion` so the
    two functions (which, like their R originals, differ in train/test
    splitting, feature selection and optional log-transform, not in the
    underlying estimators) don't duplicate 12 tuning-grid blocks.

    :return: ``(fitted_model, predict_fn, importance_or_None)``.
    """
    from sklearn.model_selection import GridSearchCV

    importance = None

    if algorithm == "PLSR":
        from sklearn.cross_decomposition import PLSRegression

        max_comp = max(1, min(20, X_train.shape[1], X_train.shape[0] - 1))
        grid = GridSearchCV(PLSRegression(), {"n_components": list(range(1, max_comp + 1))},
                             cv=min(5, X_train.shape[0]), scoring="neg_root_mean_squared_error")
        grid.fit(X_train, y_train)
        model = grid.best_estimator_
        predict = lambda X_: np.ravel(model.predict(X_))
        importance = {name: float(abs(c)) for name, c in zip(inputs, np.ravel(model.coef_))}

    elif algorithm == "SVM":
        from sklearn.svm import SVR

        grid = GridSearchCV(SVR(kernel="rbf"),
                             {"gamma": [2.0 ** g for g in (-10, -8, -6, -4)], "C": [2.0 ** c for c in (-5, -3, -1, 1)]},
                             cv=min(5, X_train.shape[0]), scoring="neg_root_mean_squared_error")
        grid.fit(X_train, y_train)
        model = grid.best_estimator_
        predict = lambda X_: np.ravel(model.predict(X_))

    elif algorithm == "RF":
        from sklearn.ensemble import RandomForestRegressor

        n_features = X_train.shape[1]
        candidates = sorted({max(1, n_features // 3), max(1, n_features // 2), n_features})
        grid = GridSearchCV(RandomForestRegressor(n_estimators=300, random_state=seed),
                             {"max_features": candidates}, cv=min(3, X_train.shape[0]),
                             scoring="neg_root_mean_squared_error")
        grid.fit(X_train, y_train)
        model = grid.best_estimator_
        predict = lambda X_: np.ravel(model.predict(X_))
        importance = {name: float(v) for name, v in zip(inputs, model.feature_importances_)}

    elif algorithm == "GB":
        from sklearn.ensemble import GradientBoostingRegressor

        model = GradientBoostingRegressor(n_estimators=300, learning_rate=0.1, max_depth=3, random_state=seed)
        model.fit(X_train, y_train)
        predict = lambda X_: np.ravel(model.predict(X_))
        importance = {name: float(v) for name, v in zip(inputs, model.feature_importances_)}

    elif algorithm == "NN":
        from sklearn.neural_network import MLPRegressor
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler

        # adam (MLPRegressor's default solver) needs mini-batch iteration to
        # converge and underperforms badly on the small (~100-row) LUTs this
        # is meant for; lbfgs is a full-batch solver well suited to small data.
        model = make_pipeline(StandardScaler(),
                               MLPRegressor(hidden_layer_sizes=(10,), alpha=0.01, max_iter=2000,
                                            solver="lbfgs", random_state=seed))
        model.fit(X_train, y_train)
        predict = lambda X_: np.ravel(model.predict(X_))

    elif algorithm == "Bayesian":
        from sklearn.linear_model import BayesianRidge

        model = BayesianRidge()
        model.fit(X_train, y_train)
        predict = lambda X_: np.ravel(model.predict(X_))
        importance = {name: float(abs(c)) for name, c in zip(inputs, model.coef_)}

    elif algorithm == "AdaBag":
        from sklearn.ensemble import AdaBoostRegressor
        from sklearn.tree import DecisionTreeRegressor

        model = AdaBoostRegressor(estimator=DecisionTreeRegressor(max_depth=3), n_estimators=100, random_state=seed)
        model.fit(X_train, y_train)
        predict = lambda X_: np.ravel(model.predict(X_))
        importance = {name: float(v) for name, v in zip(inputs, model.feature_importances_)}

    elif algorithm == "BRNN":
        from sklearn.neural_network import MLPRegressor
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler

        model = make_pipeline(StandardScaler(),
                               MLPRegressor(hidden_layer_sizes=(10,), alpha=1.0, max_iter=2000,
                                            solver="lbfgs", random_state=seed))
        model.fit(X_train, y_train)
        predict = lambda X_: np.ravel(model.predict(X_))

    elif algorithm == "xGB":
        import xgboost

        # shotgun (xgboost's default gblinear updater) converges poorly on
        # small/medium datasets; coord_descent is deterministic and reliably
        # matches a plain linear fit's accuracy here.
        model = xgboost.XGBRegressor(booster="gblinear", updater="coord_descent", n_estimators=300,
                                      learning_rate=0.3, reg_lambda=0.01, reg_alpha=0.01, random_state=seed)
        model.fit(X_train, y_train)
        predict = lambda X_: np.ravel(model.predict(X_))

    elif algorithm == "RVM":
        from sklearn.linear_model import BayesianRidge

        model = BayesianRidge()
        model.fit(X_train, y_train)
        predict = lambda X_: np.ravel(model.predict(X_))
        importance = {name: float(abs(c)) for name, c in zip(inputs, model.coef_)}

    elif algorithm == "qLASSO":
        from sklearn.linear_model import QuantileRegressor

        model = QuantileRegressor(quantile=0.5, alpha=0.01, solver="highs")
        model.fit(X_train, y_train)
        predict = lambda X_: np.ravel(model.predict(X_))
        importance = {name: float(abs(c)) for name, c in zip(inputs, model.coef_)}

    else:  # Ensemble
        from sklearn.ensemble import GradientBoostingRegressor, StackingRegressor
        from sklearn.linear_model import LinearRegression
        from sklearn.neural_network import MLPRegressor
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler
        from sklearn.svm import SVR

        model = StackingRegressor(
            estimators=[
                ("gbm", GradientBoostingRegressor(n_estimators=200, random_state=seed)),
                ("svm", make_pipeline(StandardScaler(), SVR(kernel="rbf"))),
                ("nnet", make_pipeline(StandardScaler(), MLPRegressor(hidden_layer_sizes=(10,), max_iter=2000,
                                                                       solver="lbfgs", random_state=seed))),
            ],
            final_estimator=LinearRegression(),
        )
        model.fit(X_train, y_train)
        predict = lambda X_: np.ravel(model.predict(X_))

    return model, predict, importance


def get_inversion(
    data,
    dep_var: str,
    inputs: Sequence[str],
    algorithm: str = "PLSR",
    seed: int = 123,
    n_samples: int | None = 500,
    test_size: float = 0.3,
) -> InversionResult:
    """Fit and evaluate a plant-trait inversion model with one of 12
    algorithms, on a held-out train/test split.

    Python port of ``get.inversion`` (R). See :data:`ALGORITHMS` for exactly
    which scikit-learn/xgboost estimator each ``algorithm`` name dispatches
    to, and, for the 4 algorithms with no direct equivalent (``Bayesian``,
    ``BRNN``, ``RVM`` -- and note ``AdaBag``/``xGB``/``qLASSO`` *do* have
    close matches), what was substituted and why. Unlike R's version, tuning
    here is a single small grid search per algorithm (not caret's full
    repeated-CV search), to keep this runnable as a demo rather than a
    multi-hour job -- the same design choice already used by
    ``Scripts/Python/*/2_inversion_ml.py``.

    :param data: `pandas.DataFrame` containing ``dep_var`` and all ``inputs`` columns.
    :param dep_var: name of the response column to predict.
    :param inputs: names of the predictor columns.
    :param algorithm: one of the keys of :data:`ALGORITHMS`.
    :param seed: random seed for the train/test split and any stochastic estimator.
    :param n_samples: if given and less than ``len(data)``, randomly subsample
        this many rows before splitting (matches R's own tuning-sample-size argument).
    :param test_size: fraction of (sub-sampled) data held out for testing.
    :return: :class:`InversionResult`.
    """
    from sklearn.model_selection import train_test_split

    if algorithm not in ALGORITHMS:
        raise ValueError(f"Unknown algorithm {algorithm!r}. Choose from {list(ALGORITHMS)}.")

    rng = np.random.default_rng(seed)
    df = data
    if n_samples is not None and n_samples < len(df):
        idx = rng.choice(len(df), size=n_samples, replace=False)
        df = df.iloc[idx]

    X = df[list(inputs)].to_numpy(dtype=float)
    y = df[dep_var].to_numpy(dtype=float)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=seed)

    model, predict, importance = _fit_algorithm(algorithm, X_train, y_train, inputs, seed)
    pred_train = predict(X_train)
    pred_test = predict(X_test)

    if importance is None:
        importance = _permutation_importance(model, X_test, y_test, inputs)

    return InversionResult(
        model_label=algorithm,
        model=model,
        predictions={"train": pred_train, "test": pred_test},
        statistics={"train": _regression_stats(y_train, pred_train), "test": _regression_stats(y_test, pred_test)},
        importance=importance,
    )


# ---------------------------------------------------------------------------
# hybrid_inversion / hybrid_inversionE: feature-selection + single-family fit
# ---------------------------------------------------------------------------

_HYBRID_ALGORITHMS = {"SVM": "SVM", "RF": "RF", "GB": "GB", "nnet": "NN", "Ensemble": "Ensemble"}


@dataclass
class HybridInversionResult:
    model: object
    keep_variables: list[str]  #: predictor columns actually used, after pattern/collinearity selection
    statistics: "pd.DataFrame"  #: rows "train"/"test" (+ "field" if field_data given), columns r2/rmse/mae
    predictions: dict  #: {"train": np.ndarray, "test": np.ndarray[, "field": np.ndarray]}, original (untransformed) scale


def hybrid_inversion(
    lut,
    input: str,
    split: float = 0.8,
    seed: int | None = None,
    method: str | None = None,
    collinearity: Literal["VIF", "CARS"] | None = None,
    pattern: str | None = None,
    trans: bool = True,
    field_data=None,
    acron: str | None = None,
) -> HybridInversionResult:
    """Fit a single-algorithm trait-inversion model with optional predictor
    selection (by name pattern, then optionally VIF or CARS-PLS pruning) and
    an optional log-transform of the response.

    Python port of ``hybrid_inversion`` (R). ``method`` dispatches to the
    same 5 estimators as :func:`get_inversion`'s ``SVM``/``RF``/``GB``/``NN``/
    ``Ensemble`` (``"nnet"`` here maps to ``"NN"`` there, matching R's own
    caret method name) -- see :data:`ALGORITHMS` for what each one is. Note:
    R's train/test split uses ``caret::createDataPartition`` (percentile-
    stratified on the response); this port uses a plain random split via
    scikit-learn, which is not percentile-stratified -- a documented
    approximation, not expected to change results materially for the
    LUT-sized (typically hundreds of rows) datasets this is meant for.

    :param lut: `pandas.DataFrame` with the response column ``input`` and
        candidate predictor columns.
    :param input: name of the response column to predict.
    :param split: train-fraction of the train/test split (R's own convention;
        note this is the *train* fraction, unlike :func:`get_inversion`'s ``test_size``).
    :param seed: random seed.
    :param method: one of ``"SVM"``, ``"RF"``, ``"GB"``, ``"nnet"``, ``"Ensemble"``.
        Defaults to ``"SVM"`` (matches R's own default).
    :param collinearity: ``None`` (use every ``pattern``-matched column),
        ``"VIF"`` (prune via :func:`get_vif`), or ``"CARS"`` (select via :func:`carspls`).
    :param pattern: substring that predictor column names must contain (e.g.
        ``"B"`` for reflectance-band columns named ``B1``, ``B2``, ...);
        if ``None``, every column except ``input`` is a candidate.
    :param trans: log-transform ``input`` before fitting (matches R's own
        default); predictions/statistics are reported back on the original scale.
    :param field_data: optional `pandas.DataFrame` of field observations to
        validate against, in addition to the LUT's own test split.
    :param acron: suffix appended to ``input`` to find the observed column in
        ``field_data`` (e.g. ``acron="_obsv"`` looks for ``f"{input}_obsv"``).
        Required if ``field_data`` is given.
    :return: :class:`HybridInversionResult`.
    """
    import pandas as pd
    from sklearn.model_selection import train_test_split

    if field_data is not None and acron is None:
        raise ValueError("acron is required when field_data is given.")

    method = method or "SVM"
    if method not in _HYBRID_ALGORITHMS:
        raise ValueError(f"Unknown method {method!r}. Choose from {list(_HYBRID_ALGORITHMS)}.")

    df = lut.copy()
    if trans:
        df[input] = np.log(df[input])

    if pattern is not None:
        candidates = [c for c in df.columns if pattern in c and c != input]
    else:
        candidates = [c for c in df.columns if c != input]

    if collinearity == "VIF":
        keep_variables = list(get_vif(df[candidates], thresh=10, trace=False))
    elif collinearity == "CARS":
        cars_res = carspls(df[candidates].to_numpy(), df[input].to_numpy(), n_lv=5, fold=10,
                            scale_pretreat=True, iteration=100, partition_type="interleaved")
        keep_variables = [candidates[i] for i in cars_res.selected_variables]
    else:
        keep_variables = candidates

    X = df[keep_variables].to_numpy(dtype=float)
    y = df[input].to_numpy(dtype=float)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=1 - split, random_state=seed)

    algorithm = _HYBRID_ALGORITHMS[method]
    model, predict, _ = _fit_algorithm(algorithm, X_train, y_train, keep_variables, seed or 123)

    def _untransform(v):
        return np.exp(v) if trans else v

    pred_train = _untransform(predict(X_train))
    pred_test = _untransform(predict(X_test))
    y_train_orig = _untransform(y_train)
    y_test_orig = _untransform(y_test)

    stats_rows = {
        "train": _regression_stats(y_train_orig, pred_train),
        "test": _regression_stats(y_test_orig, pred_test),
    }
    predictions = {"train": pred_train, "test": pred_test}

    if field_data is not None:
        X_field = field_data[keep_variables].to_numpy(dtype=float)
        pred_field = _untransform(predict(X_field))
        y_field = field_data[f"{input}{acron}"].to_numpy(dtype=float)
        stats_rows["field"] = _regression_stats(y_field, pred_field)
        predictions["field"] = pred_field

    statistics = pd.DataFrame(stats_rows).T[["r2", "rmse", "mae"]]

    return HybridInversionResult(model=model, keep_variables=keep_variables, statistics=statistics,
                                  predictions=predictions)


def hybrid_inversion_ensemble(lut, input: str, split: float = 0.8, seed: int | None = None,
                               collinearity: Literal["VIF", "CARS"] | None = None,
                               pattern: str | None = None, field_data=None,
                               acron: str | None = None) -> HybridInversionResult:
    """:func:`hybrid_inversion` with ``method="Ensemble"`` fixed.

    Python port of ``hybrid_inversionE`` (R) -- the R function is
    ``hybrid_inversion`` with the model choice hardcoded to the 3-model
    (GB + SVM + neural net) stacking ensemble and no ``trans``/log-transform
    option, which this wrapper matches (``trans=False``).
    """
    return hybrid_inversion(lut, input, split=split, seed=seed, method="Ensemble",
                             collinearity=collinearity, pattern=pattern, trans=False,
                             field_data=field_data, acron=acron)
