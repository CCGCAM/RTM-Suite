"""Invert leaf/canopy traits (Cab, LAI, EWT) from spectral indices using
classical (non-deep-learning) machine learning.

Python equivalent of the R course pipeline's caret-based *-inversion_ML.R
scripts (Scripts/ForPROSAIL/4-inversion_ML.R etc.). ``caret`` itself
(cross-validated tuning/resampling plumbing across 12 named algorithms) has
no meaningful 1:1 Python translation, so this uses scikit-learn's native
equivalents directly instead: a comparable spread of algorithm families
(partial least squares, tree ensembles, kernel methods, linear, a small
neural net), each fit once per trait on the indices computed in
2_spectral_indices.py and evaluated on a held-out test split. See
4_inversion_dl.py for the deep-learning counterpart (TensorFlow and
PyTorch), matching the R side's keras-based *-inversion_DL.R.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from metrics import all_metrics
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.compose import TransformedTargetRegressor

OUT_DIR = Path(__file__).resolve().parents[3] / "outs" / "Python" / "ForPROSAIL"

TARGET_TRAITS = ["Cab", "LAI", "EWT"]


def _scaled(estimator):
    """Standardize both X and y around an estimator that's sensitive to
    input/output scale (SVR, MLP) -- unlike PLS/RandomForest/GradientBoosting,
    which are scale-invariant or handle raw index/trait magnitudes fine.
    Matters especially for a trait like EWT (values ~0.005-0.03): without
    also rescaling y, MLPRegressor's default weight init assumes an O(1)
    target and its predictions diverge. Matches R's caret
    preProc=c("center","scale") for its own scale-sensitive models."""
    return TransformedTargetRegressor(
        regressor=make_pipeline(StandardScaler(), estimator), transformer=StandardScaler(),
    )


MODELS = {
    "PLS": PLSRegression(n_components=8),
    "RandomForest": RandomForestRegressor(n_estimators=300, random_state=42),
    "GradientBoosting": GradientBoostingRegressor(random_state=42),
    "SVR": _scaled(SVR(kernel="rbf", C=10, gamma="scale")),
    "Ridge": _scaled(Ridge(alpha=1.0)),
    "MLP": _scaled(MLPRegressor(hidden_layer_sizes=(64, 32), activation="relu", max_iter=2000, random_state=42)),
}


def main():
    lut = pd.read_csv(OUT_DIR / "1_LUT.csv")
    indices = pd.read_csv(OUT_DIR / "2_indices.csv")
    # Drop predictor columns with any non-finite values (some indices are
    # NaN by construction -- see toolsrtm/indices.py -- or can blow up for
    # specific samples; not usable as ML predictors either way).
    X_full = indices.replace([np.inf, -np.inf], np.nan).dropna(axis=1)
    print(f"Using {X_full.shape[1]} of {indices.shape[1]} indices as predictors "
          f"(dropped {indices.shape[1] - X_full.shape[1]} with NaN/inf values).")

    results = []
    ensemble_preds = {trait: {} for trait in TARGET_TRAITS}
    n_cols = len(MODELS) + 1  # +1 for the ensemble panel
    fig, axes = plt.subplots(len(TARGET_TRAITS), n_cols, figsize=(4 * n_cols, 4 * len(TARGET_TRAITS)))

    for i, trait in enumerate(TARGET_TRAITS):
        y = lut[trait].to_numpy()
        X_train, X_test, y_train, y_test = train_test_split(X_full, y, test_size=0.3, random_state=42)

        for j, (name, model) in enumerate(MODELS.items()):
            model.fit(X_train, y_train)
            y_pred = np.ravel(model.predict(X_test))
            ensemble_preds[trait][name] = y_pred

            m = all_metrics(y_test, y_pred)
            results.append({"trait": trait, "model": name, "n_test": len(y_test), **m})

            ax = axes[i, j]
            ax.scatter(y_test, y_pred, alpha=0.7, edgecolor="k", linewidth=0.3)
            lims = [min(y_test.min(), y_pred.min()), max(y_test.max(), y_pred.max())]
            ax.plot(lims, lims, "k--", linewidth=1)
            ax.set_title(f"{trait} - {name}\nR2={m['R2']:.2f}, RMSE={m['RMSE']:.3f}")
            ax.set_xlabel("Observed")
            ax.set_ylabel("Predicted")

        # Ensemble: simple average across every classical model's prediction.
        y_ens = np.mean(np.column_stack(list(ensemble_preds[trait].values())), axis=1)
        m_ens = all_metrics(y_test, y_ens)
        results.append({"trait": trait, "model": "Ensemble", "n_test": len(y_test), **m_ens})

        ax = axes[i, -1]
        ax.scatter(y_test, y_ens, alpha=0.7, edgecolor="k", linewidth=0.3, color="darkorange")
        lims = [min(y_test.min(), y_ens.min()), max(y_test.max(), y_ens.max())]
        ax.plot(lims, lims, "k--", linewidth=1)
        ax.set_title(f"{trait} - Ensemble\nR2={m_ens['R2']:.2f}, RMSE={m_ens['RMSE']:.3f}")
        ax.set_xlabel("Observed")
        ax.set_ylabel("Predicted")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "3_inversion_scatter.png", dpi=150)
    plt.close(fig)

    metrics_df = pd.DataFrame(results)
    metrics_df.to_csv(OUT_DIR / "3_inversion_metrics.csv", index=False)

    print(metrics_df.to_string(index=False))
    print(f"\nMetrics -> {OUT_DIR / '3_inversion_metrics.csv'}")
    print(f"Plot -> {OUT_DIR / '3_inversion_scatter.png'}")


if __name__ == "__main__":
    main()
