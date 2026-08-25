"""Invert SMC (gravimetric soil moisture content) from spectral indices
using classical ML.

Python equivalent of Scripts/R/ForMARMIT/2-inversion_ML.R. Same
scikit-learn-native approach as ForPROSAIL_fourSAIL/3_inversion_ml.py --
see that script's own docstring for why (no 1:1 caret port).
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

OUT_DIR = Path(__file__).resolve().parents[3] / "outs" / "Python" / "ForMARMIT"
TARGET_TRAITS = ["SMC"]
LUT_COLUMNS = ["L", "eps", "SMC"]


def _scaled(estimator):
    return TransformedTargetRegressor(regressor=make_pipeline(StandardScaler(), estimator),
                                       transformer=StandardScaler())


MODELS = {
    "PLS": PLSRegression(n_components=8),
    "RandomForest": RandomForestRegressor(n_estimators=300, random_state=42),
    "GradientBoosting": GradientBoostingRegressor(random_state=42),
    "SVR": _scaled(SVR(kernel="rbf", C=10, gamma="scale")),
    "Ridge": _scaled(Ridge(alpha=1.0)),
    "MLP": _scaled(MLPRegressor(hidden_layer_sizes=(64, 32), activation="relu", max_iter=2000, random_state=42)),
}


def run_for_dataset(dataset_path: Path, label: str, results: list):
    dataset = pd.read_csv(dataset_path)
    predictor_cols = [c for c in dataset.columns if c not in LUT_COLUMNS]
    X_full = dataset[predictor_cols].replace([np.inf, -np.inf], np.nan).dropna(axis=1)

    fig, axes = plt.subplots(len(TARGET_TRAITS), len(MODELS) + 1, squeeze=False,
                              figsize=(4 * (len(MODELS) + 1), 4 * len(TARGET_TRAITS)))
    for i, trait in enumerate(TARGET_TRAITS):
        y = dataset[trait].to_numpy()
        X_train, X_test, y_train, y_test = train_test_split(X_full, y, test_size=0.3, random_state=42)
        ens_preds = {}
        for j, (name, model) in enumerate(MODELS.items()):
            model.fit(X_train, y_train)
            y_pred = np.ravel(model.predict(X_test))
            ens_preds[name] = y_pred
            m = all_metrics(y_test, y_pred)
            results.append({"dataset": label, "trait": trait, "model": name, "n_test": len(y_test), **m})
            ax = axes[i, j]
            ax.scatter(y_test, y_pred, alpha=0.7, edgecolor="k", linewidth=0.3)
            lims = [min(y_test.min(), y_pred.min()), max(y_test.max(), y_pred.max())]
            ax.plot(lims, lims, "k--", linewidth=1)
            ax.set_title(f"{trait} - {name}\nR2={m['R2']:.2f}")
        y_ens = np.mean(np.column_stack(list(ens_preds.values())), axis=1)
        m_ens = all_metrics(y_test, y_ens)
        results.append({"dataset": label, "trait": trait, "model": "Ensemble", "n_test": len(y_test), **m_ens})
        ax = axes[i, -1]
        ax.scatter(y_test, y_ens, alpha=0.7, edgecolor="k", linewidth=0.3, color="darkorange")
        lims = [min(y_test.min(), y_ens.min()), max(y_test.max(), y_ens.max())]
        ax.plot(lims, lims, "k--", linewidth=1)
        ax.set_title(f"{trait} - Ensemble\nR2={m_ens['R2']:.2f}")

    fig.tight_layout()
    fig.savefig(OUT_DIR / f"2_inversion_scatter_{label}.png", dpi=150)
    plt.close(fig)


def main():
    results = []
    for label in ["native", "sentinel2a", "sentinel2b", "prisma"]:
        path = OUT_DIR / f"1_dataset_{label}.csv"
        if not path.exists():
            print(f"Skipping {label} (run 1_simulate_lut.py first)")
            continue
        print(f"=== ML inversion: {label} ===")
        run_for_dataset(path, label, results)

    metrics_df = pd.DataFrame(results)
    metrics_df.to_csv(OUT_DIR / "2_inversion_metrics.csv", index=False)
    print(metrics_df.groupby(["dataset", "trait", "model"])["R2"].mean().to_string())
    print(f"\nMetrics -> {OUT_DIR / '2_inversion_metrics.csv'}")


if __name__ == "__main__":
    main()
