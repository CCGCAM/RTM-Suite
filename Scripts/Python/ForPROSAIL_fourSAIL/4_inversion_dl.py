"""Invert leaf/canopy traits (Cab, LAI, EWT) from spectral indices using
deep learning, in **both** TensorFlow/Keras and PyTorch.

Python equivalent of the R course pipeline's keras-based *-inversion_DL.R
scripts (ToolsRTM::getMLmodel, model='Hidden-layers'): a 3-hidden-layer
dense network (64 -> 32 (dropout 0.1) -> 16 -> 1, ReLU activations, Adam
optimizer, early stopping on validation loss) trained per trait on the
indices computed in 2_spectral_indices.py -- built twice, once per
framework, so the two are directly comparable on the same data/split.

Run after 2_spectral_indices.py (and, for the combined comparison table at
the end, after 3_inversion_ml.py).
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from metrics import all_metrics
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

OUT_DIR = Path(__file__).resolve().parents[3] / "outs" / "Python" / "ForPROSAIL"

TARGET_TRAITS = ["Cab", "LAI", "EWT"]
N_EPOCHS = 200
BATCH_SIZE = 16
PATIENCE = 15
RANDOM_STATE = 42


def _load_predictors():
    lut = pd.read_csv(OUT_DIR / "1_LUT.csv")
    indices = pd.read_csv(OUT_DIR / "2_indices.csv")
    X_full = indices.replace([np.inf, -np.inf], np.nan).dropna(axis=1)
    return lut, X_full


def _train_test_split_scaled(X_full: pd.DataFrame, y: np.ndarray):
    """Standardize both X and y before handing them to a neural net.
    Without also rescaling y, a trait with tiny raw values (EWT, ~0.005-0.03)
    makes training diverge -- same fix as 3_inversion_ml.py's `_scaled()`
    for MLP/SVR, needed here for the exact same reason. Returns the fitted
    y-scaler too, so predictions can be mapped back to physical units."""
    X_train, X_test, y_train, y_test = train_test_split(X_full, y, test_size=0.3, random_state=RANDOM_STATE)
    x_scaler = StandardScaler().fit(X_train)
    y_scaler = StandardScaler().fit(y_train.reshape(-1, 1))
    X_train_s = x_scaler.transform(X_train)
    X_test_s = x_scaler.transform(X_test)
    y_train_s = y_scaler.transform(y_train.reshape(-1, 1)).ravel()
    return X_train_s, X_test_s, y_train_s, y_test, y_scaler


def train_tensorflow(X_train, y_train, X_test, n_features: int) -> np.ndarray:
    """3-hidden-layer MLP, matching R's getMLmodel(model='Hidden-layers')
    architecture exactly (64 -> 32(dropout 0.1) -> 16 -> 1, relu, adam)."""
    import tensorflow as tf
    from tensorflow import keras

    tf.random.set_seed(RANDOM_STATE)
    model = keras.Sequential([
        keras.layers.Input(shape=(n_features,)),
        keras.layers.Dense(64, activation="relu"),
        keras.layers.Dense(32, activation="relu"),
        keras.layers.Dropout(0.1),
        keras.layers.Dense(16, activation="relu"),
        keras.layers.Dense(1, activation="linear"),
    ])
    model.compile(optimizer=keras.optimizers.Adam(learning_rate=1e-3), loss="mse", metrics=["mae"])
    early_stop = keras.callbacks.EarlyStopping(monitor="val_loss", patience=PATIENCE, restore_best_weights=True)
    model.fit(
        X_train, y_train, epochs=N_EPOCHS, batch_size=BATCH_SIZE, validation_split=0.2,
        callbacks=[early_stop], verbose=0,
    )
    return model.predict(X_test, verbose=0).ravel()


def train_pytorch(X_train, y_train, X_test, n_features: int) -> np.ndarray:
    """Same 3-hidden-layer MLP architecture as train_tensorflow(), hand-rolled
    training loop (PyTorch has no built-in .fit()/EarlyStopping)."""
    import torch
    from torch import nn

    torch.manual_seed(RANDOM_STATE)

    class MLP(nn.Module):
        def __init__(self, n_in: int):
            super().__init__()
            self.net = nn.Sequential(
                nn.Linear(n_in, 64), nn.ReLU(),
                nn.Linear(64, 32), nn.ReLU(), nn.Dropout(0.1),
                nn.Linear(32, 16), nn.ReLU(),
                nn.Linear(16, 1),
            )

        def forward(self, x):
            return self.net(x)

    # Hold out a validation split from the training data for early stopping,
    # matching R's/TensorFlow's validation_split=0.2 above.
    n_val = max(1, int(0.2 * len(X_train)))
    rng = np.random.default_rng(RANDOM_STATE)
    idx = rng.permutation(len(X_train))
    val_idx, fit_idx = idx[:n_val], idx[n_val:]

    X_fit = torch.tensor(X_train[fit_idx], dtype=torch.float32)
    y_fit = torch.tensor(y_train[fit_idx], dtype=torch.float32).unsqueeze(1)
    X_val = torch.tensor(X_train[val_idx], dtype=torch.float32)
    y_val = torch.tensor(y_train[val_idx], dtype=torch.float32).unsqueeze(1)
    X_test_t = torch.tensor(X_test, dtype=torch.float32)

    model = MLP(n_features)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
    loss_fn = nn.MSELoss()

    best_val_loss = float("inf")
    best_state = None
    epochs_since_improve = 0
    n_train = len(X_fit)

    for epoch in range(N_EPOCHS):
        model.train()
        perm = torch.randperm(n_train)
        for start in range(0, n_train, BATCH_SIZE):
            batch_idx = perm[start:start + BATCH_SIZE]
            optimizer.zero_grad()
            pred = model(X_fit[batch_idx])
            loss = loss_fn(pred, y_fit[batch_idx])
            loss.backward()
            optimizer.step()

        model.eval()
        with torch.no_grad():
            val_loss = loss_fn(model(X_val), y_val).item()
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
            epochs_since_improve = 0
        else:
            epochs_since_improve += 1
            if epochs_since_improve >= PATIENCE:
                break

    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        return model(X_test_t).numpy().ravel()


def main():
    lut, X_full = _load_predictors()
    print(f"Using {X_full.shape[1]} predictors, {len(lut)} samples.")

    results = []
    fig, axes = plt.subplots(len(TARGET_TRAITS), 3, figsize=(12, 4 * len(TARGET_TRAITS)))
    preds_for_ensemble = {}

    for i, trait in enumerate(TARGET_TRAITS):
        y = lut[trait].to_numpy()
        X_train, X_test, y_train_s, y_test, y_scaler = _train_test_split_scaled(X_full, y)
        n_features = X_train.shape[1]

        # Models are trained on standardized y; map predictions back to
        # physical units before computing metrics/plotting.
        y_pred_tf = y_scaler.inverse_transform(
            train_tensorflow(X_train, y_train_s, X_test, n_features).reshape(-1, 1)
        ).ravel()
        y_pred_pt = y_scaler.inverse_transform(
            train_pytorch(X_train, y_train_s, X_test, n_features).reshape(-1, 1)
        ).ravel()
        y_pred_ens = 0.5 * (y_pred_tf + y_pred_pt)
        preds_for_ensemble[trait] = {"y_test": y_test, "TensorFlow": y_pred_tf, "PyTorch": y_pred_pt}

        for name, y_pred, ax in [
            ("TensorFlow", y_pred_tf, axes[i, 0]),
            ("PyTorch", y_pred_pt, axes[i, 1]),
            ("DL-Ensemble", y_pred_ens, axes[i, 2]),
        ]:
            m = all_metrics(y_test, y_pred)
            results.append({"trait": trait, "model": name, "n_test": len(y_test), **m})
            ax.scatter(y_test, y_pred, alpha=0.7, edgecolor="k", linewidth=0.3)
            lims = [min(y_test.min(), y_pred.min()), max(y_test.max(), y_pred.max())]
            ax.plot(lims, lims, "k--", linewidth=1)
            ax.set_title(f"{trait} - {name}\nR2={m['R2']:.2f}, RMSE={m['RMSE']:.3f}")
            ax.set_xlabel("Observed")
            ax.set_ylabel("Predicted")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "4_inversion_dl_scatter.png", dpi=150)
    plt.close(fig)

    metrics_df = pd.DataFrame(results)
    metrics_df.to_csv(OUT_DIR / "4_inversion_dl_metrics.csv", index=False)
    print(metrics_df.to_string(index=False))
    print(f"\nMetrics -> {OUT_DIR / '4_inversion_dl_metrics.csv'}")
    print(f"Plot -> {OUT_DIR / '4_inversion_dl_scatter.png'}")

    # Combined classical + deep-learning comparison, if 3_inversion_ml.py has
    # already been run.
    classical_path = OUT_DIR / "3_inversion_metrics.csv"
    if classical_path.exists():
        classical_df = pd.read_csv(classical_path)
        combined = pd.concat([classical_df, metrics_df], ignore_index=True)
        combined.to_csv(OUT_DIR / "5_inversion_combined_metrics.csv", index=False)
        best = combined.loc[combined.groupby("trait")["RMSE"].idxmin()]
        print("\nBest model per trait (lowest RMSE, classical + deep learning):")
        print(best[["trait", "model", "R2", "RMSE"]].to_string(index=False))
        print(f"\nCombined metrics -> {OUT_DIR / '5_inversion_combined_metrics.csv'}")


if __name__ == "__main__":
    main()
