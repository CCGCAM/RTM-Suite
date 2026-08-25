"""Invert Cab, LAI, EWT from native-resolution spectral indices using deep
learning, in **both** TensorFlow/Keras and PyTorch.

Python equivalent of Scripts/R/ForFoursail2/3-inversion_DL.R
(ToolsRTM::getMLmodel, model='Hidden-layers'): same 3-hidden-layer dense
network (64 -> 32(dropout 0.1) -> 16 -> 1, ReLU, Adam, early stopping) as
ForPROSAIL_fourSAIL/4_inversion_dl.py -- see that script's own docstring
for the full architecture rationale. Scoped to the native-resolution
dataset only (not all 4 sensor datasets like 2_inversion_ml.py) to keep
runtime reasonable -- swap OUT_DIR / "1_dataset_native.csv" for a sensor
dataset below if you want the sensor-band DL comparison too.
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

OUT_DIR = Path(__file__).resolve().parents[3] / "outs" / "Python" / "ForFoursail2"
TARGET_TRAITS = ["Cab", "LAI", "EWT"]
LUT_COLUMNS = ["N", "Cab", "Car", "Anth", "Cbrown", "EWT", "LMA", "alpha", "Prot", "CBC",
               "LIDFa", "LIDFb", "TypeLidf", "LAI", "hspot", "tts", "tto", "psi",
               "fraction_brown", "diss", "Cv", "Zeta"]
N_EPOCHS = 200
BATCH_SIZE = 16
PATIENCE = 15
RANDOM_STATE = 42


def _load_predictors():
    dataset = pd.read_csv(OUT_DIR / "1_dataset_native.csv")
    predictor_cols = [c for c in dataset.columns if c not in LUT_COLUMNS]
    X_full = dataset[predictor_cols].replace([np.inf, -np.inf], np.nan).dropna(axis=1)
    return dataset, X_full


def _train_test_split_scaled(X_full, y):
    X_train, X_test, y_train, y_test = train_test_split(X_full, y, test_size=0.3, random_state=RANDOM_STATE)
    x_scaler = StandardScaler().fit(X_train)
    y_scaler = StandardScaler().fit(y_train.reshape(-1, 1))
    return (x_scaler.transform(X_train), x_scaler.transform(X_test),
            y_scaler.transform(y_train.reshape(-1, 1)).ravel(), y_test, y_scaler)


def train_tensorflow(X_train, y_train, X_test, n_features):
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
    model.fit(X_train, y_train, epochs=N_EPOCHS, batch_size=BATCH_SIZE, validation_split=0.2,
              callbacks=[early_stop], verbose=0)
    return model.predict(X_test, verbose=0).ravel()


def train_pytorch(X_train, y_train, X_test, n_features):
    import torch
    from torch import nn

    torch.manual_seed(RANDOM_STATE)

    class MLP(nn.Module):
        def __init__(self, n_in):
            super().__init__()
            self.net = nn.Sequential(
                nn.Linear(n_in, 64), nn.ReLU(),
                nn.Linear(64, 32), nn.ReLU(), nn.Dropout(0.1),
                nn.Linear(32, 16), nn.ReLU(),
                nn.Linear(16, 1),
            )

        def forward(self, x):
            return self.net(x)

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
    best_val_loss, best_state, since_improve = float("inf"), None, 0
    n_train = len(X_fit)

    for _ in range(N_EPOCHS):
        model.train()
        perm = torch.randperm(n_train)
        for start in range(0, n_train, BATCH_SIZE):
            batch_idx = perm[start:start + BATCH_SIZE]
            optimizer.zero_grad()
            loss = loss_fn(model(X_fit[batch_idx]), y_fit[batch_idx])
            loss.backward()
            optimizer.step()
        model.eval()
        with torch.no_grad():
            val_loss = loss_fn(model(X_val), y_val).item()
        if val_loss < best_val_loss:
            best_val_loss, best_state, since_improve = val_loss, {k: v.clone() for k, v in model.state_dict().items()}, 0
        else:
            since_improve += 1
            if since_improve >= PATIENCE:
                break

    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        return model(X_test_t).numpy().ravel()


def main():
    dataset, X_full = _load_predictors()
    print(f"Using {X_full.shape[1]} predictors, {len(dataset)} samples (native resolution).")

    results = []
    fig, axes = plt.subplots(len(TARGET_TRAITS), 3, figsize=(12, 4 * len(TARGET_TRAITS)))
    for i, trait in enumerate(TARGET_TRAITS):
        y = dataset[trait].to_numpy()
        X_train, X_test, y_train_s, y_test, y_scaler = _train_test_split_scaled(X_full, y)
        n_features = X_train.shape[1]

        y_pred_tf = y_scaler.inverse_transform(train_tensorflow(X_train, y_train_s, X_test, n_features).reshape(-1, 1)).ravel()
        y_pred_pt = y_scaler.inverse_transform(train_pytorch(X_train, y_train_s, X_test, n_features).reshape(-1, 1)).ravel()
        y_pred_ens = 0.5 * (y_pred_tf + y_pred_pt)

        for name, y_pred, ax in [("TensorFlow", y_pred_tf, axes[i, 0]), ("PyTorch", y_pred_pt, axes[i, 1]),
                                  ("DL-Ensemble", y_pred_ens, axes[i, 2])]:
            m = all_metrics(y_test, y_pred)
            results.append({"trait": trait, "model": name, "n_test": len(y_test), **m})
            ax.scatter(y_test, y_pred, alpha=0.7, edgecolor="k", linewidth=0.3)
            lims = [min(y_test.min(), y_pred.min()), max(y_test.max(), y_pred.max())]
            ax.plot(lims, lims, "k--", linewidth=1)
            ax.set_title(f"{trait} - {name}\nR2={m['R2']:.2f}")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "3_inversion_dl_scatter.png", dpi=150)
    plt.close(fig)

    metrics_df = pd.DataFrame(results)
    metrics_df.to_csv(OUT_DIR / "3_inversion_dl_metrics.csv", index=False)
    print(metrics_df.to_string(index=False))
    print(f"\nMetrics -> {OUT_DIR / '3_inversion_dl_metrics.csv'}")


if __name__ == "__main__":
    main()
