"""Compare TRADITIONAL inversion (LUT nearest-neighbour matching by a merit
function -- RMSE, FGE) against MACHINE-LEARNING inversion (Random Forest),
on the exact same held-out test set, convolved to real Sentinel-2A bands --
so all three approaches see literally the same observations, just invert
them differently:

  - Traditional (RMSE / FGE): toolsrtm.inversion.get_inversion_opt() ranks
    every TRAIN spectrum by how well it matches each TEST spectrum under the
    chosen merit function and averages the n_opt best matches' trait values.
    No model fitting at all -- pure LUT search, one call per merit function.
  - ML (Random Forest): toolsrtm.inversion.get_inversion() fits a model on
    TRAIN (band reflectance -> trait), then predicts on the SAME TEST set
    (via the returned fitted estimator, bypassing get_inversion's own
    internal random sub-split so both methods are scored on an identical
    held-out set).

Python equivalent of Scripts/R/Comparison/compare_inversion_traditional_vs_ML.R
-- same LUT size, same 70/30 split logic, same Sentinel-2A convolution,
directly comparable R2/RMSE.

Needs the `ml` extra: pip install "toolsrtm[ml]"
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from toolsrtm import foursail
from toolsrtm.srf import srf_sentinel2a, spectral_convolution_srf
from toolsrtm.inversion import get_inversion_opt, get_inversion

OUT_DIR = Path(__file__).resolve().parents[3] / "outs" / "Python" / "Comparison"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_SAMPLES = 300
WAVELENGTHS = np.arange(400, 2501)


def r2(obs, pred):
    obs, pred = np.asarray(obs, dtype=float), np.asarray(pred, dtype=float)
    ss_res = np.sum((obs - pred) ** 2)
    ss_tot = np.sum((obs - obs.mean()) ** 2)
    return float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")


def rmse(obs, pred):
    obs, pred = np.asarray(obs, dtype=float), np.asarray(pred, dtype=float)
    return float(np.sqrt(np.mean((pred - obs) ** 2)))


def sample_lut(n: int, rng: np.random.Generator) -> pd.DataFrame:
    return pd.DataFrame({
        "N": rng.uniform(1.2, 2.5, n),
        "Cab": rng.uniform(10, 70, n),
        "Car": rng.uniform(5, 15, n),
        "Anth": rng.uniform(0, 2, n),
        "Cbrown": rng.uniform(0, 0.3, n),
        "EWT": rng.uniform(0.005, 0.03, n),
        "LMA": rng.uniform(0.003, 0.015, n),
        "alpha": np.full(n, 40.0),
        "LIDFa": rng.uniform(-0.5, 0.5, n),
        "LIDFb": rng.uniform(-0.3, 0.3, n),
        "TypeLidf": np.full(n, 1, dtype=int),
        "LAI": rng.uniform(0.5, 7.0, n),
        "hspot": np.full(n, 0.01),
        "tts": np.full(n, 30.0),
        "tto": np.full(n, 0.0),
        "psi": np.full(n, 0.0),
    })


def main():
    rng = np.random.default_rng(1)

    # 1. Simulate a LUT (PROSPECT-D + fourSAIL) and convolve to Sentinel-2A ---
    print(f"Simulating {N_SAMPLES} PROSPECT-D + fourSAIL spectra ...")
    lut = sample_lut(N_SAMPLES, rng)
    rsoil = np.full(2101, 0.15)
    reflectance = np.empty((N_SAMPLES, len(WAVELENGTHS)))
    for i, row in lut.iterrows():
        sail = foursail(row.to_dict(), rsoil, leaf_model="PROSPECT-D", spectrum_all=True)
        reflectance[i] = sail.rsot

    print("Convolving to Sentinel-2A (13 bands) ...")
    srf = srf_sentinel2a()
    first = spectral_convolution_srf(WAVELENGTHS, reflectance[0], srf)
    band_names = first.band_names
    wl_bands = first.wl  # actual SRF-weighted band-center wavelengths
    se2a = np.empty((N_SAMPLES, len(band_names)))
    se2a[0] = first.rfl
    for i in range(1, N_SAMPLES):
        se2a[i] = spectral_convolution_srf(WAVELENGTHS, reflectance[i], srf).rfl
    se2a_df = pd.DataFrame(se2a, columns=band_names)

    # 2. Train/test split (70/30) -- SAME split feeds all three methods ---
    idx = rng.permutation(N_SAMPLES)
    n_train = round(0.7 * N_SAMPLES)
    train_idx, test_idx = idx[:n_train], idx[n_train:]
    print(f"Train: {len(train_idx)} spectra, Test (held out): {len(test_idx)} spectra")

    lut_train = lut.iloc[train_idx].reset_index(drop=True)
    lut_test = lut.iloc[test_idx].reset_index(drop=True)
    se2a_train = se2a[train_idx]
    se2a_test = se2a[test_idx]
    train_df = pd.concat([lut_train, se2a_df.iloc[train_idx].reset_index(drop=True)], axis=1)

    # 3. Traditional inversion -- LUT nearest-neighbour matching (RMSE, FGE) ---
    print("Traditional inversion (merit-RMSE) ...")
    opt_rmse = get_inversion_opt(se2a_test, se2a_train, lut_train, wave=wl_bands,
                                  method="merit-RMSE", n_opt=5)

    print("Traditional inversion (merit-FGE) ...")
    opt_fge = get_inversion_opt(se2a_test, se2a_train, lut_train, wave=wl_bands,
                                 method="merit-FGE", n_opt=5)

    # 4. ML inversion -- Random Forest, same train/test split ---
    print("ML inversion (Random Forest) ...")
    ml_res = get_inversion(data=train_df, dep_var="Cab", inputs=band_names, algorithm="RF",
                            seed=1, n_samples=None)
    pred_ml_test = ml_res.model.predict(se2a_test)

    # 5. Compare: R2 and RMSE on the held-out test set, for all three methods ---
    obs_cab = lut_test["Cab"].to_numpy()
    results = pd.DataFrame({
        "method": ["Traditional (merit-RMSE)", "Traditional (merit-FGE)", "ML (Random Forest)"],
        "R2": [r2(obs_cab, opt_rmse.lut_best["Cab"]), r2(obs_cab, opt_fge.lut_best["Cab"]), r2(obs_cab, pred_ml_test)],
        "RMSE": [rmse(obs_cab, opt_rmse.lut_best["Cab"]), rmse(obs_cab, opt_fge.lut_best["Cab"]), rmse(obs_cab, pred_ml_test)],
    })
    print(f"\n=== Cab retrieval on held-out Sentinel-2A test set (n={len(test_idx)}) ===")
    print(results.to_string(index=False))
    results.to_csv(OUT_DIR / "inversion_traditional_vs_ML_stats.csv", index=False)

    # 6. Plot: predicted vs observed Cab, all three methods ---
    preds = {
        "Traditional (merit-RMSE)": opt_rmse.lut_best["Cab"].to_numpy(),
        "Traditional (merit-FGE)": opt_fge.lut_best["Cab"].to_numpy(),
        "ML (Random Forest)": pred_ml_test,
    }
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.3), sharex=True, sharey=True)
    lims = (obs_cab.min() - 5, obs_cab.max() + 5)
    for ax, (name, pred) in zip(axes, preds.items()):
        ax.plot(lims, lims, "--", color="grey")
        ax.scatter(obs_cab, pred, alpha=0.7)
        ax.set_title(name, fontsize=10)
        ax.set_xlabel("Observed Cab (ug/cm2)")
        ax.set_xlim(lims); ax.set_ylim(lims)
    axes[0].set_ylabel("Predicted Cab (ug/cm2)")
    fig.suptitle(f"Traditional (LUT merit-function) vs ML inversion  --  "
                 f"PROSPECT-D + fourSAIL -> Sentinel-2A, n.train={len(train_idx)}, n.test={len(test_idx)}",
                 fontsize=10)
    fig.tight_layout()
    plot_path = OUT_DIR / "inversion_traditional_vs_ML.png"
    fig.savefig(plot_path, dpi=150)
    print(f"\nSaved plot to '{plot_path}' and stats to 'inversion_traditional_vs_ML_stats.csv'")


if __name__ == "__main__":
    main()
