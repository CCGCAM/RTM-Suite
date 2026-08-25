"""Cross-language robustness comparison: ToolsRTM (R) vs toolsrtm (Python).

Reads the shared 600-row LUT (and R's reflectance) written by
compare_R_Python_toolsrtm_600sims.R, runs the same PROSPECT-D + fourSAIL
simulation in Python for every row, and reports the R-vs-Python error
(RMSE, R^2, max abs diff) plus a comparison figure.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from toolsrtm import prospect_d, foursail

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "_out")

lut = pd.read_csv(os.path.join(OUTDIR, "shared_LUT_600.csv"))
refl_R = pd.read_csv(os.path.join(OUTDIR, "refl_R_600.csv"))
wl_cols = [c for c in refl_R.columns if c.startswith("wl")]
wl = np.array([int(c[2:]) for c in wl_cols])
R_matrix = refl_R[wl_cols].to_numpy()

n = len(lut)
rsoil = np.full(2101, 0.15)
refl_py = np.full((n, len(wl)), np.nan)

for i in range(n):
    row = lut.iloc[i]
    inputLUT = dict(
        N=row["N"], Cab=row["Cab"], Car=row["Car"], Anth=row["Anth"], Cbrown=row["Cbrown"],
        EWT=row["EWT"], LMA=row["LMA"], alpha=row["alpha"],
        LIDFa=row["LIDFa"], LIDFb=row["LIDFb"], TypeLidf=row["TypeLidf"],
        LAI=row["LAI"], hspot=row["hspot"], tts=row["tts"], tto=row["tto"], psi=row["psi"],
    )
    try:
        sail = foursail(inputLUT, rsoil, leaf_model="PROSPECT-D", spectrum_all=True)
        refl_py[i, :] = sail.rsot
    except Exception:
        pass

ok = ~np.isnan(refl_py).any(axis=1) & ~np.isnan(R_matrix).any(axis=1)
print(f"Python side: {ok.sum()}/{n} simulations valid and comparable to R")

R_ok = R_matrix[ok]
Py_ok = refl_py[ok]
diff = Py_ok - R_ok

rmse = float(np.sqrt(np.mean(diff ** 2)))
max_abs = float(np.max(np.abs(diff)))
ss_res = np.sum((Py_ok - R_ok) ** 2)
ss_tot = np.sum((R_ok - np.mean(R_ok)) ** 2)
r2 = float(1 - ss_res / ss_tot)
bias = float(np.mean(diff))

print(f"RMSE (R vs Python, all wavelengths x {ok.sum()} sims): {rmse:.3e}")
print(f"Max abs diff: {max_abs:.3e}")
print(f"Bias (mean Python - R): {bias:.3e}")
print(f"R^2 (Python predicting R): {r2:.10f}")

# a per-simulation summary (mean reflectance across the spectrum) for the scatter plot
mean_R = R_ok.mean(axis=1)
mean_Py = Py_ok.mean(axis=1)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

axes[0].scatter(mean_R, mean_Py, s=8, alpha=0.5, color="#2166AC")
lims = [min(mean_R.min(), mean_Py.min()), max(mean_R.max(), mean_Py.max())]
axes[0].plot(lims, lims, color="grey", linestyle="--", linewidth=1)
axes[0].set_xlabel("R (ToolsRTM) mean reflectance")
axes[0].set_ylabel("Python (toolsrtm) mean reflectance")
axes[0].set_title(f"{ok.sum()} simulations, mean reflectance\nR2={r2:.8f}")

example_idx = np.where(ok)[0][0]
axes[1].plot(wl, R_ok[0], label="R (ToolsRTM)", color="black", linewidth=1.2)
axes[1].plot(wl, Py_ok[0], label="Python (toolsrtm)", color="#B2182B", linewidth=1, linestyle="--")
axes[1].set_xlabel("Wavelength (nm)")
axes[1].set_ylabel("Reflectance")
axes[1].set_title(f"Example simulation (row {example_idx + 1})")
axes[1].legend()

plt.tight_layout()
fig_path = os.path.join(OUTDIR, "comparison_toolsrtm_R_vs_Python_600sims.png")
plt.savefig(fig_path, dpi=150)
print(f"Wrote {fig_path}")

summary = pd.DataFrame({
    "metric": ["n_simulations", "rmse", "max_abs_diff", "bias", "r2"],
    "value": [int(ok.sum()), rmse, max_abs, bias, r2],
})
summary.to_csv(os.path.join(OUTDIR, "comparison_toolsrtm_summary.csv"), index=False)
print("Wrote comparison_toolsrtm_summary.csv")
