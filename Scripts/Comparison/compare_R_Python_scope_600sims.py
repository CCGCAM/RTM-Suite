"""Cross-language robustness comparison: SCOPEinR (R) vs scopeinpython (Python).

Reads the shared 600-row LUT (and R's TOC reflectance) written by
compare_R_Python_scope_600sims.R, runs the same default (lite) SCOPE
simulation in Python for every row, and reports the R-vs-Python error
(RMSE, R^2, max abs diff) plus a comparison figure. Wavelengths that come
back NA from R (a structural, non-contiguous set matching atmospheric
water-vapor absorption windows) are excluded from the comparison rather
than dropping otherwise-valid rows.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scopeinpython import ScopeOptions, get_scope

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "_out")

lut = pd.read_csv(os.path.join(OUTDIR, "shared_LUT_scope_600.csv"))
refl_R_df = pd.read_csv(os.path.join(OUTDIR, "refl_R_scope_600.csv"))
wl_grid = pd.read_csv(os.path.join(OUTDIR, "scope_wl_grid.csv"))["wl"].to_numpy()
wl_cols = [c for c in refl_R_df.columns if c.startswith("wl")]
R_matrix = refl_R_df[wl_cols].to_numpy()
row_ok_R = refl_R_df["ok"].to_numpy().astype(bool)

n = len(lut)
refl_py = np.full((n, len(wl_grid)), np.nan)
Actot_py = np.full(n, np.nan)
options = ScopeOptions()

for i in range(n):
    if not row_ok_R[i]:
        continue
    row = lut.iloc[i].to_dict()
    try:
        result = get_scope(row, options=options)
        refl_i = np.asarray(result.rtmo.refl)
        if refl_i.shape[0] == len(wl_grid):
            refl_py[i, :] = refl_i
        Actot_py[i] = result.biochemical.Actot if hasattr(result, "biochemical") else np.nan
    except Exception:
        pass
    if (i + 1) % 100 == 0:
        print(f"Python side: {i + 1}/{n} done")

row_ok_py = ~np.isnan(refl_py).all(axis=1)
both_ok = row_ok_R & row_ok_py
print(f"Python side final: {row_ok_py.sum()}/{n} valid; {both_ok.sum()} rows comparable to R")

# mask wavelength columns that are NA in either language (the structural
# water-vapor-window gap, or any row-specific hole) rather than requiring
# full-row completeness
col_valid = ~np.isnan(R_matrix[both_ok]).any(axis=0) & ~np.isnan(refl_py[both_ok]).any(axis=0)
print(f"Comparable wavelengths: {col_valid.sum()}/{len(wl_grid)}")

R_ok = R_matrix[np.ix_(both_ok, col_valid)]
Py_ok = refl_py[np.ix_(both_ok, col_valid)]
wl_ok = wl_grid[col_valid]
diff = Py_ok - R_ok

rmse = float(np.sqrt(np.mean(diff ** 2)))
max_abs = float(np.max(np.abs(diff)))
ss_res = np.sum((Py_ok - R_ok) ** 2)
ss_tot = np.sum((R_ok - np.mean(R_ok)) ** 2)
r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
bias = float(np.mean(diff))

print(f"RMSE (R vs Python, {col_valid.sum()} wavelengths x {both_ok.sum()} sims): {rmse:.3e}")
print(f"Max abs diff: {max_abs:.3e}")
print(f"Bias (mean Python - R): {bias:.3e}")
print(f"R^2 (Python predicting R): {r2:.8f}")

mean_R = R_ok.mean(axis=1)
mean_Py = Py_ok.mean(axis=1)

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

axes[0].scatter(mean_R, mean_Py, s=8, alpha=0.5, color="#2E8B57")
lims = [min(mean_R.min(), mean_Py.min()), max(mean_R.max(), mean_Py.max())]
axes[0].plot(lims, lims, color="grey", linestyle="--", linewidth=1)
axes[0].set_xlabel("R (SCOPEinR) mean TOC reflectance")
axes[0].set_ylabel("Python (scopeinpython) mean TOC reflectance")
axes[0].set_title(f"{both_ok.sum()} SCOPE simulations, mean reflectance\nR2={r2:.6f}")

# wl_grid spans SCOPE's full radiative-transfer domain: 400-2400nm optical
# (1nm step) *plus* 2500-50000nm thermal (coarser step) -- data.rad$refl's
# thermal-domain entries aren't shortwave reflectance in the same physical
# sense, so restrict the illustrative example plot to the true optical
# range only (matching this page's own "TOC reflectance" framing/axis
# label), or a stray near-zero thermal point at the 2400/2500nm boundary
# reads as a fake cliff to zero. Also pick a row with a real, closed
# vegetation canopy (not just "first valid row") -- SCOPE's random LUT
# includes very sparse canopies (LAI down to ~0) whose TOC reflectance is
# soil-dominated and doesn't show the chlorophyll/red-edge features a
# reader expects from "an example SCOPE spectrum".
optical_mask = wl_ok <= 2400
lai_col = "LAI"
veg_candidates = np.where(both_ok)[0]
veg_row = None
for idx in veg_candidates:
    if lut.iloc[idx][lai_col] >= 2.0 and lut.iloc[idx]["Cab"] >= 30:
        veg_row = idx
        break
if veg_row is None:
    veg_row = veg_candidates[0]
example_idx = int(np.where(np.where(both_ok)[0] == veg_row)[0][0])

axes[1].plot(wl_ok[optical_mask], R_ok[example_idx][optical_mask], label="R (SCOPEinR)", color="black", linewidth=1.2)
axes[1].plot(wl_ok[optical_mask], Py_ok[example_idx][optical_mask], label="Python (scopeinpython)", color="#2E8B57", linewidth=1, linestyle="--")
axes[1].set_xlabel("Wavelength (nm)")
axes[1].set_ylabel("TOC reflectance")
axes[1].set_xlim(400, 2400)
axes[1].set_title(f"Example simulation (row {veg_row + 1}, LAI={lut.iloc[veg_row][lai_col]:.1f}, Cab={lut.iloc[veg_row]['Cab']:.0f}), optical range")
axes[1].legend()

plt.tight_layout()
fig_path = os.path.join(OUTDIR, "comparison_scope_R_vs_Python_600sims.png")
plt.savefig(fig_path, dpi=150)
print(f"Wrote {fig_path}")

summary = pd.DataFrame({
    "metric": ["n_simulations", "n_wavelengths_compared", "rmse", "max_abs_diff", "bias", "r2"],
    "value": [int(both_ok.sum()), int(col_valid.sum()), rmse, max_abs, bias, r2],
})
summary.to_csv(os.path.join(OUTDIR, "comparison_scope_summary.csv"), index=False)
print("Wrote comparison_scope_summary.csv")
