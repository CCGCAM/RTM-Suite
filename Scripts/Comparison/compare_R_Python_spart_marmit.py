"""Cross-language robustness comparison: SPART and MARMIT (both versions),
ToolsRTM (R) vs toolsrtm (Python).

Reads the LUTs/reflectance written by compare_R_Python_spart_marmit.R,
runs the matching Python simulation for every row, reports RMSE/R^2/bias
per model, and produces comparison figures.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from toolsrtm import spart_toa, sentinel2a_msi, get_marmit_rsoil

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "_out")


def stats(R, Py):
    diff = Py - R
    rmse = float(np.sqrt(np.mean(diff ** 2)))
    max_abs = float(np.max(np.abs(diff)))
    bias = float(np.mean(diff))
    ss_res = np.sum(diff ** 2)
    ss_tot = np.sum((R - np.mean(R)) ** 2)
    r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else np.nan
    return dict(rmse=rmse, max_abs=max_abs, bias=bias, r2=r2)


# --- 1. SPART ---
lut = pd.read_csv(os.path.join(OUTDIR, "lut_SPART.csv"))
toc_R = pd.read_csv(os.path.join(OUTDIR, "refl_R_SPART_toc.csv"))
toa_R = pd.read_csv(os.path.join(OUTDIR, "refl_R_SPART_toa.csv"))
bands_df = pd.read_csv(os.path.join(OUTDIR, "SPART_bands.csv"))
band_cols = [c for c in toc_R.columns if c.startswith("band")]
n = len(lut)

sensor = sentinel2a_msi()
toc_py = np.full((n, len(band_cols)), np.nan)
toa_py = np.full((n, len(band_cols)), np.nan)
for i in range(n):
    row = lut.iloc[i]
    inputLUT = {k: row[k] for k in lut.columns}
    try:
        res = spart_toa(inputLUT, sensor=sensor, leaf_model="PROSPECT-PRO",
                         BSMBrightness=row["BSMBrightness"], BSMlat=row["BSMlat"],
                         BSMlon=row["BSMlon"], SMp=row["SMp"])
        if len(res.rfl_toc_brdf) == len(band_cols) and np.all(np.isfinite(res.rfl_toc_brdf)):
            toc_py[i, :] = res.rfl_toc_brdf
            toa_py[i, :] = res.rfl_toa
    except Exception:
        pass

R_toc = toc_R[band_cols].to_numpy()
R_toa = toa_R[band_cols].to_numpy()
ok = (~np.isnan(toc_py).any(axis=1)) & (~np.isnan(R_toc).any(axis=1)) & \
     (~np.isnan(toa_py).any(axis=1)) & (~np.isnan(R_toa).any(axis=1))
n_ok = int(ok.sum())
toc_stats = stats(R_toc[ok], toc_py[ok])
toa_stats = stats(R_toa[ok], toa_py[ok])
print(f"SPART TOC-BRDF: {n_ok}/{n} valid, R2={toc_stats['r2']:.8f} RMSE={toc_stats['rmse']:.3e}")
print(f"SPART TOA     : {n_ok}/{n} valid, R2={toa_stats['r2']:.8f} RMSE={toa_stats['rmse']:.3e}")

spart_summary = pd.DataFrame([
    dict(product="TOC-BRDF", n_valid=n_ok, **toc_stats),
    dict(product="TOA", n_valid=n_ok, **toa_stats),
])
spart_summary.to_csv(os.path.join(OUTDIR, "comparison_spart_summary.csv"), index=False)

wave = bands_df["wave_nm"].to_numpy()
fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
axes[0].scatter(R_toc[ok].mean(axis=1), toc_py[ok].mean(axis=1), s=10, alpha=0.6, color="#2166AC")
lims = [min(R_toc[ok].min(), toc_py[ok].min()), max(R_toc[ok].max(), toc_py[ok].max())]
axes[0].plot(lims, lims, "--", color="grey")
axes[0].set_xlabel("R (ToolsRTM) TOC BRDF"); axes[0].set_ylabel("Python (toolsrtm) TOC BRDF")
axes[0].set_title(f"TOC, mean/sim\nR2={toc_stats['r2']:.6f}")

axes[1].scatter(R_toa[ok].mean(axis=1), toa_py[ok].mean(axis=1), s=10, alpha=0.6, color="#B2182B")
lims2 = [min(R_toa[ok].min(), toa_py[ok].min()), max(R_toa[ok].max(), toa_py[ok].max())]
axes[1].plot(lims2, lims2, "--", color="grey")
axes[1].set_xlabel("R (ToolsRTM) TOA"); axes[1].set_ylabel("Python (toolsrtm) TOA")
axes[1].set_title(f"TOA, mean/sim\nR2={toa_stats['r2']:.6f}")

idx0 = np.where(ok)[0][0]
axes[2].plot(wave, R_toc[idx0], "o-", color="black", label="R TOC (ToolsRTM)")
axes[2].plot(wave, toc_py[idx0], "x--", color="#2166AC", label="Python TOC")
axes[2].plot(wave, R_toa[idx0], "o-", color="grey", label="R TOA (ToolsRTM)")
axes[2].plot(wave, toa_py[idx0], "x--", color="#B2182B", label="Python TOA")
axes[2].set_xlabel("Wavelength (nm)"); axes[2].set_ylabel("Reflectance")
axes[2].set_title(f"Example sim (row {idx0 + 1}), Sentinel-2A bands")
axes[2].legend(fontsize=7)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "comparison_spart.png"), dpi=150)
plt.close(fig)
print("Wrote comparison_spart.png")

# --- 2. MARMIT (both versions) ---
marmit_grid = pd.read_csv(os.path.join(OUTDIR, "lut_MARMIT.csv"))
wl_out = np.arange(400, 2501)
n_i, k_i, d_i = 1.53, 0.001, 0.0005

marmit_results = []
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
for ax, version in zip(axes, ["marmit1", "marmit2"]):
    refl_R = pd.read_csv(os.path.join(OUTDIR, f"refl_R_MARMIT_{version}_wet.csv"))
    wl_cols = [c for c in refl_R.columns if c.startswith("wl")]
    R_mat = refl_R[wl_cols].to_numpy()
    nrows = len(marmit_grid)
    py_mat = np.full((nrows, len(wl_cols)), np.nan)
    for i in range(nrows):
        row = marmit_grid.iloc[i]
        try:
            soil = get_marmit_rsoil(soil_id=int(row["soil_id"]), L=row["L"], eps=row["eps"],
                                     version=version, n_i=n_i, k_i=k_i, d_i=d_i, wl_out=wl_out)
            if len(soil.rsoil_wet) == len(wl_cols) and np.all(np.isfinite(soil.rsoil_wet)):
                py_mat[i, :] = soil.rsoil_wet
        except Exception:
            pass
    ok_m = (~np.isnan(py_mat).any(axis=1)) & (~np.isnan(R_mat).any(axis=1))
    n_ok_m = int(ok_m.sum())
    s = stats(R_mat[ok_m], py_mat[ok_m])
    print(f"MARMIT ({version}): {n_ok_m}/{nrows} valid, R2={s['r2']:.8f} RMSE={s['rmse']:.3e}")
    marmit_results.append(dict(version=version, n_valid=n_ok_m, **s))

    idx0 = np.where(ok_m)[0][0]
    ax.plot(wl_out, R_mat[idx0], color="black", linewidth=1.1, label="R (ToolsRTM)")
    ax.plot(wl_out, py_mat[idx0], color="#B2182B", linewidth=1, linestyle="--", label="Python (toolsrtm)")
    ax.set_title(f"{version}, example wet-soil spectrum (row {idx0+1})\nR2={s['r2']:.6f}, n={n_ok_m}/{nrows}")
    ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Wet-soil reflectance")
    ax.legend(fontsize=8)

marmit_summary = pd.DataFrame(marmit_results)
marmit_summary.to_csv(os.path.join(OUTDIR, "comparison_marmit_summary.csv"), index=False)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "comparison_marmit.png"), dpi=150)
plt.close(fig)
print("Wrote comparison_marmit.png")
print(marmit_summary.to_string(index=False))
