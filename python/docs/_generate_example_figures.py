"""Generate real figures for python/docs/examples.rst by actually running
each example's code. Run from repo root with python/toolsrtm and
python/scopeinpython importable (e.g. from within those venvs/paths)."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "_figures")
os.makedirs(OUTDIR, exist_ok=True)

palette = {"toolsrtm": "#2E8B57", "scope": "#2166AC", "s2": "#B2182B", "accent": "#6A3D9A"}

# ---------------------------------------------------------------- 1. leaf + canopy
from toolsrtm import prospect_d, foursail

leaf = prospect_d(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40)
inputLUT = dict(
    N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
    Prot=0.002, CBC=0.007, LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
    LAI=3, hspot=0.01, tts=30, tto=0, psi=0,
)
rsoil = np.full(2101, 0.15)
sail = foursail(inputLUT, rsoil, leaf_model="PROSPECT-D", spectrum_all=True)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(leaf.lambda_, leaf.refl, color=palette["toolsrtm"], label="reflectance")
axes[0].plot(leaf.lambda_, 1 - leaf.tran, color="grey", linestyle="--", label="1 - transmittance")
axes[0].set_xlabel("Wavelength (nm)"); axes[0].set_ylabel("Leaf refl / 1-tran"); axes[0].set_title("prospect_d() leaf optics")
axes[0].legend(fontsize=8)
axes[1].plot(np.arange(400, 2501), sail.rsot, color=palette["toolsrtm"])
axes[1].set_xlabel("Wavelength (nm)"); axes[1].set_ylabel("TOC reflectance"); axes[1].set_title("foursail() canopy BRDF")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "leaf_canopy.png"), dpi=140)
plt.close(fig)
print("wrote leaf_canopy.png")

# ---------------------------------------------------------------- 3. sensor convolution + indices
from toolsrtm.srf import srf_sentinel2a, spectral_convolution_srf
from toolsrtm.indices import get_indices

wave = np.arange(400, 2501)
s2a = srf_sentinel2a()
conv = spectral_convolution_srf(wave, sail.rsot, s2a)
indices = get_indices(conv.wl, conv.rfl, spectral_domain="VNIR")

fig, ax = plt.subplots(figsize=(8, 4.2))
ax.plot(wave, sail.rsot, color="grey", label="Native 1nm", linewidth=1)
ax.scatter(conv.wl, conv.rfl, color=palette["s2"], zorder=5, label="Sentinel-2A bands", s=35)
ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Reflectance")
ax.set_title(f"Sensor convolution -- NDVI={float(indices['NDVI'][0]):.3f}")
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "sensor_convolution.png"), dpi=140)
plt.close(fig)
print("wrote sensor_convolution.png")

# ---------------------------------------------------------------- 4. ML inversion
import pandas as pd
from toolsrtm.inversion import get_inversion

rng = np.random.default_rng(1)
rows = []
for _ in range(200):
    Cab, LAI = rng.uniform(10, 80), rng.uniform(0.5, 6)
    row_lut = dict(
        N=1.5, Cab=Cab, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
        LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
        LAI=LAI, hspot=0.01, tts=30, tto=0, psi=0,
    )
    sail_i = foursail(row_lut, np.full(2101, 0.15), leaf_model="PROSPECT-D", spectrum_all=True)
    row = {"Cab": Cab, "LAI": LAI}
    for wl in (490, 560, 665, 705, 740, 783, 842, 865, 1610, 2190):
        row[f"R{wl}"] = sail_i.rsot[wl - 400]
    rows.append(row)

df = pd.DataFrame(rows)
band_cols = [c for c in df.columns if c.startswith("R")]
result = get_inversion(df, dep_var="Cab", inputs=band_cols, algorithm="PLSR", n_samples=200, seed=1)

# get_inversion()'s predictions/statistics are keyed by split ("train"/"test")
# holding the predicted array directly; reproduce the same
# train_test_split(seed=1, test_size=0.3 default) externally to recover the
# matching observed y_test for the scatter plot.
from sklearn.model_selection import train_test_split as _tts
_X = df[band_cols].to_numpy(dtype=float)
_y = df["Cab"].to_numpy(dtype=float)
_, _, _, obs = _tts(_X, _y, test_size=0.3, random_state=1)
pred = result.predictions["test"]

fig, ax = plt.subplots(figsize=(5, 5))
ax.scatter(obs, pred, s=18, alpha=0.6, color=palette["toolsrtm"])
lims = [min(obs.min(), pred.min()), max(obs.max(), pred.max())]
ax.plot(lims, lims, color="grey", linestyle="--", linewidth=1)
ax.set_xlabel("Observed Cab"); ax.set_ylabel("Predicted Cab (PLSR)")
ax.set_title(f"R2={result.statistics['test']['r2']:.3f}  RMSE={result.statistics['test']['rmse']:.2f}")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "ml_inversion.png"), dpi=140)
plt.close(fig)
print("wrote ml_inversion.png")

# ---------------------------------------------------------------- 5. MARMIT
from toolsrtm.marmit import get_marmit_rsoil

soil = get_marmit_rsoil(soil_id=3, L=0.05, eps=0.4, version="marmit1")
wl_soil = np.arange(400, 2501)

fig, ax = plt.subplots(figsize=(8, 4.2))
ax.plot(wl_soil, soil.rsoil_dry, color="#8B5A2B", label="Dry soil")
ax.plot(wl_soil, soil.rsoil_wet, color=palette["scope"], label=f"Wet soil (SMC={float(soil.smc):.3f})")
ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Soil reflectance")
ax.set_title("MARMIT: dry vs. wet soil reflectance")
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "marmit_soil.png"), dpi=140)
plt.close(fig)
print("wrote marmit_soil.png")

# ---------------------------------------------------------------- 6. full SCOPE
import csv
from scopeinpython import ScopeOptions, get_scope

with open("SCOPEinR/inst/input/LUT_input.csv", newline="") as f:
    row_scope = next(csv.DictReader(f))

res = get_scope(row_scope, options=ScopeOptions(k_maxit=100, maxEBer=1.0))
wl_scope = np.arange(400, 400 + len(res.rtmo.refl))

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(wl_scope, res.rtmo.refl, color=palette["scope"])
axes[0].set_xlim(400, 2500)
axes[0].set_xlabel("Wavelength (nm)"); axes[0].set_ylabel("TOC reflectance")
axes[0].set_title(f"get_scope() -- Actot={float(res.ebal.Actot):.2f} umol/m2/s")
if res.rtmf is not None:
    wl_f = np.arange(640, 640 + len(res.rtmf.EoutF_)) if hasattr(res.rtmf, "EoutF_") else None
    if wl_f is not None:
        axes[1].plot(wl_f, res.rtmf.EoutF_, color=palette["accent"])
        axes[1].set_xlabel("Wavelength (nm)"); axes[1].set_ylabel("Fluorescence (W/m2/sr/nm)")
        axes[1].set_title(f"Emitted SIF spectrum -- EoutF={float(res.rtmf.EoutF):.4f} W/m2/sr")
    else:
        axes[1].axis("off")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "scope_full.png"), dpi=140)
plt.close(fig)
print("wrote scope_full.png")

# ---------------------------------------------------------------- 7. Liberty (conifer leaf model)
from toolsrtm import liberty

lib = liberty(cell_d=40, inter_c=0.045, baseline_abs=0.0006, leaf_thick=1.6,
              albino_abs=0, Cab=40, EWT=0.01, lign_cell=2, Nitrogen=1)

fig, ax = plt.subplots(figsize=(8, 4.2))
ax.plot(lib.lambda_, lib.refl, color=palette["toolsrtm"], label="Reflectance")
ax.plot(lib.lambda_, 1 - lib.tran, color="grey", linestyle="--", label="1 - Transmittance")
ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Leaf refl / 1-tran")
ax.set_title("liberty() -- conifer-needle leaf optics (Dawson et al. 1998)")
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "liberty_leaf.png"), dpi=140)
plt.close(fig)
print("wrote liberty_leaf.png")

# ---------------------------------------------------------------- 8. Fluspect-B (leaf optics + fluorescence excitation)
from toolsrtm.fluspect import fluspect_b

flu = fluspect_b(Cab=40, Car=8, EWT=0.01, LMA=0.009, Cs=0, N=1.5, fqe=0.01, Cx=0)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
axes[0].plot(flu.lambda_, flu.refl, color=palette["toolsrtm"], label="Reflectance")
axes[0].plot(flu.lambda_, 1 - flu.tran, color="grey", linestyle="--", label="1 - Transmittance")
axes[0].set_xlabel("Wavelength (nm)"); axes[0].set_ylabel("Leaf refl / 1-tran")
axes[0].set_title("fluspect_b() leaf optics")
axes[0].legend(fontsize=8)
Mb_total = flu.MbI + flu.MbII  # (211, 351) = (emission 640-850nm, excitation 400-750nm)
im = axes[1].imshow(Mb_total, aspect="auto", origin="lower", cmap="inferno",
                     extent=[400, 750, 640, 850])
axes[1].set_xlabel("Excitation wavelength (nm)"); axes[1].set_ylabel("Emission wavelength (nm)")
axes[1].set_title("Fluorescence excitation-emission matrix (backward, PSI+PSII)")
fig.colorbar(im, ax=axes[1], label="Fluorescence")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "fluspect_leaf.png"), dpi=140)
plt.close(fig)
print("wrote fluspect_leaf.png")

# ---------------------------------------------------------------- 9. INFORM (forest canopy)
from toolsrtm import inform

inputLUT_inform = dict(
    N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
    Prot=0.002, CBC=0.007, LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
    LAI=3, hspot=0.01, tts=30, tto=0, psi=0,
    LAIu=0.5, sd=650, cd=4.5, h=20, skyl=0.1,
)
rsoil_inform = np.full(2101, 0.15)
r_forest = inform(inputLUT_inform, rsoil_inform, leaf_model="PROSPECT-D")
sail_same = foursail(inputLUT_inform, rsoil_inform, leaf_model="PROSPECT-D", spectrum_all=True)

fig, ax = plt.subplots(figsize=(8, 4.2))
ax.plot(np.arange(400, 2501), sail_same.rsot, color="grey", label="fourSAIL (single-layer canopy)")
ax.plot(np.arange(400, 400 + len(r_forest)), r_forest, color=palette["toolsrtm"], label="INFORM (forest stand)")
ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("TOC reflectance")
ax.set_title("inform() vs. foursail() -- same leaf/LAI, explicit tree-crown geometry")
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "inform_forest.png"), dpi=140)
plt.close(fig)
print("wrote inform_forest.png")

# ---------------------------------------------------------------- 10. SPART (soil-plant-atmosphere, sensor bands)
from toolsrtm import spart_toa, sentinel2a_msi

inputLUT_spart = dict(
    N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
    Prot=0.002, CBC=0.007, LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
    LAI=3, hspot=0.01, tts=30, tto=0, psi=0,
    Pa=1000, aot550=0.3246, uo3=0.3480, uh2o=1.4116,
)
spart_res = spart_toa(inputLUT_spart, sensor=sentinel2a_msi(), leaf_model="PROSPECT-PRO",
                       BSMBrightness=0.5, BSMlat=25, BSMlon=45, SMp=15)

fig, ax = plt.subplots(figsize=(8, 4.2))
ax.plot(spart_res.wl_smac, spart_res.rfl_toc_brdf, "o-", color=palette["toolsrtm"], label="TOC (canopy BRDF)")
ax.plot(spart_res.wl_smac, spart_res.rfl_toa, "o-", color=palette["s2"], label="TOA (after SMAC)")
ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Reflectance")
ax.set_title("spart_toa() -- Sentinel-2A bands, TOC vs. TOA")
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "spart_toc_toa.png"), dpi=140)
plt.close(fig)
print("wrote spart_toc_toa.png")

# ---------------------------------------------------------------- 11. Deep-learning inversion (Keras)
from toolsrtm.deep_learning import get_ml_model

rng2 = np.random.default_rng(2)
rows_dl = []
for _ in range(600):
    Cab, LAI = rng2.uniform(10, 80), rng2.uniform(0.5, 6)
    row_lut = dict(
        N=1.5, Cab=Cab, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
        LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
        LAI=LAI, hspot=0.01, tts=30, tto=0, psi=0,
    )
    sail_i = foursail(row_lut, np.full(2101, 0.15), leaf_model="PROSPECT-D", spectrum_all=True)
    row = {"Cab": Cab, "LAI": LAI}
    for wl in (490, 560, 665, 705, 740, 783, 842, 865, 1610, 2190):
        row[f"R{wl}"] = sail_i.rsot[wl - 400]
    rows_dl.append(row)

df_dl = pd.DataFrame(rows_dl)
band_cols_dl = [c for c in df_dl.columns if c.startswith("R")]
# adam at R's own conservative default learning rate (1e-4) converges slowly --
# generous epoch budget + multiple restarts, matching what the test suite
# itself needs to reach a decent held-out R^2 (see test_deep_learning.py).
dl_result = get_ml_model(df_dl, dep_var="Cab", model="Hidden-layers",
                          n_epochs=500, n_times=3, seed=2, verbose=0)

fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
axes[0].plot(dl_result.history["loss"], color=palette["toolsrtm"], label="train")
axes[0].plot(dl_result.history["val_loss"], color=palette["s2"], label="val")
axes[0].set_xlabel("Epoch"); axes[0].set_ylabel("Loss (MSE)"); axes[0].set_title("Keras dense-network training")
axes[0].legend()
obs_dl, pred_dl = dl_result.predictions["y_true"], dl_result.predictions["y_pred"]
axes[1].scatter(obs_dl, pred_dl, s=14, alpha=0.6, color=palette["accent"])
lims_dl = [min(obs_dl.min(), pred_dl.min()), max(obs_dl.max(), pred_dl.max())]
axes[1].plot(lims_dl, lims_dl, color="grey", linestyle="--", linewidth=1)
axes[1].set_xlabel("Observed Cab"); axes[1].set_ylabel("Predicted Cab (Keras)")
axes[1].set_title(f"R2={dl_result.stats['r2']:.3f}  RMSE={dl_result.stats['rmse']:.2f}")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "deep_learning_inversion.png"), dpi=140)
plt.close(fig)
print("wrote deep_learning_inversion.png")

# ---------------------------------------------------------------- 12. Global sensitivity analysis (Johnson index)
from toolsrtm.sensitivity import spectral_sensitivity

sens = spectral_sensitivity(n_samples=500, distribution="Uniform",
                             traits=("N", "Cab", "EWT", "LMA", "LIDFa", "LAI"), wl_step=5, seed=11)
trait_order = ["N", "Cab", "EWT", "LMA", "LIDFa", "LAI", "SoilCoef"]
trait_colors = {"N": "#999999", "Cab": "#2E8B57", "EWT": "#2166AC", "LMA": "#8B5A2B",
                "LIDFa": "#B2182B", "LAI": "#6A3D9A", "SoilCoef": "#E69F00"}
wl_unique = np.unique(sens.wavelength)
stack = np.zeros((len(trait_order), len(wl_unique)))
for ti, t in enumerate(trait_order):
    mask = sens.trait == t
    order = np.argsort(sens.wavelength[mask])
    stack[ti, :] = sens.sti_pct[mask][order]

fig, ax = plt.subplots(figsize=(9, 4.5))
ax.stackplot(wl_unique, stack, labels=trait_order, colors=[trait_colors[t] for t in trait_order])
ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Relative contribution (%)")
ax.set_title("Global sensitivity analysis (Johnson index) -- foursail() + PROSPECT-D, 500 simulations")
ax.set_xlim(400, 2500); ax.set_ylim(0, 100)
ax.legend(loc="upper center", ncol=7, fontsize=7, bbox_to_anchor=(0.5, -0.15))
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "spectral_sensitivity.png"), dpi=140, bbox_inches="tight")
plt.close(fig)
print("wrote spectral_sensitivity.png")

# ---------------------------------------------------------------- 13. Real Sentinel-2 capstone (Loobos, t18 port)
from toolsrtm.satellite import get_satellite_collection, get_sentinel2_cube

rng2 = np.random.default_rng(1)
n_cap = 500
LAI_c, tts_c, soil_c = rng2.uniform(0.3, 5, n_cap), rng2.uniform(25, 45, n_cap), rng2.uniform(0.05, 0.30, n_cap)
Cab_c, Car_c, Anth_c = rng2.uniform(5, 75, n_cap), rng2.uniform(0, 20, n_cap), rng2.uniform(0, 4.5, n_cap)
EWT_c, LMA_c, N_c = rng2.uniform(0.001, 0.035, n_cap), rng2.uniform(0.001, 0.035, n_cap), rng2.uniform(1.5, 2.5, n_cap)
LIDFa_c = rng2.uniform(30, 70, n_cap)
hspot_c, tto_c, psi_c = rng2.uniform(0, 1, n_cap), rng2.uniform(15, 30, n_cap), rng2.uniform(0, 180, n_cap)

refl_c = np.stack([
    foursail(dict(N=N_c[i], Cab=Cab_c[i], Car=Car_c[i], Anth=Anth_c[i], Cbrown=0.0, EWT=EWT_c[i], LMA=LMA_c[i],
                  alpha=40.0, LIDFa=LIDFa_c[i], LIDFb=0.0, TypeLidf=1.0, LAI=LAI_c[i], hspot=hspot_c[i],
                  tts=tts_c[i], tto=tto_c[i], psi=psi_c[i]),
             np.full(2101, soil_c[i]), leaf_model="PROSPECT-D", spectrum_all=True).rsot
    for i in range(n_cap)
])

s2a_cap = srf_sentinel2a()
conv0_cap = spectral_convolution_srf(wave, refl_c[0], s2a_cap)
keep_cap = ["B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B11", "B12"]
real_names_cap = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]
keep_idx_cap = [conv0_cap.band_names.index(k) for k in keep_cap]
band_refl_cap = np.stack([spectral_convolution_srf(wave, refl_c[i], s2a_cap).rfl[keep_idx_cap] for i in range(n_cap)])

df_cap = pd.DataFrame(band_refl_cap, columns=real_names_cap); df_cap["Cab"] = Cab_c
fit_cap = get_inversion(df_cap, dep_var="Cab", inputs=real_names_cap, algorithm="RF", n_samples=n_cap, seed=42)
print("Capstone Cab R2:", fit_cap.statistics["test"]["r2"])

lat_cap, lon_cap, d_cap = 52.166447, 5.74355, 0.006
bbox_cap = (lon_cap - d_cap, lat_cap - d_cap, lon_cap + d_cap, lat_cap + d_cap)
coll_cap = get_satellite_collection(bbox_cap, collection="sentinel-2-l2a", date_range=("2024-07-01", "2024-07-31"),
                                     cloud_server="microsoft", n_limit=20, cloud_threshold=40)
cube_cap = get_sentinel2_cube(coll_cap, bbox_cap, resolution=10.0, crs="EPSG:32631", aggregation_method="mean")
r_cap = {b: cube_cap[b].values.astype(float) / 10000 for b in real_names_cap}

rep_map = 700 + 40 * (((r_cap["B04"] + r_cap["B07"]) / 2 - r_cap["B05"]) / (r_cap["B06"] - r_cap["B05"]))
pix_cap = np.column_stack([r_cap[b].ravel() for b in real_names_cap])
ok_cap = np.all(np.isfinite(pix_cap), axis=1)
cab_pixels_cap = np.full(pix_cap.shape[0], np.nan)
cab_pixels_cap[ok_cap] = fit_cap.model.predict(pix_cap[ok_cap])
cab_map = cab_pixels_cap.reshape(r_cap["B04"].shape)

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
im0 = axes[0].imshow(rep_map, cmap="viridis")
axes[0].set_title("REP (winning index, red-edge position)")
fig.colorbar(im0, ax=axes[0], fraction=0.046, label="nm")
im1 = axes[1].imshow(cab_map, cmap="viridis")
axes[1].set_title("Retrieved Cab")
fig.colorbar(im1, ax=axes[1], fraction=0.046, label="ug/cm2")
fig.suptitle("Loobos forest (NL-Loo), real Sentinel-2, July 2024")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "t18_python_capstone.png"), dpi=140)
plt.close(fig)
print("wrote t18_python_capstone.png")

# ---------------------------------------------------------------- 14. Real Sentinel-2 capstone (Speulderbos, SCOPEinR t11 port)
import csv as _csv
from scopeinpython import ScopeOptions, get_scope
from toolsrtm.sensitivity import get_cor as _get_cor, gauss_by_min_max as _gauss_by_min_max
from sklearn.ensemble import RandomForestRegressor


def _build_scope_lut(csv_path, n, seed):
    rng = np.random.default_rng(seed)
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        rows = list(_csv.DictReader(f))
    lut = {}
    for row in rows:
        trait, dist = row["variable"], row["Distribution"]
        if trait in ("startDate", "endDate"):
            continue
        if trait == "Type":
            lut[trait] = np.array([f"C{row['default']}"] * n, dtype=object)
            continue
        lo, hi = float(row["lower"]), float(row["upper"])
        if dist == "Uniform":
            lut[trait] = rng.uniform(lo, hi, size=n)
        elif dist == "Fixed":
            lut[trait] = np.full(n, float(row["default"]))
        else:
            lut[trait] = _gauss_by_min_max(n, float(row["Mean_D"]), float(row["Std_D"]), lo, hi, n * 3, rng=rng)
    return lut


n_samples_t11 = 250
lut_t11 = _build_scope_lut("SCOPEinR/inst/input/inputs_SCOPE.csv", n_samples_t11, seed=1)
cor_res_t11 = _get_cor(n_inputs=2, n_lut=n_samples_t11, distribution="Uniform", rho=0.85, seed=3,
                        var_names=["Cab", "Vcmax25"], min_range=[5, 5], max_range=[90, 250])
lut_t11["Cab"], lut_t11["Vcmax25"] = cor_res_t11.lut["Cab"], cor_res_t11.lut["Vcmax25"]

opts_t11 = ScopeOptions(k_maxit=100, maxEBer=1.0)
wl_optical_t11 = np.arange(400, 2401)
s2a_t11 = srf_sentinel2a()
real_names_t11 = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]
keep_t11 = ["B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B11", "B12"]

Actot_t11, band_refl_t11 = [], []
for i in range(n_samples_t11):
    row = {k: v[i] for k, v in lut_t11.items()}
    res = get_scope(row, options=opts_t11)
    refl = np.asarray(res.rtmo.refl)[: len(wl_optical_t11)]
    bad = ~np.isfinite(refl)
    if bad.any():
        refl[bad] = np.interp(wl_optical_t11[bad], wl_optical_t11[~bad], refl[~bad])
    conv = spectral_convolution_srf(wl_optical_t11, refl, s2a_t11)
    keep_idx = [conv.band_names.index(k) for k in keep_t11]
    Actot_t11.append(res.ebal.Actot)
    band_refl_t11.append(conv.rfl[keep_idx])
Actot_t11, band_refl_t11 = np.array(Actot_t11), np.array(band_refl_t11)
print("t11 capstone: SCOPE sims done, Actot range", Actot_t11.min(), Actot_t11.max())

df_t11 = pd.DataFrame(band_refl_t11, columns=real_names_t11)
rf_reflonly_t11 = RandomForestRegressor(n_estimators=300, random_state=1)
rf_reflonly_t11.fit(df_t11[real_names_t11], Actot_t11)

lat_t11, lon_t11, d_t11 = 52.2500, 5.6900, 0.003
bbox_t11 = (lon_t11 - d_t11, lat_t11 - d_t11, lon_t11 + d_t11, lat_t11 + d_t11)
windows_t11 = [("2024-03-01", "2024-03-31"), ("2024-05-01", "2024-05-31"),
               ("2024-07-01", "2024-07-31"), ("2024-09-01", "2024-09-30"), ("2024-11-01", "2024-11-30")]
ts_t11, cubes_t11 = [], {}
for w in windows_t11:
    coll = get_satellite_collection(bbox_t11, collection="sentinel-2-l2a", date_range=w,
                                     cloud_server="microsoft", n_limit=20, cloud_threshold=40)
    ds = get_sentinel2_cube(coll, bbox_t11, resolution=10.0, crs="EPSG:32631", aggregation_method="mean")
    cubes_t11[w[0]] = ds
    means = np.array([float(np.nanmean(ds[b].values)) / 10000 for b in real_names_t11])
    ndvi = (means[6] - means[2]) / (means[6] + means[2])
    actot_pred = float(rf_reflonly_t11.predict(means.reshape(1, -1))[0])
    ts_t11.append(dict(date=w[0], ndvi=ndvi, actot=actot_pred))
ts_df_t11 = pd.DataFrame(ts_t11)
print("t11 capstone: NDVI/Actot correlation:", np.corrcoef(ts_df_t11.ndvi, ts_df_t11.actot)[0, 1])

map_cube_t11 = cubes_t11["2024-07-01"]
r_t11 = {b: map_cube_t11[b].values.astype(float) / 10000 for b in real_names_t11}
ndvi_map_t11 = (r_t11["B08"] - r_t11["B04"]) / (r_t11["B08"] + r_t11["B04"])
pix_t11 = np.column_stack([r_t11[b].ravel() for b in real_names_t11])
ok_t11 = np.all(np.isfinite(pix_t11), axis=1)
actot_pixels_t11 = np.full(pix_t11.shape[0], np.nan)
actot_pixels_t11[ok_t11] = rf_reflonly_t11.predict(pix_t11[ok_t11])
actot_map_t11 = actot_pixels_t11.reshape(r_t11["B04"].shape)

rgb_t11 = np.stack([r_t11["B04"], r_t11["B03"], r_t11["B02"]], axis=-1)
lo_t11, hi_t11 = np.nanpercentile(rgb_t11, 2), np.nanpercentile(rgb_t11, 98)
rgb_scaled_t11 = np.clip((rgb_t11 - lo_t11) / (hi_t11 - lo_t11), 0, 1)

fig = plt.figure(figsize=(13, 9))
ax1 = fig.add_subplot(2, 3, 1); ax1.imshow(rgb_scaled_t11); ax1.set_title("True color -- 2024-07-01")
ax2 = fig.add_subplot(2, 3, 2)
im2 = ax2.imshow(ndvi_map_t11, cmap="RdYlGn"); ax2.set_title("NDVI")
fig.colorbar(im2, ax=ax2, fraction=0.046)
ax3 = fig.add_subplot(2, 3, 3)
im3 = ax3.imshow(actot_map_t11, cmap="viridis"); ax3.set_title("Retrieved Actot")
fig.colorbar(im3, ax=ax3, fraction=0.046)
ax4 = fig.add_subplot(2, 1, 2)
dates = pd.to_datetime(ts_df_t11.date)
ax4b = ax4.twinx()
ax4.plot(dates, ts_df_t11.ndvi, "o-", color="#2E8B57", label="NDVI")
ax4b.plot(dates, ts_df_t11.actot, "o-", color="#B2182B", label="Actot")
ax4.set_ylabel("NDVI", color="#2E8B57"); ax4b.set_ylabel("Actot (umol/m2/s)", color="#B2182B")
ax4.set_title("Speulderbos, 2024: NDVI vs. retrieved net photosynthesis (Actot)")
fig.suptitle("Speulderbos forest (real Sentinel-2, STAC)", y=1.0)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "t11_python_capstone.png"), dpi=140, bbox_inches="tight")
plt.close(fig)
print("wrote t11_python_capstone.png")

print("All figures written to", OUTDIR)
