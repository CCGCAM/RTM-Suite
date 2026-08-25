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

print("All figures written to", OUTDIR)
