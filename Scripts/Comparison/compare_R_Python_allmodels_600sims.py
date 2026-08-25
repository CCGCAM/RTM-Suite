"""Cross-language robustness comparison: ToolsRTM (R) vs toolsrtm (Python),
ALL 5 leaf models x ALL 3 canopy models (15 combinations).

Reads each combo's shared LUT + R reflectance written by
compare_R_Python_allmodels_600sims.R, runs the matching Python simulation
for every row, and reports RMSE/R^2/bias/max-abs-diff per combination plus
a combined figure: one simulated-spectra overlay panel per combination
(R vs Python), and a summary heatmap of R^2 across the full leaf x canopy
grid.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from toolsrtm import foursail, foursail2, inform

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "_out")

LEAF_MODELS = ["PROSPECT-D", "PROSPECT-PRO", "Liberty", "Fluspect-B", "Fluspect-B-Cx"]
CANOPY_MODELS = ["fourSAIL", "foursail2", "INFORM"]

LEAF_EXTRA_COLS = {
    "Liberty": ["cell.d", "inter.c", "baseline.abs", "leaf.thick", "albino.abs", "lign.cell", "Nitrogen"],
    "Fluspect-B": ["Cs", "fqe", "Cx"],
    "Fluspect-B-Cx": ["Cs", "fqe", "Cx"],
}
BASE_COLS = ["N", "Cab", "Car", "Anth", "Cbrown", "EWT", "LMA", "alpha", "Prot", "CBC",
             "LIDFa", "LIDFb", "TypeLidf", "LAI", "hspot", "tts", "tto", "psi"]
CANOPY_EXTRA_COLS = {
    "foursail2": ["fraction_brown", "diss", "Cv", "Zeta"],
    "INFORM": ["LAIu", "sd", "cd", "h", "skyl"],
}


def build_input_row(row, leaf_model, canopy_model):
    cols = list(BASE_COLS) + LEAF_EXTRA_COLS.get(leaf_model, []) + CANOPY_EXTRA_COLS.get(canopy_model, [])
    return {c: row[c] for c in cols}


def run_combo(canopy_model, leaf_model):
    slug = f"{canopy_model.replace('foursail2', 'foursail2') if False else canopy_model}_{leaf_model}"
    slug = "".join(ch for ch in canopy_model if ch.isalnum()) + "_" + "".join(ch for ch in leaf_model if ch.isalnum())
    lut = pd.read_csv(os.path.join(OUTDIR, f"lut_{slug}.csv"))
    refl_R = pd.read_csv(os.path.join(OUTDIR, f"refl_R_{slug}.csv"))
    wl_cols = [c for c in refl_R.columns if c.startswith("wl")]
    wl = np.array([int(c[2:]) for c in wl_cols])
    R_matrix = refl_R[wl_cols].to_numpy()

    is_fluspect = leaf_model in ("Fluspect-B", "Fluspect-B-Cx")
    n = len(lut)
    rsoil_full = np.full(2101, 0.15)
    rsoil_short = np.full(2001, 0.15)
    rsoil_in = rsoil_full if canopy_model == "INFORM" else (rsoil_short if is_fluspect else rsoil_full)

    refl_py = np.full((n, len(wl)), np.nan)
    for i in range(n):
        row = lut.iloc[i]
        inputLUT = build_input_row(row, leaf_model, canopy_model)
        try:
            if canopy_model == "fourSAIL":
                out = foursail(inputLUT, rsoil_in, leaf_model=leaf_model, spectrum_all=not is_fluspect).rsot
            elif canopy_model == "foursail2":
                out = foursail2(inputLUT, rsoil_in, leaf_model=leaf_model, spectrum_all=not is_fluspect).rsot
            else:
                out = inform(inputLUT, rsoil_in, leaf_model=leaf_model)
            if len(out) == len(wl):
                refl_py[i, :] = out
        except Exception:
            pass

    ok = ~np.isnan(refl_py).any(axis=1) & ~np.isnan(R_matrix).any(axis=1)
    R_ok, Py_ok = R_matrix[ok], refl_py[ok]
    diff = Py_ok - R_ok
    n_ok = int(ok.sum())
    if n_ok == 0:
        return dict(canopy_model=canopy_model, leaf_model=leaf_model, n_valid=0,
                     rmse=np.nan, r2=np.nan, bias=np.nan, max_abs=np.nan,
                     wl=wl, R_ok=R_ok, Py_ok=Py_ok)

    rmse = float(np.sqrt(np.mean(diff ** 2)))
    max_abs = float(np.max(np.abs(diff)))
    ss_res = np.sum(diff ** 2)
    ss_tot = np.sum((R_ok - np.mean(R_ok)) ** 2)
    r2 = float(1 - ss_res / ss_tot) if ss_tot > 0 else np.nan
    bias = float(np.mean(diff))
    print(f"{leaf_model:14s} + {canopy_model:9s}: {n_ok}/{n} valid, "
          f"R2={r2:.8f}  RMSE={rmse:.3e}  bias={bias:.3e}  max|diff|={max_abs:.3e}")
    return dict(canopy_model=canopy_model, leaf_model=leaf_model, n_valid=n_ok,
                 rmse=rmse, r2=r2, bias=bias, max_abs=max_abs,
                 wl=wl, R_ok=R_ok, Py_ok=Py_ok)


results = []
for cm in CANOPY_MODELS:
    for lm in LEAF_MODELS:
        results.append(run_combo(cm, lm))

summary = pd.DataFrame([{k: v for k, v in r.items() if k not in ("wl", "R_ok", "Py_ok")} for r in results])
summary.to_csv(os.path.join(OUTDIR, "comparison_allmodels_summary.csv"), index=False)
print("\nWrote comparison_allmodels_summary.csv")
print(summary.to_string(index=False))

# --- Figure 1: R^2 heatmap across the full leaf x canopy grid ---
pivot = summary.pivot(index="canopy_model", columns="leaf_model", values="r2").reindex(
    index=CANOPY_MODELS, columns=LEAF_MODELS)
fig, ax = plt.subplots(figsize=(8, 4))
im = ax.imshow(pivot.values, cmap="RdYlGn", vmin=0.999, vmax=1.0, aspect="auto")
ax.set_xticks(range(len(LEAF_MODELS))); ax.set_xticklabels(LEAF_MODELS, rotation=30, ha="right")
ax.set_yticks(range(len(CANOPY_MODELS))); ax.set_yticklabels(CANOPY_MODELS)
for i in range(len(CANOPY_MODELS)):
    for j in range(len(LEAF_MODELS)):
        val = pivot.values[i, j]
        txt = f"{val:.6f}" if not np.isnan(val) else "n/a"
        ax.text(j, i, txt, ha="center", va="center", fontsize=8)
ax.set_title("R vs Python agreement (R^2) -- every leaf x canopy model combination")
fig.colorbar(im, ax=ax, label="R^2")
plt.tight_layout()
heatmap_path = os.path.join(OUTDIR, "comparison_allmodels_r2_heatmap.png")
plt.savefig(heatmap_path, dpi=150)
plt.close(fig)
print(f"Wrote {heatmap_path}")

# --- Figure 2: simulated-spectra overlay, one panel per combination (R vs Python) ---
fig, axes = plt.subplots(len(CANOPY_MODELS), len(LEAF_MODELS), figsize=(20, 10), sharex=True)
for r in results:
    i = CANOPY_MODELS.index(r["canopy_model"])
    j = LEAF_MODELS.index(r["leaf_model"])
    ax = axes[i, j]
    if r["n_valid"] > 0:
        ax.plot(r["wl"], r["R_ok"][0], color="black", linewidth=1.1, label="R (ToolsRTM)")
        ax.plot(r["wl"], r["Py_ok"][0], color="#B2182B", linewidth=1, linestyle="--", label="Python (toolsrtm)")
        ax.set_title(f"{r['leaf_model']}\n+ {r['canopy_model']}  (R2={r['r2']:.6f})", fontsize=8)
    else:
        ax.set_title(f"{r['leaf_model']}\n+ {r['canopy_model']}  (no valid sims)", fontsize=8)
    if i == len(CANOPY_MODELS) - 1:
        ax.set_xlabel("Wavelength (nm)", fontsize=8)
    if j == 0:
        ax.set_ylabel("TOC reflectance", fontsize=8)
    ax.tick_params(labelsize=7)
handles, labels = axes[0, 0].get_legend_handles_labels()
if handles:
    fig.legend(handles, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.02))
fig.suptitle("Simulated reflectance spectra: R (ToolsRTM) vs Python (toolsrtm), all 15 leaf x canopy combinations",
             y=1.05, fontsize=12)
plt.tight_layout()
spectra_path = os.path.join(OUTDIR, "comparison_allmodels_spectra_grid.png")
plt.savefig(spectra_path, dpi=140, bbox_inches="tight")
plt.close(fig)
print(f"Wrote {spectra_path}")
