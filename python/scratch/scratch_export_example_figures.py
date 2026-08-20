"""One-off: generate the REAL output figures shown in docs/index.html's
"Code examples" section, from the exact Python code in each card (not
generic stock images). Run from the repo root:
    python python/scratch/scratch_export_example_figures.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from toolsrtm import inform

root = Path(__file__).resolve().parents[2]
outdir = root / "assets" / "examples"
outdir.mkdir(parents=True, exist_ok=True)

wl = np.arange(400, 2501)

# ---- 1. Simulate a forest stand spectrum (PROSPECT-PRO + INFORM) ----
row = dict(
    N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0.1, EWT=0.01, LMA=0.009, alpha=40,
    Prot=0, CBC=0,
    LAI=3, hspot=0.01, LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
    tts=30, tto=0, psi=0,
    LAIu=2, sd=500, cd=6, h=12, skyl=0.1,
)
rsoil = np.full(2101, 0.15)

reflectance = inform(row, rsoil, leaf_model="PROSPECT-PRO")

fig, ax = plt.subplots(figsize=(9, 5.6), dpi=130)
ax.plot(wl, reflectance, color="#0b3d59", lw=2)
ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel("Reflectance")
ax.set_ylim(0, reflectance.max() * 1.15)
ax.grid(color="#dddddd")
fig.tight_layout()
fig.savefig(outdir / "py_sim_inform.png")
plt.close(fig)
print("Saved py_sim_inform.png")

# ---- 2. Sensitivity to chlorophyll (Cab), INFORM ----
cab_values = [10, 25, 40, 55, 70]
colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(cab_values)))

fig, ax = plt.subplots(figsize=(9, 5.6), dpi=130)
for cab, color in zip(cab_values, colors):
    row_i = dict(row)
    row_i["Cab"] = cab
    refl_i = inform(row_i, rsoil, leaf_model="PROSPECT-PRO")
    ax.plot(wl, refl_i, color=color, lw=2, label=f"Cab={cab}")
ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel("Reflectance")
ax.set_ylim(0, 0.5)
ax.grid(color="#dddddd")
ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=9, frameon=False)
fig.tight_layout()
fig.savefig(outdir / "py_sensitivity_cab.png")
plt.close(fig)
print("Saved py_sensitivity_cab.png")
