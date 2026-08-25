"""One-off: generate the REAL output figures for the SCOPE code examples in
docs/index.html -- Python side runs 20 real get_scope() simulations sweeping
Cab, showing (1) the reflectance spectra and (2) how SIF (EoutF) responds to
Cab. Run from the repo root:
    python python/scratch/scratch_export_scope_figures.py
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scopeinpython import ScopeOptions, get_scope
from scopeinpython.spectral import get_spectra_scope

root = Path(__file__).resolve().parents[2]
outdir = root / "assets" / "examples"
outdir.mkdir(parents=True, exist_ok=True)

with open(root / "SCOPEinR/inst/input/LUT_input.csv", newline="") as f:
    row = next(csv.DictReader(f))

cab_values = np.linspace(10, 80, 20)
options = ScopeOptions(calc_fluor=True, calc_xanthophyllabs=True)

wl_refl = get_spectra_scope().wlS

refl_curves = []
sif_values = []
for cab in cab_values:
    row_i = dict(row)
    row_i["Cab"] = cab
    result = get_scope(row_i, options=options)
    refl_curves.append(result.rtmo.refl)
    sif_values.append(result.rtmf.EoutF)

optical = wl_refl <= 2400  # SCOPE's wlS also carries the thermal-region tail (up to 50000nm);
                           # reflectance is only meaningful over the optical range

fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), dpi=130)

colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(cab_values)))
for cab, refl, color in zip(cab_values, refl_curves, colors):
    axes[0].plot(wl_refl[optical], np.asarray(refl)[optical], color=color, lw=1.2)
axes[0].set_xlabel("Wavelength (nm)")
axes[0].set_ylabel("Reflectance")
axes[0].set_title("20 SCOPE simulations, Cab = 10–80")
sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(cab_values.min(), cab_values.max()))
cbar = fig.colorbar(sm, ax=axes[0])
cbar.set_label("Cab")

axes[1].plot(cab_values, sif_values, "o-", color="#0b3d59")
axes[1].set_xlabel("Cab")
axes[1].set_ylabel("SIF, EoutF (W m$^{-2}$)")
axes[1].set_title("Sun-induced fluorescence vs. Cab")
axes[1].grid(color="#dddddd")

fig.tight_layout()
fig.savefig(outdir / "py_scope_cab_sif.png")
plt.close(fig)
print("Saved py_scope_cab_sif.png")
print("SIF range:", min(sif_values), "-", max(sif_values))
