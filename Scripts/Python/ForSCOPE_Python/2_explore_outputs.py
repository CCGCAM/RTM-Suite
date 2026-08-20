"""Plot the get_scope() output from 1_run_scope.py: TOC reflectance
components, fluorescence spectrum, and canopy temperature profile.

Python equivalent of Scripts/ForSCOPE/2-Explore_outputsSCOPE.R.
"""
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from scopeinpython import ScopeOptions, get_scope, get_spectra_scope

ROOT = Path(__file__).resolve().parents[3]
LUT_CSV = ROOT / "SCOPEinR" / "inst" / "input" / "LUT_input.csv"
OUT_DIR = ROOT / "outs" / "Python" / "ForSCOPE"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main():
    with open(LUT_CSV, newline="") as f:
        row = next(csv.DictReader(f))
    result = get_scope(row, options=ScopeOptions(calc_fluor=True, calc_xanthophyllabs=True))

    wl_optical = np.arange(400, 2401)  # the uniform-1nm part of rtmo's grid
    n = len(wl_optical)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    ax = axes[0, 0]
    ax.plot(wl_optical, result.rtmo.refl[:n], label="refl (apparent)", color="black")
    ax.plot(wl_optical, result.rtmo.rdd[:n], label="rdd (bi-hemispherical)", alpha=0.7)
    ax.plot(wl_optical, result.rtmo.rso[:n], label="rso (bidirectional)", alpha=0.7)
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Reflectance")
    ax.set_title("TOC reflectance components")
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    if result.rtmf is not None:
        # rtmf's LoF_/EoutF_ are on spectral.wlF (640-850nm, 1nm step, 211 pts).
        wl_f = get_spectra_scope().wlF
        ax.plot(wl_f, result.rtmf.LoF_, color="crimson")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Fluorescence radiance (mW m-2 nm-1 sr-1)")
        ax.set_title(f"TOC fluorescence spectrum (F685={result.rtmf.F685:.2f}, F740={result.rtmf.F740:.2f})")
    else:
        ax.text(0.5, 0.5, "calc_fluor=False", ha="center", va="center", transform=ax.transAxes)

    ax = axes[1, 0]
    layers = np.arange(1, result.nlayers + 1)
    ax.plot(layers, result.ebal.Tcu, label="Tcu (sunlit)", marker=".")
    ax.plot(layers, result.ebal.Tch, label="Tch (shaded)", marker=".")
    ax.axhline(result.ebal.Tsu, color="brown", linestyle="--", label="Tsu (soil, sunlit)")
    ax.axhline(result.ebal.Tsh, color="tan", linestyle="--", label="Tsh (soil, shaded)")
    ax.set_xlabel("Canopy layer (1 = top)")
    ax.set_ylabel("Temperature (degC)")
    ax.set_title(f"Converged temperature profile ({result.ebal.counter} iterations)")
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    fluxes = {"Rntot": result.ebal.Rntot, "lEtot": result.ebal.lEtot, "Htot": result.ebal.Htot,
              "Gtot": result.ebal.Gtot}
    ax.bar(fluxes.keys(), fluxes.values(), color=["darkorange", "steelblue", "firebrick", "grey"])
    ax.set_ylabel("W m-2")
    ax.set_title("Canopy energy-balance totals")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "2_scope_overview.png", dpi=150)
    plt.close(fig)
    print(f"Plot -> {OUT_DIR / '2_scope_overview.png'}")


if __name__ == "__main__":
    main()
