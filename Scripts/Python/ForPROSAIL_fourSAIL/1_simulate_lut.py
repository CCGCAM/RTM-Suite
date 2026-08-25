"""Simulate a PROSPECT-D + fourSAIL LUT using the Python toolsrtm port.

Mirrors Scripts/ForPROSAIL/3-simulate_LUT.R (the R course pipeline), but
entirely in Python: samples 100 random leaf/canopy parameter sets, runs
toolsrtm.foursail for each, and saves the LUT + full-resolution (400-2500nm)
reflectance spectra to outs/Python/ForPROSAIL/.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from toolsrtm import foursail

OUT_DIR = Path(__file__).resolve().parents[3] / "outs" / "Python" / "ForPROSAIL"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_SAMPLES = 100
WAVELENGTHS = np.arange(400, 2501)

rng = np.random.default_rng(42)


def sample_lut(n: int) -> pd.DataFrame:
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
    lut = sample_lut(N_SAMPLES)
    rsoil = np.full(2101, 0.15)  # flat soil, matching the R course pipeline's simplification

    reflectance = np.empty((N_SAMPLES, len(WAVELENGTHS)))
    for i, row in lut.iterrows():
        sail = foursail(row.to_dict(), rsoil, leaf_model="PROSPECT-D", spectrum_all=True)
        reflectance[i] = sail.rsot

    lut.to_csv(OUT_DIR / "1_LUT.csv", index=False)
    refl_df = pd.DataFrame(reflectance, columns=[f"R.{wl}" for wl in WAVELENGTHS])
    refl_df.to_csv(OUT_DIR / "1_reflectance.csv", index=False)

    print(f"Simulated {N_SAMPLES} PROSPECT-D + fourSAIL spectra.")
    print(f"LUT -> {OUT_DIR / '1_LUT.csv'}")
    print(f"Reflectance -> {OUT_DIR / '1_reflectance.csv'}")


if __name__ == "__main__":
    main()
