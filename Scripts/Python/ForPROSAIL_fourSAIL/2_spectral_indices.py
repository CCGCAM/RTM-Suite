"""Compute spectral vegetation indices from the simulated reflectance.

Uses toolsrtm.get_indices (Python port of ToolsRTM::getIndices) on the
output of 1_simulate_lut.py. Indices are computed directly from the
full-resolution simulated spectra, which is what the ML inversion step
uses as predictors -- for sensor-band convolution (Sentinel-2/PRISMA/
EnMAP/your own sensor) see 5_sensor_convolution.py, or
Scripts/Python/ForFoursail2/, ForINFORM/, ForSPART/, ForMARMIT/ for course
pipelines that convolve before computing indices, matching their R
counterparts exactly.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from toolsrtm import get_indices

OUT_DIR = Path(__file__).resolve().parents[3] / "outs" / "Python" / "ForPROSAIL"


def main():
    refl_df = pd.read_csv(OUT_DIR / "1_reflectance.csv")
    wavelengths = np.array([float(c.split(".")[1]) for c in refl_df.columns])
    reflectance = refl_df.to_numpy(dtype=float)

    vnir = get_indices(wavelengths, reflectance, spectral_domain="VNIR")
    swir = get_indices(wavelengths, reflectance, spectral_domain="SWIR")

    indices_df = pd.DataFrame({**vnir, **swir})
    indices_df.to_csv(OUT_DIR / "2_indices.csv", index=False)

    print(f"Computed {indices_df.shape[1]} spectral indices for {indices_df.shape[0]} spectra.")
    print(f"-> {OUT_DIR / '2_indices.csv'}")


if __name__ == "__main__":
    main()
