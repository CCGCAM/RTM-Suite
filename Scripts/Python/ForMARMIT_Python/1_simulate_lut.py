"""Simulate soil reflectance spectra with toolsrtm's MARMIT port across
random surface-water combinations, convolve to Sentinel-2A/PRISMA, compute
indices.

Python equivalent of Scripts/R/ForMARMIT/1-simulate_LUT.R. MARMIT is
soil-only (no vegetation) -- there's no Cab/LAI/EWT here. The target trait
is SMC (gravimetric soil moisture content, %), MARMIT's own physical
output, driven by two knobs: L (surface water film thickness, cm) and eps
(fraction of the surface that's wet). Both are randomized here; soil `id`
(which of Bablet_2016's 17 reference soils) is held fixed so the LUT
varies wetness on one consistent soil, not soil type + wetness at once.
Feeds 2_inversion_ml.py and 3_inversion_dl.py next.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from toolsrtm.marmit import get_marmit_rsoil
from toolsrtm import get_indices
from toolsrtm.srf import srf_prisma, srf_sentinel2a, spectral_convolution_srf

OUT_DIR = Path(__file__).resolve().parents[3] / "outs" / "Python" / "ForMARMIT"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_SAMPLES = 100
SOIL_ID = 1
VERSION = "marmit1"


def main():
    rng = np.random.default_rng(42)
    L = rng.uniform(0.001, 0.15, N_SAMPLES)   # cm, water film thickness
    eps = rng.uniform(0.0, 1.0, N_SAMPLES)     # fraction wet

    wavelengths = None
    reflectance = np.empty((N_SAMPLES, 2101))
    smc = np.empty(N_SAMPLES)
    for i in range(N_SAMPLES):
        res = get_marmit_rsoil(soil_id=SOIL_ID, L=float(L[i]), eps=float(eps[i]), version=VERSION)
        if wavelengths is None:
            wavelengths = res.wavelength
        reflectance[i] = res.rsoil_wet
        smc[i] = res.smc

    lut = pd.DataFrame({"L": L, "eps": eps, "SMC": smc})
    print(f"Simulated {N_SAMPLES} MARMIT soil spectra (soil id={SOIL_ID}, {VERSION}). "
          f"SMC range: {smc.min():.1f}-{smc.max():.1f}%")

    # Most vegetation indices aren't physically meaningful for bare soil,
    # but computed anyway for a consistent predictor set, matching the R
    # course pipeline's own approach.
    indices_native = get_indices(wavelengths, reflectance, spectral_domain="VNIR")
    indices_native.update(get_indices(wavelengths, reflectance, spectral_domain="SWIR"))
    indices_native = {k: v for k, v in indices_native.items() if np.isfinite(v).all()}
    dataset_native = pd.concat([lut.reset_index(drop=True), pd.DataFrame(indices_native)], axis=1)
    dataset_native.to_csv(OUT_DIR / "1_dataset_native.csv", index=False)
    print(f"Native: {dataset_native.shape[1]} columns (L, eps, SMC + reflectance-derived indices)")

    sensors = {"Sentinel2A": srf_sentinel2a(), "PRISMA": srf_prisma()}
    for sensor_name, srf in sensors.items():
        band_refl = np.array([spectral_convolution_srf(wavelengths, reflectance[i], srf).rfl
                               for i in range(N_SAMPLES)])
        band_wl = spectral_convolution_srf(wavelengths, reflectance[0], srf).wl
        indices_sensor = get_indices(band_wl, band_refl, spectral_domain="VNIR-SWIR")
        dataset_sensor = pd.concat([lut.reset_index(drop=True), pd.DataFrame(indices_sensor)], axis=1)
        dataset_sensor.to_csv(OUT_DIR / f"1_dataset_{sensor_name.lower()}.csv", index=False)
        print(f"{sensor_name}: {band_refl.shape[1]} bands, {dataset_sensor.shape[1]} columns")

    print(f"\nSaved -> {OUT_DIR}")


if __name__ == "__main__":
    main()
