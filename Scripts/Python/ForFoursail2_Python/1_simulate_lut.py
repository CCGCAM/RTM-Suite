"""Simulate a PROSPECT-PRO + foursail2 LUT (two-layer green/brown canopy),
convolve to Sentinel-2A/2B/PRISMA, compute indices.

Python equivalent of Scripts/R/ForFoursail2/1-simulate_LUT.R: samples 100
random leaf/canopy parameter sets (including foursail2's own extra
fraction_brown/diss/Cv/Zeta columns), runs toolsrtm.foursail2 for each,
convolves onto 3 sensors, and computes spectral indices per sensor +
native resolution. Feeds 2_inversion_ml.py and 3_inversion_dl.py next.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from toolsrtm import foursail2, get_indices
from toolsrtm.srf import srf_prisma, srf_sentinel2a, srf_sentinel2b, spectral_convolution_srf

OUT_DIR = Path(__file__).resolve().parents[3] / "outs" / "Python" / "ForFoursail2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_SAMPLES = 100
WAVELENGTHS = np.arange(400, 2501)


def sample_lut(n: int, rng: np.random.Generator) -> pd.DataFrame:
    cab = rng.uniform(10, 70, n)
    lut = pd.DataFrame({
        "N": rng.uniform(1.2, 2.5, n),
        "Cab": cab,
        "Car": cab / 4 + rng.normal(0, 2, n),  # co-varies with Cab, same idea as the R course pipeline
        "Anth": rng.uniform(0, 2, n),
        "Cbrown": rng.uniform(0, 0.3, n),
        "EWT": rng.uniform(0.005, 0.03, n),
        "LMA": rng.uniform(0.003, 0.015, n),
        "alpha": np.full(n, 40.0),
        "Prot": np.full(n, 0.0), "CBC": np.full(n, 0.0),  # PROSPECT-PRO needs these even at 0
        "LIDFa": rng.uniform(-0.5, 0.5, n),
        "LIDFb": rng.uniform(-0.3, 0.3, n),
        "TypeLidf": np.full(n, 1, dtype=int),
        "LAI": rng.uniform(0.5, 7.0, n),
        "hspot": np.full(n, 0.01),
        "tts": np.full(n, 30.0),
        "tto": np.full(n, 0.0),
        "psi": np.full(n, 0.0),
        # foursail2's own extra columns -- the two-layer green/brown mix:
        "fraction_brown": rng.uniform(0.0, 0.5, n),  # fraction of canopy that's senescent/brown
        "diss": rng.uniform(0.2, 0.8, n),            # green/brown layer dissociation factor
        "Cv": rng.uniform(0.5, 1.0, n),               # vertical crown cover fraction
        "Zeta": rng.uniform(0.3, 2.0, n),             # tree shape factor
    })
    lut["Car"] = lut["Car"].clip(lower=0.5)
    return lut


def main():
    rng = np.random.default_rng(42)
    lut = sample_lut(N_SAMPLES, rng)
    rsoil = np.full(2101, 0.15)

    reflectance = np.empty((N_SAMPLES, len(WAVELENGTHS)))
    for i, row in lut.iterrows():
        sail2 = foursail2(row.to_dict(), rsoil, leaf_model="PROSPECT-PRO", spectrum_all=True)
        reflectance[i] = sail2.rsot
    print(f"Simulated {N_SAMPLES} PROSPECT-PRO + foursail2 spectra.")

    # ---- Native-resolution indices ----
    indices_native = get_indices(WAVELENGTHS, reflectance, spectral_domain="VNIR")
    indices_native.update(get_indices(WAVELENGTHS, reflectance, spectral_domain="SWIR"))
    dataset_native = pd.concat([lut.reset_index(drop=True), pd.DataFrame(indices_native)], axis=1)
    dataset_native.to_csv(OUT_DIR / "1_dataset_native.csv", index=False)
    print(f"Native: {dataset_native.shape[1]} columns (traits + reflectance-derived indices)")

    # ---- Convolve to Sentinel-2A/2B/PRISMA, then indices per sensor ----
    sensors = {"Sentinel2A": srf_sentinel2a(), "Sentinel2B": srf_sentinel2b(), "PRISMA": srf_prisma()}
    for sensor_name, srf in sensors.items():
        band_refl = np.array([spectral_convolution_srf(WAVELENGTHS, reflectance[i], srf).rfl
                               for i in range(N_SAMPLES)])
        band_wl = spectral_convolution_srf(WAVELENGTHS, reflectance[0], srf).wl
        indices_sensor = get_indices(band_wl, band_refl, spectral_domain="VNIR-SWIR")
        dataset_sensor = pd.concat([lut.reset_index(drop=True), pd.DataFrame(indices_sensor)], axis=1)
        dataset_sensor.to_csv(OUT_DIR / f"1_dataset_{sensor_name.lower()}.csv", index=False)
        print(f"{sensor_name}: {band_refl.shape[1]} bands, {dataset_sensor.shape[1]} columns")

    print(f"\nSaved -> {OUT_DIR}")


if __name__ == "__main__":
    main()
