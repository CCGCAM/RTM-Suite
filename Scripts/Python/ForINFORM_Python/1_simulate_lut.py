"""Simulate a Fluspect-B-Cx + INFORM LUT (forest canopy, explicit crown
geometry/shadowing), convolve to Sentinel-2A/2B/PRISMA, compute indices.

Python equivalent of Scripts/R/ForINFORM/1-simulate_LUT.R. INFORM (unlike
fourSAIL/foursail2) returns TOC reflectance directly as an array -- no
separate compute_brf() step. Native domain is 2001 points (400-2400nm,
Fluspect's own shorter grid) rather than fourSAIL's 2101 (400-2500nm).
Feeds 2_inversion_ml.py and 3_inversion_dl.py next.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from toolsrtm import inform, get_indices
from toolsrtm.srf import srf_prisma, srf_sentinel2a, srf_sentinel2b, spectral_convolution_srf

OUT_DIR = Path(__file__).resolve().parents[3] / "outs" / "Python" / "ForINFORM"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_SAMPLES = 100
WAVELENGTHS = np.arange(400, 2401)  # Fluspect's native domain


def sample_lut(n: int, rng: np.random.Generator) -> pd.DataFrame:
    cab = rng.uniform(10, 70, n)
    lut = pd.DataFrame({
        "N": rng.uniform(1.2, 2.5, n),
        "Cab": cab,
        "Car": (cab / 4 + rng.normal(0, 2, n)).clip(min=0.5),
        "Anth": rng.uniform(0, 2, n),
        "EWT": rng.uniform(0.005, 0.03, n),
        "LMA": rng.uniform(0.003, 0.015, n),
        "alpha": np.full(n, 40.0),
        "Cs": np.full(n, 0.1), "fqe": np.full(n, 0.01), "Cx": np.full(n, 0.5),  # Fluspect-B-Cx needs these
        "Prot": np.full(n, 0.0), "CBC": np.full(n, 0.0),  # unused by Fluspect but read unconditionally
        "LIDFa": rng.uniform(-0.5, 0.5, n),
        "LIDFb": rng.uniform(-0.3, 0.3, n),
        "TypeLidf": np.full(n, 1, dtype=int),
        "LAI": rng.uniform(0.5, 7.0, n),
        "hspot": np.full(n, 0.01),
        "tts": np.full(n, 30.0),
        "tto": np.full(n, 0.0),
        "psi": np.full(n, 0.0),
        # INFORM's own extra columns -- forest stand structure:
        "LAIu": rng.uniform(0.1, 2.0, n),   # understorey LAI
        "cd": rng.uniform(3.0, 15.0, n),    # tree crown diameter, m
        "h": rng.uniform(5.0, 30.0, n),     # tree height, m
        "sd": rng.uniform(100.0, 1500.0, n),  # stem density, ha-1
        "skyl": rng.uniform(0.05, 0.3, n),  # diffuse-radiation fraction
    })
    return lut


def main():
    rng = np.random.default_rng(42)
    lut = sample_lut(N_SAMPLES, rng)
    rsoil = np.full(2101, 0.15)  # inform() truncates this to its own 2001-domain internally

    reflectance = np.empty((N_SAMPLES, len(WAVELENGTHS)))
    for i, row in lut.iterrows():
        reflectance[i] = inform(row.to_dict(), rsoil, leaf_model="Fluspect-B-Cx")
    print(f"Simulated {N_SAMPLES} Fluspect-B-Cx + INFORM spectra.")

    indices_native = get_indices(WAVELENGTHS, reflectance, spectral_domain="VNIR")
    indices_native.update(get_indices(WAVELENGTHS, reflectance, spectral_domain="SWIR"))
    dataset_native = pd.concat([lut.reset_index(drop=True), pd.DataFrame(indices_native)], axis=1)
    dataset_native.to_csv(OUT_DIR / "1_dataset_native.csv", index=False)
    print(f"Native: {dataset_native.shape[1]} columns (traits + reflectance-derived indices)")

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
