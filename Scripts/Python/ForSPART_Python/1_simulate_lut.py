"""Simulate a PROSPECT-PRO + fourSAIL SPART LUT: soil + canopy + atmosphere
radiative transfer, giving TOC *and* TOA reflectance at a real sensor's
own bands directly (no separate "native then convolve" step).

Python equivalent of Scripts/R/ForSPART/1-simulate_LUT.R. Unlike
foursail()/foursail2()/inform(), spart_toc()/spart_toa() output straight
at sensor resolution -- "convolution" here means calling spart_toa() once
per sensor instead of a separate resampling step. Atmospheric conditions
(Pa, aot550, uo3, uh2o) use standard clear-sky values (not derived from
altitude like ToolsRTM::get.smac()'s own R fallback) -- a documented
simplification, see the note below.

Feeds 2_inversion_ml.py and 3_inversion_dl.py next.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from toolsrtm import spart_toc, spart_toa, get_indices
from toolsrtm.smac import sentinel2a_msi

OUT_DIR = Path(__file__).resolve().parents[3] / "outs" / "Python" / "ForSPART"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_SAMPLES = 100
WAVELENGTHS = np.arange(400, 2401)  # spart_toc's own native domain (Fluspect-consistent)

# Standard clear-sky atmosphere -- spart_toa() needs these in inputLUT with
# no built-in default (unlike R's get.smac(), which derives Pa from
# altitude when missing); using fixed typical values is a documented
# simplification of this course pipeline, not a physical measurement.
ATMOSPHERE = dict(Pa=1013.25, aot550=0.2, uo3=0.35, uh2o=2.0)


def sample_lut(n: int, rng: np.random.Generator) -> pd.DataFrame:
    cab = rng.uniform(10, 70, n)
    lut = pd.DataFrame({
        "N": rng.uniform(1.2, 2.5, n),
        "Cab": cab,
        "Car": (cab / 4 + rng.normal(0, 2, n)).clip(min=0.5),
        "Anth": rng.uniform(0, 2, n),
        "Cbrown": rng.uniform(0, 0.3, n),
        "EWT": rng.uniform(0.005, 0.03, n),
        "LMA": rng.uniform(0.003, 0.015, n),
        "alpha": np.full(n, 40.0),
        "Prot": np.full(n, 0.0), "CBC": np.full(n, 0.0),
        "LIDFa": rng.uniform(-0.5, 0.5, n),
        "LIDFb": rng.uniform(-0.3, 0.3, n),
        "TypeLidf": np.full(n, 1, dtype=int),
        "LAI": rng.uniform(0.5, 7.0, n),
        "hspot": np.full(n, 0.01),
        "tts": np.full(n, 30.0),
        "tto": np.full(n, 0.0),
        "psi": np.full(n, 0.0),
        "Pa": np.full(n, ATMOSPHERE["Pa"]), "aot550": np.full(n, ATMOSPHERE["aot550"]),
        "uo3": np.full(n, ATMOSPHERE["uo3"]), "uh2o": np.full(n, ATMOSPHERE["uh2o"]),
    })
    return lut


def main():
    rng = np.random.default_rng(42)
    lut = sample_lut(N_SAMPLES, rng)

    # ---- TOC reflectance (BSM soil, no atmosphere) ----
    toc = np.empty((N_SAMPLES, len(WAVELENGTHS)))
    for i, row in lut.iterrows():
        toc[i] = spart_toc(row.to_dict(), leaf_model="PROSPECT-PRO", SMC=25, BSMBrightness=0.5)
    print(f"Simulated {N_SAMPLES} SPART TOC spectra.")
    indices_toc = get_indices(WAVELENGTHS, toc, spectral_domain="VNIR")
    indices_toc.update(get_indices(WAVELENGTHS, toc, spectral_domain="SWIR"))
    dataset_toc = pd.concat([lut.reset_index(drop=True), pd.DataFrame(indices_toc)], axis=1)
    dataset_toc.to_csv(OUT_DIR / "1_dataset_toc.csv", index=False)
    print(f"TOC: {dataset_toc.shape[1]} columns (traits + reflectance-derived indices)")

    # ---- TOA reflectance at Sentinel-2A's own bands directly ----
    sensor = sentinel2a_msi()
    toa_bands = None
    band_wl = None
    for i, row in lut.iterrows():
        res = spart_toa(row.to_dict(), sensor=sensor, leaf_model="PROSPECT-PRO", SMC=25, BSMBrightness=0.5)
        if toa_bands is None:
            band_wl = np.asarray(res.wl_smac)
            toa_bands = np.empty((N_SAMPLES, len(band_wl)))
        toa_bands[i] = res.rfl_toa
    print(f"Simulated {N_SAMPLES} SPART TOA (Sentinel-2A) spectra, {toa_bands.shape[1]} bands.")
    indices_toa = get_indices(band_wl, toa_bands, spectral_domain="VNIR-SWIR")
    dataset_toa = pd.concat([lut.reset_index(drop=True), pd.DataFrame(indices_toa)], axis=1)
    dataset_toa.to_csv(OUT_DIR / "1_dataset_toa_sentinel2a.csv", index=False)
    print(f"TOA: {dataset_toa.shape[1]} columns")

    print(f"\nSaved -> {OUT_DIR}")


if __name__ == "__main__":
    main()
