"""Convolve the LUT from 1_simulate_lut.py onto real sensor bands, using all
three convolution functions toolsrtm provides -- covers the same ground as
the R course pipeline's own sensor-convolution step (e.g.
Scripts/R/ForFoursail2/1-simulate_LUT.R's "Convolving to Sentinel2a/2b/PRISMA"
section), so Python users have the exact same options as R users:

1. spectral_convolution_srf() -- a real, MEASURED per-nm SRF, no
   atmospheric-correction coefficients needed. Covers PRISMA and
   Sentinel-2A/Sentinel-2B.
2. toolsrtm.smac.spectral_convolution() -- a measured SRF bundled together
   WITH SMAC atmospheric-correction coefficients (only Sentinel-2A ships
   this way in the Python port so far).
3. spectral_convolution_gaussian() -- no measured SRF at all, only nominal
   band characteristics (center + FWHM, or center + published band edges),
   approximated as a Gaussian response. Covers EnMAP, ALI, Hyperion, MODIS
   (19-band nominal set), Quickbird, RapidEye, WorldView-2, Landsat, and any
   custom sensor/camera you supply your own band centers (+ optionally FWHM)
   for.

Run 1_simulate_lut.py first.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from toolsrtm.smac import sentinel2a_msi, spectral_convolution
from toolsrtm.srf import (
    srf_prisma,
    srf_sentinel2a,
    srf_sentinel2b,
    spectral_convolution_srf,
    spectral_convolution_gaussian,
)

OUT_DIR = Path(__file__).resolve().parents[3] / "outs" / "Python" / "ForPROSAIL"


def load_lut_reflectance():
    lut = pd.read_csv(OUT_DIR / "1_LUT.csv")
    refl_df = pd.read_csv(OUT_DIR / "1_reflectance.csv")
    wavelengths = np.array([int(c.split(".")[1]) for c in refl_df.columns])
    return lut, wavelengths, refl_df.to_numpy()


def convolve_lut_srf(wavelengths, reflectance, srf_table, fwhm=None):
    """Convolve every row of `reflectance` using spectral_convolution_srf(),
    returning a (n_samples, n_bands) array."""
    first = spectral_convolution_srf(wavelengths, reflectance[0], srf_table, fwhm=fwhm)
    out = np.empty((reflectance.shape[0], len(first.wl)))
    out[0] = first.rfl
    for i in range(1, reflectance.shape[0]):
        out[i] = spectral_convolution_srf(wavelengths, reflectance[i], srf_table, fwhm=fwhm).rfl
    return out, first.band_names


def main():
    lut, wavelengths, reflectance = load_lut_reflectance()
    n = reflectance.shape[0]
    print(f"Loaded {n} spectra ({wavelengths[0]}-{wavelengths[-1]}nm) from 1_simulate_lut.py.")

    # 1. Measured SRF, no atmospheric correction: Sentinel-2A, Sentinel-2B, PRISMA.
    s2a_refl, s2a_bands = convolve_lut_srf(wavelengths, reflectance, srf_sentinel2a())
    s2b_refl, s2b_bands = convolve_lut_srf(wavelengths, reflectance, srf_sentinel2b())
    prisma_refl, prisma_bands = convolve_lut_srf(wavelengths, reflectance, srf_prisma(), fwhm=None)
    print(f"Sentinel-2A: {s2a_refl.shape[1]} bands  |  Sentinel-2B: {s2b_refl.shape[1]} bands"
          f"  |  PRISMA: {prisma_refl.shape[1]} bands")

    # 2. Measured SRF + SMAC atmospheric-correction coefficients: Sentinel-2A only (Python port so far).
    sensor = sentinel2a_msi()
    smac_refl = np.array([spectral_convolution(wavelengths, reflectance[i], sensor) for i in range(n)])
    print(f"Sentinel-2A via SMAC-bundled SRF: {smac_refl.shape[1]} bands "
          f"(TOC reflectance values match the plain-SRF route above; SMAC path additionally "
          f"carries the atmospheric-correction coefficients needed for a TOC->TOA step, see "
          f"toolsrtm.smac.get_smac()).")

    # 3. No measured SRF -- nominal characteristics only, Gaussian-approximated.
    enmap_refl = np.array([spectral_convolution_gaussian(wavelengths, reflectance[i], sensor="EnMAP").rfl
                            for i in range(n)])
    modis_refl = np.array([spectral_convolution_gaussian(wavelengths, reflectance[i], sensor="MODIS").rfl
                            for i in range(n)])
    print(f"EnMAP (Gaussian from nominal center+FWHM): {enmap_refl.shape[1]} bands")
    print(f"MODIS (Gaussian from nominal center+edges): {modis_refl.shape[1]} bands")

    # Your own sensor/camera -- same function as #3, just pass your own centers/fwhm
    # instead of a bundled `sensor` name:
    own_centers = [444, 475, 502, 531, 550, 560, 570, 650, 668, 678, 705, 717, 740, 754, 842]
    own_fwhm = [28, 32, 18, 14, 12, 27, 14, 16, 14, 14, 10, 12, 18, 10, 57]
    own_refl = np.array([
        spectral_convolution_gaussian(wavelengths, reflectance[i], centers=own_centers, fwhm=own_fwhm).rfl
        for i in range(n)
    ])
    print(f"Your own 15-band sensor: {own_refl.shape[1]} bands")

    pd.DataFrame(s2a_refl, columns=s2a_bands).to_csv(OUT_DIR / "5_sentinel2a_bands.csv", index=False)
    pd.DataFrame(s2b_refl, columns=s2b_bands).to_csv(OUT_DIR / "5_sentinel2b_bands.csv", index=False)
    pd.DataFrame(enmap_refl).to_csv(OUT_DIR / "5_enmap_bands.csv", index=False)
    pd.DataFrame(own_refl).to_csv(OUT_DIR / "5_own_sensor_bands.csv", index=False)
    print(f"Saved convolved band tables -> {OUT_DIR}")


if __name__ == "__main__":
    main()
