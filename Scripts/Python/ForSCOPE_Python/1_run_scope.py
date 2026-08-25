"""Run one full SCOPE simulation via the Python port's get_scope().

Python equivalent of Scripts/ForSCOPE/1-getSCOPE.R: loads SCOPEinR's own
bundled example LUT row (one leaf/canopy/meteo/soil parameter set) and
runs the full leaf optics -> soil -> optical BRDF -> energy balance ->
fluorescence -> zeaxanthin pipeline in a single get_scope() call.
"""
import csv
from pathlib import Path

import pandas as pd

from scopeinpython import ScopeOptions, get_scope, get_spectra_scope

ROOT = Path(__file__).resolve().parents[3]
LUT_CSV = ROOT / "SCOPEinR" / "inst" / "input" / "LUT_input.csv"
OUT_DIR = ROOT / "outs" / "Python" / "ForSCOPE"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main():
    with open(LUT_CSV, newline="") as f:
        row = next(csv.DictReader(f))

    result = get_scope(row, options=ScopeOptions(calc_fluor=True, calc_xanthophyllabs=True))

    print(f"Canopy layers: {result.nlayers}")
    print(f"TOC reflectance at 550nm: {result.rtmo.refl[550 - 400]:.4f}")
    print(f"TOC reflectance at 800nm: {result.rtmo.refl[800 - 400]:.4f}")
    print(f"Canopy-average leaf temperature (Tcave): {result.ebal.Tcave:.2f} degC")
    print(f"Soil temperature, sunlit/shaded: {result.ebal.Tsu:.2f} / {result.ebal.Tsh:.2f} degC")
    print(f"Net radiation, total (Rntot): {result.ebal.Rntot:.1f} W/m2")
    print(f"Latent heat flux, total (lEtot): {result.ebal.lEtot:.1f} W/m2")
    print(f"Sensible heat flux, total (Htot): {result.ebal.Htot:.1f} W/m2")
    print(f"Canopy photosynthesis (Actot): {result.ebal.Actot:.2f} umol m-2 s-1")
    print(f"Energy-balance convergence: {result.ebal.counter} iterations")
    if result.rtmf is not None:
        print(f"Canopy fluorescence flux (EoutF): {result.rtmf.EoutF:.4f} W/m2")
        print(f"F685/F740: {result.rtmf.F685:.3f} / {result.rtmf.F740:.3f}")

    # rtmo's spectral fields (refl/rdd/rsd/rdo/rso) are on the full SCOPE
    # wlS grid: 400-2400nm (1nm step), then 2500-15000nm (100nm step), then
    # 16000-50000nm (1000nm step) -- NOT a uniform 1nm grid throughout.
    wlS = get_spectra_scope().wlS
    spectral_df = pd.DataFrame({
        "wave": wlS,
        "refl": result.rtmo.refl,
        "rdd": result.rtmo.rdd,
        "rsd": result.rtmo.rsd,
        "rdo": result.rtmo.rdo,
        "rso": result.rtmo.rso,
    })
    spectral_df.to_csv(OUT_DIR / "1_scope_spectral.csv", index=False)

    scalars = pd.DataFrame([{
        "nlayers": result.nlayers, "Tcave": result.ebal.Tcave, "Tsu": result.ebal.Tsu, "Tsh": result.ebal.Tsh,
        "Rntot": result.ebal.Rntot, "lEtot": result.ebal.lEtot, "Htot": result.ebal.Htot,
        "Actot": result.ebal.Actot, "counter": result.ebal.counter,
        "Pntot_Cab": result.Pntot_Cab, "Ja": result.Ja, "PNPQ": result.PNPQ, "fqe": result.fqe,
        "EoutF": result.rtmf.EoutF if result.rtmf else None,
    }])
    scalars.to_csv(OUT_DIR / "1_scope_scalars.csv", index=False)

    print(f"\nSpectral output -> {OUT_DIR / '1_scope_spectral.csv'}")
    print(f"Scalar output -> {OUT_DIR / '1_scope_scalars.csv'}")


if __name__ == "__main__":
    main()
