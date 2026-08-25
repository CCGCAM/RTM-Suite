"""Run get_scope() over a small LUT of varying leaf/canopy traits.

Python equivalent of Scripts/ForSCOPE/3-simulate_LUT.R: perturbs Cab, LAI
and EWT around SCOPEinR's bundled example row and records the resulting
canopy reflectance (at 800nm, NIR plateau), fluorescence flux, and energy
balance for each -- a small sensitivity/LUT demonstration, not a full
25,000-sample training LUT (get_scope() is far more expensive per call
than toolsrtm.foursail: ~0.25s vs ~2ms, dominated by the ebal nonlinear
solve -- see the timing note in Scripts/Python/README.md).
"""
import csv
import time
from pathlib import Path

import numpy as np
import pandas as pd

from scopeinpython import ScopeOptions, get_scope

ROOT = Path(__file__).resolve().parents[3]
LUT_CSV = ROOT / "SCOPEinR" / "inst" / "input" / "LUT_input.csv"
OUT_DIR = ROOT / "outs" / "Python" / "ForSCOPE"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_SAMPLES = 30


def main():
    with open(LUT_CSV, newline="") as f:
        base_row = next(csv.DictReader(f))

    rng = np.random.default_rng(42)
    options = ScopeOptions(calc_fluor=True, calc_xanthophyllabs=False)  # xanthophyll off: keep the loop fast

    rows = []
    t0 = time.time()
    for i in range(N_SAMPLES):
        row = dict(base_row)
        row["Cab"] = float(rng.uniform(15, 65))
        row["LAI"] = float(rng.uniform(0.5, 6.0))
        row["EWT"] = float(rng.uniform(0.005, 0.025))

        result = get_scope(row, options=options)
        rows.append({
            "Cab": row["Cab"], "LAI": row["LAI"], "EWT": row["EWT"],
            "refl_800nm": result.rtmo.refl[800 - 400],
            "refl_1600nm": result.rtmo.refl[1600 - 400],
            "Tcave": result.ebal.Tcave, "Rntot": result.ebal.Rntot, "Actot": result.ebal.Actot,
            "EoutF": result.rtmf.EoutF if result.rtmf else np.nan,
        })

    elapsed = time.time() - t0
    print(f"Ran get_scope() {N_SAMPLES} times in {elapsed:.1f}s ({elapsed / N_SAMPLES * 1000:.0f}ms/call).")

    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "3_scope_lut.csv", index=False)
    print(df.describe().loc[["min", "mean", "max"]].to_string())
    print(f"\nLUT -> {OUT_DIR / '3_scope_lut.csv'}")


if __name__ == "__main__":
    main()
