"""Numerical regression tests for the thermal/aerodynamic building blocks
(Monin-Obukhov, resistances, heatfluxes), vs. reference values generated
by running the original R SCOPEinR package (see
python/scratch/scratch_thermal_export.R). These are the first ported pieces of
the SCOPE thermal energy-balance chain (ebal.R and friends).
"""
from pathlib import Path

import numpy as np
import pandas as pd

from scopeinpython import ResistanceParams, get_heatfluxes, get_resistances, monin_obukhov

REFDATA = Path(__file__).parent / "refdata"


def test_monin_obukhov_matches_r():
    ref = pd.read_csv(REFDATA / "ref_monin_obukhov.csv")
    L = monin_obukhov(ref["ustar"].to_numpy(), ref["Ta"].to_numpy(), ref["H"].to_numpy())
    np.testing.assert_allclose(L, ref["L"].to_numpy(), rtol=1e-6, atol=1e-9)


def test_resistances_matches_r():
    ref = pd.read_csv(REFDATA / "ref_resistances.csv")
    cases = dict(
        unstable=dict(L=-50, z=40, u=3, d=12, hc=15, z0m=1.5, LAI=3, Cd=0.3, rwc=100, leafwidth=0.1, rbs=50),
        stable=dict(L=50, z=40, u=2, d=12, hc=15, z0m=1.5, LAI=3, Cd=0.3, rwc=100, leafwidth=0.1, rbs=50),
        neutral=dict(L=-1e6, z=40, u=3, d=12, hc=15, z0m=1.5, LAI=3, Cd=0.3, rwc=100, leafwidth=0.1, rbs=50),
    )
    for case_name, params in cases.items():
        row = ref[ref["case"] == case_name].iloc[0]
        res = get_resistances(ResistanceParams(
            rbs=params["rbs"], Cd=params["Cd"], LAI=params["LAI"], rwc=params["rwc"],
            z0m=params["z0m"], d=params["d"], hc=params["hc"], leafwidth=params["leafwidth"],
            z=params["z"], u=params["u"], L=params["L"],
        ))
        for field in ["ustar", "uz0", "Kh", "rai", "rar", "rac", "rws", "raa", "rawc", "raws"]:
            assert np.isclose(getattr(res, field), row[field], rtol=1e-6, atol=1e-9), \
                f"{case_name}.{field}: {getattr(res, field)} != {row[field]}"


def test_heatfluxes_matches_r():
    ref = pd.read_csv(REFDATA / "ref_heatfluxes.csv")
    e_to_q = 0.622 / 1013
    hf = get_heatfluxes(
        ra=ref["ra"].to_numpy(), rs=ref["rs"].to_numpy(), Tc=ref["Tc"].to_numpy(),
        ea=ref["ea"].to_numpy(), Ta=ref["Ta"].to_numpy(), e_to_q=e_to_q,
        Ca=ref["Ca"].to_numpy(), Ci=ref["Ci"].to_numpy(),
    )
    np.testing.assert_allclose(hf.lambda_, ref["lambda_"].to_numpy(), rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(hf.s, ref["s"].to_numpy(), rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(hf.lE, ref["lE"].to_numpy(), rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(hf.H, ref["H"].to_numpy(), rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(hf.ec, ref["ec"].to_numpy(), rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(hf.Cc, ref["Cc"].to_numpy(), rtol=1e-6, atol=1e-9)
