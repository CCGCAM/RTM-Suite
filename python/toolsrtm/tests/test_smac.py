from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from toolsrtm.smac import SENSORS, get_sensor, get_smac

REFDATA = Path(__file__).parent / "refdata"

ATMO_FIELDS = ["Ta_ss", "Ta_sd", "Ta_oo", "Ta_do", "Ta_s", "Ta_o", "Tg", "Ra_dd", "Ra_so"]

_EXPECTED_NBANDS = {
    "Sentinel2A.MSI": 13, "Sentinel2B.MSI": 13,
    "Sentinel3A.OLCI": 21, "Sentinel3B.OLCI": 21,
    "LANDSAT4.TM": 6, "LANDSAT5.TM": 6, "LANDSAT7.ETM": 6, "LANDSAT8.OLI": 9,
    "TerraAqua.MODIS": 20,
}


def test_all_9_sensors_load_with_correct_band_counts():
    assert set(SENSORS) == set(_EXPECTED_NBANDS)
    for name, nbands in _EXPECTED_NBANDS.items():
        sensor = get_sensor(name)
        assert len(sensor.wl_smac) == nbands
        for coef in sensor.coef.values():
            assert len(coef) == nbands


def test_get_sensor_rejects_unknown_name():
    with pytest.raises(ValueError):
        get_sensor("not-a-real-sensor")


@pytest.mark.parametrize("r_name,slug", [
    ("LANDSAT8.OLI", "landsat8"),
    ("Sentinel3A.OLCI", "sentinel3a"),
    ("TerraAqua.MODIS", "modis"),
    ("Sentinel2B.MSI", "sentinel2b"),
])
def test_get_smac_matches_r_reference(r_name, slug):
    # Cross-checked against a real, unmodified ToolsRTM::get.smac() call
    # (python/scratch/scratch_verify_smac_sensors.R), same geometry/atmosphere
    # inputs for all 4: tts=30, tto=0, psi=0, Pa=1013.25, aot550=0.2, uo3=0.35, uh2o=2.0.
    sensor = get_sensor(r_name)
    atmo = get_smac(sensor, tts=30, tto=0, psi=0, Pa=1013.25, taup550=0.2, uo3=0.35, uh2o=2.0)
    ref = pd.read_csv(REFDATA / f"ref_smac_{slug}.csv")
    for field in ATMO_FIELDS:
        np.testing.assert_allclose(getattr(atmo, field), ref[field].to_numpy(), rtol=1e-9, atol=1e-12)
