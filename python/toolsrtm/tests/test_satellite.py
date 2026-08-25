import pytest

pystac_client = pytest.importorskip("pystac_client")

from toolsrtm.satellite import DEFAULT_ASSETS, get_satellite_collection, get_sentinel2_cube

# Small bbox around Wageningen, NL -- a real, known area with reliable
# Sentinel-2 coverage, used only to confirm the live STAC query/cube-build
# path actually works end-to-end, not to check specific pixel values.
WAGENINGEN_BBOX = (5.60, 51.95, 5.70, 52.00)


def _skip_if_offline():
    import urllib.request
    try:
        urllib.request.urlopen("https://planetarycomputer.microsoft.com/api/stac/v1", timeout=8)
    except Exception as exc:
        pytest.skip(f"Planetary Computer STAC API not reachable: {exc}")


def test_get_satellite_collection_finds_real_sentinel2_scenes():
    _skip_if_offline()
    coll = get_satellite_collection(
        WAGENINGEN_BBOX, collection="sentinel-2-l2a", date_range=("2023-06-01", "2023-08-31"),
        cloud_server="microsoft", n_limit=5, cloud_threshold=30,
    )
    assert len(coll.items) > 0
    assert set(coll.metadata.columns) >= {"id", "datetime", "platform", "cloud_cover"}
    assert (coll.metadata["cloud_cover"] < 30).all()
    assert coll.asset_names == DEFAULT_ASSETS


def test_get_satellite_collection_rejects_unknown_cloud_server():
    with pytest.raises(ValueError):
        get_satellite_collection(WAGENINGEN_BBOX, cloud_server="not-a-real-provider")


def test_get_sentinel2_cube_builds_real_data():
    _skip_if_offline()
    coll = get_satellite_collection(
        WAGENINGEN_BBOX, collection="sentinel-2-l2a", date_range=("2023-06-01", "2023-08-31"),
        cloud_server="microsoft", n_limit=2, cloud_threshold=30,
    )
    ds = get_sentinel2_cube(coll, WAGENINGEN_BBOX, resolution=100.0, crs="EPSG:32631",
                             aggregation_method="median")
    assert "B04" in ds.data_vars
    assert "time" not in ds.dims  # collapsed by aggregation_method="median"
    b04 = ds["B04"].values
    assert b04.size > 0
    assert (b04[~__import__("numpy").isnan(b04)] >= 0).all()


def test_get_sentinel2_cube_rejects_empty_collection():
    from toolsrtm.satellite import SatelliteCollection
    import pandas as pd

    empty = SatelliteCollection(items=[], metadata=pd.DataFrame(), cloud_server="microsoft",
                                 collection="sentinel-2-l2a", asset_names=DEFAULT_ASSETS)
    with pytest.raises(ValueError):
        get_sentinel2_cube(empty, WAGENINGEN_BBOX)
