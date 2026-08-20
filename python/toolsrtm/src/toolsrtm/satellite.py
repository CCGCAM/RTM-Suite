"""Satellite image retrieval via STAC (SpatioTemporal Asset Catalog).

Python port of ``ToolsRTM::get.satellite_collection``/``get.sentinel2_cube``
(R, built on ``rstac``/``gdalcubes``). Needs the optional ``stac`` extra
(``pip install toolsrtm[stac]``: pystac-client, planetary-computer,
odc-stac, rioxarray, rasterio) and a live network connection to the chosen
STAC API -- neither is imported by ``toolsrtm/__init__.py``, so a plain
``import toolsrtm`` never requires them.

Covers the same 9 collections R's own per-collection ``switch()`` does --
Sentinel-2 L2A, Landsat Collection 2 L2, and 6 MODIS products -- via
:data:`COLLECTION_ASSETS` (default asset/band names looked up automatically
from ``collection``, same as R's ``switch()``; pass ``asset_names``
yourself to override or point this at any other STAC collection). Verified
live against Microsoft Planetary Computer for Sentinel-2 L2A and Landsat
C2 L2 (see ``tests/test_satellite.py``); the MODIS collections use the same
query/asset-selection code path but weren't each individually re-verified
live (no ``eo:cloud_cover`` property on MODIS items, so cloud filtering is
skipped for those collections, matching R).

Uses a plain ``(minx, miny, maxx, maxy)`` WGS84 bounding box rather than R's
``sf`` polygon + centroid-buffer workflow (``get_bounding_box``) -- pass your
own already-buffered box (e.g. from ``shapely``/``geopandas``).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

DEFAULT_ASSETS = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12", "SCL"]

_STAC_APIS = {
    "microsoft": "https://planetarycomputer.microsoft.com/api/stac/v1",
    "amazon": "https://earth-search.aws.element84.com/v1",
}

#: Default asset/band names per STAC collection, matching R's own
#: ``get.satellite_collection``'s per-collection ``switch()``.
COLLECTION_ASSETS = {
    "sentinel-2-l2a": DEFAULT_ASSETS,
    "sentinel-s2-l2a": DEFAULT_ASSETS,
    "sentinel-s2-l2a-cogs": DEFAULT_ASSETS,
    "landsat-c2-l2": ["blue", "green", "red", "nir08", "swir16", "swir22"],
    "modis-17A2HGF-061": ["Gpp_500m", "PsnNet_500m"],
    "modis-11A2-061": ["LST_Day_1km"],
    "modis-09A1-061": ["sur_refl_b01", "sur_refl_b02", "sur_refl_b03", "sur_refl_b04",
                        "sur_refl_b05", "sur_refl_b06", "sur_refl_b07"],
    "modis-09Q1-061": ["sur_refl_b01", "sur_refl_b02", "sur_refl_qc_250m", "sur_refl_state_250m"],
    "modis-15A2H-061": ["Lai_500m", "Fpar_500m", "LaiStdDev_500m", "FparStdDev_500m", "FparLai_QC"],
    "modis-15A3H-061": ["Lai_500m", "Fpar_500m", "LaiStdDev_500m", "FparStdDev_500m", "FparLai_QC"],
}

#: Collections that carry an ``eo:cloud_cover`` STAC property (Sentinel-2, Landsat).
#: MODIS collections don't -- cloud filtering is skipped for those, matching R.
_HAS_CLOUD_COVER = {"sentinel-2-l2a", "sentinel-s2-l2a", "sentinel-s2-l2a-cogs", "landsat-c2-l2"}


@dataclass
class SatelliteCollection:
    items: list  #: signed `pystac.Item` objects matching the search
    metadata: "pd.DataFrame"  #: one row per item: id, datetime, platform, cloud_cover (where available)
    cloud_server: str
    collection: str
    asset_names: list[str]


def get_satellite_collection(
    bbox: Sequence[float],
    collection: str = "sentinel-2-l2a",
    date_range: Sequence[str] = ("2023-06-01", "2023-09-30"),
    cloud_server: str = "microsoft",
    n_limit: int = 500,
    cloud_threshold: float = 20.0,
    asset_names: Sequence[str] | None = None,
) -> SatelliteCollection:
    """Search a STAC API for satellite scenes over a bounding box/date range.

    Python port of ``get.satellite_collection`` (R). Filters by cloud cover
    when the collection has an ``eo:cloud_cover`` property (true for
    Sentinel-2/Landsat; skipped otherwise). Items are signed
    (`planetary_computer.sign`) when ``cloud_server="microsoft"`` so their
    asset URLs are directly readable.

    :param bbox: ``(minx, miny, maxx, maxy)`` in WGS84 (EPSG:4326).
    :param collection: STAC collection id (e.g. ``"sentinel-2-l2a"``).
    :param date_range: ``(start, end)`` ISO dates, e.g. ``("2023-06-01", "2023-09-30")``.
    :param cloud_server: ``"microsoft"`` (Planetary Computer) or ``"amazon"`` (Earth Search v1).
    :param n_limit: maximum number of items to return.
    :param cloud_threshold: maximum allowed ``eo:cloud_cover`` percentage.
    :param asset_names: band/asset names to keep; defaults to the entry for
        ``collection`` in :data:`COLLECTION_ASSETS` (falling back to
        :data:`DEFAULT_ASSETS`, Sentinel-2 L2A's bands, for an unlisted collection).
    :return: :class:`SatelliteCollection`.
    """
    import pandas as pd
    import pystac_client

    if cloud_server not in _STAC_APIS:
        raise ValueError(f"cloud_server must be 'microsoft' or 'amazon', got {cloud_server!r}.")

    modifier = None
    if cloud_server == "microsoft":
        import planetary_computer
        modifier = planetary_computer.sign_inplace

    apply_cloud_filter = cloud_threshold is not None and collection in _HAS_CLOUD_COVER

    catalog = pystac_client.Client.open(_STAC_APIS[cloud_server], modifier=modifier)
    search = catalog.search(
        collections=[collection],
        bbox=list(bbox),
        datetime=f"{date_range[0]}/{date_range[1]}",
        max_items=n_limit,
        query={"eo:cloud_cover": {"lt": cloud_threshold}} if apply_cloud_filter else None,
    )
    items = list(search.items())

    rows = []
    for item in items:
        props = item.properties
        rows.append({
            "id": item.id,
            "datetime": props.get("datetime"),
            "platform": props.get("platform"),
            "cloud_cover": props.get("eo:cloud_cover"),
        })
    metadata = pd.DataFrame(rows)

    return SatelliteCollection(
        items=items,
        metadata=metadata,
        cloud_server=cloud_server,
        collection=collection,
        asset_names=list(asset_names) if asset_names is not None
        else list(COLLECTION_ASSETS.get(collection, DEFAULT_ASSETS)),
    )


def get_sentinel2_cube(
    satellite_collection: SatelliteCollection,
    bbox: Sequence[float],
    resolution: float = 10.0,
    crs: str = "EPSG:32632",
    aggregation_method: str | None = None,
    resampling_method: str = "cubic",
) -> "xr.Dataset":
    """Build a Sentinel-2 (or any STAC collection's) data cube from a
    :func:`get_satellite_collection` result, clipped to a bounding box.

    Python port of ``get.sentinel2_cube`` (R, ``gdalcubes``-based) using
    `odc.stac.load` instead -- a different resampling/mosaicking engine than
    R's ``gdalcubes``, so pixel values won't match R exactly even for the
    same scenes/bbox/resolution; both are standard, correct STAC-cube
    builders, not a numerical port.

    :param satellite_collection: result of :func:`get_satellite_collection`.
    :param bbox: ``(minx, miny, maxx, maxy)`` in WGS84 (EPSG:4326) -- the crop extent.
    :param resolution: output pixel size in ``crs`` units (metres, for a UTM ``crs``).
    :param crs: output CRS (a projected CRS -- pick the correct UTM zone for your area).
    :param aggregation_method: if given, one of ``"mean"``, ``"median"``,
        ``"min"``, ``"max"``, ``"first"`` -- collapses the time dimension
        (matching R's own ``aggregation_method``). If ``None``, every scene's
        own timestamp is kept.
    :param resampling_method: reprojection resampling: ``"nearest"``,
        ``"bilinear"``, or ``"cubic"`` (matches ``odc.stac.load``'s ``resampling``).
    :return: an `xarray.Dataset` with one data variable per asset in
        ``satellite_collection.asset_names``.
    """
    import odc.stac

    if not satellite_collection.items:
        raise ValueError("satellite_collection has no items to build a cube from.")

    ds = odc.stac.load(
        satellite_collection.items,
        bands=satellite_collection.asset_names,
        bbox=list(bbox),
        crs=crs,
        resolution=resolution,
        resampling=resampling_method,
        chunks={},
    )

    if aggregation_method is not None:
        if aggregation_method == "mean":
            ds = ds.mean(dim="time", keep_attrs=True)
        elif aggregation_method == "median":
            ds = ds.median(dim="time", keep_attrs=True)
        elif aggregation_method == "min":
            ds = ds.min(dim="time", keep_attrs=True)
        elif aggregation_method == "max":
            ds = ds.max(dim="time", keep_attrs=True)
        elif aggregation_method == "first":
            ds = ds.isel(time=0)
        else:
            raise ValueError(
                f"aggregation_method must be one of 'mean', 'median', 'min', 'max', 'first', "
                f"got {aggregation_method!r}."
            )

    return ds
