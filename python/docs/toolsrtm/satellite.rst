toolsrtm.satellite
====================

Satellite image retrieval via STAC (SpatioTemporal Asset Catalog). Port of
``ToolsRTM::get.satellite_collection``/``get.sentinel2_cube``, scoped to
Sentinel-2 L2A on Microsoft Planetary Computer / AWS Earth Search -- see the
module docstring for exactly what's in and out of scope versus the R source.

.. note::
   Needs the optional ``stac`` extra: ``pip install "toolsrtm[stac]"``
   (pystac-client, planetary-computer, odc-stac, rioxarray, rasterio) and a
   live network connection to the chosen STAC API.

Quick example
-------------

.. code-block:: python

   from toolsrtm.satellite import get_satellite_collection, get_sentinel2_cube

   lat, lon, d = 52.166447, 5.74355, 0.006   # Loobos forest, NL
   bbox = (lon - d, lat - d, lon + d, lat + d)
   coll = get_satellite_collection(bbox, collection="sentinel-2-l2a",
                                    date_range=("2024-07-01", "2024-07-31"),
                                    cloud_server="microsoft", n_limit=20, cloud_threshold=40)
   cube = get_sentinel2_cube(coll, bbox, resolution=10.0, crs="EPSG:32631", aggregation_method="mean")
   print(cube["B08"].shape)   # real Sentinel-2 NIR band, real pixel grid

.. code-block:: text

   Input                              get_sentinel2_cube()      Output
   ---------------------------        ----------------------    ---------------------------
   bbox (lon/lat box)                                            cube : xarray Dataset
   collection = "sentinel-2-l2a"      -------------------->       -- one DataArray per band
   date_range, cloud_threshold                                    -- real pixel reflectance,
   resolution, crs, aggregation                                      ready for indices/ML models

.. automodule:: toolsrtm.satellite
   :members:
   :undoc-members:
   :show-inheritance:
