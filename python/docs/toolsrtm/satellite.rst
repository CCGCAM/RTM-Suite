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

.. automodule:: toolsrtm.satellite
   :members:
   :undoc-members:
   :show-inheritance:
