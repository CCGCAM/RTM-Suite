Earth Observation
=====================

The last step: applying a model trained on simulated data
(:doc:`trait_inversion`) to a **real** satellite scene, and mapping the
result spatially. :func:`toolsrtm.satellite`'s STAC (SpatioTemporal Asset
Catalog) functions retrieve real Sentinel-2 L2A imagery -- no separate
download step, no manual band-file management.

.. note::
   Needs the optional ``stac`` extra: ``pip install "toolsrtm[stac]"``
   (pystac-client, planetary-computer, odc-stac, rioxarray, rasterio) and
   a live network connection.

1. Retrieving a real scene
-------------------------------

.. code-block:: python

   from toolsrtm.satellite import get_satellite_collection, get_sentinel2_cube

   lat, lon, d = 52.166447, 5.74355, 0.006   # Loobos forest, NL (ICOS Scots pine site)
   bbox = (lon - d, lat - d, lon + d, lat + d)
   coll = get_satellite_collection(bbox, collection="sentinel-2-l2a",
                                    date_range=("2024-07-01", "2024-07-31"),
                                    cloud_server="microsoft", n_limit=20, cloud_threshold=40)
   cube = get_sentinel2_cube(coll, bbox, resolution=10.0, crs="EPSG:32631", aggregation_method="mean")
   print(cube["B08"].shape)   # real Sentinel-2 NIR band, real pixel grid

``get_satellite_collection`` searches and cloud-filters candidate scenes
over a date range; ``get_sentinel2_cube`` retrieves and aggregates them
(``"mean"`` composites the whole window into one cloud-free-ish image --
useful over a month where any single date might be cloudy).

2. Image preparation: reflectance scaling and band naming
------------------------------------------------------------------

Sentinel-2 L2A reflectance from STAC is scaled by 10000 (integer storage);
rescale to `[0, 1]` before feeding it to anything trained on simulated
reflectance (:doc:`sensor_simulation`'s ``spectral_convolution_srf``
already returns `[0, 1]`, so both sides need to match):

.. code-block:: python

   real_names = ["B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]
   r = {b: cube[b].values.astype(float) / 10000 for b in real_names}

3. Applying a trained inversion model spatially
-----------------------------------------------------

Every pixel, through the same model trained in :doc:`trait_inversion`:

.. code-block:: python

   import numpy as np
   pix = np.column_stack([r[b].ravel() for b in real_names])
   ok = np.all(np.isfinite(pix), axis=1)   # cloud/edge pixels can be NaN
   cab_pixels = np.full(pix.shape[0], np.nan)
   cab_pixels[ok] = fit.model.predict(pix[ok])   # `fit` from get_inversion(), Section 2 of Trait Inversion
   cab_map = cab_pixels.reshape(r["B04"].shape)

4. Trait maps and uncertainty
-----------------------------------

A single ``.predict()`` call gives a point estimate per pixel, not
uncertainty -- two practical ways to get a sense of it without a fully
Bayesian model: compare the map against an independent spectral index
computed on the same real bands (:doc:`spectral_indices`; the two should
agree spatially if the retrieval is trustworthy), or, for
``algorithm="RF"``, use the per-tree prediction spread
(``fit.model.estimators_``) as a rough per-pixel confidence proxy.

Two full, real, end-to-end capstones (real Sentinel-2 data, both
directions checked against an independent signal) live on
:doc:`examples`:

- :ref:`examples:Real Sentinel-2 capstone: retrieving net photosynthesis (scopeinpython)`
  -- Speulderbos forest, SCOPE-trained Actot retrieval, checked against
  the NDVI seasonal cycle at the same site.
- :ref:`examples:Real Sentinel-2 capstone: data-driven spatial index + Cab mapping (toolsrtm)`
  -- Loobos forest, RF-retrieved Cab mapped alongside the independently
  winning spectral index (REP), both showing the same forest gap.

What's next
-----------------

- :doc:`trait_inversion` -- training the model applied spatially above.
- :doc:`sensor_simulation` -- matching simulated-spectrum band
  definitions to what STAC actually returns.
