18. Applying an Inversion Model Spatially
==============================================

What you will learn
------------------------

- How to run a model trained entirely on simulated data
  (:doc:`t13-machine-learning-inversion`) over every pixel of a real,
  prepared satellite scene (:doc:`t17-preparing-eo-observations`).
- Why this needs a real, wide training domain to work at all on real
  data.

Concept
-----------

Nothing new algorithmically -- ``.predict()`` on a prepared pixel array,
same as any other array of observations. What's new is what the model
was trained on: a LUT (:doc:`t11-lut-generation`) whose trait ranges
must be wide enough, and whose sun/soil/geometry assumptions must be
realistic enough, that real pixel spectra actually fall inside what the
model learned. A LUT built with a narrow, unrealistic domain will
"apply" without error and still produce meaningless output.

.. code-block:: text

   Prepared pixel array (n_pixels, n_bands)
                 |
                 v
        fit.model.predict(pixels)   (fit = trained in t13-machine-learning-inversion)
                 |
                 v
        Per-pixel trait predictions (n_pixels,)
                 |
                 v
        Reshape back to the scene's 2D pixel grid  ->  a trait map

Python tools used
----------------------

Just ``fit.model.predict`` (``fit`` from :func:`~toolsrtm.inversion.get_inversion`,
:doc:`t13-machine-learning-inversion`) -- the only new consideration is
handling invalid pixels correctly.

Run the example
--------------------

.. code-block:: python

   import numpy as np

   # pix: (n_pixels, n_bands) from t17-preparing-eo-observations
   ok = np.all(np.isfinite(pix), axis=1)   # cloud/edge pixels can be NaN
   trait_pixels = np.full(pix.shape[0], np.nan)
   trait_pixels[ok] = fit.model.predict(pix[ok])   # fit: a real trait model, trained on a wide LUT
   trait_map = trait_pixels.reshape(r["B04"].shape)
   print("Valid pixels predicted:", ok.sum(), "/", pix.shape[0])
   print("Trait map range:", np.nanmin(trait_map), "-", np.nanmax(trait_map))

Result
----------

.. figure:: _figures/t18_python_capstone.png
   :alt: Retrieved Cab map alongside an independent spectral index (REP), same real Loobos Sentinel-2 scene
   :width: 100%

   Real output: a Random Forest Cab model -- trained on a 500-row LUT with a deliberately wide domain (sparse-to-dense LAI, variable soil brightness, realistic geometry) -- applied to every pixel of a real Sentinel-2 scene over Loobos forest, NL. The retrieved Cab map (right) and an independently-computed spectral index (REP, left, no ML involved at all) both show the exact same forest gap, entirely independently.

Interpretation
-------------------

The two panels agreeing spatially -- the same dark gap visible in both
the ML-retrieved Cab map and the independently-computed REP index -- is
a real, meaningful consistency check: REP was never used to train the
Cab model, so this agreement isn't circular. It's evidence the retrieval
is picking up a genuine vegetation signal (a real canopy gap) rather
than an artifact of the model or the scene. This kind of independent
cross-check matters more once you're on real data than it did in
Chapters 12-15, where held-out R2 alone was a sufficient (and
achievable) standard -- real scenes don't come with ground-truth Cab to
check against directly.

Try it yourself
--------------------

- Compare the Cab map's value range against the LUT's own training
  ``Cab`` range (:doc:`t11-lut-generation`) -- values pushed to the very
  edge of (or beyond) the training range are a warning sign of
  extrapolation, not interpolation.
- Re-run with a model trained on a *narrow* LUT (e.g. ``LAI`` restricted
  to 2-4) and see the retrieved map degrade or flatten in areas the
  narrow LUT didn't anticipate.
- Mask out the cloud-flagged pixels from :doc:`t17-preparing-eo-observations`
  before predicting, and compare the resulting map's cloud-region values
  against the unmasked version.

Common mistakes
--------------------

- Predicting on cloud-contaminated pixels (:doc:`t17-preparing-eo-observations`)
  produces a confident-looking but meaningless trait value -- always mask
  first, or explicitly interpret cloudy regions as unreliable.
- A model trained on a LUT with too narrow a trait/geometry range will
  silently extrapolate on real pixels outside it, rather than erroring --
  the map will look plausible without being trustworthy.
- Reshaping predictions back to the scene's 2D grid requires the exact
  same pixel order the array was raveled in -- don't reshape into a
  transposed or differently-sized grid by mistake.

Next
--------

:doc:`t19-trait-maps-uncertainty` -- reading this map properly, including
a sense of how much to trust it.
