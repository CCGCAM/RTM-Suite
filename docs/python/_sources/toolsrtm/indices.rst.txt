toolsrtm.indices
=================

~75 VNIR + ~18 SWIR spectral vegetation indices from a reflectance
spectrum. Direct port of ``ToolsRTM/R/getIndices.R``.

Quick example
-------------

.. code-block:: python

   import numpy as np
   from toolsrtm import prospect_d, foursail, get_indices

   inputLUT = dict(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
                    LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
                    LAI=3, hspot=0.01, tts=30, tto=0, psi=0)
   sail = foursail(inputLUT, np.full(2101, 0.15), leaf_model="PROSPECT-D", spectrum_all=True)
   wl = np.arange(400, 2501)

   idx = get_indices(wl, sail.rsot, spectral_domain="VNIR")
   print(idx["NDVI"], idx["MSAVI"])

.. code-block:: text

   Input                              get_indices()             Output
   ---------------------------        ----------------------    ---------------------------
   wl   [n]  wavelengths (nm)                                   dict of ~75 VNIR (or ~18 SWIR)
   refl [n]  reflectance spectrum     -------------------->      index name -> value, e.g.
   spectral_domain = "VNIR"/"SWIR"/                               idx["NDVI"], idx["REP"], ...
                      "VNIR-SWIR"

.. note::
   ``spectral_domain="VNIR"`` is required (strict match, not a substring
   check) to get VNIR-only indices like NDVI -- ``"VNIR-SWIR"`` alone
   does **not** include them, matching the R original's exact gating
   (``ToolsRTM/R/getIndices.R``).

.. automodule:: toolsrtm.indices
   :members:
   :undoc-members:
   :show-inheritance:
