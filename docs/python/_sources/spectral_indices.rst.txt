Spectral Indices
====================

A spectral index is a simple, fixed formula over a handful of bands
(a ratio, a normalized difference, ...) built to track one trait while
staying relatively insensitive to everything else -- the oldest and still
most widely used way to go from reflectance to a trait estimate, and the
baseline :doc:`trait_inversion`'s ML/DL/LUT methods are compared against.
:func:`~toolsrtm.indices.get_indices` computes ~75 VNIR + ~18 SWIR
indices from one reflectance spectrum in a single call.

.. code-block:: python

   import numpy as np
   from toolsrtm import prospect_d, foursail, get_indices

   inputLUT = dict(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
                    LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1, LAI=3, hspot=0.01, tts=30, tto=0, psi=0)
   sail = foursail(inputLUT, np.full(2101, 0.15), leaf_model="PROSPECT-D", spectrum_all=True)
   wl = np.arange(400, 2501)

   idx_vnir = get_indices(wl, sail.rsot, spectral_domain="VNIR")   # 75 indices
   idx_swir = get_indices(wl, sail.rsot, spectral_domain="SWIR")   # 18 indices
   print("NDVI:", idx_vnir["NDVI"], " MSAVI:", idx_vnir["MSAVI"])

VNIR vs. SWIR domain
------------------------

``spectral_domain`` isn't just a filter on which formulas to compute --
it's a real, faithfully-ported behaviour worth knowing:
``spectral_domain="VNIR"`` (an exact match, not a substring check) is
required to get VNIR-only indices like NDVI; ``spectral_domain="VNIR-SWIR"``
alone does **not** include them, matching the R original's own exact
gating logic (see :func:`toolsrtm.indices.get_indices`'s reference page
for the full detail). If a formula needs a wavelength outside the
spectrum you actually have (e.g. a SWIR index needing 1510nm on a
VNIR-only sensor), it comes back ``NaN`` rather than a silently wrong
number -- always worth checking for after a real convolution
(:doc:`sensor_simulation`) rather than assuming every index is available
on every sensor.

Broad categories (not an exhaustive list -- see the reference page for
all ~93 formulas and their exact citations):

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Category
     - Example indices
     - What they track
   * - Vegetation/structure
     - ``NDVI``, ``RDVI``, ``SR``, ``MSR``, ``OSAVI``, ``MSAVI``
     - Green biomass / LAI -- the oldest, most saturating-at-high-LAI
       family.
   * - Pigment/chlorophyll
     - ``MCARI``, ``MCARI1``, ``MCARI2``, ``VOG``, ``VOG2``, ``VOG3``,
       ``CI1``, ``CI2``
     - Chlorophyll content, using the red-edge (~700-750nm) where
       chlorophyll absorption saturates less than in the red.
   * - Red-edge position
     - ``REP``
     - The wavelength of steepest reflectance rise (~700-740nm) --
       shifts with chlorophyll content, used in the
       :ref:`examples:Real Sentinel-2 capstone: data-driven spatial index + Cab mapping (toolsrtm)`
       Loobos example.
   * - SWIR / dry matter, water, nitrogen
     - ``NDNI``, ``S1080``, ``S1260``, ``N1645``, ``N870``
     - Formulas needing wavelengths only a SWIR-capable sensor
       (Sentinel-2's B11/B12, PRISMA, EnMAP) resolves --
       ``spectral_domain="SWIR"`` or ``"VNIR-SWIR"``.

R <-> Python compatibility
-------------------------------

Every index name and formula is a direct, numerically-verified port of
``ToolsRTM::getIndices()`` (see :doc:`verification`) -- the same
``spectral_domain`` argument, the same ~93 index names, the same known
asymmetry above. A LUT of indices computed in R and in Python from the
same reflectance spectrum matches to floating-point noise; nothing here
is a re-derivation or approximation of the R formulas.

What's next
-----------------

- :doc:`trait_inversion` -- ranking indices by correlation with a trait
  and using the winner as a simple, interpretable retrieval, or as one
  input among many to an ML model.
- :doc:`earth_observation` -- mapping an index spatially over a real
  Sentinel-2 scene.
