toolsrtm.spart
===============

SPART top-of-canopy (``spart_toc``) and top-of-atmosphere (``spart_toa``)
reflectance/radiance. Direct port of ``ToolsRTM/R/spart.R``, verified
against a real, unmodified ``ToolsRTM::SPART()`` call.

.. note::
   ``spart_toa``'s atmospheric correction (SMAC) is implemented in
   :doc:`smac` and is fully sensor-agnostic, but only **Sentinel-2A**'s
   sensor coefficients are bundled as package data so far, out of the 9
   the R package ships — see :doc:`smac`'s module docstring for exactly
   what adding another sensor takes (a data-export exercise, not a code
   change), and :doc:`../not_ported`.

Quick example
-------------

.. code-block:: python

   from toolsrtm import spart_toa, sentinel2a_msi

   inputLUT = dict(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
                    LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
                    LAI=3, hspot=0.01, tts=30, tto=0, psi=0,
                    Pa=1000, aot550=0.3246, uo3=0.348, uh2o=1.4116)   # atmosphere
   result = spart_toa(inputLUT, sensor=sentinel2a_msi(), leaf_model="PROSPECT-PRO",
                       BSMBrightness=0.5, BSMlat=25, BSMlon=45, SMp=15)
   print(result.rfl_toa.round(4))   # TOA reflectance, per Sentinel-2A band, after SMAC correction

.. code-block:: text

   Input                              spart_toa()               Output
   ---------------------------        ----------------------    ---------------------------
   leaf + canopy inputs (as foursail)                            result.rfl_toc_brdf  TOC refl.
   BSM soil params (Brightness,        -------------------->      result.rfl_toa       TOA refl.
     lat, lon, SMp)                                                (already resampled to the
   sensor = sentinel2a_msi()                                       sensor's own bands)
   atmosphere (Pa, aot550, uo3, uh2o)

.. automodule:: toolsrtm.spart
   :members:
   :undoc-members:
   :show-inheritance:
