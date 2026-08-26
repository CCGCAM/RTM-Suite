toolsrtm.smac
==============

SMAC (Simplified Method for Atmospheric Correction, Rahman & Dedieu 1994)
atmospheric radiative transfer, and the per-sensor coefficient/spectral-
response data it needs. Direct port of ``ToolsRTM::get.smac``/
``get.coef.SMAC``/``get.spectral.convolution``.

.. note::
   The atmospheric-correction physics (:func:`get_smac`,
   :func:`spectral_convolution`) is fully general and sensor-agnostic
   given a :class:`SmacSensor` -- only **Sentinel-2A (MSI)**'s data is
   bundled as package data so far, out of the 9 sensors the R package
   ships. See the module docstring for the exact R-side export recipe to
   add another sensor (a data-export exercise, not a code change).

Quick example
-------------

.. code-block:: python

   from toolsrtm import sentinel2a_msi

   s2a = sentinel2a_msi()          # bundled sensor coefficients (SmacSensor object)
   print(s2a.band_names)           # 13 Sentinel-2A MSI band names

.. code-block:: text

   Input                              sentinel2a_msi()          Output
   ---------------------------        ----------------------    ---------------------------
   (no arguments -- reads bundled                                s2a : SmacSensor
    package data)                     -------------------->       -- SMAC atmospheric coeffs.
                                                                   -- band spectral response fns.
                                                                   -- feeds get_smac()/spart_toa()

.. automodule:: toolsrtm.smac
   :members:
   :undoc-members:
   :show-inheritance:
