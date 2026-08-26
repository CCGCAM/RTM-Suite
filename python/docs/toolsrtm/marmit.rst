toolsrtm.marmit
================

MARMIT-1 and MARMIT-2 soil reflectance models (dry -> wet soil spectrum).
Direct port of ``ToolsRTM/R/marmit1.R``, ``marmit2.R`` / ``get.marmit.rsoil.R``.
MARMIT-2 additionally accounts for soil particle size/refractive index and
is generally more accurate for coarser soils; select it via
``get_marmit_rsoil(..., version='marmit2')``.

Quick example
-------------

.. code-block:: python

   from toolsrtm import get_marmit_rsoil

   soil = get_marmit_rsoil(soil_id=3, L=0.05, eps=0.4, version="marmit1")
   print(soil.smc)                     # derived soil moisture content
   print(soil.rsoil_dry[850 - 400], soil.rsoil_wet[850 - 400])   # dry vs. wet reflectance @850nm

.. code-block:: text

   Input                              get_marmit_rsoil()        Output
   ---------------------------        ----------------------    ---------------------------
   soil_id (dry reference spectrum)                              soil.smc   soil moisture content
   L (water-film optical thickness)   -------------------->      soil.rsoil_dry  [2101]
   eps (soil roughness parameter)                                soil.rsoil_wet  [2101]
   version = "marmit1"/"marmit2"                                  (SWIR water-absorption dips
                                                                    deepen as the soil wets)

.. automodule:: toolsrtm.marmit
   :members:
   :undoc-members:
   :show-inheritance:
