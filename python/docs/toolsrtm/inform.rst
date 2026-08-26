toolsrtm.inform
================

INFORM forest-stand reflectance (Atzberger 2000; Schlerf & Atzberger 2006):
FLIM ground-coverage geometry + fourSAIL understorey/infinite-crown +
crown transmittance. Direct port of ``ToolsRTM/R/inform.R`` and its
internal helpers (``foursail.inform.R``, ``foursail.inf.R``,
``foursail_t_s.R``, ``foursail_t_o.R``, ``Compute_BRF.R``).

Quick example
-------------

.. code-block:: python

   import numpy as np
   from toolsrtm import inform

   inputLUT = dict(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
                    LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
                    LAI=3, hspot=0.01, tts=30, tto=0, psi=0,
                    LAIu=0.5, sd=650, cd=4.5, h=20, skyl=0.1)   # INFORM-only crown/stand geometry
   r_forest = inform(inputLUT, np.full(2101, 0.15), leaf_model="PROSPECT-D")
   print(r_forest[800 - 400])   # TOC reflectance at 800nm, forest-stand geometry

.. code-block:: text

   Input                              inform()                  Output
   ---------------------------        ----------------------    ---------------------------
   same leaf + fourSAIL inputs as                                r_forest [2101]  TOC BRF,
   toolsrtm.canopy, plus:              -------------------->      400-2500nm -- explicit crown/
   LAIu (understorey LAI)                                         gap geometry lowers it vs a
   sd, cd, h (stem density,                                       homogeneous fourSAIL canopy
             crown diameter, height)
   skyl (diffuse-light fraction)

.. automodule:: toolsrtm.inform
   :members:
   :undoc-members:
   :show-inheritance:
