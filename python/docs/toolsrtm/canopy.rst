toolsrtm.canopy
================

Leaf-angle distributions (``dladgen``, ``campbell``), volume-scattering
geometry (``volscatt``), the SAIL scattering solutions (conservative and
non-conservative), and the fourSAIL / fourSAIL2 canopy BRDF models.
Direct port of ``ToolsRTM/R/dladgen.R``, ``campbell.R``, ``volscatt.R``,
``Jfunc1-4.R``, ``NonConservativeScatering.R``, ``ConservativeScattering.R``,
``foursail.R``, ``foursail2.R``.

Quick example
-------------

.. code-block:: python

   import numpy as np
   from toolsrtm import foursail

   inputLUT = dict(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
                    LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
                    LAI=3, hspot=0.01, tts=30, tto=0, psi=0)
   rsoil = np.full(2101, 0.15)   # flat soil background, 400-2500nm
   sail = foursail(inputLUT, rsoil, leaf_model="PROSPECT-D", spectrum_all=True)
   print(sail.rsot[550 - 400])   # TOC bidirectional reflectance at 550nm

.. code-block:: text

   Input                              foursail()               Output
   ---------------------------        ----------------------   ---------------------------
   leaf traits (N, Cab, Car, ...)                               sail.rsot  TOC directional refl.
   LAI, LIDFa/LIDFb (leaf angle)      -------------------->      sail.rdot  TOC hemispherical refl.
   hspot (hot-spot size)                                        sail.rsdt, sail.rddt  (diffuse terms)
   tts/tto/psi (sun-view geometry)
   rsoil (soil background spectrum)

.. automodule:: toolsrtm.canopy
   :members:
   :undoc-members:
   :show-inheritance:
