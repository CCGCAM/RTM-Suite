0-RTM-Suite — Python port
==========================

Python port of the core of the two R packages in `0-RTM-Suite
<https://gitlab.com/caminoccg>`_: **toolsrtm** (leaf + canopy radiative
transfer models) and **scopeinpython** (soil model + the optical top-of-canopy
BRDF pipeline of SCOPE). Pure NumPy/SciPy at runtime — no R dependency.

This is a **partial, deliberately scoped** port. Every function documented
here has been numerically verified against the original R package (see
:doc:`verification`) before being included — nothing here is a guess at
what the R code does.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   install
   toolsrtm/index
   scopeinpython/index
   verification
   examples
   not_ported

Quick example
-------------

.. code-block:: python

   import numpy as np
   from toolsrtm import prospect_d, foursail

   leaf = prospect_d(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40)

   inputLUT = dict(
       N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
       Prot=0.002, CBC=0.007, LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
       LAI=3, hspot=0.01, tts=30, tto=0, psi=0,
   )
   rsoil = np.full(2101, 0.15)
   sail = foursail(inputLUT, rsoil, leaf_model="PROSPECT-D", spectrum_all=True)
   print("TOC BRDF reflectance at 550 nm:", sail.rsot[550 - 400])

Indices and tables
-------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
