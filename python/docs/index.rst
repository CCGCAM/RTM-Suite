0-RTM-Suite — Python port
==========================

Python port of the core of the two R packages in `0-RTM-Suite
<https://gitlab.com/caminoccg>`_: **toolsrtm** (leaf + canopy radiative
transfer models) and **scopeinpython** (soil model + the optical top-of-canopy
BRDF pipeline of SCOPE). Pure NumPy/SciPy at runtime — no R dependency.

This is a **partial, deliberately scoped** port. Every function documented
here has been numerically verified against the original R package (see
:doc:`verification`) before being included — nothing here is a guess at
what the R code does. New to what ``N``, ``Cab``, ``LIDFa``, ``Vcmax25``,
or any other trait actually means? See the :doc:`t02-parameters-traits` first.

New here? `Python tutorials overview <../tutorials-python.html>`_ is a
five-minute visual map of all 20 chapters below (a real figure from
each, grouped by topic) — this page is the detailed version, and the
:doc:`toolsrtm/index`/:doc:`scopeinpython/index` reference pages have
exact function signatures once you know which one you need.

.. toctree::
   :maxdepth: 1
   :caption: LEARN -- I. Fundamentals & RTM Simulation

   install
   t01-getting-started
   t02-parameters-traits
   t03-leaf-models
   t04-canopy-models
   t05-soil-atmosphere
   t06-scope

.. toctree::
   :maxdepth: 1
   :caption: LEARN -- II. From Simulation to Observations

   t07-building-workflows
   t08-sensor-simulation
   t09-spectral-indices
   t10-sensitivity-analysis

.. toctree::
   :maxdepth: 1
   :caption: LEARN -- III. From Spectra to Traits

   t11-lut-generation
   t12-lut-inversion
   t13-machine-learning-inversion
   t14-deep-learning-inversion
   t15-choosing-inversion-strategy

.. toctree::
   :maxdepth: 1
   :caption: LEARN -- IV. Real Earth Observation

   t16-retrieving-eo-data
   t17-preparing-eo-observations
   t18-applying-inversion-spatially
   t19-trait-maps-uncertainty
   t20-end-to-end-workflow

.. toctree::
   :maxdepth: 1
   :caption: REFERENCE

   models
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
