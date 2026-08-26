Simulation Workflows
========================

:doc:`models` covers what each model *is* in isolation. This page is how
they chain together into the four standard simulation pipelines used
throughout this site -- each diagram is the mechanism, with a minimal
snippet showing the chain; full runnable versions with plots live on
:doc:`examples`.

1. PROSPECT -> fourSAIL -> canopy reflectance
--------------------------------------------------

The baseline pipeline: one leaf model feeds one canopy model, for a
single homogeneous crop/grassland-type canopy over a flat soil.

.. code-block:: text

   Leaf traits (N, Cab, Car, Anth, EWT, LMA, alpha)
                 |
                 v
        prospect_d() / prospect_pro()
                 |
                 v
        Leaf reflectance (rho) / transmittance (tau)
                 |
                 |    Canopy structure (LAI, LIDFa, LIDFb, hspot)
                 |    Sun-view geometry (tts, tto, psi)
                 |    Soil background (rsoil)
                 v
              foursail()
                 |
                 v
        Top-of-canopy (TOC) reflectance, 400-2500nm

.. code-block:: python

   from toolsrtm import prospect_d, foursail
   leaf = prospect_d(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40)
   sail = foursail(dict(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
                         LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1, LAI=3, hspot=0.01, tts=30, tto=0, psi=0),
                    rsoil=[0.15] * 2101, leaf_model="PROSPECT-D", spectrum_all=True)
   # sail.rsot is the resulting TOC reflectance spectrum

Full version with a plot: :ref:`examples:Leaf + canopy (toolsrtm)`.

2. PROSPECT -> INFORM -> forest reflectance
-------------------------------------------------

Same leaf model, but the canopy step adds explicit tree-crown geometry --
stem density, crown diameter, tree height, understorey LAI -- instead of
treating the whole canopy as one homogeneous layer.

.. code-block:: text

   Leaf traits  ->  prospect_d()  ->  Leaf reflectance/transmittance
                                              |
                          Understorey LAI, stem density,          |
                          crown diameter, tree height, diffuse-    |
                          light fraction (skyl)                    v
                                                                inform()
                                                                    |
                                                                    v
                                              TOC reflectance -- lower than
                                              fourSAIL at the same LAI (real
                                              crown/gap shadow geometry)

.. code-block:: python

   from toolsrtm import inform
   r_forest = inform(dict(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
                           LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1, LAI=3, hspot=0.01, tts=30, tto=0, psi=0,
                           LAIu=0.5, sd=650, cd=4.5, h=20, skyl=0.1),
                      rsoil=[0.15] * 2101, leaf_model="PROSPECT-D")

Full version comparing fourSAIL vs. INFORM side by side:
:ref:`examples:INFORM: explicit forest-canopy model (toolsrtm)`.

3. MARMIT -> fourSAIL -> SPART -> TOA
-------------------------------------------

The full soil-plant-atmosphere chain: a real, moisture-dependent soil
spectrum feeds the canopy model, and SPART's own atmosphere step (SMAC)
converts the result to what a real satellite would measure above the
atmosphere -- not just top-of-canopy.

.. code-block:: text

   Dry reference soil spectrum, water-film params (L, eps)
                 |
                 v
          get_marmit_rsoil()
                 |
                 v
          Wetted soil reflectance  ---+
                                       |
   Leaf + canopy traits                v
                 |               (used as rsoil)
                 v                     |
            foursail()  <--------------+
                 |
                 v
          TOC reflectance
                 |
                 |    Atmosphere (Pa, aot550, uo3, uh2o), sensor (Sentinel-2A)
                 v
             spart_toa()  (re-derives TOC via BSM + adds SMAC atmosphere)
                 |
                 v
          TOA reflectance, per sensor band

.. code-block:: python

   from toolsrtm import get_marmit_rsoil, spart_toa, sentinel2a_msi
   soil = get_marmit_rsoil(soil_id=3, L=0.05, eps=0.4, version="marmit1")
   result = spart_toa(dict(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
                            LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1, LAI=3, hspot=0.01, tts=30, tto=0, psi=0,
                            Pa=1000, aot550=0.3246, uo3=0.348, uh2o=1.4116),
                       sensor=sentinel2a_msi(), leaf_model="PROSPECT-PRO",
                       BSMBrightness=0.5, BSMlat=25, BSMlon=45, SMp=15)
   # result.rfl_toc_brdf (TOC) and result.rfl_toa (TOA) both per-band

Note: ``spart_toa`` builds its own soil internally via BSM parameters
(``BSMBrightness``/``BSMlat``/``BSMlon``/``SMp``), a separate soil model
from MARMIT -- see :doc:`models`'s Soil section for how the two differ.
Full versions: :ref:`examples:MARMIT soil moisture model (toolsrtm)`,
:ref:`examples:SPART: full soil-plant-atmosphere chain (toolsrtm)`.

4. SCOPE -> reflectance + SIF + energy balance
----------------------------------------------------

A structurally different pipeline: instead of one forward pass through
leaf -> canopy, SCOPE iterates leaf/soil temperature until the energy
balance closes, then derives everything else (fluorescence, carbon flux)
from that solved state. See :doc:`models`'s SCOPE section for what each
of the five internal components does.

.. code-block:: text

   One LUT row: leaf biochemistry + canopy structure + soil + meteorology
                              |
                              v
              get_fluspect_cx_scope()  (leaf optics + fluorescence matrices)
                              |
                              v
                         run_rtmo()  (optical TOC BRDF)
                              |
                              v
                    ebal()  <-----------------------+
                    (iterates leaf/soil temperature   |
                     until flux budget closes,         |
                     calling get_biochemical() at       get_biochemical()
                     each candidate temperature)  ------+  (photosynthesis,
                              |                             fluorescence yield)
                              v
                    rtmf() / rtmz()  (optional: canopy fluorescence, zeaxanthin)
                              |
                              v
          Reflectance + SIF + leaf/soil temperature + carbon/water flux

.. code-block:: python

   import csv
   from scopeinpython import ScopeOptions, get_scope
   with open("SCOPEinR/inst/input/LUT_input.csv", newline="") as f:
       row = next(csv.DictReader(f))
   res = get_scope(row, options=ScopeOptions(k_maxit=100, maxEBer=1.0))
   # res.rtmo.refl (reflectance), res.ebal.Actot (photosynthesis),
   # res.rtmf.EoutF (fluorescence, if requested)

Full version: :ref:`examples:Full SCOPE run: energy balance + fluorescence (scopeinpython)`.

What's next
-----------------

- :doc:`sensor_simulation` -- resampling any of the reflectance spectra
  above onto a real sensor's bands.
- :doc:`trait_inversion` -- going the other direction: from reflectance
  (simulated or real) back to a trait estimate.
- :doc:`earth_observation` -- applying a trained inversion model to a
  real Sentinel-2 scene.
