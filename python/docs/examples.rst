Examples
========

Every snippet on this page has actually been run against the real
``toolsrtm``/``scopeinpython`` packages -- not hand-written pseudocode.
These mirror the R tutorial series (`ToolsRTM
<https://ccgcam.github.io/RTM-Suite/toolsrtm/articles/index.html>`_,
`SCOPEinR
<https://ccgcam.github.io/RTM-Suite/scopeinr/articles/index.html>`_)
topic-for-topic where a Python port exists; see each package's own
``README.md`` for the full R-tutorial-to-Python-module bridge table,
including the gaps called out at the bottom of this page.

Leaf + canopy (toolsrtm)
-------------------------

Mirrors R Tutorials 01-02 (``prospect_d`` = leaf optics, ``foursail`` =
canopy BRDF).

.. code-block:: python

   import numpy as np
   from toolsrtm import prospect_d, foursail

   leaf = prospect_d(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40)
   print(leaf.lambda_[:3], leaf.refl[:3], leaf.tran[:3])

   inputLUT = dict(
       N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
       Prot=0.002, CBC=0.007,          # only used if leaf_model='PROSPECT-PRO'
       LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
       LAI=3, hspot=0.01, tts=30, tto=0, psi=0,
   )
   rsoil = np.full(2101, 0.15)
   sail = foursail(inputLUT, rsoil, leaf_model="PROSPECT-D", spectrum_all=True)
   print("TOC bidirectional reflectance factor (rsot) at 550 nm:", sail.rsot[550 - 400])

.. figure:: _figures/leaf_canopy.png
   :alt: PROSPECT-D leaf reflectance/transmittance and fourSAIL canopy TOC reflectance, real output of the code above
   :width: 100%

   Real output of the code above: leaf-level optics (left) and the resulting canopy-level TOC reflectance (right).

Alternative leaf models: Liberty and Fluspect-B (toolsrtm)
-------------------------------------------------------------

``prospect_d``/``prospect_pro`` aren't the only leaf models -- ``liberty``
(conifer needles, Dawson et al. 1998) and ``fluspect_b`` (PROSPECT-D
optics plus the chlorophyll-fluorescence excitation-emission matrices SCOPE
needs) both work as drop-in leaf models for ``foursail``/``foursail2``/
``inform`` via their ``leaf_model="Liberty"``/``"Fluspect-B"`` argument
(Tutorial 02's leaf-model comparison table).

.. code-block:: python

   from toolsrtm import liberty, fluspect_b

   needle = liberty(cell_d=40, inter_c=0.045, baseline_abs=0.0006, leaf_thick=1.6,
                     albino_abs=0, Cab=40, EWT=0.01, lign_cell=2, Nitrogen=1)
   print("Liberty reflectance at 800nm:", round(float(needle.refl[800 - 400]), 4))

   flu = fluspect_b(Cab=40, Car=8, EWT=0.01, LMA=0.009, Cs=0, N=1.5, fqe=0.01, Cx=0)
   print("Fluspect-B reflectance at 800nm:", round(float(flu.refl[800 - 400]), 4))
   print("Backward fluorescence matrix shape (PSI):", flu.MbI.shape)

.. figure:: _figures/liberty_leaf.png
   :alt: LIBERTY conifer-needle leaf reflectance and 1-transmittance, real output of the code above
   :width: 75%

   Real output: LIBERTY's conifer-needle optics -- flatter NIR plateau and different SWIR absorption shape than broadleaf PROSPECT, reflecting the needle-specific anatomy the model targets.

.. figure:: _figures/fluspect_leaf.png
   :alt: Fluspect-B leaf optics and fluorescence excitation-emission matrix, real output of the code above
   :width: 100%

   Real output: Fluspect-B's leaf reflectance/transmittance (left, near-identical to PROSPECT-D since it shares the same absorption physics) and its backward chlorophyll-fluorescence excitation-emission matrix (right) -- the two characteristic emission peaks near 685nm (PSII) and 740nm (PSI) are visible at both the blue (~440nm) and red (~660-680nm) chlorophyll excitation bands.

Soil (BSM) + optical canopy BRDF (scopeinpython)
-------------------------------------------------

Mirrors SCOPEinR Tutorials 01-02 (soil model + canopy optics via
``rtmo``).

.. code-block:: python

   import numpy as np
   from toolsrtm import prospect_d
   from toolsrtm.canopy import dladgen
   from scopeinpython import SoilParams, WettingParams, get_bsm, CanopyStructure, get_spectra_scope, run_rtmo

   spectral = get_spectra_scope()
   rsoil = get_bsm(SoilParams(BSMBrightness=0.5, BSMlat=25, BSMlon=45),
                    WettingParams(SMp=15, SMC=25, film=0.015))

   leaf = prospect_d(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40)
   refl_leaf, tran_leaf = leaf.refl[:2001], leaf.tran[:2001]

   lidf = dladgen(-0.35, -0.15).lidf
   canopy = CanopyStructure(LAI=3, lidf=lidf, hot=0.1 / 2.0)

   # Esun_/Esky_ must be supplied by the caller -- see python/README.md
   result = run_rtmo(
       spectral=spectral, leaf_refl=refl_leaf, leaf_tran=tran_leaf,
       rho_thermal=0.01, tau_thermal=0.01, rsoil=rsoil, canopy=canopy,
       tts=30, tto=0, psi=0, Esun_=Esun_, Esky_=Esky_,
   )
   print("TOC reflectance at 550 nm:", result.refl[550 - 400])

Sensor convolution + vegetation indices (toolsrtm)
-----------------------------------------------------

Mirrors R Tutorials 07-09: convolve a simulated hyperspectral canopy
spectrum onto real Sentinel-2A band spectral response functions, then
compute vegetation indices from the convolved bands.

.. code-block:: python

   import numpy as np
   from toolsrtm import prospect_d, foursail
   from toolsrtm.srf import srf_sentinel2a, spectral_convolution_srf
   from toolsrtm.indices import get_indices

   inputLUT = dict(
       N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
       LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
       LAI=3, hspot=0.01, tts=30, tto=0, psi=0,
   )
   rsoil = np.full(2101, 0.15)
   sail = foursail(inputLUT, rsoil, leaf_model="PROSPECT-D", spectrum_all=True)
   wave = np.arange(400, 2501)

   s2a = srf_sentinel2a()
   conv = spectral_convolution_srf(wave, sail.rsot, s2a)
   print("Sentinel-2A bands:", s2a.band_names)
   print("Convolved TOC reflectance:", np.round(conv.rfl, 4))

   indices = get_indices(conv.wl, conv.rfl, spectral_domain="VNIR")
   for name in ("NDVI", "MSAVI", "REP"):
       print(name, "=", round(float(indices[name][0]), 4))

.. figure:: _figures/sensor_convolution.png
   :alt: Native 1nm canopy spectrum with Sentinel-2A band values overlaid, real output of the code above
   :width: 75%

   Real output: the native 1nm spectrum (grey) and the same spectrum convolved onto Sentinel-2A's band spectral response functions (red points).

INFORM: explicit forest-canopy model (toolsrtm)
----------------------------------------------------

``inform`` (Atzberger, forest-stand extension of fourSAIL) adds explicit
tree-crown geometry -- stem density, crown diameter, tree height,
understorey LAI -- on top of the same leaf models. Mirrors R Tutorial 02's
``foursail``/``foursail2``/``inform`` comparison.

.. code-block:: python

   import numpy as np
   from toolsrtm import foursail, inform

   inputLUT = dict(
       N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
       Prot=0.002, CBC=0.007, LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
       LAI=3, hspot=0.01, tts=30, tto=0, psi=0,
       LAIu=0.5, sd=650, cd=4.5, h=20, skyl=0.1,   # INFORM-only: understorey LAI, stem density, crown diameter, tree height, diffuse-light fraction
   )
   rsoil = np.full(2101, 0.15)
   r_forest = inform(inputLUT, rsoil, leaf_model="PROSPECT-D")
   r_homogeneous = foursail(inputLUT, rsoil, leaf_model="PROSPECT-D", spectrum_all=True).rsot
   print("fourSAIL (homogeneous) at 800nm:", round(float(r_homogeneous[800 - 400]), 4))
   print("INFORM (forest stand) at 800nm:", round(float(r_forest[800 - 400]), 4))

.. figure:: _figures/inform_forest.png
   :alt: fourSAIL vs INFORM TOC reflectance at the same leaf and LAI, real output of the code above
   :width: 75%

   Real output: same leaf optics and LAI through both models -- INFORM's explicit crown/gap geometry produces lower reflectance than a homogeneous fourSAIL canopy, matching the expected physics of a discontinuous forest stand.

Machine-learning trait inversion (toolsrtm)
-----------------------------------------------

Mirrors R Tutorials 11-12: build a small LUT, extract Sentinel-2-like
bands, and invert Cab with a PLSR model (``get_inversion`` dispatches to
12 algorithms in total -- see :data:`toolsrtm.inversion.ALGORITHMS`).

.. code-block:: python

   import numpy as np
   import pandas as pd
   from toolsrtm import prospect_d, foursail
   from toolsrtm.inversion import get_inversion

   rng = np.random.default_rng(1)
   rows = []
   for _ in range(200):
       Cab, LAI = rng.uniform(10, 80), rng.uniform(0.5, 6)
       inputLUT = dict(
           N=1.5, Cab=Cab, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
           LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
           LAI=LAI, hspot=0.01, tts=30, tto=0, psi=0,
       )
       sail = foursail(inputLUT, np.full(2101, 0.15), leaf_model="PROSPECT-D", spectrum_all=True)
       row = {"Cab": Cab, "LAI": LAI}
       for wl in (490, 560, 665, 705, 740, 783, 842, 865, 1610, 2190):
           row[f"R{wl}"] = sail.rsot[wl - 400]
       rows.append(row)

   df = pd.DataFrame(rows)
   band_cols = [c for c in df.columns if c.startswith("R")]

   result = get_inversion(df, dep_var="Cab", inputs=band_cols, algorithm="PLSR", n_samples=200, seed=1)
   print("Test R2:", round(result.statistics["test"]["r2"], 3))
   print("Test RMSE:", round(result.statistics["test"]["rmse"], 3))

.. figure:: _figures/ml_inversion.png
   :alt: Observed vs. predicted Cab scatter plot on the held-out test set, real output of the code above
   :width: 55%

   Real output: predicted vs. observed Cab on the held-out test set (R2=0.981).

SPART: full soil-plant-atmosphere chain (toolsrtm)
-------------------------------------------------------

Mirrors R Tutorial 03: ``spart_toa`` chains ``foursail`` with the BSM soil
model and SMAC atmospheric correction, returning both top-of-canopy and
top-of-atmosphere reflectance already resampled to a real sensor's bands
(Sentinel-2A here) -- there's no separate "simulate native, then convolve"
step, unlike plain ``foursail``.

.. code-block:: python

   from toolsrtm import spart_toa, sentinel2a_msi

   inputLUT = dict(
       N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
       Prot=0.002, CBC=0.007, LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
       LAI=3, hspot=0.01, tts=30, tto=0, psi=0,
       Pa=1000, aot550=0.3246, uo3=0.3480, uh2o=1.4116,   # atmosphere: pressure, aerosol optical thickness, ozone, water vapour
   )
   result = spart_toa(inputLUT, sensor=sentinel2a_msi(), leaf_model="PROSPECT-PRO",
                       BSMBrightness=0.5, BSMlat=25, BSMlon=45, SMp=15)
   print("Sentinel-2A band centers (nm):", result.wl_smac)
   print("TOC (canopy BRDF):", result.rfl_toc_brdf.round(4))
   print("TOA (after SMAC atmospheric correction):", result.rfl_toa.round(4))

.. figure:: _figures/spart_toc_toa.png
   :alt: SPART TOC and TOA reflectance across Sentinel-2A bands, real output of the code above
   :width: 75%

   Real output: TOC and TOA reflectance diverge sharply in the water-vapour bands (~940/1370nm), where the atmosphere absorbs most of the signal before it reaches the sensor -- exactly what real atmospheric correction has to undo.

Deep-learning trait inversion (toolsrtm, optional ``dl`` extra)
------------------------------------------------------------------

Mirrors R Tutorial 13: a Keras dense network (matching R's own
``getMLmodel``/``getMLmodel.withRetrain`` layer sizes, dropout placement,
and 7-optimizer choice) trained to invert Cab from the same 10
Sentinel-2-like bands used in the PLSR example above. Needs
``pip install "toolsrtm[dl]"`` (TensorFlow); nothing else in the package
requires it. Adam at R's own conservative default learning rate (1e-4)
converges slowly, so this needs a generous epoch budget and several
random restarts (``n_times``) to land a good held-out fit -- the same
trade-off the test suite itself budgets for.

.. code-block:: python

   import numpy as np
   import pandas as pd
   from toolsrtm import foursail
   from toolsrtm.deep_learning import get_ml_model

   rng = np.random.default_rng(2)
   rows = []
   for _ in range(600):
       Cab, LAI = rng.uniform(10, 80), rng.uniform(0.5, 6)
       inputLUT = dict(
           N=1.5, Cab=Cab, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
           LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
           LAI=LAI, hspot=0.01, tts=30, tto=0, psi=0,
       )
       sail = foursail(inputLUT, np.full(2101, 0.15), leaf_model="PROSPECT-D", spectrum_all=True)
       row = {"Cab": Cab, "LAI": LAI}
       for wl in (490, 560, 665, 705, 740, 783, 842, 865, 1610, 2190):
           row[f"R{wl}"] = sail.rsot[wl - 400]
       rows.append(row)

   df = pd.DataFrame(rows)
   band_cols = [c for c in df.columns if c.startswith("R")]

   result = get_ml_model(df, dep_var="Cab", model="Hidden-layers",
                          n_epochs=500, n_times=3, seed=2)
   print("Held-out R2:", round(result.stats["r2"], 3))
   print("Held-out RMSE:", round(result.stats["rmse"], 2))

.. figure:: _figures/deep_learning_inversion.png
   :alt: Keras training loss curve and observed vs predicted Cab scatter plot, real output of the code above
   :width: 100%

   Real output: training/validation loss over 500 epochs (left) and predicted vs. observed Cab on the held-out split (right, R2=0.882, RMSE=7.00 ug/cm2).

MARMIT soil moisture model (toolsrtm)
------------------------------------------

Mirrors R Tutorial 16: build a wetted soil spectrum from a dry
reference and couple it into a canopy simulation.

.. code-block:: python

   from toolsrtm.marmit import get_marmit_rsoil
   from toolsrtm import prospect_d, foursail

   soil = get_marmit_rsoil(soil_id=3, L=0.05, eps=0.4, version="marmit1")
   print("SMC (soil moisture content):", round(float(soil.smc), 4))
   for wl in (550, 850, 1600):
       i = wl - 400
       print(f"{wl}nm: dry={soil.rsoil_dry[i]:.4f}  wet={soil.rsoil_wet[i]:.4f}")

   inputLUT = dict(
       N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
       LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
       LAI=1.5, hspot=0.01, tts=30, tto=0, psi=0,
   )
   sail = foursail(inputLUT, soil.rsoil_wet, leaf_model="PROSPECT-D", spectrum_all=True)
   print("Canopy TOC reflectance at 850nm with wet MARMIT soil:", round(float(sail.rsot[850 - 400]), 4))

.. figure:: _figures/marmit_soil.png
   :alt: Dry vs. wet soil reflectance spectrum from MARMIT, real output of the code above
   :width: 75%

   Real output: MARMIT's dry-reference vs. wetted soil reflectance -- the SWIR water-absorption dips (~1400/1900nm) deepen and overall brightness drops as the soil wets.

Full SCOPE run: energy balance + fluorescence (scopeinpython)
-------------------------------------------------------------------

Mirrors SCOPEinR Tutorials 03-04: one full ``get_scope()`` call --
optics, energy balance, photosynthesis and fluorescence together --
against the same bundled example LUT row SCOPEinR's own test suite and
R vignettes use (``SCOPEinR/inst/input/LUT_input.csv``).

.. code-block:: python

   import csv
   from pathlib import Path
   from scopeinpython import ScopeOptions, get_scope

   with open("SCOPEinR/inst/input/LUT_input.csv", newline="") as f:
       row = next(csv.DictReader(f))

   res = get_scope(row, options=ScopeOptions(k_maxit=100, maxEBer=1.0))
   print("Canopy layers:", res.nlayers)
   print("TOC reflectance at 550/700/850 nm:",
         round(float(res.rtmo.refl[550 - 400]), 4),
         round(float(res.rtmo.refl[700 - 400]), 4),
         round(float(res.rtmo.refl[850 - 400]), 4))
   print("Net radiation, total (Rntot):", round(float(res.ebal.Rntot), 2), "W/m2")
   print("Total photosynthesis (Actot):", round(float(res.ebal.Actot), 2), "umol CO2/m2/s")
   if res.rtmf is not None:
       print("Emitted fluorescence (EoutF):", round(float(res.rtmf.EoutF), 4), "W/m2/sr")

.. figure:: _figures/scope_full.png
   :alt: Full SCOPE TOC reflectance and emitted SIF spectrum, real output of the code above
   :width: 100%

   Real output of the single ``get_scope()`` call above: TOC reflectance (left; the dashed gaps are the water-vapor-absorption wavelengths SCOPE itself leaves undefined) and the emitted SIF spectrum (right).

Full simulate -> indices -> ML-invert pipeline
------------------------------------------------

A complete, actually-executed pipeline (100 samples, spectral indices,
scikit-learn trait inversion, real R² 0.7-0.9) lives outside this
package as plain scripts + a Jupyter notebook, kept separate from the R
scripts:

``Scripts/Python/ForPROSAIL_fourSAIL/`` — ``1_simulate_lut.py`` →
``2_spectral_indices.py`` → ``3_inversion_ml.py``, plus ``pipeline.ipynb``.
See ``Scripts/Python/README.md`` for how to run it.

What isn't ported yet
----------------------

Two tutorial-level gaps not covered by any example above, on top of the
finer-grained implementation gaps in :doc:`not_ported`:

- **Sobol/Johnson sensitivity analysis** (``toolsrtm``, R Tutorial 10's
  ``get.sobol.indices()``/``get.spectral.sensitivity()``) and **LUT
  correlation-distribution helpers** (``get_distributionLUT()``/``getCor()``,
  R Tutorial 05) have no Python module yet.
- **The real-Sentinel-2/STAC capstone** (SCOPEinR Tutorial 11 --
  ``Actot`` retrieved from a real satellite time series and mapped
  spatially; ToolsRTM Tutorial 18 -- a spatial index map from a real
  STAC-retrieved image) has no ``scopeinpython``/``toolsrtm`` equivalent
  yet, even though ``toolsrtm.satellite`` already has the STAC retrieval
  machinery a Python port of it would reuse.

See each package's own ``README.md`` for the full R-tutorial-to-Python-module
bridge table.
