Examples
========

Leaf + canopy (toolsrtm)
-------------------------

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

Soil (BSM) + optical canopy BRDF (scopeinpython)
-------------------------------------------------

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

Full simulate -> indices -> ML-invert pipeline
------------------------------------------------

A complete, actually-executed pipeline (100 samples, spectral indices,
scikit-learn trait inversion, real R² 0.7-0.9) lives outside this
package as plain scripts + a Jupyter notebook, kept separate from the R
scripts:

``Scripts/Python/ForPROSAIL_fourSAIL/`` — ``1_simulate_lut.py`` →
``2_spectral_indices.py`` → ``3_inversion_ml.py``, plus ``pipeline.ipynb``.
See ``Scripts/Python/README.md`` for how to run it.
