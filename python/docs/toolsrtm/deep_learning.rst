toolsrtm.deep_learning
========================

Deep-learning trait inversion: dense ("Hidden-layers") and 1D-CNN Keras
architectures with a configurable optimizer. Direct port of
``ToolsRTM::getMLmodel``/``getMLmodel.withRetrain``.

.. note::
   Optional -- needs the ``dl`` extra: ``pip install "toolsrtm[dl]"``
   (TensorFlow). Not required for the rest of the package;
   :mod:`toolsrtm.inversion`'s scikit-learn-based dispatcher covers most
   trait-inversion needs without it.

Quick example
-------------

.. code-block:: python

   import numpy as np, pandas as pd
   from toolsrtm import foursail
   from toolsrtm.deep_learning import get_ml_model

   rng = np.random.default_rng(2)
   rows = []
   for _ in range(600):
       Cab, LAI = rng.uniform(10, 80), rng.uniform(0.5, 6)
       inputLUT = dict(N=1.5, Cab=Cab, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
                        LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
                        LAI=LAI, hspot=0.01, tts=30, tto=0, psi=0)
       sail = foursail(inputLUT, np.full(2101, 0.15), leaf_model="PROSPECT-D", spectrum_all=True)
       row = {"Cab": Cab, "LAI": LAI}
       for wl in (490, 560, 665, 705, 740, 783, 842, 865, 1610, 2190):
           row[f"R{wl}"] = sail.rsot[wl - 400]
       rows.append(row)
   df = pd.DataFrame(rows)
   band_cols = [c for c in df.columns if c.startswith("R")]

   result = get_ml_model(df, dep_var="Cab", model="Hidden-layers", n_epochs=500, n_times=3, seed=2)
   print(result.stats["r2"])   # held-out R2

.. code-block:: text

   Input                              get_ml_model()            Output
   ---------------------------        ----------------------    ---------------------------
   df   [n rows]  LUT: predictor                                 result.model     fitted Keras model
                  bands + dep_var     -------------------->      result.x_scaler  fitted StandardScaler
   dep_var  = trait to invert                                     (re-apply to new X before
   model = "Hidden-layers"/"CNN"                                    predict() -- see the R
   n_epochs, n_times, seed                                          Tutorial 13 scaling-bug story)

.. note::
   ``result.x_scaler`` must be applied to any new predictor data before
   calling ``result.model.predict(...)`` -- training happens in scaled
   space, so predicting on raw reflectance directly produces silently
   wrong (often catastrophically bad) results. This is exactly the bug
   documented and fixed in `ToolsRTM Tutorial 13
   <https://ccgcam.github.io/RTM-Suite/toolsrtm/articles/t13-deep-learning-inversion.html>`_.

.. automodule:: toolsrtm.deep_learning
   :members:
   :undoc-members:
   :show-inheritance:
