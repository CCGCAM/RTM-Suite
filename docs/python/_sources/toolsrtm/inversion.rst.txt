toolsrtm.inversion
====================

Trait-inversion tools: CARS-PLS and VIF predictor selection, LUT nearest-
neighbour ("merit function") matching, and a 12-algorithm ML dispatcher
(PLSR/SVM/RF/GB/NN/Bayesian/AdaBag/BRNN/xGB/RVM/qLASSO/Ensemble) built on
scikit-learn/xgboost. Direct port of ``ToolsRTM::carspls``/``get.cars.pls``,
``getVIF``, ``get.inversionOpt``, ``get.inversion``, ``hybrid_inversion``/
``hybrid_inversionE``.

.. note::
   Needs the optional ``ml`` extra: ``pip install "toolsrtm[ml]"``. Nothing
   in this module is imported by ``toolsrtm/__init__.py``'s own import chain
   -- a plain ``import toolsrtm`` never requires scikit-learn/xgboost.

.. note::
   R's ``get.inversion``/``hybrid_inversion`` dispatch to specific ``caret``
   methods (``bartMachine``, ``rqlasso``, ``rvmLinear``, ``AdaBag``, ``brnn``,
   ...). See :data:`~toolsrtm.inversion.ALGORITHMS` for exactly which
   scikit-learn/xgboost estimator each algorithm name maps to, and, where
   there's no direct equivalent, what was substituted and why.

Quick example
-------------

.. code-block:: python

   import numpy as np, pandas as pd
   from toolsrtm import foursail
   from toolsrtm.inversion import get_inversion

   rng = np.random.default_rng(1)
   rows = []
   for _ in range(200):
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

   result = get_inversion(df, dep_var="Cab", inputs=band_cols, algorithm="PLSR", n_samples=200, seed=1)
   print(result.statistics["test"]["r2"])   # held-out test R2

.. code-block:: text

   Input                              get_inversion()           Output
   ---------------------------        ----------------------    ---------------------------
   df   [n rows]  LUT: predictor                                 result.model      fitted estimator
                  bands + dep_var     -------------------->      result.statistics  train/test R2, RMSE
   dep_var  = trait to invert (e.g. "Cab")                       result.predictions test-set predicted
   inputs   = predictor column names                                                 vs. observed
   algorithm = "PLSR"/"RF"/"SVM"/... (see ALGORITHMS)

.. automodule:: toolsrtm.inversion
   :members:
   :undoc-members:
   :show-inheritance:
