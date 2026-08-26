Trait Inversion
===================

The whole point of the models on :doc:`models`, running forward
(traits -> reflectance), is usually to go **backward**: given a real or
simulated reflectance spectrum, recover the trait that produced it. This
package offers three genuinely different ways to do that, from simplest
to most flexible:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Method
     - How it works
     - Best when
   * - LUT inversion (merit function)
     - Find the LUT row(s) whose simulated spectrum most closely matches
       the observation, average their trait values.
     - You already have (or can quickly build) a well-sampled LUT and
       want an interpretable, no-training-required retrieval.
   * - Machine learning
     - Fit a regression model (PLSR/RF/SVM/...) on (spectrum -> trait)
       pairs from a LUT, then predict on new spectra.
     - You have a reasonably-sized LUT and want a fast, reusable trained
       model -- the most common choice.
   * - Deep learning
     - A neural network (dense or 1D-CNN), same idea as ML but able to
       learn more complex, nonlinear structure, especially across many
       contiguous hyperspectral bands.
     - Hyperspectral data (many contiguous bands) and a large-enough LUT
       to train on; needs the optional ``dl`` extra.

All three consume a LUT built as in :doc:`lut_generation`.

1. LUT inversion (merit function)
--------------------------------------

:func:`~toolsrtm.inversion.get_inversion_opt` ranks every LUT row by how
well its simulated spectrum matches an observation (a merit/error
function -- RMSE by default), then averages the ``n_opt`` best matches'
trait values. No training step at all -- just a LUT and a comparison:

.. code-block:: python

   import numpy as np
   import pandas as pd
   from toolsrtm import foursail
   from toolsrtm.inversion import get_inversion_opt

   rng = np.random.default_rng(1)
   rows = []
   for _ in range(300):
       Cab, LAI = rng.uniform(10, 80), rng.uniform(0.5, 6)
       inputLUT = dict(N=1.5, Cab=Cab, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
                        LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
                        LAI=LAI, hspot=0.01, tts=30, tto=0, psi=0)
       sail = foursail(inputLUT, np.full(2101, 0.15), leaf_model="PROSPECT-D", spectrum_all=True)
       row = {"Cab": Cab, "LAI": LAI}
       for wl in (490, 560, 665, 705, 740, 783, 842, 865, 1610, 2190):
           row[f"R{wl}"] = sail.rsot[wl - 400]
       rows.append(row)
   lut_df = pd.DataFrame(rows)
   band_cols = [c for c in lut_df.columns if c.startswith("R")]

   # one held-out "observation" (in practice: a real convolved sensor spectrum)
   obs_sail = foursail(dict(N=1.5, Cab=42.0, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
                             LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
                             LAI=3.2, hspot=0.01, tts=30, tto=0, psi=0),
                        np.full(2101, 0.15), leaf_model="PROSPECT-D", spectrum_all=True)
   obs_refl = np.array([[obs_sail.rsot[wl - 400] for wl in (490, 560, 665, 705, 740, 783, 842, 865, 1610, 2190)]])

   result = get_inversion_opt(rfl_sensor=obs_refl, rfl_rtm=lut_df[band_cols].values,
                               lut=lut_df[["Cab", "LAI"]], method="merit-RMSE", n_opt=5)
   print(result.lut_best["Cab"].iloc[0], result.lut_best["LAI"].iloc[0])   # retrieved Cab, LAI

``method`` also accepts ``"merit-NRMSE"``, ``"merit-MAE"``, ``"merit-NMB"``,
``"merit-FGE"``, or a fully custom ``custom_stat=f(sim, obs) -> error``.
``n_opt`` trades off noise (``n_opt=1``, most exact match, most sensitive
to a single lucky/unlucky LUT row) against bias (larger ``n_opt``,
smoother but pulled toward the LUT's own average).

2. Machine learning
------------------------

:func:`~toolsrtm.inversion.get_inversion` dispatches to 12 algorithms
(PLSR/SVM/RF/GB/NN/Bayesian/AdaBag/BRNN/xGB/RVM/qLASSO/Ensemble -- see
:data:`~toolsrtm.inversion.ALGORITHMS` for exactly which scikit-learn/
xgboost estimator each maps to) behind one interface:

.. code-block:: python

   from toolsrtm.inversion import get_inversion

   lut_df["Cab"] = [r["Cab"] for r in rows]   # already in lut_df above
   fit = get_inversion(lut_df, dep_var="Cab", inputs=band_cols, algorithm="PLSR",
                        n_samples=len(lut_df), seed=1)
   print(fit.statistics["test"]["r2"], fit.statistics["test"]["rmse"])

Swap ``algorithm="RF"``/``"SVM"``/... for any of the 12; the held-out
train/test split and evaluation stats are computed the same way
regardless of which one is chosen.

3. Deep learning
---------------------

:func:`toolsrtm.deep_learning.get_ml_model` -- a dense network or a 1D-CNN
(needs ``pip install "toolsrtm[dl]"``). The 1D-CNN treats predictor bands
*in spectral order*, so it benefits from wide, contiguous hyperspectral
input rather than a handful of scattered multispectral bands:

.. code-block:: python

   from toolsrtm.deep_learning import get_ml_model

   result = get_ml_model(lut_df, dep_var="Cab", model="Hidden-layers",
                          n_epochs=500, n_times=3, seed=2)
   print(result.stats["r2"])

.. important::
   ``result.x_scaler`` (a fitted ``sklearn.StandardScaler``) must be
   applied to any new predictor data before calling
   ``result.model.predict(...)`` -- training happens in scaled space, so
   predicting on raw reflectance directly produces silently wrong
   results. This is a real bug found and fixed in `ToolsRTM Tutorial 13
   <https://ccgcam.github.io/RTM-Suite/toolsrtm/articles/t13-deep-learning-inversion.html>`_
   -- worth reading before relying on a DL model's predictions.

Full DL examples (dense + 1D-CNN, with training-curve and true-vs-
predicted figures): `Deep-learning trait inversion, on examples.rst
<examples.html#deep-learning-trait-inversion-toolsrtm-optional-dl-extra>`_.

Which one should I use?
----------------------------

- **Small LUT, need it fast, want it interpretable**: LUT inversion.
- **Default choice for most trait-retrieval work**: ML (``RF`` or
  ``PLSR`` are good starting points).
- **Wide hyperspectral input, large LUT, willing to tune a network**: DL,
  and only after confirming ML doesn't already do the job -- it usually
  does.

Spectral indices (:doc:`t09-spectral-indices`) are a fourth, even simpler
option worth trying first for a quick sanity check: rank indices by
correlation with the trait, and see how far a single-formula retrieval
already gets before reaching for any of the three methods above.

What's next
-----------------

- :doc:`lut_generation` -- building the LUT these methods all consume.
- :doc:`earth_observation` -- applying a trained ML/DL model to a real
  Sentinel-2 scene.
- :doc:`examples` -- every method above with full runnable code and
  result figures.
