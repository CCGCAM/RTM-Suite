LUT Generation
==================

What is a LUT?
------------------

A look-up table (LUT) is a set of RTM simulations, each row a different
combination of trait values and its resulting reflectance spectrum -- the
raw material every method on :doc:`trait_inversion` needs: rank spectral
indices by correlation with a trait, train an ML/DL model, or match a
real observation against simulated neighbours. Building one well is
mostly about *how* traits are sampled, not the RTM call itself (that part
is just :doc:`workflows`, called once per row).

Parameter distributions
----------------------------

Each trait needs its own realistic range and shape, not one arbitrary
sweep -- :func:`~toolsrtm.sensitivity.get_distribution_lut` builds a
multi-trait LUT from a per-trait distribution choice in one call:

.. code-block:: python

   from toolsrtm.sensitivity import get_distribution_lut

   lut = get_distribution_lut(
       minval=dict(Cab=10, LAI=0.5), maxval=dict(Cab=80, LAI=6),
       n_samples=500, type_distrib=dict(Cab="Gaussian", LAI="Uniform"),
       mean_gauss=dict(Cab=45), std_gauss=dict(Cab=15), seed=1,
   )
   print(lut["Cab"].mean(), lut["Cab"].std())   # ~45, ~15 (Gaussian)
   print(lut["LAI"].min(), lut["LAI"].max())    # spans [0.5, 6] (Uniform)

``"Uniform"`` samples flat across ``[minval, maxval]``; ``"Gaussian"``
samples around ``mean_gauss``/``std_gauss`` while still clipped to
``[minval, maxval]`` (:func:`~toolsrtm.sensitivity.gauss_by_min_max`) --
the same two distributions SCOPE's own ``inputs_SCOPE.csv``-driven LUT
uses (see the SCOPE-LUT sampling example later in this chapter).

Parameter constraints: real trait co-variation, not independence
------------------------------------------------------------------------

Sampling every trait fully independently is the default, but real leaves
don't vary that way -- chlorophyll and carotenoids, or ``Cab`` and
``Vcmax25``, co-vary. Two ways to impose that:

.. code-block:: python

   # 1. Built into get_distribution_lut: Car tied to Cab (r=0.8), a fixed
   #    empirical relationship from leaf pigment data.
   lut = get_distribution_lut(
       minval=dict(Cab=10, Car=2), maxval=dict(Cab=80, Car=25),
       n_samples=500, type_distrib=dict(Cab="Uniform", Car="Uniform"),
       dep_cab=True, seed=1,
   )

   # 2. General-purpose: correlate any two traits at any target correlation.
   from toolsrtm.sensitivity import get_cor

   cor_res = get_cor(n_inputs=2, n_lut=500, distribution="Uniform", rho=0.85, seed=3,
                      var_names=["Cab", "Vcmax25"], min_range=[5, 5], max_range=[90, 250])
   lut["Cab"], lut["Vcmax25"] = cor_res.lut["Cab"], cor_res.lut["Vcmax25"]
   print(cor_res.covariance)   # realized correlation matrix, close to rho

RTM simulation: from a trait table to a spectral LUT
----------------------------------------------------------

Once the trait table exists, run the chosen model (:doc:`workflows`) once
per row -- the pattern every ML/DL example on this site uses:

.. code-block:: python

   import numpy as np
   from toolsrtm import foursail

   rows = []
   for i in range(len(lut["Cab"])):
       row_lut = dict(N=1.5, Cab=lut["Cab"][i], Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009,
                       alpha=40, LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
                       LAI=lut["LAI"][i], hspot=0.01, tts=30, tto=0, psi=0)
       sail = foursail(row_lut, np.full(2101, 0.15), leaf_model="PROSPECT-D", spectrum_all=True)
       rows.append(sail.rsot)
   spectral_lut = np.stack(rows)   # (n_samples, 2101)

Adding noise
----------------

Simulated reflectance is noise-free -- real sensor data never is. Adding
realistic noise before training an inversion model (:doc:`trait_inversion`)
avoids a model that looks good on simulated test data but degrades
sharply on real observations:

.. code-block:: python

   rng = np.random.default_rng(7)
   # Additive Gaussian, ~1% reflectance -- a simple, common sensor-noise proxy.
   noisy_lut = spectral_lut + rng.normal(0, 0.01, spectral_lut.shape)
   noisy_lut = np.clip(noisy_lut, 0, 1)   # reflectance can't be negative or >1

   # Or relative (multiplicative) noise, proportionally larger at low reflectance:
   noisy_lut_rel = spectral_lut * (1 + rng.normal(0, 0.03, spectral_lut.shape))

Train / validation / test
------------------------------

Every ML/DL example on this site (:doc:`trait_inversion`) splits the LUT
before fitting, never evaluates on rows the model trained on:

.. code-block:: python

   from sklearn.model_selection import train_test_split

   idx = np.arange(spectral_lut.shape[0])
   train_idx, test_idx = train_test_split(idx, test_size=0.3, random_state=1)
   X_train, X_test = spectral_lut[train_idx], spectral_lut[test_idx]
   y_train, y_test = lut["Cab"][train_idx], lut["Cab"][test_idx]

:func:`toolsrtm.deep_learning.get_ml_model` and
:func:`toolsrtm.inversion.get_inversion` do this split internally
(``n_samples``/an internal train/val ratio) -- see :doc:`trait_inversion`
for the exact behaviour of each.

What's next
-----------------

- :doc:`trait_inversion` -- what to do with a LUT once it's built: LUT
  matching, ML, or DL.
- :doc:`t02-parameters-traits` -- realistic ranges for every trait, so
  ``minval``/``maxval`` above aren't guessed.
- SCOPEinR's `Getting LUTs for SCOPE
  <https://ccgcam.github.io/RTM-Suite/scopeinr/articles/getting-luts-scope.html>`_
  -- the R-side equivalent, including SCOPE's own CSV-driven sampling.
