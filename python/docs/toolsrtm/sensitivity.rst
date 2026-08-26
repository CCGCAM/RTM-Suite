toolsrtm.sensitivity
=====================

Global sensitivity analysis (Johnson relative-importance index, verified to
reproduce R's ``sensitivity::johnson()`` to 8 decimal places; a simplified
two-sample-split Sobol-like ``Si``/``STi`` estimator) and correlated/
multi-distribution LUT builders. Direct port of
``ToolsRTM/R/get.sobol.indices.R``, ``get.spectral.sensitivity.R``,
``Correlated_value.R``, ``Gaussian_MinMax.R``, ``get_distributionLUT.R``,
``getCor.R``.

Quick example
-------------

.. code-block:: python

   from toolsrtm.sensitivity import spectral_sensitivity

   result = spectral_sensitivity(n_samples=500, distribution="Uniform",
                                  traits=("N", "Cab", "EWT", "LMA", "LIDFa", "LAI"),
                                  wl_step=5, seed=11)
   at_700nm = result.sti_pct[result.wavelength == 700]
   print(dict(zip(result.trait[result.wavelength == 700], at_700nm.round(1))))

.. code-block:: text

   Input                              spectral_sensitivity()    Output
   ---------------------------        ----------------------    ---------------------------
   n_samples (LUT runs per trait)                                long-format arrays:
   distribution = "Uniform"/...       -------------------->       result.wavelength
   traits = tuple of trait names                                  result.trait
   wl_step (nm, output resolution)                                result.sti_pct  (sums to 100%
                                                                     per wavelength -- relative
                                                                     contribution to variance)

.. automodule:: toolsrtm.sensitivity
   :members:
   :undoc-members:
   :show-inheritance:
