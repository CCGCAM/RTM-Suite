scopeinpython.scope
=====================

The end-to-end SCOPE simulation wrapper: one LUT input row in, leaf optics
(Fluspect-Cx) → soil (BSM or the bundled reference spectra) → optical BRDF
(``run_rtmo``) → thermal energy balance (``ebal``) → fluorescence
(``rtmf``, optional) → zeaxanthin (``rtmz``, optional) out. Direct,
partial port of ``SCOPEinR::get.SCOPE``.

.. note::
   Verified against a real, **unmodified** ``SCOPEinR::get.SCOPE()`` call
   using the R package's own bundled example inputs
   (``SCOPEinR/inst/input/LUT_input.csv`` + ``setoptions.csv``) -- see
   ``python/scratch/scratch_scope_export.R``. The exact-formula outputs
   (``LAIsunlit``/``Pnsun_Cab``/``Pnsha_Cab``/``Pntot_Cab``, the TOC
   reflectance spectrum) match to floating-point noise; the
   iterative-convergence outputs (temperatures, energy-balance totals)
   match at the same ~1-2% tolerance established for :func:`scopeinpython.ebal.ebal`
   itself, for the same reason (small floating-point divergence compounding
   over ~7-10 nonlinear iterations along a per-layer-loop biochemistry
   path -- not a functional bug).

   See the module docstring below, and :doc:`../not_ported`, for the full
   list of ``options.SCOPE`` branches this wrapper does not expose at all
   (directional BRDF, ``RTMt_planck``, multi-layer mSCOPE, time-series
   mode, angle-file LIDF, measurement-file/MODTRAN irradiance) and which
   canopy-level "derived data products" beyond :class:`~scopeinpython.scope.ScopeResult`
   aren't computed yet.

Quick example
-------------

.. code-block:: python

   import csv
   from scopeinpython import ScopeOptions, get_scope

   with open("SCOPEinR/inst/input/LUT_input.csv", newline="") as f:
       row = next(csv.DictReader(f))

   res = get_scope(row, options=ScopeOptions(k_maxit=100, maxEBer=1.0))
   print(res.rtmo.refl[550 - 400])   # TOC reflectance at 550nm
   print(res.ebal.Actot)             # total photosynthesis, umol CO2/m2/s

.. code-block:: text

   Input                              get_scope()               Output
   ---------------------------        ----------------------    ---------------------------
   row : dict  (one LUT row -- leaf,                             res.rtmo   optical BRDF
     canopy, soil, meteo traits,       -------------------->      res.ebal   energy balance,
     matching inputs_SCOPE.csv)                                              photosynthesis
   options : ScopeOptions                                        res.rtmf   fluorescence (opt.)
                                                                   res.rtmz   zeaxanthin (opt.)

.. automodule:: scopeinpython.scope
   :members:
   :undoc-members:
   :show-inheritance:
