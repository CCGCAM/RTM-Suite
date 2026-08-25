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

.. automodule:: scopeinpython.scope
   :members:
   :undoc-members:
   :show-inheritance:
