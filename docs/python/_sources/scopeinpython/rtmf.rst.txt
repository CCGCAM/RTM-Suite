scopeinpython.rtmf
====================

SCOPE canopy fluorescence radiative transfer model. Direct port of
``SCOPEinR/R/RTMf.R`` (``get.RTMf``), computing the TOC fluorescence
radiance in the observation direction and the TOC hemispherical upward
fluorescence flux, given ``Mb``/``Mf`` (from
:func:`scopeinpython.fluspect_mscope.fluspect_mscope`) and per-layer
fluorescence quantum efficiencies (from
:func:`scopeinpython.biochemical.get_biochemical`) as composable inputs.

.. note::
   The final upsampling from the native 53-point fluorescence
   wavelength grid to the 211-point display grid (``spectral.wlF``) uses
   ``scipy``'s ``not-a-knot`` cubic spline, a close but not bit-identical
   approximation of R's ``fmm``-method spline
   (``signal::interp1(...,'spline')``) — see the module docstring for
   measured error bounds. This is the one place in the whole Python port
   with a deliberate, documented sub-floating-point-precision
   approximation.

.. automodule:: scopeinpython.rtmf
   :members:
   :undoc-members:
   :show-inheritance:
