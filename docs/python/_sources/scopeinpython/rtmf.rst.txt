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

Where this fits
----------------

This is a pipeline-internal component, not something you typically call
directly: :func:`scopeinpython.scope.get_scope` calls this after
:func:`scopeinpython.ebal.ebal`, when fluorescence output is requested,
using the ``Mb``/``Mf`` matrices from
:func:`scopeinpython.fluspect_mscope.fluspect_mscope` and the per-layer
``eta`` from :func:`scopeinpython.biochemical.get_biochemical`. Most
users should just call :func:`~scopeinpython.scope.get_scope` -- see
:doc:`scope` for a full end-to-end example.

.. code-block:: text

   fluspect_mscope()  ->  Mb, Mf  (excitation-emission matrices)
   ebal()              ->  eta   (per-layer fluorescence yield)
                                  \
                                   -->  rtmf()  ->  TOC fluorescence radiance/flux

.. automodule:: scopeinpython.rtmf
   :members:
   :undoc-members:
   :show-inheritance:
