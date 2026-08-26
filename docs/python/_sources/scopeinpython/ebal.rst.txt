scopeinpython.ebal
====================

The SCOPE energy-balance closure loop: iterates sunlit/shaded leaf and
soil temperature until sensible+latent heat flux matches net radiation,
coupling radiative transfer, aerodynamics, photosynthesis/fluorescence
and heat-flux partitioning. Direct, **"SCOPE-lite"-only** port of
``SCOPEinR::get.ebal`` -- the top of the whole thermal energy-balance
chain built this session.

.. note::
   Depends on :func:`scopeinpython.rtmo.net_radiation_lite` and
   :func:`scopeinpython.rtmt_sb.rtmt_sb`.

   Only the default (non-MD12) fluorescence-model branch, only the simple
   ``G = 0.35*Rn`` ground-heat-flux method, and only ``meanleaf.v2``'s
   ``'layers'`` aggregation mode are ported (matches every reference case
   built during this whole port -- see the module docstring for the full
   reasoning).

   A subtlety not present in R's own (fully vectorized) implementation:
   this Python port calls :func:`scopeinpython.biochemical.get_biochemical`
   once per canopy layer with scalar inputs (see
   :func:`~scopeinpython.ebal._biochemical_per_layer`), since that
   function's Brent Ci-solver was only proven correct for scalar leaf
   micro-environments.

Where this fits
----------------

This is a pipeline-internal component, not something you typically call
directly: :func:`scopeinpython.scope.get_scope` calls this once, after
:func:`scopeinpython.rtmo.run_rtmo`, to close the energy balance. Most
users should just call :func:`~scopeinpython.scope.get_scope` -- see
:doc:`scope` for a full end-to-end example.

.. code-block:: text

   get_scope()  ->  run_rtmo()  [optical BRDF]
                ->  ebal()      [iterates leaf/soil temperature until
                                 sensible+latent heat flux matches net
                                 radiation, calling get_biochemical() and
                                 rtmt_sb() at each candidate temperature]
                ->  rtmf()/rtmz()  [fluorescence/zeaxanthin, optional]

.. automodule:: scopeinpython.ebal
   :members:
   :undoc-members:
   :show-inheritance:
