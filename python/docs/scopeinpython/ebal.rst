scopeinpython.ebal
====================

The SCOPE energy-balance closure loop: iterates sunlit/shaded leaf and
soil temperature until sensible+latent heat flux matches net radiation,
coupling radiative transfer, aerodynamics, photosynthesis/fluorescence
and heat-flux partitioning. Direct, **"SCOPE-lite"-only** port of
``SCOPEinR::get.ebal`` -- the top of the whole thermal energy-balance
chain built this session.

.. note::
   Depends on :func:`scopeinpython.rtmo.net_radiation_lite` (a **real bug
   found and fixed in R** here: the "lite" direct-beam term used a stray
   leftover loop variable, making direct-beam-absorbed net radiation come
   out constant across all canopy layers -- see that function's
   docstring) and :func:`scopeinpython.rtmt_sb.rtmt_sb`.

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

.. automodule:: scopeinpython.ebal
   :members:
   :undoc-members:
   :show-inheritance:
