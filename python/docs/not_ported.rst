Known limitations
===================

Everything under the API reference sections is ported and tested,
including all 9 SMAC sensors (:mod:`toolsrtm.smac`) and Sentinel-2/Landsat/MODIS
satellite retrieval (:mod:`toolsrtm.satellite`). The remaining, genuine gaps
are all in R's own less-used code paths, not missing pieces of the core
workflow:

- **SCOPE thermal/fluorescence chain**: the per-wavelength thermal RTM
  alternative (``RTMt_planck.R``), the full ``(13,36,nl)`` per-leaf-angle
  branches (only "SCOPE-lite" scalar-per-layer is ported), the ``obsdir``
  (directional brightness temperature) branch, the two time-series-history
  soil-heat-flux methods, and the alternative MD12 fluorescence-model
  branch.
- **PAR / net-radiation breakdown** inside RTMo: only the 6 quantities
  :func:`scopeinpython.rtmo.net_radiation_lite` needs are ported.
- **mSCOPE per-layer wiring**: the multi-layer leaf-optics wrapper
  (:func:`scopeinpython.fluspect_mscope.fluspect_mscope`) is ported, but
  nothing yet feeds its per-layer output into ``run_rtmo``'s canopy
  radiative transfer, so :func:`scopeinpython.scope.get_scope` only
  supports a single mSCOPE profile layer.
- **MODTRAN atmospheric-file irradiance mode** for ``get.calcTOCirr``.
- ``get.SCOPE.parallel``, directional BRDF (``options.calc_directional``),
  time-series mode, and angle-file LIDF.
- General statistics/ML helpers and CARS-PLS's plotting helper
  (``get.plot.cars.pls``).
