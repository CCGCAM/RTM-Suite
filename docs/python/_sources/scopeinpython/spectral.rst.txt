scopeinpython.spectral
======================

SCOPE spectral-region/wavelength-grid definitions. Direct port of
``SCOPEinR/R/define_bands.R`` (``get.spectra.SCOPE``).

Quick example
-------------

.. code-block:: python

   from scopeinpython import get_spectra_scope

   spectral = get_spectra_scope()
   print(spectral.wlO.shape)   # optical wavelength grid, 400-2400nm
   print(spectral.wlT.shape)   # thermal wavelength grid, 2500-50000nm

.. code-block:: text

   Input                              get_spectra_scope()       Output
   ---------------------------        ----------------------    ---------------------------
   (no arguments -- fixed SCOPE                                  spectral : SpectralConfig
    wavelength-grid definitions)      -------------------->       wlS, wlO, wlT, wlP, wlPAR,
                                                                   wlE, wlF, IwlP, IwlT -- fed
                                                                   into run_rtmo()/get_scope()

.. automodule:: scopeinpython.spectral
   :members:
   :undoc-members:
   :show-inheritance:
