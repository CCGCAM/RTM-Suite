scopeinpython.rtmo
==================

Optical top-of-canopy BRDF: volume-scattering geometry, hot-spot gap
probability, multi-layer 4-stream reflectance propagation, vertical flux
profile, and the top-level ``run_rtmo`` orchestrator. Direct port of
``SCOPEinR/R/RTMo_functions.R`` and ``RTMo.R`` (``getRTMo``) — stops at
the **optical** BRDF; does not compute leaf/soil temperatures, sensible/
latent heat fluxes, or photosynthesis (see :doc:`../not_ported`).

Where this fits
----------------

``run_rtmo`` is the one function in this module callable on its own (see
the soil+canopy-optics example on the :doc:`../examples` page), but it
needs solar/sky irradiance spectra (``Esun_``/``Esky_``) as caller-supplied
inputs -- :func:`scopeinpython.scope.get_scope` is what actually loads
those (from SCOPE's bundled irradiance file) and wires everything
together. Most users should just call
:func:`~scopeinpython.scope.get_scope` -- see :doc:`scope` for a full
end-to-end example.

.. code-block:: text

   get_scope()  ->  loads Esun_, Esky_ (irradiance)
                ->  run_rtmo(spectral, leaf_refl/tran, rsoil, canopy,
                              tts, tto, psi, Esun_, Esky_)
                              -> TOC BRDF (rsot, rdot, ...), flux profile

.. automodule:: scopeinpython.rtmo
   :members:
   :undoc-members:
   :show-inheritance:
