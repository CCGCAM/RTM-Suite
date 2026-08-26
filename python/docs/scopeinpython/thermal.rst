scopeinpython.thermal
=======================

Scalar-per-timestep thermal/aerodynamic building blocks used by the
SCOPE energy-balance loop (``ebal.R``, not yet ported). Direct ports of
``SCOPEinR/R/Monin_ObuKhov.R``, ``SCOPEinR/R/resistances.R``, and
``SCOPEinR/R/heatfluxes.R`` -- the first pieces of the thermal
energy-balance chain to be ported.

.. note::
   ``get_resistances``'s ``rac``/``rws`` outputs reproduce a real R
   quirk: they use the *stability-uncorrected* eddy diffusivity (the
   local ``Kh`` variable in ``SCOPEinR::get.resistances``), not the
   stability-corrected value returned as the function's own ``Kh``
   output field. Confirmed by direct reading of the R source (the local
   ``Kh`` is never reassigned after the corrected value is written into
   ``resist_out[['Kh']]``) -- reproduced exactly rather than "fixed",
   since there's no independent way to tell whether this was intentional.

Quick example
-------------

.. code-block:: python

   import numpy as np
   from scopeinpython.thermal import stefan_boltzmann

   M = stefan_boltzmann(np.array([15.0, 25.0]))   # leaf/soil temperature, degC
   print(M)   # blackbody radiant exitance, W/m2

.. code-block:: text

   Input                              stefan_boltzmann()        Output
   ---------------------------        ----------------------    ---------------------------
   T_C  temperature(s), degC          -------------------->      M  blackbody radiant
                                                                    exitance, W/m2 (Stefan-
                                                                    Boltzmann law)

.. automodule:: scopeinpython.thermal
   :members:
   :undoc-members:
   :show-inheritance:
