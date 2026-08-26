scopeinpython.utils
===================

Small numerical/radiometric helpers used by RTMo: ``sint`` (Simpson-like
integration), photon-energy conversions, Planck's law.

Quick example
-------------

.. code-block:: python

   import numpy as np
   from scopeinpython.utils import get_planck, satvap

   L = get_planck(np.array([500.0, 800.0]), np.array([5778.0]))   # solar radiance, W/m2/sr/m
   print(L)
   es = satvap(np.array([20.0]))   # saturated vapour pressure at 20degC
   print(es)

.. code-block:: text

   Input                              get_planck()              Output
   ---------------------------        ----------------------    ---------------------------
   wl_nm  wavelength(s), nm                                      L  spectral radiance,
   Tb     blackbody temperature, K   -------------------->        W m-2 sr-1 m-1 (Planck's law)
   em     emissivity (optional)

.. automodule:: scopeinpython.utils
   :members:
   :undoc-members:
   :show-inheritance:
