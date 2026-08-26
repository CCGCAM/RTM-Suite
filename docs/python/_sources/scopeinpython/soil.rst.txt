scopeinpython.soil
==================

Brightness-Shape-Moisture (BSM) soil reflectance model. Direct port of
``SCOPEinR/R/BSM.R``.

Quick example
-------------

.. code-block:: python

   from scopeinpython import SoilParams, WettingParams, get_bsm

   rsoil = get_bsm(SoilParams(BSMBrightness=0.5, BSMlat=25, BSMlon=45),
                    WettingParams(SMp=15, SMC=25, film=0.015))
   print(rsoil[800 - 400])   # soil reflectance at 800nm

.. code-block:: text

   Input                              get_bsm()                 Output
   ---------------------------        ----------------------    ---------------------------
   SoilParams(BSMBrightness,                                     rsoil  [2001]  400-2400nm
              BSMlat, BSMlon)         -------------------->       soil reflectance spectrum,
   WettingParams(SMp, SMC, film)                                  ready as rsoil for run_rtmo()

.. automodule:: scopeinpython.soil
   :members:
   :undoc-members:
   :show-inheritance:
