toolsrtm.leaf
=============

Leaf optical property models: ``calctav``, PROSPECT-D, PROSPECT-PRO.
Direct port of ``ToolsRTM/R/calctav.R``, ``prospect_DB.R``, ``prospect_PRO.R``.

Quick example
-------------

.. code-block:: python

   from toolsrtm import prospect_d

   leaf = prospect_d(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0,
                      EWT=0.01, LMA=0.009, alpha=40)
   print(leaf.refl[550 - 400], leaf.tran[550 - 400])   # reflectance/transmittance at 550nm

.. code-block:: text

   Input                            prospect_d()              Output
   ---------------------------      ----------------------    ---------------------------
   N        = 1.5  (leaf structure)                           leaf.refl  [2101]  400-2500nm
   Cab      = 40   (ug/cm2)         -------------------->      leaf.tran  [2101]  400-2500nm
   Car      = 8    (ug/cm2)
   EWT      = 0.01 (cm, water)
   LMA      = 0.009 (g/cm2, dry matter)
   alpha    = 40   (deg, incidence)

.. automodule:: toolsrtm.leaf
   :members:
   :undoc-members:
   :show-inheritance:
