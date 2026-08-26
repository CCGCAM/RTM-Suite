toolsrtm.fluspect
===================

FLUSPECT-B / FLUSPECT-B-Cx leaf model: PROSPECT reflectance/transmittance
plus chlorophyll-fluorescence excitation-emission matrices. Direct port
of ``ToolsRTM/R/fluspect_B.R`` (``getFluspect.B``) and ``fluspect_Cx.R``
(``getFluspect.Cx``).

Quick example
-------------

.. code-block:: python

   from toolsrtm import fluspect_b

   flu = fluspect_b(N=1.5, Cab=40, Car=8, Anth=1, EWT=0.01, LMA=0.009,
                     Cs=0, Cx=0, fqe=0.01)
   print(flu.refl[800 - 400])       # reflectance at 800nm, same physics as PROSPECT-D
   print(flu.MbI.shape)             # backward fluorescence excitation-emission matrix (PSI)

.. code-block:: text

   Input                              fluspect_b()              Output
   ---------------------------        ----------------------    ---------------------------
   N, Cab, Car, Anth, EWT, LMA                                  flu.refl / flu.tran  [2101]
   Cs (senescent pigment)             -------------------->      flu.MbI / flu.MbII
   Cx (xanthophyll de-epoxidation)                                (fluorescence excitation-
   fqe (fluorescence quantum eff.)                                 emission matrices, PSI/PSII)

.. automodule:: toolsrtm.fluspect
   :members:
   :undoc-members:
   :show-inheritance:
