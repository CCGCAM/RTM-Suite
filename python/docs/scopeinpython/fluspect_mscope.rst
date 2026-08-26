scopeinpython.fluspect_mscope
==============================

Multi-layer (mSCOPE) leaf-optics wrapper around
:func:`scopeinpython.fluspect.get_fluspect_cx_scope`. Direct port of
``SCOPEinR/R/fluspect_mSCOPE.R`` (``get.fluspect_mSCOPE``).

Computes leaf optical properties (reflectance, transmittance,
fluorescence excitation-emission matrices ``Mb``/``Mf``, pigment
contribution factors) once per distinct leaf-biochemistry profile layer,
then replicates across the canopy sublayers each profile layer spans,
weighted by ``pLAI``.

.. note::
   Two things worth knowing, ported/documented rather than "fixed":

   - **R quirk reproduced exactly**: at each profile-layer boundary, the
     canopy sublayer shared between two consecutive profile layers is
     assigned twice (once per layer's replication loop) -- the later
     assignment wins. This mirrors R's own overlapping-range
     ``rho_temp[in1:in2,] <- ...`` assignment.
   - **Real R crash, not offered as a Python default**: calling
     ``get.fluspect_mSCOPE()`` in R *without* the ``step`` argument always
     crashes ("number of items to replace is not a multiple of
     replacement length"), confirmed via a standalone repro -- its
     ``Mb``/``Mf`` output array is pre-allocated assuming a fixed 1 nm
     grid, but the function actually calls
     ``getFluspect.Cx.SCOPE(..., step=5)`` internally when ``step`` is
     missing, producing an irreconcilable shape mismatch. Since this R
     code path can never succeed, ``step`` is a **required** parameter
     here rather than optional-with-a-crashing-default.

Quick example
-------------

.. code-block:: python

   import numpy as np
   from scopeinpython import get_spectra_scope
   from scopeinpython.fluspect_mscope import MultiLayerLeafBio, fluspect_mscope

   spectral = get_spectra_scope()
   mly = MultiLayerLeafBio(nly=2, pLAI=np.array([0.5, 0.5]),
                            pCab=np.array([40.0, 30.0]), pEWT=np.array([0.01, 0.01]),
                            pCar=np.array([8.0, 8.0]), pLMA=np.array([0.009, 0.009]),
                            pCs=np.array([0.0, 0.0]), pN=np.array([1.5, 1.5]))
   res = fluspect_mscope(mly, spectral, nl=10, Cx=0.0, fqe=0.01,
                          Prot=0.0, CBC=0.0, Anth=1.0, step=5.0)
   print(res.refl.shape)   # (nl, nwlP) -- one reflectance spectrum per canopy sublayer

.. code-block:: text

   Input                              fluspect_mscope()         Output
   ---------------------------        ----------------------    ---------------------------
   mly : MultiLayerLeafBio                                      res.refl / res.tran  (nl, nwlP)
     (nly biochem. profile layers,     -------------------->     res.Mb / res.Mf  (nwlf,nwle,nl)
      pLAI-weighted)                                              -- replicated across all nl
   spectral, nl (canopy sublayers)                                 canopy sublayers, feeds
   Cx, fqe, Prot, CBC, Anth, step                                  run_rtmo()/rtmf()

.. automodule:: scopeinpython.fluspect_mscope
   :members:
   :undoc-members:
   :show-inheritance:
