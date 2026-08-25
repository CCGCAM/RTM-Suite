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

.. automodule:: scopeinpython.fluspect_mscope
   :members:
   :undoc-members:
   :show-inheritance:
