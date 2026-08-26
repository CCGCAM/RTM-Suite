scopeinpython.fluspect
======================

FLUSPECT-B-Cx, the SCOPE variant (``getFluspect.Cx.SCOPE``) that SCOPE's
own leaf-optics pipeline actually calls -- **not the same function** as
``toolsrtm.fluspect_cx`` despite very similar code. Direct port of
``SCOPEinR/R/fluspect_Cx_forSCOPE.R``.

Key differences from ``toolsrtm.fluspect_cx``:

- A caller-configurable ``step`` (nm) controls the ``Mb``/``Mf``
  excitation-emission matrix resolution. SCOPE's own default ``step=5``
  gives 53x71 matrices; ``step=1`` gives 211x351 (like the non-SCOPE
  version, but see below -- still not numerically identical).
- The SIF response is scaled by ``step`` itself, not a fixed constant
  (the non-SCOPE version always uses ``int=5`` regardless of its own
  fixed 1 nm grid).
- A single combined ``Mb``/``Mf`` pair (one ``phi`` spectrum), same as
  the non-SCOPE Cx version.
- Uses ``match()`` (correct) for locating excitation/emission wavelengths
  within the full grid -- does **not** have the ``intersect()`` indexing
  bug documented in :func:`toolsrtm.fluspect.fluspect_cx`.
- Uses ``SCOPEinR::optipar2021.Pro.CX`` as its real default optical
  parameter table (``optipar2017.ProspectD`` is missing ``Kp``/``Kcbc``
  data needed for the PROSPECT-PRO/Cx ``Kall`` formula).

Quick example
-------------

.. code-block:: python

   from scopeinpython.fluspect import get_fluspect_cx_scope

   fl = get_fluspect_cx_scope(Cab=40, Car=8, EWT=0.01, LMA=0.009, Cs=0, N=1.5,
                               fqe=0.01, Cx=0, Prot=0.0, CBC=0.0, Anth=1.0, step=5.0)
   print(fl.refl.shape, fl.Mb.shape)   # (2001,) reflectance; (53, 71) fluorescence matrix

.. code-block:: text

   Input                              get_fluspect_cx_scope()   Output
   ---------------------------        ----------------------    ---------------------------
   Cab, Car, EWT, LMA, Cs, N                                     fl.refl / fl.tran  [2001]
   fqe (fluorescence quantum eff.)    -------------------->      fl.Mb / fl.Mf  [53x71 @step=5]
   Cx, Prot, CBC, Anth                                            (single-layer input to
   step (nm, EEM matrix resolution)                                fluspect_mscope())

.. automodule:: scopeinpython.fluspect
   :members:
   :undoc-members:
   :show-inheritance:
