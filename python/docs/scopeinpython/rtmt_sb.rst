scopeinpython.rtmt_sb
=======================

Total thermal-infrared outgoing radiation and net radiation per
leaf/soil component (Stefan-Boltzmann, spectrally-integrated), given
already-solved leaf/soil temperatures. Direct, **partial** port of
``SCOPEinR/R/RTMt.sb.R`` (``get.RTMt.sb``).

.. note::
   Only the "SCOPE-lite" scalar-per-layer branch is ported (matches
   every reference case built during this port). The full per-leaf-angle
   ``(13, 36, nl)`` array branch and the ``obsdir`` (observation-direction
   brightness temperature) branch are **not ported** -- the latter has an
   unresolved, likely-buggy R indexing expression
   (``data.rad$vb[1, nl]``, using the *layer count* as a *wavelength*
   index) flagged in the module docstring but not independently confirmed.

.. automodule:: scopeinpython.rtmt_sb
   :members:
   :undoc-members:
   :show-inheritance:
