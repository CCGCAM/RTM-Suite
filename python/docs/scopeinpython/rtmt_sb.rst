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
   brightness temperature) branch are **not ported**.

Where this fits
----------------

This is a pipeline-internal component, not something you typically call
directly: :func:`scopeinpython.ebal.ebal` calls this at each candidate
leaf/soil temperature during its iteration, converting temperatures (via
:func:`scopeinpython.thermal.stefan_boltzmann`) into outgoing thermal
radiation and net radiation per component. Most users should just call
:func:`~scopeinpython.scope.get_scope` -- see :doc:`scope` for a full
end-to-end example.

.. code-block:: text

   ebal()  [iterating leaf/soil temperature]
       ->  rtmt_sb(leaf/soil temperatures, ...)
             ->  stefan_boltzmann()  ->  outgoing LW, net radiation per component

.. automodule:: scopeinpython.rtmt_sb
   :members:
   :undoc-members:
   :show-inheritance:
