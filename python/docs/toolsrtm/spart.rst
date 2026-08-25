toolsrtm.spart
===============

SPART top-of-canopy (``spart_toc``) and top-of-atmosphere (``spart_toa``)
reflectance/radiance. Direct port of ``ToolsRTM/R/spart.R``, verified
against a real, unmodified ``ToolsRTM::SPART()`` call.

.. note::
   ``spart_toa``'s atmospheric correction (SMAC) is implemented in
   :doc:`smac` and is fully sensor-agnostic, but only **Sentinel-2A**'s
   sensor coefficients are bundled as package data so far, out of the 9
   the R package ships — see :doc:`smac`'s module docstring for exactly
   what adding another sensor takes (a data-export exercise, not a code
   change), and :doc:`../not_ported`.

.. automodule:: toolsrtm.spart
   :members:
   :undoc-members:
   :show-inheritance:
