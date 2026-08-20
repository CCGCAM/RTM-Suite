toolsrtm.fluspect
===================

FLUSPECT-B / FLUSPECT-B-Cx leaf model: PROSPECT reflectance/transmittance
plus chlorophyll-fluorescence excitation-emission matrices. Direct port
of ``ToolsRTM/R/fluspect_B.R`` (``getFluspect.B``) and ``fluspect_Cx.R``
(``getFluspect.Cx``).

.. note::
   A real indexing bug was found in **``ToolsRTM::getFluspect.Cx``**
   during this port: its fluorescence-matrix code used
   ``intersect(wlp, wlf)`` where it should use ``which()``/``match()``,
   silently reading Kubelka-Munk coefficients at the wrong wavelengths
   (1039-1249 nm instead of 640-850 nm) when building ``Mb``/``Mf``.
   **Since fixed in the R source**; ``fluspect_cx`` here matches the fixed
   behavior. See the function's docstring and :doc:`../verification` for
   the full history. ``fluspect_b`` never had this bug.

.. automodule:: toolsrtm.fluspect
   :members:
   :undoc-members:
   :show-inheritance:
