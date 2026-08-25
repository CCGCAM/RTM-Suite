scopeinpython.rtmz
====================

SCOPE canopy zeaxanthin RTM: the small modification of TOC outgoing
radiance due to violaxanthin-to-zeaxanthin conversion in leaves
(photoprotection, 500-600nm). Direct port of ``SCOPEinR/R/RTMz.R``
(``get.RTMz``). Structurally close to :func:`scopeinpython.rtmf.rtmf`
(same geometric-factor setup, same layer-recursion pattern) but works
directly in the wavelength domain (no excitation-emission matrix) and
returns *corrections* to be added onto an existing
:class:`~scopeinpython.rtmo.RTMoResult`, rather than a self-contained
new spectrum.

.. automodule:: scopeinpython.rtmz
   :members:
   :undoc-members:
   :show-inheritance:
