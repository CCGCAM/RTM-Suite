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

Where this fits
----------------

This is a pipeline-internal component, not something you typically call
directly: :func:`scopeinpython.scope.get_scope` calls this after
:func:`scopeinpython.rtmo.run_rtmo`, when the zeaxanthin-correction
output is requested. Most users should just call
:func:`~scopeinpython.scope.get_scope` -- see :doc:`scope` for a full
end-to-end example.

.. code-block:: text

   run_rtmo()  ->  RTMoResult  (TOC BRDF)
                        \
                         -->  rtmz(RTMoResult, ...)  ->  correction to add onto
                                                          the existing TOC spectrum

.. automodule:: scopeinpython.rtmz
   :members:
   :undoc-members:
   :show-inheritance:
