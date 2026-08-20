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

.. note::
   Ported against a **fixed** R source. Five real issues were found and
   (mostly) fixed in ``SCOPEinR::get.RTMz`` during this port:

   - The same column-recycling bug family as ``get.RTMf``.
   - Its real output used to live entirely inside an unreachable
     ``if (get.plots==TRUE)`` block *after* the loop that produced its
     inputs — normal calls silently returned the input unmodified.
   - A ``Po[1:nl+1]`` indexing bug (should be the scalar ``Po[nl+1]``).
   - A hardcoded ``dim = c(30, 13, 36)`` instead of ``c(nl, 13, 36)``.
   - A MATLAB-to-R mistranslation, ``sum(LoF_, 2)`` (R's ``sum()`` has no
     dimension argument — fixed to ``rowSums(LoF_)``).

   See :doc:`../verification` for the full history, plus two things
   flagged but **not** independently confirmed/fixed (documented in the
   module docstring): a likely 7th bug (``vfEplu_u`` using ``fsctl_nl``
   where the parallel ``vfEplu_h``/``get.RTMf`` code uses ``foctl_nl``),
   and R's *other* ``data.Knu``-is-already-an-array code path, which
   never applies the ``Kn2Cx()`` NPQ conversion at all (discovered while
   building this port's own R reference case).

.. automodule:: scopeinpython.rtmz
   :members:
   :undoc-members:
   :show-inheritance:
