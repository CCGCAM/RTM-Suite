scopeinpython.biochemical
=========================

Leaf-level Farquhar/Collatz photosynthesis + van der Tol et al. (2014)
fluorescence yield, given an assumed leaf micro-environment. Direct port
of ``SCOPEinR/R/biochemical.R`` (``get.biochemical``) and its helpers in
``Biochemical_functions.R``. This is the piece SCOPE's (unported) energy-
balance iteration calls repeatedly to get ``eta`` at each candidate leaf
temperature -- it does not itself solve for temperature.

.. warning::
   ``Type='C4'`` with ``temp_correction=False`` reproduces a **real crash
   in the R source**: ``Vcmax``/``Rd`` are never assigned in that branch
   combination of ``biochemical.R`` (only the ``tempcor==1`` C4 branch and
   a separate C3-only block set them). Not worked around here.

.. automodule:: scopeinpython.biochemical
   :members:
   :undoc-members:
   :show-inheritance:
