scopeinpython
=============

Python port of the SCOPE model: soil (BSM), leaf optics (Fluspect-Cx,
including the multi-layer mSCOPE wrapper), the optical top-of-canopy BRDF
pipeline (RTMo), leaf biochemistry/fluorescence yield, the thermal
energy-balance closure loop (``ebal``), the fluorescence (RTMf) and
zeaxanthin (RTMz) canopy radiative-transfer models, and an end-to-end
``get_scope`` wrapper tying all of the above together from a single LUT
input row. Depends on ``toolsrtm`` for leaf optics, exactly as the R
``SCOPEinR`` package depends on ``ToolsRTM``.

.. toctree::
   :maxdepth: 1

   soil
   spectral
   rtmo
   biochemical
   fluspect
   fluspect_mscope
   rtmf
   rtmz
   thermal
   rtmt_sb
   ebal
   scope
   utils
