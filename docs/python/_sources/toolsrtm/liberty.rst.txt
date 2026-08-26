toolsrtm.liberty
==================

LIBERTY conifer-needle leaf model (Dawson, Curran & Plummer 1998).
Direct port of ``ToolsRTM/R/liberty.R``.

Quick example
-------------

.. code-block:: python

   from toolsrtm import liberty

   needle = liberty(cell_d=40, inter_c=0.045, baseline_abs=0.0006, leaf_thick=1.6,
                     albino_abs=0, Cab=40, EWT=0.01, lign_cell=2, Nitrogen=1)
   print(needle.refl[800 - 400])   # reflectance at 800nm

.. code-block:: text

   Input                              liberty()                 Output
   ---------------------------        ----------------------    ---------------------------
   cell_d (avg cell diameter, um)                                needle.refl  [2101]  400-2500nm
   inter_c, baseline_abs,             -------------------->      needle.tran  [2101]  400-2500nm
   leaf_thick, albino_abs                                        (flatter NIR plateau than
   Cab, EWT, lign_cell, Nitrogen                                  broadleaf PROSPECT)

.. automodule:: toolsrtm.liberty
   :members:
   :undoc-members:
   :show-inheritance:
