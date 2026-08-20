Installation
============

Both packages require only ``numpy`` and ``scipy`` at runtime.

.. code-block:: bash

   pip install -e python/toolsrtm
   pip install -e python/scopeinpython    # depends on toolsrtm

For the test suites (numeric regression tests against R reference data):

.. code-block:: bash

   pip install -e "python/toolsrtm[test]"
   pip install -e "python/scopeinpython[test]"
   cd python/toolsrtm && python -m pytest tests -q
   cd python/scopeinpython && python -m pytest tests -q

Building these docs
--------------------

.. code-block:: bash

   pip install sphinx sphinx-rtd-theme myst-parser
   cd python/docs
   sphinx-build -b html . _build/html

The built HTML in this repo (``docs/python/``) is committed, the same way
``docs/toolsrtm/`` and ``docs/scopeinr/`` (the R pkgdown sites) are —
regenerate it with the command above and copy ``_build/html/*`` into
``docs/python/`` after any docstring or API change.
