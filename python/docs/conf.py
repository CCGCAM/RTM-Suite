"""Sphinx configuration for the 0-RTM-Suite Python port (toolsrtm + scopeinpython).

Style/tooling modeled on https://scope-model.readthedocs.io/en/master/
(Sphinx + Read the Docs theme + autodoc from docstrings), the same way the
R side's docs/toolsrtm and docs/scopeinr are pkgdown sites auto-generated
from roxygen docstrings.
"""
import os
import sys

sys.path.insert(0, os.path.abspath("../toolsrtm/src"))
sys.path.insert(0, os.path.abspath("../scopeinpython/src"))

project = "0-RTM-Suite Python port"
author = "Carlos Camino"
copyright = "Wageningen University & Research"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",  # parses the numpy-style docstrings used throughout toolsrtm/scopeinpython
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosectionlabel",  # lets :ref:`page:Section Title` link across pages (workflows.rst -> examples.rst, etc.)
    "myst_parser",
]

autosectionlabel_prefix_document = True

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_use_param = True
napoleon_use_rtype = True

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
autodoc_member_order = "bysource"
autodoc_typehints = "description"

intersphinx_mapping = {
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

# toolsrtm.srf's SrfTable.wl and SrfConvolutionResult.wl are two unrelated
# dataclasses that both, correctly, name their wavelength-array field "wl" --
# autodoc's typehints-in-description mode then tries to cross-reference the
# bare attribute name "wl" and finds both, which Sphinx warns about (not
# nitpick-only, so nitpick_ignore has no effect without nitpicky=True).
# Suppressing this one class of warning rather than renaming either field,
# since "wl" is the correct/consistent name for what both hold.
suppress_warnings = ["ref.python"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_theme_options = {
    # collapse_navigation=True: only the section containing the current page
    # expands in the sidebar (e.g. reading a LEARN chapter no longer shows
    # all 14 toolsrtm/13 scopeinpython function pages pre-expanded under
    # REFERENCE -- "toolsrtm"/"scopeinpython" show as single collapsed
    # links, matching the LEARN chapters' own general-page-first pattern).
    "collapse_navigation": True,
    "navigation_depth": 3,
}

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
