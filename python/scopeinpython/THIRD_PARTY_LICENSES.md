# Third-party model license

`scopeinpython` is a Python port of **SCOPE** (Soil Canopy Observation,
Photochemistry and Energy fluxes), a model published and maintained by
Christiaan van der Tol, Wouter Verhoef, Peiqi Yang, Egor Prikaziuk and
collaborators. SCOPE's own reference implementation,
[`Christiaanvandertol/SCOPE`](https://github.com/Christiaanvandertol/SCOPE)
(GitHub, MATLAB/R), is **GPL-3.0-licensed**.

Since this package's entire purpose is porting SCOPE, `scopeinpython` is
distributed under **GPL-3.0-only** (see [`LICENSE`](LICENSE)) -- this
matches the R sibling package `SCOPEinR`'s own `License: GPL-3` field.
`scopeinpython` depends on [`toolsrtm`](https://github.com/CCGCAM/ToolsRTMinPython)
for leaf-level optics (PROSPECT-D/-PRO, Fluspect-B/-Cx) -- see that
package's own [`THIRD_PARTY_LICENSES.md`](https://github.com/CCGCAM/ToolsRTMinPython/blob/main/THIRD_PARTY_LICENSES.md)
for the licensing of those specific leaf models.

| Model | Original developers | Reference | DOI | Source-code provenance | License |
|---|---|---|---|---|---|
| **SCOPE** | Van der Tol, Verhoef, Yang, Prikaziuk et al. | Van der Tol et al. (2009), Biogeosciences 6(12), 3109-29; Yang et al. (2021), Geosci. Model Dev. 14, 4697-4712 (SCOPE 2.0, final published version) | [10.5194/bg-6-3109-2009](https://doi.org/10.5194/bg-6-3109-2009); [10.5194/gmd-14-4697-2021](https://doi.org/10.5194/gmd-14-4697-2021) | Original authors' own repository, `Christiaanvandertol/SCOPE` (GitHub), is GPL-3.0-licensed | **GPL-3.0** |

The scientific model above remains the work of its original authors.
`scopeinpython`'s inclusion of it does not imply authorship of SCOPE by the
`scopeinpython` developers -- always cite the original publication(s) above
when using this package in scientific work, in addition to citing
RTM-Suite/SCOPEinR where appropriate.

Full write-up of what's ported, and its numerical verification against R:
[`../README.md`](../README.md). Root-repo license overview:
[`../../THIRD_PARTY_LICENSES.md`](https://github.com/CCGCAM/RTM-Suite/blob/main/THIRD_PARTY_LICENSES.md).

## Questions

If you plan to redistribute or commercially use this package, note it is
GPL-3.0 as a whole (see above) -- verify SCOPE's license terms directly
with its original authors before any use beyond GPL-3.0's own terms.
