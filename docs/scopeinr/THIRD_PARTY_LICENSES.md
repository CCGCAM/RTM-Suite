# Third-party model license

`SCOPEinR` is an R port of **SCOPE** (Soil Canopy Observation,
Photochemistry and Energy fluxes), a model published and maintained by
Christiaan van der Tol, Wouter Verhoef, Peiqi Yang, Egor Prikaziuk and
collaborators. SCOPE’s own reference implementation,
[`Christiaanvandertol/SCOPE`](https://github.com/Christiaanvandertol/SCOPE)
(GitHub, MATLAB/R), is **GPL-3.0-licensed**.

Since this package’s entire purpose is porting SCOPE, `SCOPEinR` is
distributed under **GPL-3.0** (`License: GPL-3` in `DESCRIPTION`) – this
matches the Python port `scopeinpython`’s own license. `SCOPEinR`
depends on [`ToolsRTM`](https://gitlab.com/caminoccg/toolsrtm) for
leaf-level optics (PROSPECT-D/-PRO) – see that package’s own
[`THIRD_PARTY_LICENSES.md`](https://gitlab.com/caminoccg/toolsrtm/-/blob/main/THIRD_PARTY_LICENSES.md)
for the licensing of those specific leaf models.

| Model | Original developers | Reference | DOI | Source-code provenance | License |
|----|----|----|----|----|----|
| **SCOPE** | Van der Tol, Verhoef, Yang, Prikaziuk et al. | Van der Tol et al. (2009), Biogeosciences 6(12), 3109-29; Yang et al. (2021), Geosci. Model Dev. 14, 4697-4712 (SCOPE 2.0, final published version) | [10.5194/bg-6-3109-2009](https://doi.org/10.5194/bg-6-3109-2009); [10.5194/gmd-14-4697-2021](https://doi.org/10.5194/gmd-14-4697-2021) | Original authors’ own repository, `Christiaanvandertol/SCOPE` (GitHub), is GPL-3.0-licensed | **GPL-3.0** |

The scientific model above remains the work of its original authors.
`SCOPEinR`’s inclusion of it does not imply authorship of SCOPE by the
`SCOPEinR` developers – always cite the original publication(s) above
when using this package in scientific work, in addition to citing
RTM-Suite/SCOPEinR.

Full write-up of the package: [`README.md`](README.md). Root-repo
license overview:
[`../THIRD_PARTY_LICENSES.md`](https://github.com/CCGCAM/RTM-Suite/blob/main/THIRD_PARTY_LICENSES.md).

## Questions

If you plan to redistribute or commercially use this package, note it is
GPL-3.0 as a whole (see above) – verify SCOPE’s license terms directly
with its original authors before any use beyond GPL-3.0’s own terms.
