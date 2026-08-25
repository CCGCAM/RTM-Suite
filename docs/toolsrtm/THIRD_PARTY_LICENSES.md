# Third-party model licenses

`ToolsRTM` is an R port of several radiative transfer models published
by the remote-sensing research community, integrated behind one common
interface. The models it bundles are not all under the same license, so
`ToolsRTM` as a whole is distributed under **GPL-3.0** (`License: GPL-3`
in `DESCRIPTION`): three of the bundled models
(`getFluspect.B`/`getFluspect.Cx`, `foursail2`, `SPART`) are ports of
GPL-3.0-licensed original models, and GPL-3.0 requires any combined work
incorporating GPL-3.0 code to be distributed as GPL-3.0 as a whole –
this matches the Python port `toolsrtm`’s own license.

Individually, some of `ToolsRTM`’s functions (e.g. PROSPECT-D/-PRO) are
themselves original/independent implementations that would be
MIT-licensable on their own; that does not change the GPL-3.0 status of
the combined `ToolsRTM` distribution as a whole, per GPL-3.0’s terms.

If a row below says “not independently verified”: `ToolsRTM` has not
confirmed which upstream source (if any) its specific implementation was
derived from. Treat those models as used under academic/research terms
pending verification.

| Model | R function(s) | Original developers | Reference | DOI | Source-code provenance | License |
|----|----|----|----|----|----|----|
| **PROSPECT-D** | `prospect_PRO` (leaf.model=“PROSPECT-D”) | Féret, Gitelson, Noble, Jacquemoud | Féret et al. (2017), Remote Sens. Environ. 193, 204-215 | [10.1016/j.rse.2017.03.004](https://doi.org/10.1016/j.rse.2017.03.004) | Original authors’ own reference R package (`jbferet/prospect`, GitLab) is MIT-licensed | MIT (as implemented; the combined `ToolsRTM` package is GPL-3.0, see above) |
| **PROSPECT-PRO** | `prospect_PRO` | Féret et al. | Féret et al. (2021), Remote Sens. Environ. 252, 112173 | [10.1016/j.rse.2020.112173](https://doi.org/10.1016/j.rse.2020.112173) | Same `jbferet/prospect` package (MIT) | MIT (as implemented; combined package is GPL-3.0) |
| **LIBERTY** | `liberty` | Dawson, Curran, Plummer | Dawson et al. (1998), Remote Sens. Environ. 65(1), 50-60 | [10.1016/S0034-4257(98)00007-8](https://doi.org/10.1016/S0034-4257(98)00007-8) | No public original-author reference implementation located | Not independently verified |
| **Fluspect-B** | `getFluspect.B` | Vilfan, van der Tol, Muller, Rascher, Verhoef | Vilfan et al. (2016), Remote Sens. Environ. 186, 596-615 | [10.1016/j.rse.2016.09.017](https://doi.org/10.1016/j.rse.2016.09.017) | Ships as a bundled submodel inside `Christiaanvandertol/SCOPE` (GPL-3.0) | **GPL-3.0** |
| **Fluspect-B-Cx** | `getFluspect.Cx` | Vilfan, van der Tol, Yang, Wyber, Malenovský, Robinson, Verhoef | Vilfan et al. (2018), Remote Sens. Environ. 211, 345-356 | [10.1016/j.rse.2018.04.012](https://doi.org/10.1016/j.rse.2018.04.012) | Same as Fluspect-B – bundled inside SCOPE’s GPL-3.0 codebase | Not independently verified (treat as GPL-3.0-adjacent) |
| **fourSAIL** | `foursail` | Verhoef | Verhoef (1984), Remote Sens. Environ. 16(2), 125-141; Verhoef (1998), PhD thesis, Wageningen University | [10.1016/0034-4257(84)90057-9](https://doi.org/10.1016/0034-4257(84)90057-9) | No confirmed source for this port. Note: the most prominent modern reference implementation, `jbferet/prosail` (R, bundles 4SAIL), is GPL-3 + file LICENSE | Not independently verified |
| **fourSAIL2** | `foursail2` | Verhoef & Bach | Verhoef & Bach (2007), Remote Sens. Environ. 109, 166-182 | [10.1016/j.rse.2006.12.013](https://doi.org/10.1016/j.rse.2006.12.013) | Ships as a bundled submodel inside `Christiaanvandertol/SCOPE` (GPL-3.0) | **GPL-3.0** |
| **INFORM** | `inform`, `foursail.inform`, `foursail.inf` | Atzberger | Atzberger (2000), 20th EARSeL Symposium, Dresden, 39-44 | No DOI (conference proceedings) | No public original-author reference implementation located | Not independently verified |
| **MARMIT** / **MARMIT-2** | `get.marmit1`, `get.marmit2`, `get.marmit.rsoil` | Bablet, Dupiau, Jacquemoud, Briottet et al. | Bablet et al. (2018), Remote Sens. Environ. 217, 1-17; Dupiau et al. (2022), Remote Sens. Environ. 272, 112951 | [10.1016/j.rse.2018.07.031](https://doi.org/10.1016/j.rse.2018.07.031); [10.1016/j.rse.2022.112951](https://doi.org/10.1016/j.rse.2022.112951) | Original authors’ own reference implementation (Python): [`marmit/marmit`](https://pss-gitlab.math.univ-paris-diderot.fr/marmit/marmit), maintained by Alice Dupiau (<dupiau@ipgp.fr>), Stéphane Jacquemoud (<jacquemoud@ipgp.fr>) and Xavier Briottet (<xavier.briottet@onera.fr>) | Used with the original authors’ permission (not a formal OSI license) |
| **SPART** | `SPART`, `SPART.simN` | Yang, van der Tol, Yin, Verhoef | Yang et al. (2020), Remote Sens. Environ. 247, 111870 | [10.1016/j.rse.2020.111870](https://doi.org/10.1016/j.rse.2020.111870) | Original authors’ own repository, `peiqiyang/SPART` (GitHub), is GPL-3.0-licensed | **GPL-3.0** |

**MARMIT note:** only `Bablet_2016` ships with this package (keeps
install size small). All 8 official MARMIT databases are available
directly from the monorepo’s own `databases/` folder (repo root, one
level up from `ToolsRTM/`) – pass `db_root` to
[`get.marmit.rsoil()`](reference/get.marmit.rsoil.md), see
[`?get.marmit.rsoil`](reference/get.marmit.rsoil.md).

Full write-up of the package: [`README.md`](README.md). Root-repo
license overview:
[`../THIRD_PARTY_LICENSES.md`](https://github.com/CCGCAM/RTM-Suite/blob/main/THIRD_PARTY_LICENSES.md).

## Questions

If you plan to redistribute or commercially use this package, note it is
GPL-3.0 as a whole (see above) and, for any specific model row marked
“not independently verified”, verify that model’s license terms directly
with its original authors first.
