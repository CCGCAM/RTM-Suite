# Third-party model licenses

RTM-Suite's **own** code — the `ToolsRTM`, `SCOPEinR`, `toolsrtm`, and
`scopeinpython` packages' integration layer, utilities, inversion methods
(classical ML and deep learning), sensor convolution, apps, tutorials and
workflows — is distributed under the **MIT License** (see [`LICENSE`](LICENSE)).

That MIT license covers the *framework*: the code that lets you run
simulations and inversions against these models, from R or Python, through a
common interface. It does **not** relicense the radiative transfer models
themselves. Each model below is a scientific method published by its own
authors, and RTM-Suite's implementation of it retains that model's own
license where RTM-Suite is a **port/translation of the original model's own
source code**. Where RTM-Suite's implementation is instead an **independent
re-implementation from the model's published equations** (no original source
code involved), the resulting R/Python code is RTM-Suite's own work and is
MIT-licensed like the rest of the framework — but the underlying *science*
still belongs to, and must be cited to, its original authors regardless of
which license box applies below.

If a row below says "not independently verified": RTM-Suite has not
confirmed which upstream source (if any) its specific R/Python
implementation was derived from, and therefore cannot confirm which license
applies. Treat those models as used under academic/research terms pending
verification — do not assume a permissive license, and check with the
original authors before any commercial or redistribution use.

## Model-by-model record

| Model | Domain | Implementation | Original developers | Reference | DOI | Source-code provenance | License of RTM-Suite's implementation |
|---|---|---|---|---|---|---|---|
| **PROSPECT-D** | Leaf optical properties | R / Python | Féret, Gitelson, Noble, Jacquemoud | Féret et al. (2017), Remote Sens. Environ. 193, 204-215 | [10.1016/j.rse.2017.03.004](https://doi.org/10.1016/j.rse.2017.03.004) | Original authors' own reference R package (`jbferet/prospect`, GitLab) is MIT-licensed | MIT |
| **PROSPECT-PRO** | Leaf optical properties | R / Python | Féret et al. | Féret et al. (2021), Remote Sens. Environ. 252, 112173 | [10.1016/j.rse.2020.112173](https://doi.org/10.1016/j.rse.2020.112173) | Same `jbferet/prospect` package (MIT) implements PROSPECT-PRO alongside -D | MIT |
| **LIBERTY** | Leaf optical properties | R / Python | Dawson, Curran, Plummer | Dawson et al. (1998), Remote Sens. Environ. 65(1), 50-60 | [10.1016/S0034-4257(98)00007-8](https://doi.org/10.1016/S0034-4257(98)00007-8) | No public original-author reference implementation located | Not independently verified |
| **Fluspect-B** | Leaf optics & fluorescence | R / Python | Vilfan, van der Tol, Muller, Rascher, Verhoef | Vilfan et al. (2016), Remote Sens. Environ. 186, 596-615 | [10.1016/j.rse.2016.09.017](https://doi.org/10.1016/j.rse.2016.09.017) | Ships as a bundled submodel inside `Christiaanvandertol/SCOPE` (GPL-3.0); the standalone `Christiaanvandertol/Fluspect` repo carries no separate LICENSE file | GPL-3.0 (port of the model as bundled in the GPL-3.0 SCOPE codebase) |
| **Fluspect-B-Cx** | Leaf optics & fluorescence | R / Python | Vilfan, van der Tol, Yang, Wyber, Malenovský, Robinson, Verhoef | Vilfan et al. (2018), Remote Sens. Environ. 211, 345-356 | [10.1016/j.rse.2018.04.012](https://doi.org/10.1016/j.rse.2018.04.012) | Same as Fluspect-B — bundled inside SCOPE's GPL-3.0 codebase | Not independently verified (treat as GPL-3.0-adjacent pending confirmation) |
| **fourSAIL** | Canopy radiative transfer | R / Python | Verhoef | Verhoef (1984), Remote Sens. Environ. 16(2), 125-141; Verhoef (1998), PhD thesis, Wageningen University | [10.1016/0034-4257(84)90057-9](https://doi.org/10.1016/0034-4257(84)90057-9) (1998 thesis has no DOI) | No confirmed source for RTM-Suite's specific port. Note: the most prominent modern reference implementation, `jbferet/prosail` (R, bundles 4SAIL), is licensed **GPL-3 + file LICENSE**, not MIT | Not independently verified |
| **fourSAIL2** | Two-layer canopy radiative transfer | R / Python | Verhoef & Bach | Verhoef & Bach (2007), Remote Sens. Environ. 109, 166-182 | [10.1016/j.rse.2006.12.013](https://doi.org/10.1016/j.rse.2006.12.013) | Ships as a bundled submodel inside `Christiaanvandertol/SCOPE` (GPL-3.0) | **GPL-3.0** (port of the model as bundled in the GPL-3.0 SCOPE codebase) |
| **INFORM** | Forest canopy radiative transfer | R / Python | Atzberger | Atzberger (2000), 20th EARSeL Symposium, Dresden, 39-44 | No DOI (conference proceedings) | No public original-author reference implementation located | Not independently verified |
| **MARMIT** | Soil reflectance & moisture | R / Python | Bablet, Dupiau, Jacquemoud, Briottet et al. | Bablet et al. (2018), Remote Sens. Environ. 217, 1-17; Dupiau et al. (2022), MARMIT-2, Remote Sens. Environ. 272, 112951 | [10.1016/j.rse.2018.07.031](https://doi.org/10.1016/j.rse.2018.07.031); [10.1016/j.rse.2022.112951](https://doi.org/10.1016/j.rse.2022.112951) | Original authors' own reference implementation (Python): [`marmit/marmit`](https://pss-gitlab.math.univ-paris-diderot.fr/marmit/marmit), maintained by Alice Dupiau (dupiau@ipgp.fr), Stéphane Jacquemoud (jacquemoud@ipgp.fr) and Xavier Briottet (xavier.briottet@onera.fr) | Used with the original authors' permission (not a formal OSI license) |
| **SPART** | Soil-plant-atmosphere RT | R | Yang, van der Tol, Yin, Verhoef | Yang et al. (2020), Remote Sens. Environ. 247, 111870 | [10.1016/j.rse.2020.111870](https://doi.org/10.1016/j.rse.2020.111870) | Original authors' own repository, `peiqiyang/SPART` (GitHub), is GPL-3.0-licensed | **GPL-3.0** (port of the original GPL-3.0 model) |
| **SCOPE** | RT, photosynthesis, SIF & energy balance | SCOPEinR / SCOPEinPython | Van der Tol, Verhoef, Yang, Prikaziuk et al. | Van der Tol et al. (2009), Biogeosciences 6(12), 3109-29; Yang et al. (2021), Geosci. Model Dev. 14, 4697-4712 (SCOPE 2.0, final published version) | [10.5194/bg-6-3109-2009](https://doi.org/10.5194/bg-6-3109-2009); [10.5194/gmd-14-4697-2021](https://doi.org/10.5194/gmd-14-4697-2021) | Original authors' own repository, `Christiaanvandertol/SCOPE` (GitHub), is GPL-3.0-licensed | **GPL-3.0** (port of the original GPL-3.0 model) |

**MARMIT note:** only the `Bablet_2016` database ships inside `ToolsRTM`/`toolsrtm` (keeps install size small). All 8 official MARMIT databases (Bablet 2016, Dupiau 2020, Humper 2015, Lesaignoux 2008, Liu 2002, Lobell 2002, Marcq 2012, Philpot 2014) are available directly from this monorepo's own [`databases/`](databases/) folder (repo root, ~200MB) -- pass `db_root` to `get.marmit.rsoil()`/`get_marmit_rsoil()` to use any of them, no download needed.

The scientific models listed above remain the work of their respective
authors. Their inclusion in RTM-Suite does not imply authorship of the
underlying models by the RTM-Suite developers, and RTM-Suite's own MIT
license does not extend to the GPL-3.0 rows above — code that imports or
links against those specific implementations (`Fluspect-B`, `fourSAIL2`,
`SPART`, `SCOPE`/`SCOPEinR`/`SCOPEinPython`) is itself subject to GPL-3.0's
copyleft terms if redistributed.

## Questions

If you plan to redistribute or commercially use any specific model's
implementation from RTM-Suite, verify that model's license terms directly
with its original authors first — the "not independently verified" rows
above are exactly that: unverified, not confirmed-permissive.
