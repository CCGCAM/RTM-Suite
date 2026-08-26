### ToolsRTM package

[![R \>=
4.3](https://img.shields.io/badge/R-%3E%3D%204.3-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![ToolsRTM on
GitLab](https://img.shields.io/badge/GitLab-ToolsRTM-FC6D26?logo=gitlab&logoColor=white)](https://gitlab.com/caminoccg/toolsrtm)
[![RTM–Suite on
GitHub](https://img.shields.io/badge/GitHub-RTM--Suite-181717?logo=github&logoColor=white)](https://github.com/CCGCAM/RTM-Suite)
[![toolsrtm (Python port) on
GitHub](https://img.shields.io/badge/GitHub-toolsrtm%20(Python)-181717?logo=github&logoColor=white)](https://github.com/CCGCAM/RTM-Suite/tree/main/python/toolsrtm)

The **ToolsRTM** package provides a comprehensive suite of tools for
simulating canopy reflectance using various radiative transfer (RT)
models at multiple satellite resolutions. Currently in the testing
phase, this package is designed to facilitate detailed simulations,
enabling versatile and accurate analyses of canopy reflectance
characteristics.

[**ToolsRTM documentation — A Package for Simulating RT
Models**](https://ccgcam.github.io/RTM-Suite/toolsrtm/index.html)

For canopy-level simulations, the package features models such as
**INFORM**, **fourSAIL**, and **fourSAIL2**. When it comes to leaf-level
simulations, it includes the **PROSPECT** model (with D and PRO
variants), **Liberty**, and **FLUSPECT** (B-Cx). These models empower
users to conduct sophisticated simulations that capture the intricate
dynamics of reflectance behavior in both leaf and canopy contexts.

![](reference/figures/rtm_sims.png)

**Fig. 1.** Simulations performed with ToolsRTM Package for several
radiative transfer models .

Additionally, the **SPART** model (Soil-Plant-Atmosphere Radiative
Transfer model) is tailored for satellite measurements in the solar
spectrum. It integrates three computationally efficient RT models: the
**BSM** model for soil, **PROSAIL** for vegetation canopies, and
**SMAC** for the atmosphere. These components are interconnected using
the four-stream theory and the adding method, allowing SPART to simulate
directional top-of-atmosphere (TOA) spectral observations. This approach
accounts for significant effects, such as sun-observer geometries and
the non-Lambertian reflectance of the land surface.

**ToolsRTM** includes **MARMIT** (Multilayer rAdiative tRansfer Model of
soIl reflecTance), a radiative transfer model for simulating the
spectral reflectance of bare soil as a function of surface water
content. MARMIT simulates soil reflectance across the **400–2500 nm**
spectral range at **1 nm resolution**. The high-level
[`get.marmit.rsoil()`](reference/get.marmit.rsoil.md) wrapper, built on
[`get.marmit1()`](reference/get.marmit1.md) and
[`get.marmit2()`](reference/get.marmit2.md), loads a dry-soil reference
spectrum and runs the MARMIT simulation directly.

The package includes the **Bablet 2016** soil database by default to
keep the installation lightweight. Seven additional MARMIT soil
databases are available in the RTM-Suite [`databases/`](../databases/)
directory:

- Dupiau 2020
- Humper 2015
- Lesaignoux 2008
- Liu 2002
- Lobell 2002
- Marcq 2012
- Philpot 2014

Together, these databases represent approximately **200 MB** of
additional soil spectral data. They can be used directly with
[`get.marmit.rsoil()`](reference/get.marmit.rsoil.md) by specifying
their location through the `db_root` argument, without copying the files
into the R package.

``` r

soil <- ToolsRTM::get.marmit.rsoil(database = "Liu_2002", id = 1, db_root = "databases")
```

### RTM-Suite: A Unified R and Python Framework

**ToolsRTM** is one library within **RTM-Suite**, which links both the R
packages (`ToolsRTM`, `SCOPEinR`) and their Python ports (`toolsrtm`,
`scopeinpython`) behind one common site — with reference manuals, worked
tutorials, and runnable example pipelines for both languages side by
side.

![RTM-Suite website](reference/figures/Website_rtm-suite.png)

RTM-Suite website

**Fig. 2.** The [RTM-Suite website](https://ccgcam.github.io/RTM-Suite/)
— see **Documentation** for R/Python reference manuals, **Tutorials**
for step-by-step walkthroughs (R and Python side by side), and
**Examples** for copy-paste runnable code with real generated figures.

### Installation in R

**ToolsRTM** requires R 4.3 or newer. Install the development version
directly from GitLab:

``` r

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

if (!requireNamespace("ToolsRTM", quietly = TRUE)) {
  remotes::install_gitlab("caminoccg/toolsrtm", upgrade = "never")
}
```

Alternatively, install a downloaded source archive:

``` r

install.packages(
  "path/to/toolsrtm-main.tar.gz",
  repos = NULL,
  type = "source"
)
```

Check the installed version:

``` r

packageVersion("ToolsRTM")
```

The version described by this README is 0.62.5.

### Quick example

The following example creates a small lookup table and runs PROSPECT-PRO
coupled to fourSAIL. Both [`getLUT()`](reference/getLUT.md) and
[`simulate_RTM()`](reference/simulate_RTM.md) are exported by the
current package.

``` r

inputs <- ToolsRTM::inputsPROSAIL
lut <- as.data.frame(
  ToolsRTM::getLUT(inputs = inputs, nLUT = 30, setseed = 1234)
)
rsoil <- rep(0.2, 2101)

simulation <- ToolsRTM::simulate_RTM(
  inputLUT = lut[1, , drop = FALSE],
  rsoil = rsoil,
  leaf.model = "PROSPECT-PRO",
  canopy.model = "fourSAIL"
)
```

### SCOPEinR

Install [SCOPEinR](https://gitlab.com/caminoccg/scopeinr) when the
workflow requires the SCOPE model. ToolsRTM is installed first because
SCOPEinR imports it.

``` r

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_gitlab("caminoccg/toolsrtm", upgrade = "never")
remotes::install_gitlab("caminoccg/scopeinr", upgrade = "never")
packageVersion("SCOPEinR")
```

For an offline installation, replace the second `install_gitlab()` call
with:

``` r

install.packages(
  "path/to/scopeinr-main.tar.gz",
  repos = NULL,
  type = "source"
)
```

### Manuals and RTM-Suite resources

| Resource | Description |
|----|----|
| [RTM-Suite documentation](https://ccgcam.github.io/RTM-Suite/) | Entry point for the complete R suite |
| [ToolsRTM reference](https://ccgcam.github.io/RTM-Suite/toolsrtm/index.html) | Function reference and package articles |
| [Tutorials 01-18](https://ccgcam.github.io/RTM-Suite/toolsrtm/articles/index.html) | The full numbered tutorial series: leaf/canopy RT, SPART, model comparison, LUTs, parallel simulation, hyperspectral sensors, indices, sensitivity, hybrid/ML/DL inversion, end-to-end pipeline, real Sentinel-2 (STAC) applications, MARMIT+SPART soil integration, and a forest time-series capstone |
| [Real EO application](https://ccgcam.github.io/RTM-Suite/toolsrtm/articles/t15-real-eo-application.html) | Real Sentinel-2 retrieval (STAC) and a genuine spatial trait map |
| [Reference manuals](https://ccgcam.github.io/RTM-Suite/toolsrtm/articles/index.html) | The older comprehensive `ToolsRTM`/`Getting-LUTs`/`InversionOpt` manuals, kept alongside the tutorial series |
| [SCOPEinR reference](https://ccgcam.github.io/RTM-Suite/scopeinr/index.html) | SCOPE model documentation |
| [SCOPEinR tutorials 01-11](https://ccgcam.github.io/RTM-Suite/scopeinr/articles/index.html) | Energy balance, fluorescence, hybrid inversion, and a real-Sentinel-2 photosynthesis capstone |
| [Tutorials](../Tutorials/) | Reproducible R tutorials included in RTM-Suite |
| [Pipeline scripts](../Scripts/R/Pipeline/) | Adaptable simulate-to-invert workflows |

### Interactive application

Use
[`Apps/RTMs`](https://github.com/CCGCAM/RTM-Suite/tree/main/Apps/RTMs)
(run locally via `shiny::runApp("Apps/RTMs")`) to configure models and
inspect simulations interactively, no code required. Its source lives
alongside the rest of RTM-Suite; the Shiny application is not launched
through the core ToolsRTM API.

![](reference/figures/shiny.png)

**Fig. 3.** Interactive reflectance simulator using PROSAIL model based
on shiny app.

### Citation

If you use **ToolsRTM** or **SCOPEinR**, please consider citing:

1.  Camino et al. (2024). **RT-Simulator: An Online Platform to Simulate
    Canopy Reflectance from Biochemical and Structural Plant Properties
    Using Radiative Transfer Models**. *IGARSS 2024*, Athens, Greece,
    pp. 2811-2814. [doi:
    10.1109/IGARSS53475.2024.10642442](https://doi.org/10.1109/IGARSS53475.2024.10642442)

2.  Arano et al. (2024). **Enhancing Chlorophyll Content Estimation with
    Sentinel-2 Imagery: A Fusion of Deep Learning and Biophysical
    Models**. *IGARSS 2024*, Athens, Greece, pp. 4486-4489. doi:
    [10.1109/IGARSS53475.2024.10641613](https://doi.org/10.1109/IGARSS53475.2024.10641613)

3.  Camino et al. (in preparation). **Integrating Physiological Plant
    Traits with Sentinel-2 Imagery for Monitoring Gross Primary
    Production and Detecting Forest Disturbances**.

#### **References**

**Leaf models**

- Jacquemoud, S., Baret, F. (1990). *PROSPECT: A model of leaf optical
  properties spectra.* Remote Sensing of Environment, 34(2), 75-91.
  [10.1016/0034-4257(90)90100-Z](https://doi.org/10.1016/0034-4257(90)90100-Z)

- Féret, J.-B. et al. (2017). *PROSPECT-D: Towards modeling leaf optical
  properties through a complete lifecycle.* Remote Sensing of
  Environment, 193, 204-215.
  [10.1016/j.rse.2017.03.004](https://doi.org/10.1016/j.rse.2017.03.004)
  (PROSPECT-D, adds anthocyanins)

- Féret, J.-B. et al. (2021). *PROSPECT-PRO for estimating content of
  nitrogen-containing leaf proteins and other carbon-based
  constituents.* Remote Sensing of Environment, 252, 112173.
  [10.1016/j.rse.2020.112173](https://doi.org/10.1016/j.rse.2020.112173)
  (PROSPECT-PRO, splits dry matter into protein + carbon-based
  constituents)

- Dawson, T.P., Curran, P.J., Plummer, S.E. (1998). *LIBERTY — Modelling
  the effects of leaf biochemical concentration on reflectance spectra.*
  Remote Sensing of Environment, 65(1), 50-60.
  [10.1016/S0034-4257(98)00007-8](https://doi.org/10.1016/S0034-4257(98)00007-8)

- Vilfan, N., van der Tol, C., Muller, O., Rascher, U., Verhoef, W.
  (2016). *Fluspect-B: A model for leaf fluorescence, reflectance and
  transmittance spectra.* Remote Sensing of Environment, 186, 596-615.
  [10.1016/j.rse.2016.09.017](https://doi.org/10.1016/j.rse.2016.09.017)

- Vilfan, N., Van der Tol, C., Yang, P., Wyber, R., Malenovský, Z.,
  Robinson, S.A., Verhoef, W. (2018). *Extending Fluspect to simulate
  xanthophyll driven leaf reflectance dynamics.* Remote Sensing of
  Environment, 211, 345-356.
  [10.1016/j.rse.2018.04.012](https://doi.org/10.1016/j.rse.2018.04.012)
  (Fluspect-B-Cx, adds the xanthophyll/Cx de-epoxidation state)

**Canopy models**

- Verhoef, W. (1984). *Light scattering by leaf layers with application
  to canopy reflectance modeling: The SAIL model.* Remote Sensing of
  Environment, 16(2), 125-141.
  [10.1016/0034-4257(84)90057-9](https://doi.org/10.1016/0034-4257(84)90057-9)

- Verhoef, W. (1998). *Theory of radiative transfer models applied in
  optical remote sensing of vegetation canopies.* PhD thesis, Wageningen
  University. (4SAIL, the extended/corrected SAIL version this suite’s
  `fourSAIL` implements)

- Verhoef, W., Bach, H. (2007). *Coupled soil-leaf-canopy and atmosphere
  radiative transfer modeling to simulate hyperspectral multi-angular
  surface reflectance and TOA radiance data.* Remote Sensing of
  Environment, 109, 166-182.
  [10.1016/j.rse.2006.12.013](https://doi.org/10.1016/j.rse.2006.12.013)
  (introduces 4SAIL2, this suite’s `fourSAIL2`)

- Atzberger, C. (2000). *Development of an invertible forest reflectance
  model: The INFOR-model.* In: *A Decade of Trans-European Remote
  Sensing Cooperation*, Proceedings of the 20th EARSeL Symposium,
  Dresden, Germany, 39-44. (no DOI, conference proceedings)

**Soil, atmosphere & SCOPE**

- Bablet, A., Vu, P.V.H., Jacquemoud, S., Viallefont-Robinet, F., Fabre,
  S., Briottet, X., Sadeghi, M., Whiting, M.L., Baret, F., Tian, J.
  (2018). *MARMIT: a multilayer radiative transfer model of soil
  reflectance to estimate surface soil moisture content in the solar
  domain (400-2500 nm).* Remote Sensing of Environment, 217:1-17.
  [10.1016/j.rse.2018.07.031](https://doi.org/10.1016/j.rse.2018.07.031)

- Dupiau, A., Jacquemoud, S., Briottet, X., Fabre, S.,
  Viallefont-Robinet, F., Philpot, W., Di Biagio, C., Hébert, H.,
  Formenti, P. (2022). *MARMIT-2: an improved version of the MARMIT
  model to predict soil reflectance as a function of surface water
  content in the solar domain.* Remote Sensing of Environment,
  272:112951.
  [10.1016/j.rse.2022.112951](https://doi.org/10.1016/j.rse.2022.112951)

- Rahman, H., Dedieu, G. (1994). *SMAC: a simplified method for the
  atmospheric correction of satellite measurements in the solar
  spectrum.* International Journal of Remote Sensing, 15(1), 123-143.

- Yang, P., van der Tol, C., Yin, T., Verhoef, W. (2020). *The SPART
  model: A soil-plant-atmosphere radiative transfer model for satellite
  measurements in the solar spectrum.* Remote Sensing of Environment,
  247, 111870.
  [10.1016/j.rse.2020.111870](https://doi.org/10.1016/j.rse.2020.111870)

- Van der Tol, C., Verhoef, W., Timmermans, J., Verhoef, A., Su, Z.
  (2009). *An integrated model of soil-canopy spectral radiances,
  photosynthesis, fluorescence, temperature and energy balance.*
  Biogeosciences 6(12), 3109-29.
  [10.5194/bg-6-3109-2009](https://doi.org/10.5194/bg-6-3109-2009)

- Yang, P., Prikaziuk, E., Verhoef, W., van der Tol, C. (2021). *SCOPE
  2.0: A model to simulate vegetated land surface fluxes and satellite
  signals.* Geoscientific Model Development, 14, 4697-4712.
  [10.5194/gmd-14-4697-2021](https://doi.org/10.5194/gmd-14-4697-2021)

### License

[![RTM-Suite code:
MIT](https://img.shields.io/badge/RTM--Suite%20code-MIT-yellow.svg)](#id_0)
[![Ported GPL models:
GPL--3.0](https://img.shields.io/badge/Ported%20GPL%20models-GPL--3.0-blue.svg)](#id_0)
[![Model
licenses](https://img.shields.io/badge/Per--model%20licenses-THIRD__PARTY__LICENSES.md-informational.svg)](#id_0)

**ToolsRTM** is a port of several radiative transfer models bundled
behind one common R interface, and not all of them carry the same
license. Three of the bundled models – **Fluspect-B**, **fourSAIL2**,
and **SPART** – are ports of GPL-3.0-licensed original models, and
GPL-3.0 requires any combined work incorporating GPL-3.0 code to be
distributed as GPL-3.0 as a whole. `ToolsRTM` is therefore distributed
under **GPL-3.0** (`License: GPL-3` in `DESCRIPTION`), matching its
Python port `toolsrtm`’s own license.

Within that GPL-3.0 distribution, two kinds of code coexist:

- **The ported radiative transfer models themselves** (leaf and canopy
  functions listed in
  [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md)) – GPL-3.0 for
  the three GPL-derived ports above, and independently MIT-licensable
  for the rest (e.g. PROSPECT-D/-PRO).
- **Original utilities developed in this package**– LUT generation
  (`getLUT`, `getLUTs`, `getLUT_liberty`), spectral indices
  (`getIndices`, `getIndices_SE2*`, `getSpectraIndices`), sensor
  convolution wrappers (`Spectral.convolution`,
  `get.spectral.convolution.*`, `get.smac`, `get.coef.SMAC`),
  sensitivity analysis (`get.sobol.indices`), and the trait-inversion
  tooling (`get.inversion`, `get.inversionOpt`, `hybrid_inversion`,
  `hybrid_inversionE`, `carspls`, `get.cars.pls`, `getVIF`) – are
  original, independent work and are **MIT** individually. Because
  GPL-3.0 requires the combined, distributed package to be GPL-3.0 as a
  whole, the package you install is still `License: GPL-3` end to end;
  the MIT notice above is about authorship/reuse of those specific
  original files on their own, not a separate installable subset.

See [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) for the license
and source-code provenance of every individual model this package
implements. Always cite the original publication(s) of each model you
use, in addition to citing RTM-Suite/ToolsRTM.
