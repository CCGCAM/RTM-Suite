# ToolsRTM in Python

[![ToolsRTM (R)](https://img.shields.io/badge/ToolsRTM-GitLab-FC6D26?logo=gitlab&logoColor=white)](https://gitlab.com/caminoccg/toolsrtm) [![RTM-Suite monorepo](https://img.shields.io/badge/python%2Ftoolsrtm-RTM--Suite-181717?logo=github&logoColor=white)](https://github.com/CCGCAM/RTM-Suite/tree/main/python/toolsrtm)

Python port of [`ToolsRTM`](https://gitlab.com/caminoccg/toolsrtm) (R) — **leaf and canopy radiative transfer models for optical remote sensing of vegetation**. Simulate top-of-canopy reflectance and fluorescence from leaf biochemistry and canopy structure, convolve simulated spectra to satellite and airborne sensor bands, and invert observations to retrieve biophysical traits.

Every function is verified against a real, unmodified call to the original R package — outputs match to floating-point precision unless noted otherwise. See [`verification.rst`](https://github.com/CCGCAM/RTM-Suite/blob/main/python/docs/verification.rst) for full per-function verification details.

**How this fits together:** this package is not an independent reimplementation of the R code. It is developed, tested, and verified within the same [`RTM-Suite`](https://github.com/CCGCAM/RTM-Suite) monorepo as `ToolsRTM` (R), helping both implementations remain consistent.

It complements [`scopeinpython`](https://github.com/CCGCAM/scopeinpython), the Python port of `SCOPEinR`. `scopeinpython` uses `toolsrtm` for leaf optical models such as PROSPECT-D and PROSPECT-PRO, mirroring the relationship between `SCOPEinR` and `ToolsRTM` in R.

Together, the libraries provide an end-to-end workflow for vegetation radiative transfer modelling: **leaf traits → leaf optics → canopy reflectance and fluorescence → sensor bands → SCOPE energy balance and fluorescence → retrieval of biophysical traits using classical machine learning or deep learning.**

## What's in it

| Component | What it does |
|------------------------------------|------------------------------------|
| **Leaf optics** | PROSPECT-D, PROSPECT-PRO, Fluspect-B, Fluspect-B-Cx, Liberty -- reflectance, transmittance, and (Fluspect) sun-induced fluorescence (SIF) |
| **Canopy models** | fourSAIL, fourSAIL2 (two-layer green/brown), INFORM (forest understorey) -- bidirectional canopy reflectance from LAI, leaf angle distribution, soil, and sun/view geometry |
| **Sensor convolution** | Three functions covering every real-world case: a measured SRF + atmospheric-correction coefficients (`smac.py`), a measured SRF alone (`srf.py`, PRISMA/Sentinel-2A/2B), or just nominal center+FWHM -- including **your own sensor or camera** (`srf.py`'s Gaussian convolution, also covers EnMAP, Landsat, MODIS, Hyperion, WorldView-2, ...) |
| **Spectral indices** | The common vegetation/water/pigment indices computed from a simulated or convolved spectrum |
| **Atmospheric correction** | SPART's TOC -\> TOA path via SMAC -- all 9 sensors R ships are bundled (Landsat 4/5/7/8, Sentinel-2A/B, Sentinel-3A/B, Terra/Aqua MODIS) |
| **Trait inversion** (`toolsrtm.inversion`) | CARS-PLS and VIF-based predictor selection, LUT nearest-neighbour ("merit function") matching, and a 12-algorithm ML dispatcher (PLSR/SVM/RF/GB/NN/Bayesian/AdaBag/BRNN/xGB/RVM/qLASSO/Ensemble via scikit-learn/xgboost) plus feature-selection wrappers (`hybrid_inversion`). Needs the optional `ml` extra: `pip install "toolsrtm[ml]"` |
| **Deep-learning inversion** (`toolsrtm.deep_learning`, optional) | Dense and 1D-CNN Keras architectures for trait inversion, matching R's `getMLmodel`. Not required for the rest of the package -- needs the optional `dl` extra: `pip install "toolsrtm[dl]"` (TensorFlow). scikit-learn's own estimators (above) cover most trait-inversion needs without this. |
| **Satellite retrieval** (`toolsrtm.satellite`, optional) | Search and download real scenes via STAC (Microsoft Planetary Computer or AWS Earth Search) for a bounding box/date range, and build a cropped multi-band data cube -- Sentinel-2 L2A, Landsat C2 L2, and 6 MODIS products. Needs the optional `stac` extra: `pip install "toolsrtm[stac]"` and live network access. |

`toolsrtm` is one library within [**RTM-Suite**](https://ccgcam.github.io/RTM-Suite/), which links both the R packages (`ToolsRTM`, `SCOPEinR`) and their Python ports (`toolsrtm`, `scopeinpython`) behind one common site -- with reference manuals, worked tutorials, and runnable example pipelines for both languages side by side.

![RTM-Suite website](docs/images/Website_rtm-suite.png)

The [RTM-Suite website](https://ccgcam.github.io/RTM-Suite/) -- see **Documentation** for R/Python reference manuals, **Tutorials** for step-by-step walkthroughs (R and Python side by side), and **Examples** for copy-paste runnable code with real generated figures.

## Install

``` bash
pip install git+https://github.com/CCGCAM/ToolsRTMinPython.git
```

or, editable, from a local clone:

``` bash
git clone https://github.com/CCGCAM/ToolsRTMinPython.git
cd ToolsRTMinPython
pip install -e ".[test]"
pytest tests -q
```

## Quick example

``` python
import numpy as np
from toolsrtm import foursail, compute_brf

row = dict(
    N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0.1, EWT=0.01, LMA=0.009, alpha=40,
    LAI=3, hspot=0.01, LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,
    tts=30, tto=0, psi=0,
)
rsoil = np.full(2101, 0.15)

sail = foursail(row, rsoil, leaf_model="PROSPECT-D")
reflectance = compute_brf(sail.rdot, sail.rsot, row["tts"])
```

<p align="center">

<img src="docs/images/simulated_spectrum.png" alt="Simulated canopy reflectance spectrum" width="70%"/>

</p>

Sweep a trait (e.g. chlorophyll content) and see how the whole spectrum responds:

<p align="center">

<img src="docs/images/sensitivity_sweep.png" alt="Sensitivity sweep across Cab" width="70%"/>

</p>

Convolve onto a real sensor's bands (here: Sentinel-2A) -- coarser and fewer bands than the native 1nm simulation, exactly what a satellite actually observes:

<p align="center">

<img src="docs/images/sensor_convolution.png" alt="Convolved onto Sentinel-2A bands" width="70%"/>

</p>

## Documentation & tutorials

This repo is deliberately just the installable package -- everything else (manuals, worked tutorials, the course materials this was built for) lives in [`CCGCAM/RTM-Suite`](https://github.com/CCGCAM/RTM-Suite), the monorepo this package is developed and verified in and kept in sync with:

- **Manual** (every function's full docstring, browsable): [`docs/python/index.html`](https://github.com/CCGCAM/RTM-Suite/blob/main/docs/python/index.html)
- **The R tutorial series, 01-18** ([`docs/toolsrtm/articles/index.html`](https://ccgcam.github.io/RTM-Suite/toolsrtm/articles/index.html)) is this package's fullest worked-example coverage -- topic-to-module bridge, since this package doesn't mirror every R tutorial 1:1 (deliberate model differences are part of this suite's own design, not a gap):
  - Getting started / leaf-to-canopy / model comparison (R 01-02, 04) -\> `toolsrtm.leaf`, `toolsrtm.fluspect`, `toolsrtm.liberty`, `toolsrtm.canopy`, `toolsrtm.inform`
  - SPART, soil-plant-atmosphere (R 03) -\> `toolsrtm.spart`, `toolsrtm.smac`
  - Sensor convolution incl. hyperspectral/VNIR (R 07-08) -\> `toolsrtm.srf` (measured-SRF, SMAC-bundled, and Gaussian-from-nominal-characteristics convolution, matching PRISMA/Sentinel-2/EnMAP/custom-sensor coverage)
  - Vegetation indices (R 09) -\> `toolsrtm.indices`
  - Hybrid/ML inversion (R 11-12) -\> `toolsrtm.inversion`
  - Deep learning (R 13) -\> `toolsrtm.deep_learning` (TensorFlow + PyTorch, not R's single-backend `getMLmodel()`)
  - Real EO application via STAC (R 15, 17-18) -\> `toolsrtm.satellite`
  - MARMIT soil integration (R 16) -\> `toolsrtm.marmit`
  - **Not currently ported**: R's LUT-distribution helpers (`get_distributionLUT()`/`getCor()`, R Tutorial 05) and its Sobol/Johnson sensitivity tooling (R Tutorial 09-10's ToolsRTM equivalent) have no Python module yet -- a real gap, not something this package silently works around.
- **Tutorials, step by step** (simulate -\> sweep a trait -\> convolve onto a sensor, incl. your own sensor/camera -\> invert with ML, in R side-by-side with Python): [`Tutorials/How-in-Python.ipynb`](https://github.com/CCGCAM/RTM-Suite/blob/main/Tutorials/How-in-Python.ipynb) (R version: `How-in-R.Rmd`) -- matches the `Apps/RTMs` Shiny app's own **"How in Python"** tab
- **Complete R reference manual** (every leaf/canopy model, trait sampling, all 12 inversion algorithms via `caret`, TensorFlow/Keras deep learning): [`Tutorials/ToolsRTM_PROSAIL_tutorial.Rmd`](https://github.com/CCGCAM/RTM-Suite/blob/main/Tutorials/ToolsRTM_PROSAIL_tutorial.Rmd) -- the Python equivalents of its inversion sections are `toolsrtm.inversion`/`toolsrtm.deep_learning` (this package) plus the runnable pipeline scripts linked below
- **A real, runnable pipeline script** (simulate a LUT -\> compute indices -\> train an ML inversion model): [`Scripts/Python/README.md`](https://github.com/CCGCAM/RTM-Suite/blob/main/Scripts/Python/README.md)
- **Full writeup** -- every function ported, with its numerical verification status: [`python/README.md`](https://github.com/CCGCAM/RTM-Suite/blob/main/python/README.md)

## License

[![RTM-Suite code: MIT](https://img.shields.io/badge/RTM--Suite%20code-MIT-yellow.svg)](#0) [![Ported GPL models: GPL--3.0](https://img.shields.io/badge/Ported%20GPL%20models-GPL--3.0-blue.svg)](#0) [![Model licenses](https://img.shields.io/badge/Per--model%20licenses-THIRD__PARTY__LICENSES.md-informational.svg)](#0)

`toolsrtm` is a Python port of several radiative transfer models bundled behind one common interface, and not all of them carry the same license. Three of the bundled models -- **Fluspect-B**, **fourSAIL2**, and **SPART** -- are ports of GPL-3.0-licensed original models, and GPL-3.0 requires any combined work incorporating GPL-3.0 code to be distributed as GPL-3.0 as a whole. `toolsrtm` is therefore distributed under **GPL-3.0-only** (see [`LICENSE`](LICENSE)), matching its R sibling package `ToolsRTM`'s own `License: GPL-3` field.

Within that GPL-3.0 distribution, two kinds of code coexist:

- **The ported radiative transfer models themselves** (`leaf`, `liberty`, `fluspect`, `canopy`, `inform`, listed individually in [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md)) -- GPL-3.0 for the three GPL-derived ports above, and independently MIT-licensable for the rest (e.g. PROSPECT-D/-PRO in `leaf`).
- **Original utilities Carlos Camino wrote on top of those models** -- sensor convolution (`srf`, `smac`), spectral indices (`indices`), satellite/STAC retrieval (`satellite`), and the trait-inversion tooling (`inversion`, `deep_learning`) -- are original, independent work and are **MIT** individually. Because GPL-3.0 requires the combined, distributed package to be GPL-3.0 as a whole, the package you `pip install` is still GPL-3.0-only end to end; the MIT notice above is about authorship/reuse of those specific original modules on their own, not a separate installable subset.

See [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md) for the license and source-code provenance of every individual model this package implements. Always cite the original publication(s) of each model you use, in addition to citing RTM-Suite/ToolsRTM.
