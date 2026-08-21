# ToolsRTM 0.62.5

## Documentation

* Replaced the single comprehensive `ToolsRTM` vignette and a handful of
  standalone articles with a full numbered tutorial series (18 vignettes,
  `01-getting-started` through `18-spatial-index-mapping`), covering leaf/
  canopy radiative transfer, SPART, model comparison, LUT construction,
  parallel simulation, hyperspectral/VNIR sensor convolution, vegetation
  indices, sensitivity analysis, hybrid/ML/deep-learning inversion, a full
  end-to-end pipeline, and closing with real-data applications: live
  Sentinel-2 retrieval via STAC with genuine spatial trait/index maps over
  two real Dutch forest sites (Speulderbos, Loobos), and a real multi-date
  forest time series.
* The pre-existing comprehensive manuals (`ToolsRTM`, `Getting-LUTs`,
  `InversionOpt`, `course-pipeline`, `model-comparison-and-sensitivity`)
  are kept alongside the new series as deeper reference material.
* Network-dependent (live STAC/Sentinel-2) and TensorFlow/PyTorch-dependent
  tutorials moved to `vignettes/articles/` (pkgdown-only articles, not
  built or checked as package vignettes) so `R CMD check`/CRAN builds stay
  fast and network-independent.

## Bug fixes

Real bugs found and fixed while building/verifying the tutorial series
above (each confirmed against known physics or by direct numerical
comparison, not just code inspection):

* `get.plots` argument handling in `SPART()` ignored the caller's actual
  value and forced `TRUE` whenever the argument was supplied at all.
* `Compute_BRF()`'s `short.waves` handling silently recycled solar-
  irradiance vectors instead of truncating them for Fluspect's shorter
  spectral domain when the argument was omitted.
* `foursail2()` never truncated `rsoil` for Fluspect leaf models, silently
  recycling a full-length soil spectrum against the shorter leaf-optics
  domain throughout its SAIL2 core.
* `getFluspect.B()`/`getFluspect.Cx()` returned `NULL` (rather than valid
  zero-fluorescence output) whenever `fqe = 0`, crashing every downstream
  caller with a zero-length reflectance/transmittance vector.
* `get.spectral.convolution.rfl()` mislabelled band-center wavelengths for
  6 of 9 bundled sensors (Landsat 4/5/7, Sentinel-3A/B OLCI, Terra/Aqua
  MODIS) due to a column-ordering assumption that only happened to hold
  for Sentinel-2.
* `get.sentinel2_cube()` and `get.satellite_collection()` called `rstac`/
  `gdalcubes`/`terra` functions without namespace-qualifying them, so a
  fresh session with only `library(ToolsRTM)` attached failed with
  "could not find function" errors; both are now fully namespace-qualified
  and work standalone.
* `get.inversion(algorithm = "Ensemble")` referenced an undefined `fmla.n`
  variable and passed the wrong argument shape to `caret::createFolds()`
  in its internal `caretEnsemble` stacking step.
* `get.sobol.indices()`'s `STi` (total Sobol index) column had a formula
  bug (never referenced the actual input variable), producing a near-
  uniform, physically meaningless split across traits; `get.spectral.
  sensitivity()` now uses the function's separately-correct
  `I.Johnson_norm` column instead.
* `get.sentinel2_cube()`'s NDVI-band-arithmetic output was checked against
  the wrong class name (`"RasterLayer"` instead of `"SpatRaster"`) after
  the raster-to-terra migration in downstream Shiny apps, silently
  disabling the NDVI-included GeoTIFF download path.

## Internal

* Migrated every remaining use of the archived/legacy `raster` package to
  `terra` across the package source and the bundled Shiny apps (`raster`
  removed from `Imports`); repo-wide `grep -rl "raster::"` (excluding
  generated documentation) now returns nothing.
* Removed two redundant staging copies from `TOcheck/` that were already
  properly ported into the package source.

# ToolsRTM 0.1.0 - 0.62.5 (prior to this changelog)

* Initial development: PROSPECT-D/-PRO, Liberty, Fluspect-B/-Cx leaf
  models; fourSAIL, fourSAIL2, INFORM canopy models; SPART soil-plant-
  atmosphere model; MARMIT soil reflectance model; LUT generation and
  correlated sampling; sensor convolution (measured SRF, SMAC-bundled,
  and Gaussian-from-nominal-characteristics); vegetation indices; hybrid,
  classical-ML (`caret`, 12 algorithms), and deep-learning (TensorFlow/
  Keras) trait inversion; real Sentinel-2 retrieval via STAC. See git
  history for the detailed development record predating this file.
