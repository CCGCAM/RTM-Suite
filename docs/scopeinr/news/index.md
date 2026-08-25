# Changelog

## SCOPEinR 0.46

### Documentation

- Added a full numbered tutorial series (11 vignettes,
  `01-getting-started` through `11-photosynthesis-capstone`), covering
  the SCOPE input structure, soil/canopy BRDF, energy balance,
  fluorescence (SIF), LUT construction, parallel simulation, sensitivity
  (direct vs. indirect trait effects), hybrid ML inversion, and whether
  SIF adds retrievable information about photosynthesis beyond
  reflectance-based greenness.
- Closes with a capstone applying the trained model to a real Sentinel-2
  time series and spatial map (via STAC) over a real forest site
  (Speulderbos, NL), including an explicit accuracy comparison between a
  SIF-inclusive model (not valid on real Sentinel-2 data, since
  Sentinel-2 cannot observe SIF) and the reflectance-only model actually
  applied to the real retrieval.
- The pre-existing comprehensive manuals (`SCOPEinR`,
  `getting-luts-scope`, `scope-pipeline`, `sensitivity-scope`,
  `inversion-scope`, `sif-photosynthesis-proxy`) are kept alongside the
  new series as deeper reference material.
- The live-STAC capstone tutorial is a pkgdown-only article (not built
  or checked as a package vignette), so `R CMD check`/CRAN builds stay
  fast and network-independent.

## SCOPEinR 0.1.0 - 0.46 (prior to this changelog)

- Initial development: full R port of SCOPE v2.1 (Van der Tol et
  al. 2009; Yang et al. 2021) – Fluspect-Cx leaf optics, RTMo canopy
  optical BRDF, Farquhar/Collatz photosynthesis and fluorescence-yield
  biochemistry, and the full canopy energy-balance solve
  ([`get.SCOPE()`](../reference/get.SCOPE.md)/
  [`get.SCOPE.parallel()`](../reference/get.SCOPE.parallel.md)), coupled
  with `ToolsRTM` for LUT generation and sensor convolution. See git
  history for the detailed development record predating this file.
