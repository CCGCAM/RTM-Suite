# RTM-Suite verification status

Rolling checklist, updated as each item is actually run and confirmed
(not just read). ✅ = executed and confirmed correct. 🔄 = in progress.
⬜ = not yet checked this pass.

## Cross-language robustness

- ✅ ToolsRTM (R) vs toolsrtm (Python), **all 15 leaf×canopy combinations** (PROSPECT-D/PROSPECT-PRO/Liberty/Fluspect-B/Fluspect-B-Cx × fourSAIL/foursail2/INFORM), 100 sims each, independent shared LUT per combination — 15/15 combos R²=1.00000000; RMSE ≈1e-16 for PROSPECT-D/PRO/Fluspect, ≈1e-12 for Liberty (its own iterative-solver precision, not a discrepancy)
  - Found+fixed during this pass: Python comparison harness was missing `Prot`/`CBC` columns, silently zeroing out PROSPECT-PRO (0/100) and degrading Fluspect-B/-B-Cx (Fluspect-B was giving R²≈0.17-0.53 before the fix — real bug in the *comparison script*, not the ported model; `foursail()`/`foursail2()`/`inform()` themselves were correct all along)
- ✅ ToolsRTM (R) vs toolsrtm (Python), PROSPECT-D+fourSAIL, 600 sims, 1 shared LUT — R²=1.00000000, RMSE=8.1e-16
- ✅ SCOPEinR (R) vs scopeinpython (Python), full SCOPE 2.1, 600 sims, 1 shared LUT — R²=1.00000000, RMSE=2.5e-16
- ✅ SPART (R) vs spart_toa/spart_toc (Python), 100 sims, Sentinel-2A bands — TOC R²=1.00000000 RMSE=4.2e-16, TOA R²=1.00000000 RMSE=6.4e-15
- ✅ MARMIT-1 and MARMIT-2 (both versions, R vs Python), 60 sims each (Bablet_2016 soil DB, random id/L/eps) — both R²=1.00000000, RMSE=3.5e-10 (marmit1) / 1.8e-8 (marmit2)
- ✅ **Real-data validation (not R vs Python)**: PROSPECT-D (R) vs the bundled Angers database, 308 real measured leaves (measured Cab/Car/Anth/EWT/LMA in, N fit per leaf) — median R²=0.9912, RMSE=1.47% reflectance units, best-fit N range 1.0-2.7 matching published literature values
- Published: `docs/comparison.html`, linked from main nav + Documentation section, reorganized into clearly-labeled sections (leaf×canopy, soil/atmosphere, full SCOPE, real-data validation, deep-dive, reproduce)

## SCOPEinR — R package

- ✅ Default (`lite=1`) baseline — Actot=6.95, unchanged across every fix this session
- ✅ `lite=0` (full SCOPE) — fixed (RTMt_planck.R dimension bug), Actot=6.76
- ✅ `calc_vert_profiles=1` — fixed (RTMo.R missing lite-mode branch), real 24-layer profile confirmed
- ✅ `irradiance=2` (MODTRAN) — fixed (4 chained bugs in get.calcTOCirr), Actot=7.11
- ✅ `Fluorescence_model=2` (Magnani 2012) — fixed (was wired to value 1 instead of 2), Actot=17.17
- ✅ `mSCOPE=1`, `calc_rss_rbs=1` — fixed and verified
- ⬜ `Fluorescence_model=1` ("sigmoid Kn") — confirmed no separate implementation exists; documented honestly
- ⬜ `irradiance=1` (WithE) — runs, flagged for caution (possible Esun/Esky swap), not independently resolved
- ⬜ `simulation=1` (time-series mode) — needs different input format, out of scope for LUT-based tutorials

### SCOPEinR tutorials — ALL 12 CONFIRMED RENDERING CLEAN ✅
t00 (scope-options, rewritten+republished), t01 (getting-started), t02
(soil-canopy-brdf), t03 (energy-balance), t04 (fluorescence), t05
(building-scope-luts), t06 (parallel-scope-runs), t07 (sensitivity), t08
(hybrid-inversion), t09 (sif-photosynthesis), t10 (end-to-end-pipeline),
t11 (photosynthesis-capstone, incl. live STAC pull) — every one
regression-checked this session, 0 render errors.

## ToolsRTM — R package

- ✅ Core PROSPECT-D + fourSAIL simulation — smoke-tested directly
- ✅ SPART `psoil`/BSM real-parameter fix — verified via diff review
- ✅ All comment-cleanup edits this session — diffed individually, confirmed no logic change
- ✅ `get.spectral.sensitivity()` Johnson/Sobol fix — confirmed physically sane wavelength pattern

### ToolsRTM tutorials — ALL 19 CONFIRMED RENDERING CLEAN ✅
t01 (getting-started) through t19 (FLEX Cal/Val ESU heterogeneity, incl.
live STAC pull + RF trait retrieval at FLEX 300m scale) — every one
regression-checked this session (t13 deep-learning, t14 end-to-end, t15
real-EO, t17 forest-time-series all needed a package reinstall for their
`doParallel` workers to find `ToolsRTM`, then confirmed clean), 0 render
errors. t19 added this pass, based on user-supplied
`1FLORA_ESU_generation_from_S2.R`, fully ported from `raster::` to
`terra::`.

- ⬜ The detailed 9-tutorial scientific-review requests (03, 06, 11-17,
  content-quality asks like "add more figures", "use a larger LUT",
  "add a domain-gap check") from earlier in this session — tutorials
  all RUN correctly, but that separate, deeper content-quality pass has
  not been actioned in this direct-work phase.

## Python

- ✅ `toolsrtm` test suite — 100/100 passed
- ✅ `scopeinpython` test suite — 18/18 passed
- ✅ `python/docs/examples.rst` — 5 real generated figures added (leaf+canopy,
  sensor convolution, ML inversion, MARMIT soil, full SCOPE+SIF), Sphinx
  rebuilt and republished

## Scripts (sampled, not exhaustive)

- ✅ `Scripts/R/ForSCOPE/1-getSCOPE.R` — runs clean
- ✅ `Scripts/R/Sensibility/1-Sobol_spectral_sensitivity.R` — runs clean, produces output
- ⬜ Remaining `Scripts/R/*`, `Scripts/Python/*` — not individually re-run this session (most don't touch the R functions fixed tonight)
- ✅ `Scripts/Comparison/*` (new) — all 5 R+Python pairs run clean end-to-end (single-combo ToolsRTM, single-combo SCOPEinR, all-15-combos ToolsRTM, SPART+MARMIT, Angers real-data validation)

## Documentation

- ✅ `docs/comparison.html` — new, published
- ✅ `docs/python/` — rebuilt with example figures
- ⬜ Cross-references from tutorials/READMEs to the new comparison page
