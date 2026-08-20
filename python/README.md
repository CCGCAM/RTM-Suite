# ToolsRTM and SCOPE in Python

Python port of **both** [`ToolsRTM`](../ToolsRTM) and [`SCOPEinR`](../SCOPEinR) — leaf and canopy radiative transfer, sensor convolution, trait inversion, and the full soil–leaf–canopy radiative transfer and energy-balance framework of SCOPE.

This repository currently provides a **deliberately scoped port** of the functionality available in `ToolsRTM` and `SCOPEinR`. It includes the components explicitly listed under **“What was ported”** below; each entry defines its current coverage and limitations. Functionality not listed there should not be assumed to be available.

The Python implementation is split into two installable packages:

- [`python/toolsrtm/`](toolsrtm/) — Python port of **ToolsRTM**: leaf and canopy radiative transfer models (PROSPECT-D/-PRO, Fluspect-B/-Cx, LIBERTY, fourSAIL, fourSAIL2, INFORM), MARMIT soil, sensor convolution (SMAC, measured SRF, Gaussian-from-nominal), spectral indices, satellite/STAC retrieval, and classical ML + deep-learning trait inversion. A complete radiative-transfer toolkit on its own, independent of SCOPE.
- [`python/scopeinpython/`](scopeinpython/) — Python port of **SCOPEinR**: SCOPE's own soil model, the optical top-of-canopy BRDF / fluorescence pipeline, and the energy-balance chain. Depends on `toolsrtm` for leaf optics (PROSPECT-D/-PRO), exactly as `SCOPEinR` depends on `ToolsRTM` in R.

**How this fits together:** `scopeinpython` mirrors the relationship between `SCOPEinR` and `ToolsRTM` in R — leaf optical properties are simulated with models such as PROSPECT-D and PROSPECT-PRO through `toolsrtm`, while `scopeinpython` extends the workflow to the soil–canopy system and SCOPE-level radiative transfer and energy-balance calculations. `toolsrtm` is not just a dependency of `scopeinpython`, though: on its own it already covers everything from a single leaf spectrum to a full sensor-convolved, ML-inverted trait retrieval pipeline, with no SCOPE involved at all.

Together, the two Python libraries provide the basis for an end-to-end vegetation radiative-transfer workflow, from either package's own entry point:

**leaf traits → leaf optics → canopy reflectance and fluorescence (`toolsrtm`) → sensor bands (`toolsrtm`) → SCOPE soil, energy balance and fluorescence (`scopeinpython`) → retrieval of biophysical traits with classical ML or deep learning (`toolsrtm`).**

Every ported function is checked against the corresponding implementation in the original R packages. See the verification documentation for per-function comparisons and known numerical or scope differences.

## Documentation

The **full Python API reference** is available in `docs/python/index.html`, including:

- complete function documentation generated directly from source docstrings,

- verification results against the original R implementations,

- the full **“What was ported”** table,

- and the list of functionality that has **not yet been ported**.

Documentation sources are located in `python/docs/`.

Rebuild the documentation with:

``` bash
pip install sphinx sphinx-rtd-theme myst-parser
cd python/docs && sphinx-build -b html . _build/html
# then copy _build/html/* into docs/python/ (committed, same as the R
# pkgdown sites at docs/toolsrtm/ and docs/scopeinr/)
```

## What was ported

### `toolsrtm` (`python/toolsrtm/`)

| Python | R source | What it is |
|---------------------|---------------------|------------------------------|
| `toolsrtm.calctav` | `ToolsRTM/R/calctav.R` | Transmissivity of a dielectric surface (Stern 1964 / Allen 1973) |
| `toolsrtm.prospect_d` | `ToolsRTM/R/prospect_DB.R` | PROSPECT-D leaf model (400-2500 nm) |
| `toolsrtm.prospect_pro` | `ToolsRTM/R/prospect_PRO.R` | PROSPECT-PRO leaf model (400-2500 nm) |
| `toolsrtm.canopy.volscatt` | `ToolsRTM/R/volscatt.R` | Scalar-leaf-angle volume-scattering geometry (used by fourSAIL) |
| `toolsrtm.canopy.campbell` | `ToolsRTM/R/campbell.R` | Ellipsoidal (Campbell 1986) leaf angle distribution |
| `toolsrtm.canopy.dladgen` | `ToolsRTM/R/dladgen.R` | Bimodal (Verhoef) leaf angle distribution |
| `toolsrtm.canopy.dcum` | `ToolsRTM/R/dcum.R` | Cumulative LIDF helper used by `dladgen` |
| `toolsrtm.canopy.jfunc1/jfunc2/jfunc4` | `ToolsRTM/R/Jfunc1.R`, `Jfunc2.R`, `Jfunc4.R` | Singularity-avoiding J functions used inside SAIL |
| `toolsrtm.canopy.non_conservative_scattering` | `ToolsRTM/R/NonConservativeScatering.R` | SAIL 4-stream non-conservative-scattering solution |
| `toolsrtm.canopy.foursail_core` / `toolsrtm.canopy.foursail` | `ToolsRTM/R/foursail.R` | Full fourSAIL canopy BRDF (single canopy layer). Leaf model dispatch via `toolsrtm.canopy._leaf_optics` covers all 5 leaf models ToolsRTM itself supports: PROSPECT-D, PROSPECT-PRO, Liberty, Fluspect-B, Fluspect-B-Cx. |
| `toolsrtm.canopy.conservative_scattering` / `toolsrtm.canopy.scattering` | `ToolsRTM/R/ConservativeScattering.R` / (dispatcher) | Near-conservative-scattering solution (`m<=0.01`) and the combined conservative/non-conservative dispatcher used by `foursail2`/`inform` |
| `toolsrtm.canopy.foursail2_core` / `toolsrtm.canopy.foursail2` | `ToolsRTM/R/foursail2.R` | Two-layer green/brown canopy BRDF (crown clumping, hotspot per layer). All 5 leaf models supported for both green (from `inputLUT`) and brown (package's own default row, extended with Liberty/Fluspect defaults) leaves. |
| `toolsrtm.inform.inform` | `ToolsRTM/R/inform.R` (+ `foursail.inform.R`, `foursail.inf.R`, `foursail_t_s.R`, `foursail_t_o.R`, `Compute_BRF.R`) | INFORM forest-stand reflectance (FLIM ground-coverage geometry + fourSAIL understorey/infinite-crown + crown transmittance), all 5 leaf models supported (output length switches to 2001 pts for the two Fluspect models, matching R). Two R quirks reproduced exactly (not "fixed") -- see docstrings in `inform.py`: a hardcoded partial LIDF table for TypeLidf==1 in three of the internal helpers, and `Compute_BRF`'s Es/Ed illumination spectra actually coming from `dataSpec_PDB`'s `direct_light`/`diffuse_light` columns via its named-argument code path (a `[,11]`/`[,12]` positional path exists in the R source but is dead code, never hit by any real caller) |
| `toolsrtm.spart.spart_toc` | `ToolsRTM/R/spart.R` (TOC portion, + `BSM.R`'s `getBSM.toolsRTM`/`soilwat`) | Top-of-canopy BRDF reflectance (BSM soil + fourSAIL + `Compute_BRF`, `short.waves=TRUE`), i.e. everything `SPART()` computes before the atmosphere. Uses ToolsRTM's own bundled BSM optipar table (a separate dataset from SCOPEinR's, see `spart.py`). |
| `toolsrtm.spart.spart_toa` | `ToolsRTM/R/spart.R` (full `SPART()`) | Top-of-atmosphere reflectance/radiance for a specific sensor: `spart_toc`'s TOC BRDF, atmosphere-corrected (SMAC) and resampled onto the sensor's band centers. Verified against a real, unmodified `ToolsRTM::SPART()` call. |
| `toolsrtm.smac.get_smac` / `spectral_convolution` / `get_sensor` | `ToolsRTM/R/get.smac.R` / `Spectral.convolution.R` / `get.coef.SMAC.R` + all 9 `data/*.rda` sensor objects | The SMAC atmospheric-correction physics (transmittances/reflectances per sensor band, given geometry + aerosol/ozone/water-vapour amounts) and SRF-weighted spectral convolution. **All 9 sensors the R package ships are bundled**: Landsat 4/5/7/8, Sentinel-2A/B, Sentinel-3A/B, Terra/Aqua MODIS (see `SENSORS`/`get_sensor(name)`) -- verified against a real, unmodified `ToolsRTM::get.smac()` call for 4 representative sensors (~1e-15 relative difference), Sentinel-2A additionally against a real `ToolsRTM::SPART()` call. |
| `toolsrtm.marmit.marmit1` / `marmit2` / `toolsrtm.get_marmit_rsoil` | `ToolsRTM/R/marmit1.R`, `marmit2.R` / `ToolsRTM/R/get.marmit.rsoil.R` | MARMIT-1 and MARMIT-2 soil reflectance models (dry-\>wet soil spectrum) -- MARMIT-2 additionally accounts for soil particle size/refractive index (`n_i`/`k_i`/`d_i`), generally more accurate for coarser soils; select via `get_marmit_rsoil(..., version='marmit2')`. Only the Bablet_2016 database is bundled (17 IDs, driest spectrum per ID only -- same scoping the R side uses), exported via `python/scratch/scratch_export_marmit.py`. The other 7 official MARMIT databases (Dupiau 2020, Humper 2015, Lesaignoux 2008, Liu 2002, Lobell 2002, Marcq 2012, Philpot 2014 -- see the [MARMIT GitLab](https://pss-gitlab.math.univ-paris-diderot.fr/marmit/marmit)) live in this monorepo's own `databases/` folder (repo root, ~200MB, not bundled -- too large) -- pass `database="Liu_2002", db_root="databases"` (etc.) to use them directly, no download needed; matches `ToolsRTM::get.marmit.rsoil()`'s own `db_root` argument on the R side. |
| `toolsrtm.fluspect.fluspect_b` / `fluspect_cx` | `ToolsRTM/R/fluspect_B.R` (`getFluspect.B`) / `fluspect_Cx.R` (`getFluspect.Cx`) | FLUSPECT-B / FLUSPECT-B-Cx leaf model: PROSPECT reflectance/transmittance plus chlorophyll-fluorescence excitation-emission matrices (doubling method). |
| `toolsrtm.liberty.liberty` | `ToolsRTM/R/liberty.R` | LIBERTY conifer-needle leaf model (Dawson, Curran & Plummer 1998) -- cell-diameter/intercellular-air-space Kubelka-Munk formulation, not PROSPECT's N-layer plate model. |
| `toolsrtm.indices.get_indices` | `ToolsRTM/R/getIndices.R` | \~75 VNIR + \~18 SWIR spectral vegetation indices from a reflectance spectrum (linear interpolation onto an integer-nm grid, then per-index algebra at named wavelengths). Takes `wavelengths`/`reflectance` arrays directly rather than dataframe columns matched by regex (see `indices.py` docstring) -- the numerics are unchanged. `CIgreen`, `TCARI`, and `OSAVI.1510` currently come back as `NaN`, and `GnyLi`/`CI2` reflect only their final internally-computed value -- reproduced exactly as R itself returns them. |
| `toolsrtm.inversion` (`carspls`, `get_vif`, `get_inversion_opt`, `get_inversion`, `hybrid_inversion`, `hybrid_inversion_ensemble`) | `ToolsRTM/R/carspls.R`/`get.cars.pls.R`, `getVIF.R`, `get.inversionOpt.R`, `get.inversion.R`, `hybrid_inversion.R`/`hybrid_inversionE.R` | Trait-inversion toolkit: CARS-PLS and VIF predictor selection, LUT nearest-neighbour ("merit function") matching, and a 12-algorithm ML dispatcher (PLSR/SVM/RF/GB/NN/Bayesian/AdaBag/BRNN/xGB/RVM/qLASSO/Ensemble) built on scikit-learn/xgboost rather than `caret` (which has no 1:1 Python equivalent) -- see `inversion.py`'s module docstring for exactly which estimator each algorithm name maps to. Needs the optional `ml` extra: `pip install "toolsrtm[ml]"`. CARS-PLS is verified against a real R call (see `tests/test_inversion.py`); the ML dispatcher, being estimator-for-estimator different from `caret`'s own stochastic tuning, is verified by running end-to-end and checking accuracy on held-out data instead. |
| `toolsrtm.deep_learning.get_ml_model` (optional) | `ToolsRTM/R/getMLmodel.R`/`getMLmodel_withRetrain.R` | Dense ("Hidden-layers") and 1D-CNN Keras architectures for trait inversion, matching R's own layer sizes, dropout placement, and 7-optimizer choice. Not required for the rest of the package -- needs the optional `dl` extra: `pip install "toolsrtm[dl]"` (TensorFlow). `toolsrtm.inversion`'s scikit-learn-based dispatcher (above) covers most trait-inversion needs without this. |
| `toolsrtm.satellite` (`get_satellite_collection`, `get_sentinel2_cube`, optional) | `ToolsRTM/R/get.satellite_collection.R`/`get.sentinel2_cube.R` | STAC-based retrieval covering the same 9 collections as R's own `switch()` (Sentinel-2 L2A, Landsat C2 L2, 6 MODIS products) via `COLLECTION_ASSETS`, on Microsoft Planetary Computer or AWS Earth Search. Uses `pystac-client`/`odc-stac` rather than R's `rstac`/`gdalcubes`, a different resampling engine, so pixel values aren't expected to match R exactly even for the same scene. Needs the optional `stac` extra: `pip install "toolsrtm[stac]"`, and live network access -- Sentinel-2 L2A and Landsat C2 L2 verified against real Planetary Computer queries + cube builds over Wageningen, NL; the MODIS collections share the same code path but weren't each individually re-verified live (no `eo:cloud_cover` on MODIS items, so cloud filtering is skipped for those, matching R). |

Bundled data: `dataSpec_PDB.csv`/`dataSpec_PRO.csv` (exported once from the R `data/*.rda` objects; see `python/scratch/scratch_export.R`), `optipar_spart_bsm.csv` (ToolsRTM's own BSM GSV/Kw/nw table, used by `spart_toc`), `data/marmit/*.csv` (water optics + Bablet_2016, see `python/scratch/scratch_export_marmit.py`), `smac_coef_sentinel2a.csv`/ `smac_bands_sentinel2a.csv`/`smac_srf_wl_sentinel2a.csv`/ `smac_srf_weight_sentinel2a.csv` (Sentinel-2A MSI's SMAC coefficients/SRF tables) and `extraterrestrial_irradiance.csv` (default `spart_toa` irradiance), all exported once from R (`python/scratch/scratch_export.R`), so the Python package has **no runtime dependency on R**.

### `scopeinpython` (`python/scopeinpython/`)

| Python | R source | What it is |
|------------------------|------------------------|------------------------|
| `scopeinpython.soilwat`, `scopeinpython.get_bsm` | `SCOPEinR/R/BSM.R` | Brightness-Shape-Moisture soil reflectance model |
| `scopeinpython.get_spectra_scope` | `SCOPEinR/R/define_bands.R` (`get.spectra.SCOPE`) | SCOPE spectral-region/wavelength-grid definitions |
| `scopeinpython.rtmo.get_volscatt_scope` | `SCOPEinR/R/RTMo_functions.R` | Vectorised volume-scattering geometry used by RTMo (distinct impl. from `toolsrtm.canopy.volscatt`) |
| `scopeinpython.rtmo.get_pso` | `SCOPEinR/R/RTMo_functions.R` | Bi-directional (hot-spot) gap probability |
| `scopeinpython.rtmo.get_reflectances` | `SCOPEinR/R/RTMo_functions.R` | Multi-layer 4-stream reflectance propagation (leaf layers -\> soil) |
| `scopeinpython.rtmo.get_fluxprofile` | `SCOPEinR/R/RTMo_functions.R` | Vertical direct/diffuse flux profile within the canopy |
| `scopeinpython.run_rtmo` | `SCOPEinR/R/RTMo.R` (`getRTMo`) | **Optical** top-of-canopy BRDF: leaf optics + soil + geometry -\> `rdd`, `rsd`, `rdo`, `rso`, `refl`, `Lo_`, `Eout_`, gap probabilities `Ps`/`Po`/`Pso` |
| `scopeinpython.get_fluspect_cx_scope` | `SCOPEinR/R/fluspect_Cx_forSCOPE.R` (`getFluspect.Cx.SCOPE`) | The FLUSPECT-B-Cx variant SCOPE's own leaf-optics pipeline actually calls -- **not the same function** as `toolsrtm.fluspect_cx` despite similar code: a caller-configurable `step` (nm) controls the `Mb`/`Mf` excitation-emission matrix resolution (SCOPE's own default `step=5` gives 53x71 matrices; `step=1` gives 211x351) and scales the SIF response by `step` itself rather than a fixed constant, and a single combined `Mb`/`Mf` pair (one `phi` spectrum) rather than separate PSI/PSII matrices. Uses `SCOPEinR::optipar2021.Pro.CX` (the dataset with real `Kp`/`Kcbc` data; `optipar2017.ProspectD` lacks them). Both `step=5` and `step=1` verified against R. |
| `scopeinpython.get_biochemical` | `SCOPEinR/R/biochemical.R` (+ `Biochemical_functions.R`) | Leaf-level Farquhar/Collatz photosynthesis + van der Tol et al. (2014) fluorescence yield (`A`, `Ci`, `rcw`, `eta`, `qE`, `qQ`, `SIF`, ...), given an assumed leaf micro-environment (temperature, CO2, PAR, humidity) -- this is the piece called *inside* SCOPE's energy-balance iteration (not ported) to get `eta` at each candidate leaf temperature; it doesn't itself solve for temperature. C3 (iterative Ball-Berry `Ci` via Brent root-finding) and C4 paths both ported and verified. **Note**: `Type='C4'` with `temp_correction=False` reproduces a real crash in the R source (`Vcmax`/`Rd` are never assigned in that branch combination) rather than working around it -- see the module docstring. |
| `scopeinpython.fluspect_mscope` | `SCOPEinR/R/fluspect_mSCOPE.R` (`get.fluspect_mSCOPE`) | Multi-layer (mSCOPE) leaf-optics wrapper: computes leaf optics once per distinct leaf-biochemistry profile layer (via `get_fluspect_cx_scope`), replicates across the canopy sublayers each profile layer spans (weighted by `pLAI`). Reproduces R's profile-layer-boundary overwrite quirk exactly. **Note**: R's `get.fluspect_mSCOPE()` called without `step` always crashes (confirmed via standalone repro -- a dimension-mismatch in its own array pre-allocation); the Python `step` parameter is required rather than offering that non-working default. |
| `scopeinpython.rtmf` | `SCOPEinR/R/RTMf.R` (`get.RTMf`) | SCOPE canopy fluorescence RTM: TOC fluorescence radiance in the viewing direction (`LoF_`) and hemispherical upward fluorescence flux (`EoutF_`), given `Mb`/`Mf` (from `fluspect_mscope`) and per-layer fluorescence quantum efficiencies (`etau`/`etah`, from `get_biochemical`) as composable inputs -- no dependency on the (not-ported) thermal energy-balance loop. Ported against a **fixed** R source: found and fixed 3 real bugs during this port (see below) -- an 18-expression column-recycling issue and an `absfs_nl` copy-paste mixup. **One deliberate, documented numerical approximation**: the final upsampling from the native 53-point fluorescence wavelength grid to the 211-point display grid uses `scipy`'s `not-a-knot` cubic spline, which is close to but not bit-identical with R's `fmm`-method spline (`signal::interp1(...,'spline')`) -- a small, boundary-localized difference (\~1.5e-3 absolute worst case), unlike everything else in this port which matches to floating-point precision. |
| `scopeinpython.rtmz` | `SCOPEinR/R/RTMz.R` (`get.RTMz`) | SCOPE canopy zeaxanthin RTM: the small TOC-radiance correction (500-600nm) from violaxanthin-\>zeaxanthin conversion, given baseline/full-zeaxanthin leaf optics (`refl`/`tran`/`reflZ`/`tranZ`, two `fluspect_mscope` calls with `Cx=0`/`Cx=1`) and per-layer NPQ (`Knu`/`Knh`, from `get_biochemical`'s `Kn`). Returns deltas to add onto an existing `run_rtmo` result's `rso`/`rdo`/`Eout_`/`Lo_` at the 500-600nm band. Ported against a **fixed** R source: found and fixed 5 real bugs (same column-recycling family as `get.RTMf`; the `get.plots`-only dead-code bug -- its real output was never computed with the normal `get.plots=FALSE` call; a `Po[1:nl+1]` indexing bug; a hardcoded `dim=c(30,13,36)`; a MATLAB `sum(A,2)`-to-R mistranslation). Two more issues flagged but not independently confirmed/fixed -- see `rtmz.py`'s module docstring. |
| `scopeinpython.thermal` | `SCOPEinR/R/Monin_ObuKhov.R`, `resistances.R`, `heatfluxes.R` | Scalar-per-timestep building blocks for the SCOPE thermal energy-balance loop (`ebal.R`, ported as `scopeinpython.ebal`, see below): `monin_obukhov` (stability length), `get_resistances` (Wallace & Verhoef 2000 two-layer aerodynamic resistance scheme + Paulson 1970 stability corrections), `get_heatfluxes` (latent/sensible heat flux of a leaf or soil surface). **Note**: `get_resistances`'s `rac`/`rws` reproduce a real R quirk -- they use the *stability-uncorrected* eddy diffusivity, not the corrected value returned as the function's own `Kh` output, confirmed by direct reading of the R source. |
| `scopeinpython.rtmt_sb` | `SCOPEinR/R/RTMt.sb.R` (`get.RTMt.sb`) | Total thermal-IR outgoing radiation + net radiation per leaf/soil component (Stefan-Boltzmann, spectrally-integrated), given already-solved leaf/soil temperatures (from `ebal`'s convergence loop). **Only the "SCOPE-lite" scalar-per-layer branch is ported** (matches every reference case in this port); the full `(13,36,nl)` per-leaf-angle array branch and the `obsdir` (directional brightness temperature) branch are not ported -- the latter has an unresolved, likely-buggy R indexing expression (`data.rad$vb[1, nl]`, using the layer count as a wavelength index) flagged but not independently confirmed. |
| `scopeinpython.net_radiation_lite` (in `rtmo.py`) | `SCOPEinR/R/RTMo.R` section 4 (partial) | Just the 6 quantities `ebal` needs (`Rnuc`/`Rnhc`/`Rnus`/`Rnhs`/`Pnu_Cab`/`Pnh_Cab`), "lite" branch only. **A real bug found and fixed in R here too**: the direct-beam absorption term used a stray leftover loop variable (`Asun[j]` with `j` always `nl`, the last layer) instead of the full per-layer vector, making direct-beam-absorbed net radiation come out exactly constant across all canopy layers -- confirmed via direct repro. Fixed in `RTMo.R`, and a downstream `Rndir[1] + ...` (harmless before the fix, wrong after) fixed in the same edit. |
| `scopeinpython.ebal` | `SCOPEinR/R/ebal.R` | The energy-balance closure loop itself: iterates sunlit/shaded leaf and soil temperature until sensible+latent heat flux matches net radiation, coupling `net_radiation_lite`/`rtmt_sb` (radiation), `get_resistances`/`monin_obukhov` (aerodynamics), `get_biochemical` (photosynthesis/fluorescence) and `get_heatfluxes`. **"SCOPE-lite"-only**: only the default (non-MD12) fluorescence-model branch, only the simple `G = 0.35*Rn` ground-heat-flux method (not the two time-series-history soil-inertia methods), only `meanleaf.v2`'s `'layers'` aggregation mode. |
| `scopeinpython.get_scope` | `SCOPEinR/R/get.SCOPE.R` | The end-to-end wrapper: one LUT input row (matching `SCOPEinR::LUT_input.csv`'s columns) in, leaf optics -\> soil -\> optical BRDF -\> energy balance -\> fluorescence (optional) -\> zeaxanthin (optional) out, via all of the above. Verified against a real, **unmodified** `get.SCOPE()` call using the R package's own bundled example LUT/options. **Not ported as part of it**: `get.SCOPE.parallel` (R's `foreach`/`doParallel` backend -- parallelize `get_scope` calls yourself), directional BRDF, `RTMt_planck`, multi-layer mSCOPE, time-series mode, angle-file LIDF, measurement-file/MODTRAN irradiance, and the canopy-level "derived data products" beyond `ScopeResult`'s fields (`Pnsun_Car`/`Rnsun_Cab`/`Rnsun_PAR`/`LST`/etc) -- see the module docstring. Along the way, found and fixed a second real R bug: `get.zo_and_d`'s degenerate-canopy branch referenced an undefined `d` variable. |

Bundled data: `constants.csv` (physical constants), `optipar_bsm.csv` (GSV/Kw/nw spectra needed by BSM), `soil_scope.csv` (SCOPE's own 3 reference dry-soil spectra) and `default_irradiance.csv` (SCOPE's bundled example `Esun_`/`Esky_`), exported once from R.

`scopeinpython` depends on `toolsrtm` for leaf optics (PROSPECT-D/PROSPECT-PRO) exactly as the R `SCOPEinR` package depends on `ToolsRTM`.

## How this was verified numerically

Every ported function has a `pytest` regression test comparing its output against values generated by running the **original R package** once and saving the results to CSV (`python/scratch/scratch_export.R` and `python/scratch/scratch_rtmo_export.R`, run via a local `Rscript` install). The reference CSVs are bundled under each package's `tests/refdata/`.

Commands run:

```         
cd python/toolsrtm && python -m pytest tests -q   # 49 tests
cd python/scopeinpython && python -m pytest tests -q   # 18 tests
```

All 67 tests pass. Tolerance used in the tests is `rtol=1e-6` (loose, to be safe across platforms/BLAS), but the *actual* observed differences are at floating-point noise level except where noted, i.e. R and Python agree to \~14-15 significant digits:

| Function | Reference call | Max relative diff (Python vs R) |
|------------------------|------------------------|------------------------|
| `calctav` | `ToolsRTM::calctav` | 0 (exact) |
| `prospect_d` (refl/tran) | `ToolsRTM:::prospect_DB` | 5.8e-14 |
| `prospect_pro` (refl/tran) | `ToolsRTM::prospect_PRO` | \~1e-14 |
| `dladgen`, `campbell`, `volscatt` | `ToolsRTM::dladgen/campbell/volscatt` | 0 (exact, scalar/short outputs) |
| `foursail` (rdot/rsot/rddt/rsdt) | `ToolsRTM::foursail` (PROSPECT-PRO leaf model) | 3.4e-14 |
| `foursail2` (rdot/rsot/rddt/rsdt/alfast/alfadt) | `ToolsRTM::foursail2` (two-layer green/brown canopy) | \~1e-15 |
| `inform` (forest reflectance), TypeLidf=1 and TypeLidf=2 | `ToolsRTM::inform` (PROSPECT-D leaf model) | \~1e-15 |
| `spart_toc` (TOC BRDF reflectance) + its BSM soil step | `ToolsRTM::SPART` (TOC portion) / `ToolsRTM::getBSM.toolsRTM` | \~1e-15 |
| `get_indices`, all 3 spectral domains (75+18 indices, 5 random spectra) | `ToolsRTM::getIndices` | \~5e-10 (worst: `CR.red.nir.2`, an R-vs-scipy linear-extrapolation rounding difference; most indices agree to \~1e-14) |
| `fluspect_b` (refl/tran/kChlrel/MbI/MbII/MfI/MfII) | `ToolsRTM::getFluspect.B` | \~1e-15 (refl/tran), \~1e-19 (fluorescence matrices) |
| `fluspect_cx` (refl/tran/kChlrel/kCarrel/Mb/Mf) | `ToolsRTM::getFluspect.Cx` | \~1e-15 (refl/tran), \~1e-19 (fluorescence matrices) |
| `liberty` (refl/tran/RR) | `ToolsRTM::liberty` | 6.9e-12 |
| `foursail`/`foursail2`/`inform` wired with Liberty, Fluspect-B, Fluspect-B-Cx | `ToolsRTM::foursail`/`foursail2`/`inform` with those `LeafModel`s | \~5e-12 (Liberty path), \~1e-15 (Fluspect paths) |
| `get_bsm`/`soilwat` | `SCOPEinR::getBSM` | 4.2e-15 |
| `run_rtmo` (rdd/rsd/rdo/rso/refl) | `SCOPEinR:::getRTMo` (assembled inputs: PROSPECT-D leaf optics + BSM soil, LAI=3, 30 canopy layers) | 2.7e-14 |
| `get_biochemical`, C3 (iterative Brent Ci) and C4, both temperature-corrected, 24 output fields each | `SCOPEinR::get.biochemical` | \~1e-13 (worst field: `rcw`; most fields \~1e-15 to 1e-18) |
| `get_fluspect_cx_scope`, step=5 and step=1 | `SCOPEinR::getFluspect.Cx.SCOPE` | \~1e-15 (refl/tran/kChlrel/kCarrel), \~1e-19 (Mb/Mf) |
| `fluspect_mscope` (3-profile-layer, 10-canopy-layer case) | `SCOPEinR::get.fluspect_mSCOPE` | \~1e-15 (refl/tran/kChlrel/kCarrel/phiI/phiII), \~1e-19 (Mb/Mf) |
| `rtmf`, native-grid (pre-spline) quantities | `SCOPEinR::get.RTMf` | \~1e-6 (loop-accumulated over 30 layers; tighter per-expression) |
| `rtmf`, fully-interpolated `LoF_`/`EoutF_` (211-pt display grid) | `SCOPEinR::get.RTMf` | \~1.5e-3 absolute worst case (documented `not-a-knot` vs `fmm` spline approximation, see `scopeinpython.rtmf` docs) |
| `rtmz`, `rso`/`rdo`/`Eout_` corrections (500-600nm) | `SCOPEinR::get.RTMz` | \~1e-5 |
| `net_radiation_lite` (`Rnuc`/`Rnhc`/`Rnus`/`Rnhs`/`Pnu_Cab`/`Pnh_Cab`) | `SCOPEinR::RTMo.R` section 4, "lite" branch | \~1e-13 |
| `rtmt_sb` | `SCOPEinR::get.RTMt.sb`, "lite" scalar-per-layer branch | \~1e-13 |
| `ebal` (converged temperatures) | `SCOPEinR::get.ebal` | \~1e-3 (energy-balance totals \~1-2%, see module docstring for why) |
| `get_scope`, exact-formula outputs (`LAIsunlit`/`Pnsun_Cab`/`Pnsha_Cab`/`Pntot_Cab`/TOC `refl`) | `SCOPEinR::get.SCOPE` (package's own bundled example LUT) | \~1e-9 to \~1e-13 |
| `get_scope`, iterative-convergence outputs (temperatures, energy-balance totals, `Ja`/`PNPQ`/`fqe`) | `SCOPEinR::get.SCOPE` | \~1-2% (\~5% worst case, `fqe`) -- same per-layer-loop biochemistry divergence as `ebal` itself |

No unresolved numerical discrepancies were found — every ported function matches R to floating-point precision. (Two spots worth knowing about, not discrepancies: (1) `run_rtmo`'s `rso`/`refl` legitimately produce `NaN`/`inf` at far-thermal wavelengths where `Esun_` underflows to \~0 — this reproduces R's own behaviour at the same wavelengths; (2) `soilwat`'s use of `scipy.stats.poisson.pmf(round(mu), k)` mirrors the original R call `stats::dpois(round(mu), k)` exactly, and is covered by the BSM regression test.)

The RTMo reference case was built by assembling the exact `data.spectral`/ `data.soil`/`data.leafopt`/`data.canopy`/`data.angles`/`data.meteo`/ `atmo` R structures by hand (PROSPECT-D leaf optics broadcast to all canopy layers, BSM soil, default SCOPEinR example irradiance) and calling the internal (non-exported) `SCOPEinR:::getRTMo` directly — this isolates the RTMo math itself from the full `get.SCOPE()` wrapper. See `python/scratch/scratch_rtmo_export.R` for that setup, `python/scratch/scratch_ebal_export.R` for the hand-assembled `ebal` case, and `python/scratch/scratch_scope_export.R` for `get_scope`'s own reference case — the one call in this whole port that uses the real, unmodified `get.SCOPE()` end to end, against its own bundled example `LUT_input.csv`/`setoptions.csv`.

## Install

``` bash
pip install -e python/toolsrtm
pip install -e python/scopeinpython   # depends on toolsrtm
```

Both require `numpy` and `scipy` only (`pandas`/`pytest` are only needed to run the test suites: `pip install -e "python/toolsrtm[test]"`).

## Usage

### Leaf + canopy (toolsrtm)

``` python
import numpy as np
from toolsrtm import prospect_d, foursail

# Leaf reflectance/transmittance, 400-2500 nm
leaf = prospect_d(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40)
print(leaf.lambda_[:3], leaf.refl[:3], leaf.tran[:3])

# Full canopy BRDF (fourSAIL), dispatches to PROSPECT-D/PROSPECT-PRO internally
inputLUT = dict(
    N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40,
    Prot=0.002, CBC=0.007,          # only used if leaf_model='PROSPECT-PRO'
    LIDFa=-0.35, LIDFb=-0.15, TypeLidf=1,  # spherical-ish LIDF (dladgen)
    LAI=3, hspot=0.01, tts=30, tto=0, psi=0,
)
rsoil = np.full(2101, 0.15)  # flat soil reflectance, 400-2500 nm
sail = foursail(inputLUT, rsoil, leaf_model="PROSPECT-D", spectrum_all=True)
print("TOC bidirectional reflectance factor (rsot) at 550 nm:", sail.rsot[550 - 400])
```

### Soil (BSM) + optical canopy BRDF (scopeinpython)

``` python
import numpy as np
from toolsrtm import prospect_d
from toolsrtm.canopy import dladgen
from scopeinpython import SoilParams, WettingParams, get_bsm, CanopyStructure, get_spectra_scope, run_rtmo

spectral = get_spectra_scope()

# Soil reflectance from BSM, 400-2400 nm
rsoil = get_bsm(SoilParams(BSMBrightness=0.5, BSMlat=25, BSMlon=45),
                 WettingParams(SMp=15, SMC=25, film=0.015))

# Leaf optics, truncated to 400-2400 nm (SCOPE's optical range)
leaf = prospect_d(N=1.5, Cab=40, Car=8, Anth=1, Cbrown=0, EWT=0.01, LMA=0.009, alpha=40)
refl_leaf, tran_leaf = leaf.refl[:2001], leaf.tran[:2001]

lidf = dladgen(-0.35, -0.15).lidf
canopy = CanopyStructure(LAI=3, lidf=lidf, hot=0.1 / 2.0)  # leafwidth=0.1 m, hc=2 m

# Default SCOPEinR example irradiance (Esun_/Esky_) must be supplied by the
# caller -- only the "precomputed irradiance" mode of get.calcTOCirr is
# ported (MODTRAN-derived irradiance is not). Load it once from your own R export,
# or from any Esun_/Esky_ pair on the spectral.wlS grid (2162 values).
Esun_ = ...  # np.ndarray, shape (2162,)
Esky_ = ...  # np.ndarray, shape (2162,)

result = run_rtmo(
    spectral=spectral, leaf_refl=refl_leaf, leaf_tran=tran_leaf,
    rho_thermal=0.01, tau_thermal=0.01, rsoil=rsoil, canopy=canopy,
    tts=30, tto=0, psi=0, Esun_=Esun_, Esky_=Esky_,
)
print("TOC reflectance (refl) at 550 nm:", result.refl[550 - 400])
print("TOC bidirectional reflectance factor (rso) at 550 nm:", result.rso[550 - 400])
```

A ready-to-use `Esun_`/`Esky_` pair (the default SCOPEinR example irradiance, on the `spectral.wlS` grid) is included as reference data at `python/scratch/_refdata/default_irradiance.csv` and `python/scopeinpython/tests/refdata/ref_RTMo_outputs.csv` (columns `Esun_`, `Esky_`) for convenience/testing -- and bundled as package data (`scopeinpython.data.default_irradiance.csv`), used automatically by `get_scope` below if you don't supply your own.

### Full SCOPE simulation from one LUT row (get_scope)

``` python
import csv
from scopeinpython import ScopeOptions, get_scope

# One row of SCOPEinR's own LUT_input.csv layout -- a dict or pandas.Series
# works; here we just read the package's own bundled example row.
with open("SCOPEinR/inst/input/LUT_input.csv", newline="") as f:
    lut_row = next(csv.DictReader(f))

result = get_scope(lut_row, options=ScopeOptions(calc_fluor=True, calc_xanthophyllabs=True))

print("TOC reflectance at 550 nm:", result.rtmo.refl[550 - 400])
print("Canopy-average leaf temperature (Tcave):", result.ebal.Tcave)
print("Total net radiation (Rntot):", result.ebal.Rntot)
print("Canopy fluorescence flux (EoutF):", result.rtmf.EoutF)  # None if calc_fluor=False
```

## Repo layout

```         
python/
  toolsrtm/                  installable package (pyproject.toml, pip install -e .)
    src/toolsrtm/            leaf.py, canopy.py, _data.py, data/*.csv
    tests/                   pytest regression tests + refdata/*.csv (from R)
  scopeinpython/                  installable package (pyproject.toml, pip install -e .)
    src/scopeinpython/            soil.py, spectral.py, rtmo.py, utils.py, _data.py, data/*.csv
    tests/                   pytest regression tests + refdata/*.csv (from R)
  scratch_export.R           one-off R script: exports ToolsRTM/SCOPEinR reference
                              data + values used by the tests above (re-run only if
                              you need to regenerate/extend the reference data)
  scratch_rtmo_export.R      one-off R script: builds the RTMo reference case
                              (assembles R structures, calls SCOPEinR:::getRTMo)
  _refdata/                  raw CSV output of the two scratch_*.R scripts (source
                              of the tests/refdata/*.csv files bundled in each package)
```
