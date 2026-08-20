# ForSCOPE_Python

Python equivalent of `Scripts/R/ForSCOPE/` (the R course pipeline for the
full SCOPE model), built entirely on `scopeinpython.get_scope` -- the
end-to-end wrapper (leaf optics -> soil -> optical BRDF -> energy balance
-> fluorescence -> zeaxanthin) added to the Python port.

## Setup

```bash
pip install -e python/toolsrtm
pip install -e python/scopeinpython
pip install pandas matplotlib
```

## Scripts

| Step | Script | What it does |
|---|---|---|
| 1 | `1_run_scope.py` | Runs one full `get_scope()` simulation on SCOPEinR's own bundled example LUT row, prints the key outputs (reflectance, temperatures, fluxes, fluorescence), saves spectral and scalar CSVs. |
| 2 | `2_explore_outputs.py` | Plots the TOC reflectance components, the fluorescence spectrum (640-850nm, showing the classic 685/740nm double peak), the converged leaf/soil temperature profile, and the energy-balance totals. |
| 3 | `3_simulate_lut.py` | Runs `get_scope()` 30 times with Cab/LAI/EWT perturbed around the example row, records canopy reflectance/temperature/photosynthesis/fluorescence per sample, reports timing. |

Run in order:

```bash
cd Scripts/Python/ForSCOPE_Python
python 1_run_scope.py
python 2_explore_outputs.py
python 3_simulate_lut.py
```

Real results from the last verified run (SCOPEinR's bundled example row,
LAI=3, Cab=40, C3 vegetation): 30 canopy layers, TOC reflectance 0.044 at
550nm / 0.36 at 800nm, canopy-average leaf temperature 22.0°C (soil
23.8-39.9°C sunlit/shaded), net radiation 495 W/m², latent/sensible heat
181/275 W/m², canopy photosynthesis 19.9 µmol m⁻² s⁻¹, converges in 7
iterations, canopy fluorescence flux 0.39 W/m² with the expected
double-peaked spectrum (F685=0.31, F740=1.81).

**Timing**: `get_scope()` takes ~0.1-0.25s per call (dominated by the
`ebal` nonlinear leaf/soil-temperature solve), vs. ~2ms for
`toolsrtm.foursail` -- so `3_simulate_lut.py` uses 30 samples, not the
25,000-sample LUTs a real training run would use (that scale is a
minutes-not-seconds job, left for a real training pipeline rather than
this demo).

## What's NOT here (and why)

- **Trait inversion (ML/DL) from the simulated LUT** — the technique is
  identical to `ForPROSAIL_fourSAIL/3_inversion_ml.py`/`4_inversion_dl.py`
  (same scikit-learn/TensorFlow/PyTorch models, same metrics), just applied
  to `get_scope()`'s outputs instead of `foursail`'s -- not duplicated here
  to avoid repeating the same code twice; see the AEO-Course material for a
  worked example that goes through SCOPE (see the repo's course-template
  work).
- **`get.SCOPE.parallel` / batch LUT runs at scale** — `get_scope` is the
  composable per-row building block (see `python/README.md`); parallelize
  calls to it yourself (`multiprocessing`, `joblib`, ...) for a large LUT.
- **`6-validate_ebal_convergence.R`'s convergence-diagnostics deep dive** —
  `get_scope()`'s `EbalResult` already exposes `counter`/`maxEBercu`/
  `maxEBerch`/`maxEBers`; a dedicated convergence-analysis script wasn't
  built here, only used inline in `1_run_scope.py`'s printed summary.
