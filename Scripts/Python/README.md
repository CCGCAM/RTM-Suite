# Scripts/Python

Python demonstration of the `toolsrtm`/`scopeinpython` Python port (`python/`
at the repo root) -- kept entirely separate from the R scripts in the rest
of `Scripts/`. Nothing here reads or writes anything under `ToolsRTM/`,
`SCOPEinR/`, or the R `Scripts/*` folders; outputs go to `outs/Python/`,
not `outs/<RFolderName>/`.

## Setup

The Python port must be installed (editable) first:

```bash
pip install -e python/toolsrtm
pip install -e python/scopeinpython          # depends on toolsrtm
pip install scikit-learn pandas matplotlib   # classical ML + plotting
pip install tensorflow torch xgboost         # deep-learning inversion (4_inversion_dl.py)
```

## ForPROSAIL_fourSAIL/

Python equivalent of the R course pipeline in `Scripts/R/ForPROSAIL/`
(simulate -> spectral indices -> ML inversion, classical **and** deep
learning), built entirely on the Python port:

| Step | Script | What it does |
|---|---|---|
| 1 | `1_simulate_lut.py` | Samples 100 random PROSPECT-D + fourSAIL parameter sets, simulates 400-2500 nm TOC reflectance via `toolsrtm.foursail`, saves the LUT and spectra. |
| 2 | `2_spectral_indices.py` | Computes ~75 VNIR + ~18 SWIR spectral indices from the simulated spectra via `toolsrtm.get_indices`. |
| 3 | `3_inversion_ml.py` | Trains 6 scikit-learn regressors (PLS, RandomForest, GradientBoosting, SVR, Ridge, MLP) plus an average ensemble, per trait (Cab, LAI, EWT), on the indices; evaluates on a held-out test split; saves R²/RMSE/NRMSE/NMB/FGE metrics and predicted-vs-observed scatter plots. |
| 4 | `4_inversion_dl.py` | Trains a 3-hidden-layer dense network in **both TensorFlow/Keras and PyTorch** (same architecture, same train/test split, directly comparable), per trait; also writes a combined classical+deep-learning "best model per trait" table. |
| 5 | `5_sensor_convolution.py` | Convolves the step-1 LUT's spectra onto real sensor bands using all 3 of `toolsrtm.srf`'s convolution functions: measured SRF (Sentinel-2A/2B/PRISMA), measured SRF + SMAC atmospheric coefficients (Sentinel-2A), and Gaussian-from-nominal-characteristics (EnMAP/MODIS/... and your own sensor/camera). |

`metrics.py` holds the shared evaluation functions (R², RMSE, NRMSE, NMB,
FGE) used by both inversion scripts.

`pipeline.ipynb` runs all 4 steps inline with markdown explanations and
embedded plots -- already executed once (outputs are saved in the
notebook) so it can be read without re-running, and can be re-run with
`jupyter nbconvert --to notebook --execute --inplace pipeline.ipynb` (or
opened directly in Jupyter/VS Code).

Run the scripts in order; each reads the previous step's output from
`outs/Python/ForPROSAIL/`:

```bash
cd Scripts/Python/ForPROSAIL_fourSAIL
python 1_simulate_lut.py
python 2_spectral_indices.py
python 3_inversion_ml.py
python 4_inversion_dl.py
python 5_sensor_convolution.py
```

### Real results from the last verified run

100 samples, 70/30 train/test split, `random_state=42`/seeded throughout
for reproducibility. SVR/Ridge/MLP (both the scikit-learn one and the two
deep-learning ones) standardize predictors **and** target internally --
without that, a small-magnitude trait like EWT (~0.005-0.03) makes neural
net training diverge (confirmed empirically while building this pipeline;
see the scripts' own comments).

**Classical ML (R²):**

| Trait | PLS | RandomForest | GradientBoosting | SVR | Ridge | MLP | **Ensemble** |
|---|---|---|---|---|---|---|---|
| Cab | 0.87 | 0.83 | 0.82 | 0.77 | 0.90 | 0.87 | **0.91** |
| LAI | 0.85 | 0.69 | 0.72 | 0.69 | 0.81 | 0.74 | **0.82** |
| EWT | 0.90 | 0.82 | 0.82 | 0.88 | 0.87 | 0.92 | **0.93** |

**Deep learning (R²):**

| Trait | TensorFlow | PyTorch | DL-Ensemble |
|---|---|---|---|
| Cab | 0.79 | 0.87 | 0.84 |
| LAI | 0.75 | 0.74 | 0.78 |
| EWT | 0.83 | 0.87 | 0.87 |

The classical Ensemble is the single best model for Cab and EWT here; PLS
is best for LAI. This is expected for a 100-sample LUT: deep nets need
more data to show their usual advantage over classical methods, and both
frameworks reach comparable accuracy on the same architecture/data (as
they should -- they're two implementations of the same optimization
problem). Re-run with a larger `N_SAMPLES` in `1_simulate_lut.py` to see
the deep-learning models pull ahead.

## ForFoursail2_Python/, ForINFORM_Python/, ForSPART_Python/, ForMARMIT_Python/

Python equivalents of `Scripts/R/ForFoursail2/`, `ForINFORM/`, `ForSPART/`,
`ForMARMIT/` -- same 3-step shape as those R course pipelines (simulate +
convolve + indices -> classical ML inversion -> deep learning inversion),
each verified running end-to-end:

| Folder | Canopy/soil model | Target trait(s) | Sensors |
|---|---|---|---|
| `ForFoursail2_Python/` | PROSPECT-PRO + foursail2 (two-layer green/brown canopy) | Cab, LAI, EWT | native, Sentinel-2A, Sentinel-2B, PRISMA |
| `ForINFORM_Python/` | Fluspect-B-Cx + INFORM (forest, explicit crown geometry) | Cab, LAI, EWT | native, Sentinel-2A, Sentinel-2B, PRISMA |
| `ForSPART_Python/` | PROSPECT-PRO + fourSAIL + BSM soil + SMAC atmosphere | Cab, LAI, EWT | TOC (native) + TOA (Sentinel-2A) |
| `ForMARMIT_Python/` | Soil only, no vegetation | SMC (soil moisture, %) | native, Sentinel-2A, PRISMA |

Each folder: `1_simulate_lut.py` (100 samples -> LUT + reflectance + spectral
indices, per dataset), `2_inversion_ml.py` (6 scikit-learn regressors +
ensemble, per trait per dataset -- same approach as `ForPROSAIL_fourSAIL/
3_inversion_ml.py`), `3_inversion_dl.py` (TensorFlow + PyTorch, same
3-hidden-layer architecture as `ForPROSAIL_fourSAIL/4_inversion_dl.py`,
scoped to one dataset to keep runtime reasonable). Run in order:

```bash
cd Scripts/Python/ForFoursail2_Python   # or ForINFORM_Python / ForSPART_Python / ForMARMIT_Python
python 1_simulate_lut.py
python 2_inversion_ml.py
python 3_inversion_dl.py
```

**ForSPART_Python's atmosphere is a documented simplification**: `spart_toa()`
needs `Pa`/`aot550`/`uo3`/`uh2o` in every row with no built-in default
(unlike R's `get.smac()`, which derives `Pa` from altitude when missing) --
this pipeline uses fixed standard clear-sky values (`Pa=1013.25`,
`aot550=0.2`, `uo3=0.35`, `uh2o=2.0`), not a physical measurement or R's
own altitude-derived fallback.

## Comparison/

`compare_inversion_traditional_vs_ML.py` -- traditional LUT nearest-neighbour
inversion (`toolsrtm.inversion.get_inversion_opt`, merit-RMSE and merit-FGE)
against ML inversion (`toolsrtm.inversion.get_inversion`, Random Forest), on
the exact same held-out Sentinel-2A-convolved test set, so R²/RMSE are
directly comparable. Python equivalent of
`Scripts/R/Comparison/compare_inversion_traditional_vs_ML.R` (same LUT size,
split logic, and sensor). Needs the `ml` extra (`pip install "toolsrtm[ml]"`).

```bash
cd Scripts/Python/Comparison
python compare_inversion_traditional_vs_ML.py
```

Outputs -> `outs/Python/Comparison/` (predicted-vs-observed plot + R²/RMSE
table). Last verified run: ML (Random Forest) R²=0.76 > merit-FGE R²=0.72 >
merit-RMSE R²=0.66 for Cab retrieval -- ML wins, but the traditional
merit-function approach (no training data, no model fitting) is
surprisingly competitive, and FGE noticeably beats plain RMSE as a merit
function here.

## Where the tutorials are (R and Python, side by side)

- **`Tutorials/How-in-R.Rmd` / `How-in-Python.ipynb`** — the shorter
  on-ramp: one simulation -> 500 simulations -> sensitivity -> sensor
  convolution (all 3 convolution functions) -> ML inversion, matching the
  `Apps/RTMs` Shiny app's own tutorial tabs.
- **`Tutorials/ToolsRTM_PROSAIL_tutorial.Rmd`** — the comprehensive R
  reference manual (every leaf/canopy model, trait sampling, all 12 `caret`
  algorithms, TensorFlow/Keras deep learning) -- no direct Python
  equivalent yet (this `Scripts/Python/` folder's own scripts are the
  closest Python counterpart to its inversion sections).
- **`Tutorials/How-in-R-SCOPEinR.Rmd` / `How-in-Python-SCOPEinR.ipynb`** and
  **`Tutorials/SCOPEinR_tutorial.Rmd`** — the SCOPE equivalents of the two
  pairs above.

See the repo root [`README.md`](../../README.md)'s "Tutorials & manuals"
section for the full picture across both languages.

## What's NOT here (and why)

- **`caret`'s own cross-validated tuning/resampling plumbing**
  (`createMultiFolds`, `trainControl`, `tuneLength`, parallel workers) --
  `caret` has no 1:1 Python equivalent; the scripts here (`3_inversion_ml.py`
  etc.) use scikit-learn's native equivalents directly instead (a comparable
  spread of algorithm families: partial least squares, tree ensembles,
  kernel methods, linear, a small neural net), each fit once rather than
  grid-searched, to keep this a runnable demo rather than a multi-hour
  tuning job. The full 12-algorithm dispatcher (`ToolsRTM::get.inversion`,
  same algorithm names, single small grid search per algorithm) is ported as
  an installable package function, `toolsrtm.inversion.get_inversion` --
  see [`python/README.md`](../../../python/README.md) -- rather than
  duplicated here as another standalone script.

## For the ForSCOPE_Python/ pipeline (SCOPE energy balance + fluorescence)

See `Scripts/Python/ForSCOPE_Python/README.md` for the `get_scope()`
end-to-end demonstration (leaf optics -> soil -> optical BRDF -> energy
balance -> fluorescence -> zeaxanthin), the Python equivalent of
`Scripts/R/ForSCOPE/`.
