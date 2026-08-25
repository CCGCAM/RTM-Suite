# Pipeline

End-to-end simulate -> convolve -> index -> invert workflows, built on the
current `ToolsRTM`/`SCOPEinR` packages. Two tracks, same shape:

| | ToolsRTM (optical only) | SCOPEinR (full energy balance) |
|---|---|---|
| 1. Simulate + convolve + index | `1-simulate_LUT.R` | `SCOPE-1-simulate.R` |
| 2. Classic ML inversion (12 algorithms) | `2-inversion_ML.R` | `SCOPE-2-inversion_ML.R` |
| 3. Deep learning inversion (CNN / hidden-layers) | `3-inversion_deep_learning.R` | `SCOPE-3-inversion_deep_learning.R` |

Each step reads the `.rds` the previous step saved under `outs/rtm_sims/pipeline*`
-- run them in order, 1 then (2 and/or 3).

Every script starts with `use.dev.source <- TRUE/FALSE` (load the local
`ToolsRTM/R`, `SCOPEinR/R` source vs. the installed package) and writes all
outputs under `outs/rtm_sims/pipeline` (ToolsRTM) or
`outs/rtm_sims/pipeline_scope` (SCOPEinR) -- never inside `Scripts/`.

## One script, not one script per model

`1-simulate_LUT.R` is not fourSAIL-only -- it already covers all three
ToolsRTM canopy models and all three convolution sensors. There's no
`1-simulate_LUT_INFORM.R` or `1-simulate_LUT_foursail2.R` because there's
nothing those would do differently; the model/sensor/trait choice is a
variable at the top of the file, not a separate codepath:

```r
canopy.model <- "fourSAIL"     # or "foursail2" (two-layer canopy) or "INFORM" (forest, explicit crown geometry)
sensors      <- c("Sentinel2a", "Sentinel2b", "PRISMA")   # convolved every run
```

`1-datasets.rds` comes back with all three convolutions ready to invert from:
`datasets$native` (full simulated spectrum), `datasets$se2a`, `datasets$prisma`.
Steps 2/3 pick one via `dataset <- datasets$native` (swap for `$se2a`/`$prisma`).

Same story for trait choice -- `target_trait <- "Cab"` in steps 2/3 works for
any column of the LUT (`"Car"`, `"Anth"`, `"LAI"`, `"CBC"`, `"EWT"`, ...), and
for SCOPE also biochemistry parameters like `"Vcmax25"` (verified: mechanism
identical to Cab/LAI, though Vcmax's spectral signal is weaker/more indirect
so it needs more training samples for a comparable R2).

get.SCOPE() itself always runs its own integral Fluspect-Cx + RTMo (its
`leaf.model`/`canopy.model` arguments aren't functional -- see
`Scripts/compare_SCOPE_models.R`), so there's one leaf/canopy configuration
on the SCOPE side, not a model choice to expose here.

## Other traits/correlations (ToolsRTM track only)

`1-simulate_LUT.R` also exposes, at the top of the file:

```r
trait.distribution <- list(Cab = "Gaussian")                        # override a trait's sampling distribution
correlate.traits   <- list(Car = list(with = "Cab", scale = 1/4, r = 0.8))  # make one trait co-vary with another
```

## 0-*.R helper scripts

- `0-setup_python_env.R` -- run once before step 3 (deep learning): provisions
  the Python/TensorFlow environment `reticulate`/`keras` need.
- `0-integrate_MARMIT_soil.R`, `0-evaluate_MARMIT.R` -- MARMIT soil model
  exploration, independent of the simulate/invert steps above.

## Related: Scripts/ForPROSAIL

`Scripts/ForPROSAIL/` is a separate, older PROSAIL-specific pipeline (its own
LUT reduction against real field spectroradiometer data, its own deep-learning
model scripts) -- see its own `data/README.md`. Prefer this `Pipeline/` folder
for new work; `ForPROSAIL/` exists because it's tied to a specific field
campaign's private data and workflow.
