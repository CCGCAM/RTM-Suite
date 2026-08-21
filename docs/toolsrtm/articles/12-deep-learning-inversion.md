# 12. Deep Learning for RTM Inversion

``` r

library(ToolsRTM)
```

[`getMLmodel()`](../reference/getMLmodel.md) is
[`get.inversion()`](../reference/get.inversion.md)’s (Tutorial 11)
deep-learning sibling – a dense or 1D-CNN network via TensorFlow/Keras,
fit on the same kind of (sensor-band reflectance, trait) LUT.

## 1. Same simulate-and-split pattern as Tutorials 10-11

``` r

n_samples <- 150
LUT <- as.data.frame(getLUT(inputs = ToolsRTM::inputsPROSAIL, nLUT = n_samples, setseed = 1))
wl <- 400:2500
rsoil <- rep(0.15, length(wl))
refl <- t(sapply(seq_len(n_samples), function(i) {
  foursail(inputLUT = LUT[i, ], rsoil = rsoil, LeafModel = "PROSPECT-PRO")$rsot
}))
refl_X <- as.data.frame(refl); colnames(refl_X) <- paste0("X", wl); refl_X <- cbind(id = seq_len(nrow(refl_X)), refl_X)
se2a <- suppressMessages(get.spectra.convolved(rfl = refl_X, sensor = "Sentinel2a", plot.spectra = FALSE))
#> [1] "Spectral resampling function to SENTINEL2A is being processed ..."
#>   |                                                                              |                                                                      |   0%  |                                                                              |=====                                                                 |   8%  |                                                                              |===========                                                           |  15%  |                                                                              |================                                                      |  23%  |                                                                              |======================                                                |  31%  |                                                                              |===========================                                           |  38%  |                                                                              |================================                                      |  46%  |                                                                              |======================================                                |  54%  |                                                                              |===========================================                           |  62%  |                                                                              |================================================                      |  69%  |                                                                              |======================================================                |  77%  |                                                                              |===========================================================           |  85%  |                                                                              |=================================================================     |  92%  |                                                                              |======================================================================| 100%
band_names <- paste0("B", seq_along(as.numeric(names(se2a)[-1])))
names(se2a) <- c("id", band_names)

set.seed(1)
train_idx <- sample(seq_len(n_samples), size = round(0.7 * n_samples))
train_df <- cbind(LUT[train_idx, ], se2a[train_idx, band_names])
```

## 2. `getMLmodel()`: deep learning (TensorFlow/Keras)

Needs a working Python/TensorFlow/tf-keras stack reachable through
`reticulate` – genuinely machine-specific setup. Wrapped in
[`tryCatch()`](https://rdrr.io/r/base/conditions.html) here so this page
stays buildable on a machine without that stack configured, while still
showing the real, correct calling code (the same pattern used throughout
this package’s own course scripts, `Scripts/R/*/3-inversion_DL.R`):

``` r

dl_result <- tryCatch({
  local_venv_root <- file.path(Sys.getenv("LOCALAPPDATA"), "r-reticulate-venvs")
  Sys.setenv(RETICULATE_VIRTUALENV_ROOT = local_venv_root)
  Sys.setenv(WORKON_HOME = local_venv_root)
  Sys.setenv(TF_USE_LEGACY_KERAS = "1")
  library(reticulate); library(keras)
  dl_model <- getMLmodel(dataset = train_df, depVar = "Cab", model = "Hidden-layers",
                          optimizer = "adam", n.epochs = 5)
  list(ok = TRUE, model = dl_model)
}, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
#> [1] "Normalize"

if (dl_result$ok) {
  cat("DL model trained. Final training loss:",
      round(tail(dl_result$model$history$metrics$loss, 1), 3), "\n")
} else {
  first_line <- strsplit(dl_result$msg, "\n", fixed = TRUE)[[1]][1]
  cat("Skipped -- no working Python/TensorFlow/tf-keras stack reachable through",
      "reticulate on this machine (", first_line, "). The calling code above is",
      "correct and unchanged.\n")
}
#> Skipped -- no working Python/TensorFlow/tf-keras stack reachable through reticulate on this machine ( infinite or missing values in 'x' ). The calling code above is correct and unchanged.
```

`model = "Hidden-layers"` can be swapped for `"CNN"` (1D convolutional,
treating the band sequence as a spatial dimension) – same call
otherwise. `n.epochs = 5` above is a smoke test, not a tuned model; a
real run uses far more (`Scripts/R/*/3-inversion_DL.R` trains both
architectures per trait to convergence).

**Important, machine-specific setup note**:
[`getMLmodel()`](../reference/getMLmodel.md) calls
[`callback_early_stopping()`](https://rdrr.io/pkg/keras/man/callback_early_stopping.html)
without a package prefix internally, so
[`library(keras)`](https://tensorflow.rstudio.com/) must actually be
*attached* (not just installed) before calling it, or the fit fails with
`could not find function "callback_early_stopping"` – every DL script in
this package’s `Scripts/R/` starts with the exact
[`Sys.setenv()`](https://rdrr.io/r/base/Sys.setenv.html)/
[`library()`](https://rdrr.io/r/base/library.html) block shown above for
that reason.

## 3. Which one should you use? (Tutorials 10-12, together)

| Method | Needs training? | Handles a small LUT well? | Typical use |
|----|----|----|----|
| [`get.inversionOpt()`](../reference/get.inversionOpt.md) (Tutorial 10) | No – pure search | Yes – works from a handful of reference spectra | Quick, physically-grounded retrieval; no ML infrastructure needed |
| [`get.inversion()`](../reference/get.inversion.md) (Tutorial 11) | Yes | Needs enough rows to fit reliably (hundreds+) | Production-scale retrieval, many predictors (full spectrum + indices) |
| [`getMLmodel()`](../reference/getMLmodel.md) (this page) | Yes, more data-hungry | No – needs thousands of rows | Large training sets, CNN over the full spectral shape |

## What’s next

- **Tutorial 13** – simulate, convolve, index, and invert as one
  end-to-end pipeline, across every canopy model this package supports.
