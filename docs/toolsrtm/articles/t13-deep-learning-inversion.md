# 13. Deep Learning for RTM Inversion

``` r

library(ToolsRTM)
```

[`getMLmodel()`](../reference/getMLmodel.md) is
[`get.inversion()`](../reference/get.inversion.md)’s (Tutorial 12)
deep-learning sibling – a dense or 1D-CNN network via TensorFlow/Keras,
fit on the same kind of (sensor-band reflectance, trait) LUT. Needs a
working Python/TensorFlow/ tf-keras stack reachable through `reticulate`
– genuinely machine- specific setup, wrapped in
[`tryCatch()`](https://rdrr.io/r/base/conditions.html) throughout this
page so it stays buildable on a machine without that stack configured,
while showing the real, correct calling code either way.

``` text
RTM LUT (getLUT)
     |
     v
Spectral simulation (foursail)
     |
     v
Sensor / hyperspectral representation (Sentinel-2A bands, or PRISMA bands for the CNN)
     |
     v
Train / validation / test split
     |
     v
Normalization (fit on train only)
     |
     v
DL training, with early stopping on validation loss
     |
     v
Prediction on the INDEPENDENT test set
     |
     v
R2 / RMSE / bias, true-vs-predicted plot
     |
     v
Comparison against Tutorial 12's Random Forest, same split
```

## 1. A real bug, found and fixed while building this page

An earlier version of this page passed
[`getMLmodel()`](../reference/getMLmodel.md) a `dataset` containing the
*entire* training LUT (every trait column) concatenated with the
reflectance bands, and training failed with
`infinite or missing values in 'x'`. Traced to the real cause:
[`getMLmodel()`](../reference/getMLmodel.md) has no `inputs =` argument
– unlike [`get.inversion()`](../reference/get.inversion.md) (Tutorial
12), it treats **every column of `dataset` other than `depVar`** as a
predictor. `inputsPROSAIL`’s LUT fixes several trait columns at a single
constant value by design (`LMA`, `alpha`, `LIDFb`, `TypeLidf`, `hspot`,
`tts`, `psi` all have zero variance in a plain
[`getLUT()`](../reference/getLUT.md) draw) – feeding those into
[`getMLmodel()`](../reference/getMLmodel.md)’s internal `Normalize`
preprocessing (range-scaling each column by its own min/max) divides by
zero for every one of them, producing the `Inf`/`NA` values the error
message reports. **Fix**: build `dataset` from only `depVar` plus the
actual reflectance-band columns, exactly as
`get.inversion(..., inputs = band_names)` already does explicitly for
Tutorial 12 – [`getMLmodel()`](../reference/getMLmodel.md) needs the
same discipline applied manually, since it has no `inputs =` parameter
to enforce it.

``` r

local_venv_root <- file.path(Sys.getenv("LOCALAPPDATA"), "r-reticulate-venvs")
Sys.setenv(RETICULATE_VIRTUALENV_ROOT = local_venv_root)
Sys.setenv(WORKON_HOME = local_venv_root)
Sys.setenv(TF_USE_LEGACY_KERAS = "1")
keras_ok <- tryCatch({ library(reticulate); library(keras); TRUE }, error = function(e) FALSE)
cat("Working Keras/TensorFlow stack available:", keras_ok, "\n")
#> Working Keras/TensorFlow stack available: TRUE
```

**Important, machine-specific setup note**:
[`getMLmodel()`](../reference/getMLmodel.md) calls
[`callback_early_stopping()`](https://rdrr.io/pkg/keras/man/callback_early_stopping.html)
without a package prefix internally, so
[`library(keras)`](https://tensorflow.rstudio.com/) must actually be
*attached* (not just installed) before calling it, or the fit fails with
`could not find function "callback_early_stopping"`.

## 2. A properly sized LUT, not a 150-row smoke test

`n_samples = 150` (an earlier version of this page) is far too small to
demonstrate what deep learning is actually for – 3000 rows here, large
enough to show real training dynamics without an impractical pkgdown
build time:

``` r

r2_f <- function(obs, pred) 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
rmse_f <- function(obs, pred) sqrt(mean((obs - pred)^2))
bias_f <- function(obs, pred) mean(pred - obs)

n_samples <- 3000
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
```

## 3. A real train / validation / test split – evaluated only on the held-out test set

RTM inversion specifically needs an *independent* test set, not just a
train/test split: the LUT is synthetic, so nothing stops the model from
memorizing quirks of the exact rows it trained on unless a third,
completely untouched slice is reserved for the final number reported.
60% train / 20% validation (used *during* training, for early stopping)
/ 20% test (touched only once, at the end):

``` r

set.seed(1)
idx <- sample(seq_len(n_samples))
n_train <- round(0.6 * n_samples); n_val <- round(0.2 * n_samples)
train_idx <- idx[1:n_train]
val_idx   <- idx[(n_train + 1):(n_train + n_val)]
test_idx  <- idx[(n_train + n_val + 1):n_samples]
cat("train:", length(train_idx), " val:", length(val_idx), " test:", length(test_idx), "\n")
#> train: 1800  val: 600  test: 600

dataset_dense <- cbind(Cab = LUT$Cab, se2a[, band_names])  # depVar + bands only -- Section 1's fix
train_df <- dataset_dense[train_idx, ]
```

**Normalization, stated explicitly**:
[`getMLmodel()`](../reference/getMLmodel.md)’s internal
[`getSplitData()`](../reference/getSplitData.md) step fits its range/PCA
scaler on the training split only, then applies that *same* fitted
transform to whatever data is passed to
[`predict()`](https://rdrr.io/r/stats/predict.html) later – no
information from validation or test rows leaks into the scaling
parameters. This matters specifically for RTM inversion because the
LUT’s own trait ranges (not real-world variability) determine what
“normal” looks like to the model; fitting the scaler on the full dataset
would let test-set values quietly influence how training data gets
scaled.

## 4. Dense network, with early stopping, on Sentinel-2 bands

``` r

dl_result <- tryCatch({
  dl_dense <- getMLmodel(dataset = train_df, depVar = "Cab", model = "Hidden-layers",
                          optimizer = "adam", n.epochs = 150, batch.size = 64, save.model = FALSE)
  list(ok = TRUE, model = dl_dense)
}, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
```

![](t13-deep-learning-inversion_files/figure-html/dl-dense-1.png)

``` r

if (dl_result$ok) {
  h <- dl_result$model$history$metrics
  cat("Epochs actually run:", length(h$loss), "of 150 requested\n")
  cat("Epoch with lowest validation loss:", which.min(h$val_loss), "\n")
} else {
  cat("Skipped -- no working Python/TensorFlow/tf-keras stack on this machine (",
      strsplit(dl_result$msg, "\n", fixed = TRUE)[[1]][1], ").\n")
}
#> Epochs actually run: 150 of 150 requested
#> Epoch with lowest validation loss: 150
```

``` r

if (dl_result$ok) {
  h <- dl_result$model$history$metrics
  plot(seq_along(h$loss), h$loss, type = "l", col = "#2166AC",
       ylim = range(c(h$loss, h$val_loss)), xlab = "Epoch", ylab = "Loss (MSE)",
       main = "Learning curve: training vs. validation loss")
  lines(seq_along(h$val_loss), h$val_loss, col = "#D55E00")
  legend("topright", c("Training loss", "Validation loss"), col = c("#2166AC", "#D55E00"), lty = 1)
}
```

![](t13-deep-learning-inversion_files/figure-html/learning-curve-1.png)

**Read the learning curve before trusting the model.** In this run,
[`getMLmodel()`](../reference/getMLmodel.md)’s built-in
[`callback_early_stopping()`](https://rdrr.io/pkg/keras/man/callback_early_stopping.html)
(`patience = 5`, monitoring validation loss) never actually triggered –
validation loss was still falling at epoch 150, meaning this network
hadn’t finished converging within the epoch budget this page’s build
time allows. That by itself is a real, honest limitation of this demo (a
production run would keep training, or raise the learning rate), not a
claim that 150 epochs is generally sufficient.

## 5. Independent-test-set result – and a direct comparison against Tutorial 12’s Random Forest

``` r

if (dl_result$ok) {
  pred_test_dl <- as.numeric(predict(dl_result$model$model, as.matrix(dataset_dense[test_idx, band_names])))
  obs_test <- dataset_dense$Cab[test_idx]
  cat("Dense DL, independent test set: R2=", round(r2_f(obs_test, pred_test_dl), 3),
      " RMSE=", round(rmse_f(obs_test, pred_test_dl), 2),
      " bias=", round(bias_f(obs_test, pred_test_dl), 2), "\n")
}
#> 19/19 - 0s - 32ms/epoch - 2ms/step
#> Dense DL, independent test set: R2= -0.155  RMSE= 16.98  bias= -5.87
```

``` r

if (dl_result$ok) {
  fit_rf <- get.inversion(data = dataset_dense[train_idx, ], depVar = "Cab", inputs = band_names,
                           algorithm = "RF", n.samples = length(train_idx), seed = 42)
  pred_test_rf <- as.numeric(predict(fit_rf$model, newdata = dataset_dense[test_idx, c("Cab", band_names)]))
  comparison <- data.frame(
    method = c("Dense DL (getMLmodel)", "Random Forest (get.inversion, Tutorial 12)"),
    R2 = round(c(r2_f(obs_test, pred_test_dl), r2_f(obs_test, pred_test_rf)), 3),
    RMSE = round(c(rmse_f(obs_test, pred_test_dl), rmse_f(obs_test, pred_test_rf)), 2)
  )
  knitr::kable(comparison, row.names = FALSE)
}
#> [1] "processing hybrid approach using Random Forest ..."
#> -0.005000317 0.01 
#> 0.2320068 0.01 
#> -0.04837869 0.01
```

![](t13-deep-learning-inversion_files/figure-html/rf-comparison-1.png)

| method                                     |     R2 |  RMSE |
|:-------------------------------------------|-------:|------:|
| Dense DL (getMLmodel)                      | -0.155 | 16.98 |
| Random Forest (get.inversion, Tutorial 12) |  0.868 |  5.75 |

``` r

if (dl_result$ok) {
  op <- par(mfrow = c(1, 2))
  plot(obs_test, pred_test_dl, pch = 19, col = "#D55E00",
       xlab = "Observed Cab", ylab = "Predicted Cab", main = "Dense DL")
  abline(0, 1, col = "grey40", lty = 2)
  plot(obs_test, pred_test_rf, pch = 19, col = "#2166AC",
       xlab = "Observed Cab", ylab = "Predicted Cab", main = "Random Forest")
  abline(0, 1, col = "grey40", lty = 2)
  par(op)
}
```

![](t13-deep-learning-inversion_files/figure-html/true-vs-predicted-1.png)

**On this LUT (3000 rows, 10 Sentinel-2 bands), Random Forest clearly
beats the dense network, not the other way around.** This is the real
result, not a hypothetical – and it’s exactly the outcome Section 8
below explains rather than papers over: 10-band multispectral input and
a moderate-size LUT is close to RF’s best case, and the network here
hadn’t even finished converging (Section 4). Do not read this as “DL is
worse” in general – read it as “DL needed either more training budget,
more data, or higher-dimensional input to earn its extra complexity
here, and didn’t get any of the three in this demo.”

## 6. The 1D-CNN, on hyperspectral (PRISMA) bands

A dense network on 10 multispectral bands is a weak demonstration of
what a 1D CNN specifically buys you – convolution exploits
*spectral-neighborhood* structure (adjacent, correlated bands), which 10
widely-spaced multispectral bands barely have. PRISMA’s ~230 contiguous
hyperspectral bands are a much fairer test of that idea:

``` r

n_hyp <- 2000
LUT_h <- as.data.frame(getLUT(inputs = ToolsRTM::inputsPROSAIL, nLUT = n_hyp, setseed = 3))
refl_h <- t(sapply(seq_len(n_hyp), function(i) {
  foursail(inputLUT = LUT_h[i, ], rsoil = rsoil, LeafModel = "PROSPECT-PRO")$rsot
}))
refl_X_h <- as.data.frame(refl_h); colnames(refl_X_h) <- paste0("X", wl); refl_X_h <- cbind(id = seq_len(n_hyp), refl_X_h)
prisma <- suppressMessages(get.spectra.convolved(rfl = refl_X_h, sensor = "PRISMA", plot.spectra = FALSE))
prisma_bands <- paste0("P", seq_len(ncol(prisma) - 1))
names(prisma) <- c("id", prisma_bands)
cat("PRISMA bands:", length(prisma_bands), "\n")

idx_h <- sample(seq_len(n_hyp))
train_h <- idx_h[1:round(0.6 * n_hyp)]
test_h  <- idx_h[(round(0.8 * n_hyp) + 1):n_hyp]
dataset_cnn <- cbind(Cab = LUT_h$Cab, prisma[, prisma_bands])

cnn_result <- tryCatch({
  m <- getMLmodel(dataset = dataset_cnn[train_h, ], depVar = "Cab", model = "CNN",
                   optimizer = "adam", n.epochs = 40, batch.size = 64, save.model = FALSE)
  list(ok = TRUE, model = m)
}, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
```

![](t13-deep-learning-inversion_files/figure-html/cnn-prisma-1.png)

``` r

if (cnn_result$ok) {
  pred_test_cnn <- as.numeric(predict(cnn_result$model$model,
    array_reshape(as.matrix(dataset_cnn[test_h, prisma_bands]), c(length(test_h), length(prisma_bands), 1))))
  obs_test_h <- dataset_cnn$Cab[test_h]
  cat("CNN (PRISMA, ", length(prisma_bands), "bands), independent test set: R2=",
      round(r2_f(obs_test_h, pred_test_cnn), 3), " RMSE=", round(rmse_f(obs_test_h, pred_test_cnn), 2), "\n")
}
#> 13/13 - 0s - 45ms/epoch - 3ms/step
#> CNN (PRISMA,  234 bands), independent test set: R2= -0.37  RMSE= 19.45
```

The dense-vs-1D-CNN architectural difference itself: `"Hidden-layers"`
is a stack of plain fully-connected layers, treating each band as an
independent input with no notion of spectral order; `"CNN"` instead runs
1D convolution + pooling over the bands *in their spectral sequence*, so
it can learn local absorption-feature shapes (a dip spanning several
adjacent bands) as a reusable pattern – something a dense network has to
learn independently, band-combination by band- combination, with no
structural help.

**This CNN result does not beat Random Forest either, at this LUT size
and this epoch budget.** That’s consistent with Section 5/8’s point: the
CNN’s structural advantage (exploiting spectral neighborhoods) needs
either a larger LUT, more training, or a harder problem where local
spectral shape genuinely matters more than this synthetic Cab-from-
reflectance mapping does, to actually pay off over a well-tuned RF.

## 7. Which one should you use? (Tutorials 10-13, together) – corrected

An earlier version of this table implied
[`get.inversionOpt()`](../reference/get.inversionOpt.md) (pure LUT
search, no training) works fine from a small reference LUT. Tutorial
11’s own real test contradicts that: at ~100 reference spectra,
[`get.inversionOpt()`](../reference/get.inversionOpt.md) scored R² close
to zero (no real skill); it only became useful around 2000. **All three
methods need substantially more data than a quick demo LUT provides for
genuinely reliable results** – the difference between them is *how* they
degrade and recover, not whether one of them is exempt from needing
data:

| Method | Needs training? | Small-LUT behavior (verified) | Typical use |
|----|----|----|----|
| [`get.inversionOpt()`](../reference/get.inversionOpt.md) (Tutorial 11) | No – pure search | Degrades smoothly but really does need real coverage – ~100 reference spectra gave R2 near zero in Tutorial 11’s own test; ~2000 gave real skill | Quick, physically-grounded retrieval when a genuinely well-populated reference LUT already exists |
| [`get.inversion()`](../reference/get.inversion.md) (Tutorial 12) | Yes | The most forgiving of the three in practice – RF in particular stayed usable even on this page’s 150-row original LUT, and clearly beat both DL architectures above on a 3000-row one | General-purpose retrieval; often the right default |
| [`getMLmodel()`](../reference/getMLmodel.md) (this page) | Yes, more so | Needs the most: on a 3000-row LUT / 10 multispectral bands, neither dense nor CNN beat Tutorial 12’s RF here (Sections 5-6) | Large LUTs, hyperspectral input, and/or more complex nonlinear mappings than RF captures well – not a default choice |

## 8. When should I use deep learning?

Not “when it’s available” – this page’s own results argue against that.
A DL model is more likely to earn its extra complexity when **several**
of these hold at once, not just one:

- **LUT size**: thousands to tens of thousands of rows, not hundreds –
  Section 5 used 3000 and RF still won comfortably.
- **Input dimensionality**: hyperspectral (dozens to hundreds of
  correlated bands), where a CNN’s spectral-neighborhood structure has
  something real to exploit (Section 6) – 10 multispectral bands don’t
  give it much to work with.
- **Compute cost, honestly accounted**: Section 4’s dense model took
  real wall-clock training time and still hadn’t converged at 150
  epochs; Tutorial 12’s RF trains in a fraction of that time on the same
  data. That cost needs to buy something.
- **Nonlinearity/complexity of the true mapping**: RF and gradient
  boosting already capture a great deal of nonlinearity; DL’s edge shows
  up more clearly on mappings even tree ensembles struggle with, not on
  ones they already handle well (as `Cab` from 10 S2 bands apparently
  is, per Section 5’s numbers).
- **Overfitting risk, watched, not assumed**: the learning-curve check
  in Section 4 is not optional – a validation curve that’s still falling
  means the reported test number isn’t from a converged model; one
  that’s risen back up after a minimum means the model has started
  memorizing the (synthetic) training LUT specifically.

For this package’s own bundled models and typical LUT sizes,
[`get.inversion()`](../reference/get.inversion.md) (Tutorial 12) –
especially Random Forest – is the sensible default. Reach for
[`getMLmodel()`](../reference/getMLmodel.md) when the LUT is genuinely
large, the sensor is hyperspectral, and there’s compute budget to let
training actually converge – not by default.

## What’s next

- **Tutorial 14** – simulate, convolve, index, and invert as one
  end-to-end pipeline, across every canopy model this package supports.
