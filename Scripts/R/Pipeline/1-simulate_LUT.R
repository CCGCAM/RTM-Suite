# ==============================================================================
# Pipeline step 1: simulate a LUT of leaf+canopy reflectance spectra, convolve
# it to real sensor bands, and compute vegetation indices -- all in one script,
# since these three things are tightly coupled (same LUT/simulation feeds all
# of it) and splitting them across files just meant re-reading/re-writing the
# same intermediate data. Steps 2 (classic ML inversion) and 3 (deep-learning
# inversion) stay separate scripts, since either can be re-run independently
# against the dataset this script produces.
#
# User controls: n.samples, leaf/canopy model, soil source, AND (new):
#   - trait.distribution: override the sampling distribution (Uniform/
#     Gaussian) for individual traits, drawn AFTER the main getLUT() call via
#     ToolsRTM::gauss_byMin_Max() (truncated-normal rejection sampling, still
#     respects each trait's normal min/max bounds from inputsPROSAIL).
#   - correlate.traits: make one trait co-vary with another (e.g. Car with
#     Cab -- real leaf pigments are correlated, not independent) via
#     ToolsRTM::correlatedValue(), applied AFTER getLUT()/distribution
#     overrides so it always uses each trait's final sampled values.
# ==============================================================================

rm(list = ls())  # avoid leftover objects from a previous run/session leaking in

if (!requireNamespace("ToolsRTM", quietly = TRUE)) {
  remotes::install_gitlab("caminoccg/toolsrtm")
}
if (!requireNamespace("SCOPEinR", quietly = TRUE)) {
  remotes::install_gitlab("caminoccg/scopeinr")
}
library(ToolsRTM)
library(SCOPEinR)
library(ggplot2)  # ToolsRTM doesn't re-export it; needed for the diagnostic plot below

## ----------------------------------------------------------------------------
## User-configurable settings
## ----------------------------------------------------------------------------
##
## This ONE script already covers fourSAIL, foursail2, and INFORM -- there is
## no separate script per canopy model. Switch `canopy.model` below and re-run;
## everything downstream (soil, sensor convolution, indices, saved datasets,
## and steps 2/3's inversion) works unchanged for all three:
##
##   canopy.model <- "fourSAIL"    # classic turbid-medium canopy (default)
##   canopy.model <- "foursail2"   # two-layer canopy (understory + overstory)
##   canopy.model <- "INFORM"      # forest canopy with explicit crown geometry
##
## Sensor choice is likewise not a separate script: `sensors` below convolves
## to all three every run, and `datasets$native` / `$se2a` / `$prisma` (saved
## to 1-datasets.rds) select which one steps 2/3 invert from -- see those
## scripts' `dataset <- datasets$native` line.
n.samples    <- 500
leaf.model   <- "PROSPECT-D"     # one of ToolsRTM's 5 leaf models
canopy.model <- "fourSAIL"       # "fourSAIL", "foursail2", or "INFORM" -- see note above
soil.source  <- "marmit"         # "flat" (rep(0.15, ...)) or "marmit" (ToolsRTM::get.marmit.rsoil())
seed         <- 1
sensors      <- c("Sentinel2a", "Sentinel2b", "PRISMA")

# Override the sampling distribution for specific traits (default for every
# trait comes from inputsPROSAIL$Distribution -- mostly "Uniform", Cab is
# "Gaussian"). Each entry re-draws that column via a truncated normal
# (mean/sd from inputsPROSAIL if set there, else midpoint/range-based).
# Set to list() to disable and use inputsPROSAIL's defaults untouched.
trait.distribution <- list(Cab = "Gaussian")

# Make one trait co-vary with another (real leaf pigments/traits aren't
# independent -- e.g. Car tracks Cab). Each entry: redraw `name` as
# correlatedValue(LUT[[with]] * scale, r). Set to list() to disable.
correlate.traits <- list(Car = list(with = "Cab", scale = 1/4, r = 0.8))

# Own subfolder per canopy model, so running fourSAIL/foursail2/INFORM in
# turn doesn't overwrite each other's saved dataset/plot.
out_dir <- file.path("outs/rtm_sims/pipeline", canopy.model)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("=== Simulating %d spectra (%s + %s, soil = %s) ===\n",
            n.samples, leaf.model, canopy.model, soil.source))

## ----------------------------------------------------------------------------
## 1. Build the LUT (trait values)
## ----------------------------------------------------------------------------

inputs <- ToolsRTM::inputsPROSAIL
LUT <- as.data.frame(getLUT(inputs = inputs, nLUT = n.samples, setseed = seed))

# --- distribution overrides ---
for (trait in names(trait.distribution)) {
  dist_i <- trait.distribution[[trait]]
  row_i <- inputs[inputs$variable == trait, ]
  if (nrow(row_i) == 0) { warning("trait.distribution: '", trait, "' not in inputsPROSAIL, skipped"); next }
  lwr <- as.numeric(row_i$lower); upr <- as.numeric(row_i$upper)
  if (dist_i == "Gaussian") {
    m <- suppressWarnings(as.numeric(row_i$Mean_D))
    s <- suppressWarnings(as.numeric(row_i$Std_D))
    if (is.na(m)) m <- (lwr + upr) / 2   # inputsPROSAIL stores "-" (non-numeric) when unset
    if (is.na(s)) s <- (upr - lwr) / 4
    LUT[[trait]] <- gauss_byMin_Max(n = n.samples, m = m, s = s, lwr = lwr, upr = upr, nnorm = n.samples * 4)
    cat(sprintf("Distribution override: %s -> Gaussian(mean=%.3g, sd=%.3g, bounds=[%.3g, %.3g])\n", trait, m, s, lwr, upr))
  } else if (dist_i == "Uniform") {
    LUT[[trait]] <- stats::runif(n.samples, lwr, upr)
    cat(sprintf("Distribution override: %s -> Uniform(%.3g, %.3g)\n", trait, lwr, upr))
  } else {
    warning("trait.distribution: unknown distribution '", dist_i, "' for '", trait, "', skipped")
  }
}

# --- trait correlations (applied after distribution overrides, so it uses final values) ---
for (trait in names(correlate.traits)) {
  spec <- correlate.traits[[trait]]
  if (!spec$with %in% names(LUT)) { warning("correlate.traits: '", spec$with, "' not in LUT, skipped"); next }
  scale_ <- if (is.null(spec$scale)) 1 else spec$scale
  LUT[[trait]] <- correlatedValue(x = LUT[[spec$with]] * scale_, r = spec$r)
  cat(sprintf("Correlation: %s ~ %s (scale=%.3g, target r=%.2f, realized r=%.2f)\n",
              trait, spec$with, scale_, spec$r, cor(LUT[[trait]], LUT[[spec$with]])))
}

# inputsPROSAIL only ships PROSPECT-PRO/-D columns -- add what Fluspect-B/-Cx
# and Liberty need too (same requirement as foursail() itself, see
# Scripts/compare_RTM_models.R).
LUT$Cs <- 0; LUT$fqe <- 0.01; LUT$Cx <- 0
LUT$cell.d <- 40; LUT$inter.c <- 0.045; LUT$baseline.abs <- 0.0006
LUT$leaf.thick <- 1.6; LUT$albino.abs <- 0; LUT$lign.cell <- 2; LUT$Nitrogen <- 1
# NOTE: `is.null(LUT[["x"]])` (exact match), NOT `is.null(LUT$x)` -- data.frame's
# `$` does partial name matching, so e.g. `LUT$h` silently resolved to the
# existing `hspot` column (no exact `h` column ever existed, but is.null(LUT$h)
# came back FALSE), so this never actually added the `h` column INFORM needs.
if (is.null(LUT[["fraction_brown"]])) LUT$fraction_brown <- 0.1
if (is.null(LUT[["diss"]])) LUT$diss <- 0.5
if (is.null(LUT[["Cv"]])) LUT$Cv <- 1
if (is.null(LUT[["Zeta"]])) LUT$Zeta <- 0
if (is.null(LUT[["LAIu"]])) LUT$LAIu <- 0.5
if (is.null(LUT[["sd"]])) LUT$sd <- 650
if (is.null(LUT[["cd"]])) LUT$cd <- 4.5
if (is.null(LUT[["h"]])) LUT$h <- 20
if (is.null(LUT[["skyl"]])) LUT$skyl <- 0.1

## ----------------------------------------------------------------------------
## 2. Soil spectrum
## ----------------------------------------------------------------------------

short_domain <- leaf.model %in% c("Fluspect-B", "Fluspect-B-Cx")
wl <- if (short_domain) 400:2400 else 400:2500

if (soil.source == "marmit") {
  soil <- get.marmit.rsoil(database = "Bablet_2016", id = 1, L = 0.05, eps = 0.3, wl.out = wl)
  rsoil <- soil$rsoil.wet
  cat(sprintf("Soil: MARMIT (Bablet_2016 id=1, SMC = %.1f%%)\n", soil$SMC))
} else {
  rsoil <- rep(0.15, length(wl))
  cat("Soil: flat 0.15\n")
}

## ----------------------------------------------------------------------------
## 3. Run the canopy model for every LUT row
## ----------------------------------------------------------------------------

refl <- t(sapply(seq_len(n.samples), function(i) {
  sim <- suppressMessages(simulate_RTM(inputLUT = LUT[i, ], rsoil = rsoil,
                                        leaf.model = leaf.model, canopy.model = canopy.model))
  sim$rsot
}))
colnames(refl) <- paste0("R.", wl)
refl <- as.data.frame(refl)

cat(sprintf("Simulated: %d x %d reflectance matrix\n", nrow(refl), ncol(refl)))

## ----------------------------------------------------------------------------
## 4. Convolve to sensors
## ----------------------------------------------------------------------------

refl_X <- refl
colnames(refl_X) <- paste0("X", wl)
refl_X <- cbind(id = seq_len(nrow(refl_X)), refl_X)

convolved <- list()
for (s in sensors) {
  cat(sprintf("Convolving to %s ... ", s))
  convolved[[s]] <- suppressMessages(get.spectra.convolved(rfl = refl_X, sensor = s, plot.spectra = FALSE))
  cat(sprintf("%d bands\n", ncol(convolved[[s]]) - 1))
}

## ----------------------------------------------------------------------------
## 5. Vegetation indices per sensor + assembled modeling datasets
## ----------------------------------------------------------------------------

cat("Computing indices: native ... ")
indices_native <- suppressMessages(getIndices(refl, pattern.rfl = "R.", spectral.domain = "VNIR"))
indices_native <- indices_native[, colSums(is.na(indices_native)) == 0, drop = FALSE]
cat(sprintf("%d indices\n", ncol(indices_native)))

cat("Computing indices: Sentinel-2A ... ")
se2a_bands <- convolved$Sentinel2a
names(se2a_bands) <- c("id", "B1","B2","B3","B4","B5","B6","B7","B8","B8A","B9","B10","B11","B12")
indices_se2a <- suppressMessages(getIndicesSE2(se2a_bands[, -1], sensor = "Sentinel-2a", df.data = NULL, fast.process = TRUE))
cat(sprintf("%d indices\n", ncol(indices_se2a)))

cat("Computing indices: PRISMA ... ")
prisma_bands <- convolved$PRISMA
wl_prisma <- as.numeric(names(prisma_bands)[-1])
names(prisma_bands) <- c("id", paste0("R.", round(wl_prisma)))
indices_prisma <- tryCatch(
  suppressMessages(getIndices(prisma_bands[, -1], pattern.rfl = "R.", spectral.domain = "VNIR-SWIR")),
  error = function(e) { cat("failed:", conditionMessage(e), " "); NULL }
)
if (!is.null(indices_prisma)) {
  indices_prisma <- indices_prisma[, colSums(is.na(indices_prisma)) == 0, drop = FALSE]
  cat(sprintf("%d indices\n", ncol(indices_prisma)))
}

dataset_native <- cbind(LUT, refl, indices_native)
dataset_native <- dataset_native[, !duplicated(names(dataset_native))]

dataset_se2a <- cbind(LUT, se2a_bands[, -1], indices_se2a)
dataset_se2a <- dataset_se2a[, !duplicated(names(dataset_se2a))]

dataset_prisma <- NULL
if (!is.null(indices_prisma)) {
  dataset_prisma <- cbind(LUT, prisma_bands[, -1], indices_prisma)
  dataset_prisma <- dataset_prisma[, !duplicated(names(dataset_prisma))]
}

## ----------------------------------------------------------------------------
## 6. Diagnostic plot -- quick visual check the run looks physically sensible
## ----------------------------------------------------------------------------
## print()ed so it shows up if you're running this interactively, and always
## saved to outs/ regardless (so a batch/Rscript run leaves visible proof it
## worked, not just an .rds file).

refl.cols.plot <- colnames(refl)[seq(1, ncol(refl), by = max(1, ncol(refl) %/% 400))]  # subsample for a lighter plot
plot.rfl <- cbind(LUT[, "Cab", drop = FALSE], stack(refl[, refl.cols.plot]))
colnames(plot.rfl) <- c("Cab", "RFL", "wave")
plot.rfl$wave <- as.numeric(sub("R.", "", plot.rfl$wave, fixed = TRUE))
plot.rfl$GroupID <- cut(plot.rfl$Cab, breaks = c(0, 20, 30, 40, 60, 80, Inf),
                         labels = c('0-20', '20-30', '30-40', '40-60', '60-80', '>80'))
plot.rfl.mean <- aggregate(RFL ~ GroupID + wave, data = plot.rfl, FUN = "mean")

spectral_plot <- ggplot(plot.rfl.mean, aes(x = wave, y = RFL, group = GroupID)) +
  theme_bw() +
  geom_line(aes(color = GroupID), linetype = "dashed", linewidth = 0.7) +
  labs(color = 'Cab (ug/cm2)', x = 'Wavelength (nm)', y = 'Reflectance') +
  ggtitle(sprintf("%s + %s reflectance by Cab (n=%d)", canopy.model, leaf.model, n.samples)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

print(spectral_plot)
plot_path <- file.path(out_dir, "0-reflectance_spectra_byCab.png")
ggsave(plot_path, plot = spectral_plot, width = 14, height = 9, dpi = 300, units = "cm")
cat("Saved diagnostic plot to '", plot_path, "'\n", sep = "")

## ----------------------------------------------------------------------------
## 7. Save
## ----------------------------------------------------------------------------

saveRDS(list(LUT = LUT, refl = refl, wl = wl, convolved = convolved,
             datasets = list(native = dataset_native, se2a = dataset_se2a, prisma = dataset_prisma),
             leaf.model = leaf.model, canopy.model = canopy.model,
             soil.source = soil.source, n.samples = n.samples, seed = seed),
        file.path(out_dir, "1-datasets.rds"))

cat("\nSaved to '", out_dir, "/1-datasets.rds' -- read back with readRDS() in steps 2/3 ",
    "($datasets$native / $se2a / $prisma are ready for inversion).\n", sep = "")
