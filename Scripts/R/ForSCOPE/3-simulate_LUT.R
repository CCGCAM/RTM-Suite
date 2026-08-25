# ==============================================================================
# ForSCOPE course pipeline, step 1/3: simulate N SCOPE runs, convolve the
# resulting TOC reflectance to real sensor bands, compute vegetation indices,
# and extract solar-induced fluorescence (SIF) at 687/760nm as extra
# features. The SCOPEinR equivalent of Scripts/ForPROSAIL/3-simulate_LUT.R
# etc, adapted from the verified Scripts/Pipeline/SCOPE-1-simulate.R.
#
# get.SCOPE()'s leaf.model/canopy.model arguments are NOT functional (SCOPE
# always runs its own integral multi-layer Fluspect-Cx + RTMo, see
# Scripts/compare_SCOPE_models.R) -- there is nothing to configure there.
#
# Soil: get.SCOPE() sources soil reflectance either from a bundled file
# (soil_spectra/soil_scope.txt, 3 selectable spectra via the LUT's `spectrum`
# column) or from the BSM model (options.SCOPE row 5, options.soilspectrum),
# controlled by SCOPEinR/inst/input/setoptions.csv -- NOT via a passed-in
# rsoil vector like ToolsRTM::foursail()/SPART(). This script uses SCOPE's
# own built-in soil options as-is.
#
# get.SCOPE.parallel() already has a working `parallel = TRUE/FALSE` switch
# (no ToolsRTM::get.inversion()-style bug here) -- exposed below as
# `use.parallel`.
#
# SIF (solar-induced fluorescence): calc_fluorescence=1 is on by default in
# setoptions.csv, so every run already computes fluorescence radiance
# (data.rad$LoF_ over data.spectral$wlF, ~640-850nm) -- this script pulls out
# SIF687/SIF760 (the two classic single-band SIF metrics used in remote
# sensing) as extra columns, so 4-inversion_ML.R can test whether they help
# predict Vcmax25 beyond reflectance/indices alone (fluorescence is a direct
# byproduct of photosynthetic electron transport, so there's real physical
# reason to expect it correlates with photosynthetic capacity).
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
library(ggplot2)

## ----------------------------------------------------------------------------
## User-configurable settings
## ----------------------------------------------------------------------------

n.samples    <- 100              # always 100 for this course pipeline
use.parallel <- FALSE            # TRUE = get.SCOPE.parallel(parallel = TRUE); FALSE = sequential
seed         <- 1
sensors      <- c("Sentinel2a", "Sentinel2b", "PRISMA")

out_dir <- "outs/ForSCOPE"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
pdf(file.path(out_dir, "Rplots.pdf"))  # captures any stray plot() call, see ForPROSAIL/3-simulate_LUT.R for why

cat(sprintf("=== Simulating %d SCOPE runs (parallel = %s) ===\n", n.samples, use.parallel))

## ----------------------------------------------------------------------------
## 1. Build the LUT and run SCOPE
## ----------------------------------------------------------------------------

opts <- read.table(system.file("input", "setoptions.csv", package = "SCOPEinR"), header = TRUE, sep = ",")
inputLUT <- read.table(system.file("input", "inputs_SCOPE.csv", package = "SCOPEinR"), header = TRUE, sep = ",")
LUT <- getLUT.SCOPE(inputLUT = inputLUT, nLUT = n.samples, setseed = seed)

## ----------------------------------------------------------------------------
## Diagnostic figures: trait distributions and correlations (incl. Vcmax25)
## ----------------------------------------------------------------------------

key_traits <- intersect(c("Cab", "LAI", "EWT", "Vcmax25"), names(LUT))

hist_data <- stack(LUT[, key_traits])
colnames(hist_data) <- c("value", "trait")
hist_plot <- ggplot(hist_data, aes(x = value)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  facet_wrap(~trait, scales = "free") +
  theme_bw() +
  labs(x = NULL, y = "count", title = sprintf("Sampled SCOPE LUT distributions (n=%d)", n.samples)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
print(hist_plot)
ggsave(file.path(out_dir, "0-trait_histograms.png"), plot = hist_plot, width = 16, height = 10, dpi = 300, units = "cm")

corr_mat <- cor(LUT[, key_traits])
png(file.path(out_dir, "0-trait_correlations.png"), width = 1300, height = 1300, res = 200)
corrplot::corrplot(corr_mat, method = "color", type = "upper",
                    addCoef.col = "black", tl.col = "black", tl.srt = 45, number.cex = 0.9,
                    title = "Sampled LUT trait correlations", mar = c(0, 0, 2, 0))
dev.off()
corrplot::corrplot(corr_mat, method = "color", type = "upper",
                    addCoef.col = "black", tl.col = "black", tl.srt = 45, number.cex = 0.9,
                    title = "Sampled LUT trait correlations", mar = c(0, 0, 2, 0))

# get.SCOPE.parallel() runs the whole LUT (or each parallel chunk) through a
# single get.SCOPE() call with no per-row error isolation -- one row with an
# extreme random parameter combination that makes the energy-balance solver
# blow up can crash the entire batch with no R-level error message (a native/
# numerical fault, not a clean stop()). Calling get.SCOPE() one row at a time,
# each wrapped in tryCatch, fixes that: a bad row is skipped and logged, not
# fatal to the rest. `use.parallel` still controls whether those per-row calls
# run on a cluster (foreach %dopar%) or sequentially -- robustness and
# parallelism aren't mutually exclusive here, unlike the old bulk-call approach.
run_one_row <- function(i) {
  tryCatch(
    suppressMessages(suppressWarnings(
      get.SCOPE(LUT = LUT[i, ], n.LUT = 1, options.SCOPE = opts, optipar = SCOPEinR::optipar2021.Pro.CX,
                get.outputs = "ALL", get.plots = FALSE)[[1]]
    )),
    error = function(e) { message(sprintf("Row %d FAILED: %s", i, conditionMessage(e))); NULL }
  )
}

t0 <- Sys.time()
n.cores <- if (use.parallel) max(1, parallel::detectCores() - 2) else 1

if (n.cores > 1) {
  cl <- parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl)
  i <- NULL  # avoid "no visible binding" NOTE for the foreach loop variable
  sim_list <- foreach::foreach(i = seq_len(n.samples), .packages = "SCOPEinR") %dopar% run_one_row(i)
  parallel::stopCluster(cl)
} else {
  sim_list <- vector("list", n.samples)
  for (i in seq_len(n.samples)) {
    sim_list[[i]] <- run_one_row(i)
    if (i %% 50 == 0) cat(sprintf("  ... %d/%d done\n", i, n.samples))
  }
}

ok <- !sapply(sim_list, is.null)
sims <- sim_list[ok]
LUT <- LUT[ok, , drop = FALSE]
failed_rows <- which(!ok)
cat(sprintf("Simulated %d/%d runs in %.1fs (%d row(s) failed: %s)\n",
            length(sims), n.samples, as.numeric(Sys.time() - t0, units = "secs"),
            length(failed_rows), paste(failed_rows, collapse = ", ")))

## ----------------------------------------------------------------------------
## 2. Extract TOC reflectance into a matrix (native SCOPE wavelength grid)
## ----------------------------------------------------------------------------

wl <- sims[[1]]$data.spectral$wlS  # SCOPE's native output wavelength grid (extends to thermal, up to 50000nm)
refl <- t(sapply(sims, function(s) s$data.rad$refl))
colnames(refl) <- paste0("R.", round(wl))
refl <- as.data.frame(refl)

## ----------------------------------------------------------------------------
## 2b. Extract SIF687/SIF760 -- the two classic single-band solar-induced
## fluorescence metrics, from data.rad$LoF_ (fluorescence radiance) over
## data.spectral$wlF (~640-850nm, only populated when calc_fluorescence=1,
## the default in setoptions.csv). Used in 4-inversion_ML.R to test whether
## SIF adds predictive value for Vcmax25 beyond reflectance/indices alone.
## ----------------------------------------------------------------------------

wlF <- sims[[1]]$data.spectral$wlF
i687 <- which.min(abs(wlF - 687))
i760 <- which.min(abs(wlF - 760))
LUT$SIF687 <- sapply(sims, function(s) s$data.rad$LoF_[i687])
LUT$SIF760 <- sapply(sims, function(s) s$data.rad$LoF_[i760])
cat(sprintf("SIF687 range: %.3f-%.3f, SIF760 range: %.3f-%.3f (W m-2 um-1 sr-1)\n",
            min(LUT$SIF687), max(LUT$SIF687), min(LUT$SIF760), max(LUT$SIF760)))

# A small number of bands (strong water-vapor absorption regions, e.g. ~1350-1450nm)
# come back NA from SCOPE's radiative transfer -- physically expected (near-zero
# photons to work with there), not a bug, but downstream inversion/index code
# can't handle NA columns, so drop them here rather than at every consumer.
na_cols <- colSums(is.na(refl)) > 0
if (any(na_cols)) {
  cat(sprintf("Dropping %d band(s) with NA values (e.g. strong water-vapor absorption regions): %s\n",
              sum(na_cols), paste(head(wl[na_cols], 5), collapse = ", ")))
  refl <- refl[, !na_cols, drop = FALSE]
  wl <- wl[!na_cols]
}

cat(sprintf("Reflectance matrix: %d x %d (%.0f-%.0fnm)\n", nrow(refl), ncol(refl), min(wl), max(wl)))

## ----------------------------------------------------------------------------
## 3. Convolve to sensors
## ----------------------------------------------------------------------------

refl_X <- refl
colnames(refl_X) <- paste0("X", round(wl))
refl_X <- cbind(id = seq_len(nrow(refl_X)), refl_X)

convolved <- list()
for (s in sensors) {
  cat(sprintf("Convolving to %s ... ", s))
  convolved[[s]] <- tryCatch(
    suppressMessages(get.spectra.convolved(rfl = refl_X, sensor = s, plot.spectra = FALSE)),
    error = function(e) { cat("failed:", conditionMessage(e), " "); NULL }
  )
  if (!is.null(convolved[[s]])) cat(sprintf("%d bands\n", ncol(convolved[[s]]) - 1))
}

## ----------------------------------------------------------------------------
## 4. Vegetation indices (native + per sensor)
## ----------------------------------------------------------------------------

cat("Computing indices: native ... ")
indices_native <- suppressMessages(getIndices(refl, pattern.rfl = "R.", spectral.domain = "VNIR"))
indices_native <- indices_native[, colSums(is.na(indices_native)) == 0, drop = FALSE]
cat(sprintf("%d indices\n", ncol(indices_native)))

se2a_bands <- convolved$Sentinel2a
indices_se2a <- NULL
if (!is.null(se2a_bands)) {
  cat("Computing indices: Sentinel-2A ... ")
  names(se2a_bands) <- c("id", "B1","B2","B3","B4","B5","B6","B7","B8","B8A","B9","B10","B11","B12")
  indices_se2a <- suppressMessages(getIndicesSE2(se2a_bands[, -1], sensor = "Sentinel-2a", df.data = NULL, fast.process = TRUE))
  cat(sprintf("%d indices\n", ncol(indices_se2a)))
}

## ----------------------------------------------------------------------------
## 5. Assemble modeling datasets: target traits (from LUT) + bands + indices
## ----------------------------------------------------------------------------

dataset_native <- cbind(LUT, refl, indices_native)
dataset_native <- dataset_native[, !duplicated(names(dataset_native))]

dataset_se2a <- NULL
if (!is.null(se2a_bands)) {
  dataset_se2a <- cbind(LUT, se2a_bands[, -1], indices_se2a)
  dataset_se2a <- dataset_se2a[, !duplicated(names(dataset_se2a))]
}

## ----------------------------------------------------------------------------
## 6. Diagnostic plot -- quick visual check the run looks physically sensible
## ----------------------------------------------------------------------------
## print()ed so it shows up if you're running this interactively, and always
## saved to outs/ regardless (so a batch/Rscript run leaves visible proof it
## worked, not just an .rds file).

# SCOPE's native grid extends to thermal wavelengths (up to 50000nm, see note
# at the top of this script) -- restrict to the optical range first, or the
# subsample below would spend most of its points on the physically-flat
# thermal tail instead of the actual reflectance signal.
refl.cols.optical <- colnames(refl)[wl <= 2500]
refl.cols.plot <- refl.cols.optical[seq(1, length(refl.cols.optical), by = max(1, length(refl.cols.optical) %/% 400))]
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
  ggtitle(sprintf("SCOPE TOC reflectance by Cab (n=%d)", nrow(LUT))) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

print(spectral_plot)
plot_path <- file.path(out_dir, "0-reflectance_spectra_byCab.png")
ggsave(plot_path, plot = spectral_plot, width = 14, height = 9, dpi = 300, units = "cm")
cat("Saved diagnostic plot to '", plot_path, "'\n", sep = "")

## ----------------------------------------------------------------------------
## 7. Diagnostic plot: SIF760 vs Vcmax25 -- is there a visible relationship?
## ----------------------------------------------------------------------------

sif_vcmax_plot <- ggplot(LUT, aes(x = Vcmax25, y = SIF760)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red", linetype = "dashed") +
  theme_bw() +
  labs(x = "Vcmax25 (umol/m2/s)", y = "SIF760 (W/m2/um/sr)",
       title = sprintf("SIF760 vs Vcmax25 (n=%d, r=%.2f)", nrow(LUT), cor(LUT$Vcmax25, LUT$SIF760))) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
print(sif_vcmax_plot)
ggsave(file.path(out_dir, "1-SIF760_vs_Vcmax25.png"), plot = sif_vcmax_plot, width = 14, height = 9, dpi = 300, units = "cm")

saveRDS(list(LUT = LUT, refl = refl, wl = wl, convolved = convolved,
             datasets = list(native = dataset_native, se2a = dataset_se2a),
             n.samples = n.samples, seed = seed),
        file.path(out_dir, "SCOPE-1-datasets.rds"))

dev.off()

cat("\nSaved to '", out_dir, "/SCOPE-1-datasets.rds'\n", sep = "")
