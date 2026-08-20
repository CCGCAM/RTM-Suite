# ==============================================================================
# SCOPEinR pipeline step 1: simulate N SCOPE runs, convolve the resulting TOC
# reflectance to real sensor bands, and compute vegetation indices -- the
# SCOPEinR equivalent of Scripts/Pipeline/1-simulate_LUT.R (ToolsRTM).
#
# get.SCOPE()'s leaf.model/canopy.model arguments are NOT functional (SCOPE
# always runs its own integral multi-layer Fluspect-Cx + RTMo, see
# Scripts/compare_SCOPE_models.R) -- there is nothing to configure there.
#
# Soil: get.SCOPE() sources soil reflectance either from a bundled file
# (soil_spectra/soil_scope.txt, 3 selectable spectra via the LUT's `spectrum`
# column) or from the BSM model (options.SCOPE row 5, options.soilspectrum),
# controlled by SCOPEinR/inst/input/setoptions.csv -- NOT via a passed-in
# rsoil vector like ToolsRTM::foursail()/SPART(). Feeding a MARMIT-generated
# soil into SCOPE the same way it was fed into ToolsRTM's models would need
# an analogous new `rsoil` parameter added to get.SCOPE()/get.SCOPE.parallel()
# (the same pattern already added to ToolsRTM::SPART() this session) -- not
# done here; this script uses SCOPE's own built-in soil options as-is.
#
# get.SCOPE.parallel() already has a working `parallel = TRUE/FALSE` switch
# (no ToolsRTM::get.inversion()-style bug here) -- exposed below as
# `use.parallel`.
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
library(ggplot2)  # needed for the diagnostic plot below

## ----------------------------------------------------------------------------
## User-configurable settings
## ----------------------------------------------------------------------------

n.samples    <- 100              # how many LUT rows / SCOPE simulations to run
use.parallel <- FALSE            # TRUE = get.SCOPE.parallel(parallel = TRUE); FALSE = sequential
seed         <- 1
sensors      <- c("Sentinel2a", "Sentinel2b", "PRISMA")

out_dir <- "outs/rtm_sims/pipeline_scope"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("=== Simulating %d SCOPE runs (parallel = %s) ===\n", n.samples, use.parallel))

## ----------------------------------------------------------------------------
## 1. Build the LUT and run SCOPE
## ----------------------------------------------------------------------------

opts <- read.table(system.file("input", "setoptions.csv", package = "SCOPEinR"), header = TRUE, sep = ",")
inputLUT <- read.table(system.file("input", "inputs_SCOPE.csv", package = "SCOPEinR"), header = TRUE, sep = ",")
LUT <- getLUT.SCOPE(inputLUT = inputLUT, nLUT = n.samples, setseed = seed)

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

saveRDS(list(LUT = LUT, refl = refl, wl = wl, convolved = convolved,
             datasets = list(native = dataset_native, se2a = dataset_se2a),
             n.samples = n.samples, seed = seed),
        file.path(out_dir, "SCOPE-1-datasets.rds"))

cat("\nSaved to '", out_dir, "/SCOPE-1-datasets.rds'\n", sep = "")
