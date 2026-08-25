## Cross-language robustness comparison: SPART (soil-plant-atmosphere) and
## MARMIT (wet-soil reflectance, both versions) -- ToolsRTM (R) vs
## toolsrtm (Python). Same shared-LUT approach as the other
## Scripts/Comparison/*.R scripts.

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
outdir <- file.path(root, "Scripts/Comparison/_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressMessages(devtools::load_all(file.path(root, "ToolsRTM"), quiet = TRUE))

## ---------------------------------------------------------------------
## 1. SPART: TOC + TOA reflectance, Sentinel-2A bands, realistic fixed
##    atmosphere (matching t03-spart.Rmd's own choice to avoid the known
##    unphysical-atmosphere edge case at getLUT()'s wide default ranges).
## ---------------------------------------------------------------------
n_spart <- 100
LUT_spart <- as.data.frame(getLUT(inputs = ToolsRTM::inputsSPART, nLUT = n_spart, setseed = 21))
LUT_spart$Cs <- 0; LUT_spart$fqe <- 0.01; LUT_spart$Cx <- 0
LUT_spart$Pa <- 1000; LUT_spart$aot550 <- 0.3246; LUT_spart$uo3 <- 0.3480; LUT_spart$uh2o <- 1.4116
LUT_spart$alt_m <- 0; LUT_spart$Pa0 <- 1000
LUT_spart$BSMBrightness <- 0.5; LUT_spart$BSMlat <- 25; LUT_spart$BSMlon <- 45; LUT_spart$SMp <- 15

rfl_toc <- matrix(NA_real_, n_spart, 13)
rfl_toa <- matrix(NA_real_, n_spart, 13)
wave_bands <- NULL
for (i in seq_len(n_spart)) {
  sim <- tryCatch(suppressWarnings(SPART(inputLUT = LUT_spart[i, , drop = FALSE], CanopyModel = "fourSAIL",
                                          LeafModel = "PROSPECT-PRO", sensor.i = ToolsRTM::Sentinel2A.MSI,
                                          rsoil = NULL, get.plots = FALSE)),
                  error = function(e) NULL)
  if (!is.null(sim) && all(is.finite(sim$output$rfl.toc.BRDF)) && all(is.finite(sim$output$rfl.toa))) {
    if (is.null(wave_bands)) wave_bands <- sim$output$wave
    if (length(sim$output$rfl.toc.BRDF) == 13) {
      rfl_toc[i, ] <- sim$output$rfl.toc.BRDF
      rfl_toa[i, ] <- sim$output$rfl.toa
    }
  }
}
ok_spart <- stats::complete.cases(rfl_toc) & stats::complete.cases(rfl_toa)
cat("R side, SPART:", sum(ok_spart), "/", n_spart, "valid\n")

write.csv(LUT_spart, file.path(outdir, "lut_SPART.csv"), row.names = FALSE)
out_toc <- as.data.frame(rfl_toc); colnames(out_toc) <- paste0("band", seq_len(13))
out_toc$row <- seq_len(n_spart); out_toc$ok <- ok_spart
write.csv(out_toc, file.path(outdir, "refl_R_SPART_toc.csv"), row.names = FALSE)
out_toa <- as.data.frame(rfl_toa); colnames(out_toa) <- paste0("band", seq_len(13))
out_toa$row <- seq_len(n_spart); out_toa$ok <- ok_spart
write.csv(out_toa, file.path(outdir, "refl_R_SPART_toa.csv"), row.names = FALSE)
write.csv(data.frame(band = seq_len(13), wave_nm = wave_bands), file.path(outdir, "SPART_bands.csv"), row.names = FALSE)

## ---------------------------------------------------------------------
## 2. MARMIT: both versions (marmit1, marmit2), varying soil id/L/eps,
##    against the bundled Bablet_2016 database (17 soil IDs).
## ---------------------------------------------------------------------
set.seed(22)
n_marmit <- 60
marmit_grid <- data.frame(
  soil_id = sample(1:17, n_marmit, replace = TRUE),
  L = runif(n_marmit, 0.01, 0.15),
  eps = runif(n_marmit, 0.1, 0.9)
)
n_i <- 1.53; k_i <- 0.001; d_i <- 0.0005
wl_out <- 400:2500

run_marmit_version <- function(version) {
  wet <- matrix(NA_real_, n_marmit, length(wl_out))
  dry <- matrix(NA_real_, n_marmit, length(wl_out))
  for (i in seq_len(n_marmit)) {
    soil <- tryCatch(get.marmit.rsoil(database = "Bablet_2016", id = marmit_grid$soil_id[i], version = version,
                                       L = marmit_grid$L[i], eps = marmit_grid$eps[i],
                                       n_i = n_i, k_i = k_i, d_i = d_i, wl.out = wl_out),
                      error = function(e) NULL)
    if (!is.null(soil) && length(soil$rsoil.wet) == length(wl_out)) {
      wet[i, ] <- soil$rsoil.wet
      dry[i, ] <- soil$rsoil.dry
    }
  }
  ok <- stats::complete.cases(wet)
  cat("R side, MARMIT (", version, "):", sum(ok), "/", n_marmit, "valid\n")
  list(wet = wet, dry = dry, ok = ok)
}

write.csv(marmit_grid, file.path(outdir, "lut_MARMIT.csv"), row.names = FALSE)
for (version in c("marmit1", "marmit2")) {
  res <- run_marmit_version(version)
  out_wet <- as.data.frame(res$wet); colnames(out_wet) <- paste0("wl", wl_out)
  out_wet$row <- seq_len(n_marmit); out_wet$ok <- res$ok
  write.csv(out_wet, file.path(outdir, paste0("refl_R_MARMIT_", version, "_wet.csv")), row.names = FALSE)
}

cat("Done. SPART + MARMIT (both versions) written to", outdir, "\n")
