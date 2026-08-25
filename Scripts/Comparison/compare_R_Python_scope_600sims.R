## Cross-language robustness comparison: SCOPEinR (R) vs scopeinpython (Python).
## 600 simulations, one shared LUT, default (lite=1) SCOPE options.
## This script builds the LUT and the R-side reference TOC reflectance; a
## companion Python script (compare_R_Python_scope_600sims.py) reads the
## same LUT and computes the Python-side reflectance for direct comparison.

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
outdir <- file.path(root, "Scripts/Comparison/_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressMessages({
  devtools::load_all(file.path(root, "ToolsRTM"), quiet = TRUE)
  devtools::load_all(file.path(root, "SCOPEinR"), quiet = TRUE)
})

path_input <- system.file("input", package = "SCOPEinR")
scope_options <- read.table(file.path(path_input, "setoptions.csv"), header = TRUE, sep = ",")
inputLUT <- read.table(file.path(path_input, "inputs_SCOPE.csv"), header = TRUE, sep = ",")

n_samples <- 600
set.seed(7)
LUT <- as.data.frame(SCOPEinR::getLUT.SCOPE(inputLUT = inputLUT, nLUT = n_samples))
write.csv(LUT, file.path(outdir, "shared_LUT_scope_600.csv"), row.names = FALSE)

spectral <- SCOPEinR::get.spectra.SCOPE()
wl <- spectral[['wlS']]
write.csv(data.frame(wl = wl), file.path(outdir, "scope_wl_grid.csv"), row.names = FALSE)
refl_R <- matrix(NA_real_, nrow = n_samples, ncol = length(wl))
Actot_R <- rep(NA_real_, n_samples)

t0 <- Sys.time()
for (i in seq_len(n_samples)) {
  r <- tryCatch(
    SCOPEinR::get.SCOPE(LUT = LUT[i, , drop = FALSE], options.SCOPE = scope_options,
                         optipar = SCOPEinR::optipar2021.Pro.CX, get.outputs = "ALL", get.plots = FALSE),
    error = function(e) NULL
  )
  if (!is.null(r)) {
    refl_i <- r[[1]]$data.rad$refl
    if (!is.null(refl_i) && length(refl_i) == length(wl)) refl_R[i, ] <- refl_i
    Actot_R[i] <- r[[1]]$data.fluxes$Actot
  }
  if (i %% 100 == 0) cat("R side:", i, "/", n_samples, "done,", round(as.numeric(difftime(Sys.time(), t0, units="secs"))), "s elapsed\n")
}

# A real, structural set of wavelengths (scattered, non-contiguous, ~109 of
# 2162, matching known atmospheric water-vapor absorption windows around
# 1400/1900/2500nm) come back NA from every simulation -- not a per-row
# failure. Row validity is judged on Actot (a scalar flux total) instead of
# full-spectrum completeness; the comparison script masks the NA wavelength
# columns instead of dropping otherwise-valid rows.
ok <- !is.na(Actot_R)
cat("R side final: ", sum(ok), "/", n_samples, " simulations produced valid output\n", sep = "")

out_R <- as.data.frame(refl_R)
colnames(out_R) <- paste0("wl", wl)
out_R$row <- seq_len(n_samples)
out_R$ok <- ok
out_R$Actot <- Actot_R
write.csv(out_R, file.path(outdir, "refl_R_scope_600.csv"), row.names = FALSE)
cat("Wrote shared_LUT_scope_600.csv and refl_R_scope_600.csv\n")
