## Cross-language robustness comparison: ToolsRTM (R) vs toolsrtm (Python)
## 600 simulations, one shared LUT, PROSPECT-D + fourSAIL.
## This script builds the LUT and the R-side reference reflectance; a
## companion Python script (compare_R_Python_toolsrtm_600sims.py) reads the
## same LUT and computes the Python-side reflectance for direct comparison.

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
outdir <- file.path(root, "Scripts/Comparison/_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressMessages(devtools::load_all(file.path(root, "ToolsRTM"), quiet = TRUE))

n_samples <- 600
set.seed(42)
LUT <- as.data.frame(ToolsRTM::getLUT(inputs = ToolsRTM::inputsPROSAIL, nLUT = n_samples, setseed = 42))
write.csv(LUT, file.path(outdir, "shared_LUT_600.csv"), row.names = FALSE)

rsoil <- rep(0.15, 2101)
wl <- 400:2500

refl_R <- matrix(NA_real_, nrow = n_samples, ncol = length(wl))
for (i in seq_len(n_samples)) {
  sim <- tryCatch(
    suppressMessages(ToolsRTM::simulate_RTM(inputLUT = LUT[i, , drop = FALSE], rsoil = rsoil,
                                             leaf.model = "PROSPECT-D", canopy.model = "fourSAIL")),
    error = function(e) NULL
  )
  if (!is.null(sim) && length(sim$rsot) == length(wl)) refl_R[i, ] <- sim$rsot
}

ok <- stats::complete.cases(refl_R)
cat("R side: ", sum(ok), "/", n_samples, " simulations produced a valid reflectance spectrum\n", sep = "")

out_R <- as.data.frame(refl_R)
colnames(out_R) <- paste0("wl", wl)
out_R$row <- seq_len(n_samples)
out_R$ok <- ok
write.csv(out_R, file.path(outdir, "refl_R_600.csv"), row.names = FALSE)
cat("Wrote", file.path(outdir, "shared_LUT_600.csv"), "and refl_R_600.csv\n")
