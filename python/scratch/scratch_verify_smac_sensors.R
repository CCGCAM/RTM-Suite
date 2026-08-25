## Verify the newly-exported SMAC sensors against real ToolsRTM::get.smac() calls.
## Run from the repo root: Rscript python/scratch/scratch_verify_smac_sensors.R
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
devtools::load_all(file.path(root, "ToolsRTM/R"), quiet = TRUE)

outdir <- file.path(root, "python/toolsrtm/tests/refdata")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sensors <- list(landsat8 = LANDSAT8.OLI, sentinel3a = Sentinel3A.OLCI,
                modis = TerraAqua.MODIS, sentinel2b = Sentinel2B.MSI)

inputLUT <- data.frame(tts = 30, tto = 0, psi = 0, Pa = 1013.25, aot550 = 0.2,
                        uo3 = 0.35, uh2o = 2.0, alt_m = 0, Pa0 = 1013.25, DOY = 180, FWHM = 10)

for (slug in names(sensors)) {
  sensor <- sensors[[slug]]
  atmo <- get.smac(inputLUT, sensor)
  df <- data.frame(Ta_ss = atmo$Ta_ss, Ta_sd = atmo$Ta_sd, Ta_oo = atmo$Ta_oo, Ta_do = atmo$Ta_do,
                    Ta_s = atmo$Ta_s, Ta_o = atmo$Ta_o, Tg = atmo$Tg, Ra_dd = atmo$Ra_dd, Ra_so = atmo$Ra_so)
  write.csv(df, file.path(outdir, paste0("ref_smac_", slug, ".csv")), row.names = FALSE)
  cat("Exported ref_smac_", slug, ".csv (", nrow(df), " bands)\n", sep = "")
}
