## Export SMAC coefficient/SRF data for the 8 sensors not yet bundled in the
## Python port (Sentinel-2A was already exported earlier). Closes the "only
## Sentinel-2A bundled" gap documented in toolsrtm.smac's module docstring --
## the atmospheric-correction physics itself is already sensor-agnostic.
## Run from the repo root:
##   Rscript python/scratch/scratch_export_smac_sensors.R
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
outdir <- file.path(root, "python/toolsrtm/src/toolsrtm/data")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

devtools::load_all(file.path(root, "ToolsRTM/R"), quiet = TRUE)

## slug -> R sensor object name
sensors <- list(
  landsat4    = "LANDSAT4.TM",
  landsat5    = "LANDSAT5.TM",
  landsat7    = "LANDSAT7.ETM",
  landsat8    = "LANDSAT8.OLI",
  sentinel2b  = "Sentinel2B.MSI",
  sentinel3a  = "Sentinel3A.OLCI",
  sentinel3b  = "Sentinel3B.OLCI",
  modis       = "TerraAqua.MODIS"
)

for (slug in names(sensors)) {
  obj_name <- sensors[[slug]]
  sensor <- get(obj_name, envir = asNamespace("ToolsRTM"))
  nbands <- length(sensor$wl_smac)

  ## 1. Band centers (only wl_smac is actually used by the Python loader --
  ## center_wvl is a separate, sometimes differently-lengthed field on some
  ## sensors, e.g. LANDSAT4.TM, not needed here)
  write.csv(data.frame(wl_smac = as.vector(sensor$wl_smac)),
            file.path(outdir, paste0("smac_bands_", slug, ".csv")), row.names = FALSE)

  ## 2. SMAC coefficients: one row per coefficient name, one column per band.
  ## SMAC_coef is a list of 49 named vectors (one per coefficient), each
  ## length nbands -- do.call(rbind, ...) does NOT keep the list names as
  ## row names, they must be set explicitly.
  coef_mat <- do.call(rbind, sensor$SMAC_coef)
  rownames(coef_mat) <- names(sensor$SMAC_coef)
  colnames(coef_mat) <- paste0("V", seq_len(nbands))
  write.csv(as.data.frame(coef_mat), file.path(outdir, paste0("smac_coef_", slug, ".csv")))

  ## 3. SRF wavelengths + weights, NaN-padded rectangular (nsrf x nbands)
  wl_srf <- sensor$wl_srf_smac
  p_srf <- sensor$p_srf_smac
  colnames(wl_srf) <- paste0("band", seq_len(nbands))
  colnames(p_srf) <- paste0("band", seq_len(nbands))
  write.csv(as.data.frame(wl_srf), file.path(outdir, paste0("smac_srf_wl_", slug, ".csv")), row.names = FALSE, na = "nan")
  write.csv(as.data.frame(p_srf), file.path(outdir, paste0("smac_srf_weight_", slug, ".csv")), row.names = FALSE, na = "nan")

  cat("Exported", slug, "(", obj_name, "):", nbands, "bands\n")
}
