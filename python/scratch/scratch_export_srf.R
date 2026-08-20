## Export PRISMA/Sentinel-2A/Sentinel-2B plain SRF tables (+ PRISMA's bundled
## FWHM table) as Python package data, matching the existing smac_*.csv
## convention in python/toolsrtm/src/toolsrtm/data/ -- closes the Python
## parity gap for get.spectral.convolution.srf() (ToolsRTM's new shared
## plain-SRF-table convolution function, see
## ToolsRTM/R/get.spectral.convolution.srf.R).
## Repo root, derived from this script's own location (python/scratch/<file>.R)
## rather than hardcoded, so it works regardless of who runs it or from where.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
outdir <- file.path(root, "python/toolsrtm/src/toolsrtm/data")

suppressPackageStartupMessages(library(ToolsRTM))

## ---- PRISMA: 234 hyperspectral bands, no individually meaningful names ----
prisma_band_cols <- setdiff(colnames(srf.prisma), "wavelength")
srf_prisma_out <- data.frame(wl = srf.prisma$wavelength)
for (i in seq_along(prisma_band_cols)) {
  srf_prisma_out[[paste0("band", i)]] <- srf.prisma[[prisma_band_cols[i]]]
}
write.csv(srf_prisma_out, file.path(outdir, "srf_prisma.csv"), row.names = FALSE)

write.csv(data.frame(wl = fwhm.prisma$wavelength, fwhm = fwhm.prisma$fwhm, qb = fwhm.prisma$QB),
          file.path(outdir, "fwhm_prisma.csv"), row.names = FALSE)

## ---- Sentinel-2A/B: 13 named bands (B1..B8,B8A,B9..B12) -- keep the real
## band names, they're meaningful (unlike PRISMA's 234 generic channels).
export_s2 <- function(srf_table, prefix, filename) {
  band_cols <- setdiff(colnames(srf_table), "SR_WL")
  out <- data.frame(wl = srf_table$SR_WL)
  for (col in band_cols) {
    band_name <- sub(paste0("^", prefix, "_SR_AV_"), "", col)  # e.g. "S2A_SR_AV_B8A" -> "B8A"
    out[[band_name]] <- srf_table[[col]]
  }
  write.csv(out, file.path(outdir, filename), row.names = FALSE)
}
export_s2(srf.sentinel2a, "S2A", "srf_sentinel2a.csv")
export_s2(srf.sentinel2b, "S2B", "srf_sentinel2b.csv")

## ---- Nominal sensor characteristics (no measured SRF, just center/FWHM
## or center/band-edges) for get.spectral.convolution.gaussian() /
## toolsrtm.srf.spectral_convolution_gaussian() -- EnMAP (center+fwhm
## already given) and the multi-sensor sensor.characteristics catalog
## (ALI/Hyperion/Landsat4-8/MODIS/Quickbird/RapidEye/Sentinel2a-b/
## WorldView2-4/2-8, which ship band edges lb/ub instead of fwhm).
write.csv(data.frame(sensor = "EnMAP", channel = EnMap.characteristics$channel,
                      center = EnMap.characteristics$center, fwhm = EnMap.characteristics$fwhm),
          file.path(outdir, "enmap_characteristics.csv"), row.names = FALSE)
write.csv(sensor.characteristics, file.path(outdir, "sensor_characteristics.csv"), row.names = FALSE)

cat("SRF export done:", paste(list.files(outdir, pattern = "^(srf|fwhm|enmap|sensor)_"), collapse = ", "), "\n")
