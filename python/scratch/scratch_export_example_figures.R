## One-off: generate the REAL output figures shown in docs/index.html's
## "Code examples" section, from the exact code in each card (not generic
## stock images). Run from the repo root:
##   Rscript python/scratch/scratch_export_example_figures.R
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
outdir <- file.path(root, "assets/examples")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

devtools::load_all(file.path(root, "ToolsRTM/R"), quiet = TRUE)

## ---- 1. Simulate a vegetation spectrum (PROSPECT-PRO + fourSAIL) ----
LUT <- as.data.frame(getLUT(inputs = ToolsRTM::inputsPROSAIL,
                            nLUT = 10, setseed = 123))
rsoil <- rep(0.15, 2101)

sim <- simulate_RTM(
    inputLUT    = LUT[1, ],
    rsoil       = rsoil,
    leaf.model  = "PROSPECT-PRO",
    canopy.model = "fourSAIL")

wl <- 400:2500
png(file.path(outdir, "r_sim_foursail.png"), width = 900, height = 560, res = 130)
par(mar = c(4, 4, 1, 1))
plot(wl, sim$rsot, type = "l", col = "#0b3d59", lwd = 2,
     xlab = "Wavelength (nm)", ylab = "Reflectance",
     ylim = c(0, max(sim$rsot) * 1.15))
grid(col = "grey85")
dev.off()
cat("Saved r_sim_foursail.png\n")

## ---- 2. Sensitivity to chlorophyll (Cab), fourSAIL ----
cab_values <- c(10, 25, 40, 55, 70)
palette_cab <- colorRampPalette(c("#f2e9b3", "#2bb3a3", "#0b3d59"))(length(cab_values))

png(file.path(outdir, "r_sensitivity_cab.png"), width = 900, height = 560, res = 130)
par(mar = c(4, 4, 1, 8), xpd = TRUE)
plot(NULL, xlim = c(400, 2500), ylim = c(0, 0.5),
     xlab = "Wavelength (nm)", ylab = "Reflectance")
grid(col = "grey85")
for (i in seq_along(cab_values)) {
  row <- LUT[1, ]
  row$Cab <- cab_values[i]
  sim_i <- simulate_RTM(inputLUT = row, rsoil = rsoil, leaf.model = "PROSPECT-PRO", canopy.model = "fourSAIL")
  lines(wl, sim_i$rsot, col = palette_cab[i], lwd = 2)
}
legend(2560, 0.5, legend = paste0("Cab=", cab_values), col = palette_cab, lwd = 2, bty = "n", cex = 0.8)
dev.off()
cat("Saved r_sensitivity_cab.png\n")
