## One-off: generate the REAL R output figure for the "50 simulations" code
## example in docs/index.html. Run from the repo root:
##   Rscript python/scratch/scratch_export_50sims_figures.R
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
outdir <- file.path(root, "assets/examples")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

devtools::load_all(file.path(root, "ToolsRTM/R"), quiet = TRUE)

LUT <- as.data.frame(getLUT(inputs = ToolsRTM::inputsPROSAIL, nLUT = 50, setseed = 42))
rsoil <- rep(0.15, 2101)
wl <- 400:2500

palette_lai <- colorRampPalette(c("#f2e9b3", "#2bb3a3", "#0b3d59"))(50)
lai_order <- order(LUT$LAI)

png(file.path(outdir, "r_50sims.png"), width = 900, height = 560, res = 130)
par(mar = c(4, 4, 1, 1))
plot(NULL, xlim = c(400, 2500), ylim = c(0, 0.5),
     xlab = "Wavelength (nm)", ylab = "Reflectance")
grid(col = "grey85")
for (i in seq_len(50)) {
  row <- LUT[i, ]
  sim <- simulate_RTM(inputLUT = row, rsoil = rsoil, leaf.model = "PROSPECT-PRO", canopy.model = "fourSAIL")
  color_rank <- which(lai_order == i)
  lines(wl, sim$rsot, col = palette_lai[color_rank], lwd = 1)
}
dev.off()
cat("Saved r_50sims.png\n")
