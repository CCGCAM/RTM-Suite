## One-off: generate the REAL R output figure for the SCOPE code example in
## docs/index.html -- a single get.SCOPE() run, plain reflectance (no SIF).
## Run from the repo root:
##   Rscript python/scratch/scratch_export_scope_figure_r.R
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
outdir <- file.path(root, "assets/examples")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

devtools::load_all(file.path(root, "ToolsRTM/R"), quiet = TRUE)
devtools::load_all(file.path(root, "SCOPEinR/R"), quiet = TRUE)

opts <- read.table(system.file("input", "setoptions.csv", package = "SCOPEinR"),
                    header = TRUE, sep = ",")
lut <- read.table(system.file("input", "LUT_input.csv", package = "SCOPEinR"),
                   header = TRUE, sep = ",")

scope <- get.SCOPE(
  LUT = lut[1, ], options.SCOPE = opts, optipar = optipar2021.Pro.CX,
  leaf.model = "fluspect-CX", canopy.model = "fourSAIL",
  get.outputs = "ALL", get.plots = FALSE)

res <- scope[[1]]
refl <- res$data.rad$refl
wl <- res$data.spectral$wlS
optical <- wl <= 2400

png(file.path(outdir, "r_scope_refl.png"), width = 900, height = 560, res = 130)
par(mar = c(4, 4, 1, 1))
plot(wl[optical], refl[optical], type = "l", col = "#0b3d59", lwd = 2,
     xlab = "Wavelength (nm)", ylab = "Reflectance")
grid(col = "grey85")
dev.off()
cat("Saved r_scope_refl.png\n")
