# ==============================================================================
# Johnson relative-importance (dominance analysis) indices for Sentinel-2A
# bands -- how much each trait contributes to explaining variance in each
# band's reflectance, via sensitivity::johnson().
#
# Originally read a pre-simulated LUT CSV that isn't bundled with this repo --
# rewritten to build and simulate its own LUT via ToolsRTM instead, so it's
# fully self-contained and actually runs. Unlike the Sobol script next to
# this one, johnson() is a regression-based method and doesn't need a special
# quasi-random design -- an ordinary random LUT (getLUT()) is enough.
# ==============================================================================

rm(list = ls())

if (!require("sensitivity")) { install.packages("sensitivity"); require("sensitivity") }
if (!require("ggplot2")) { install.packages("ggplot2"); require("ggplot2") }
if (!require("dplyr")) { install.packages("dplyr"); require("dplyr") }

if (!requireNamespace("ToolsRTM", quietly = TRUE)) {
  remotes::install_gitlab("caminoccg/toolsrtm")
}
library(ToolsRTM)

out_dir <- "outs/Sensibility"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
pdf(file.path(out_dir, "Rplots.pdf"))  # catches any stray plot()/print() call, keeps it out of Scripts/

## ----------------------------------------------------------------------------
## 1. Build a random LUT (always 100 -- johnson() doesn't need a Sobol design)
## ----------------------------------------------------------------------------

n.samples <- 100
inputs <- ToolsRTM::inputsPROSAIL
LUT <- as.data.frame(getLUT(inputs = inputs, nLUT = n.samples, setseed = 1234))

LUT$Cs <- 0; LUT$fqe <- 0.01; LUT$Cx <- 0
LUT$cell.d <- 40; LUT$inter.c <- 0.045; LUT$baseline.abs <- 0.0006
LUT$leaf.thick <- 1.6; LUT$albino.abs <- 0; LUT$lign.cell <- 2; LUT$Nitrogen <- 1

cat(sprintf("=== Johnson relative importance: %d simulations ===\n", n.samples))

## ----------------------------------------------------------------------------
## 2. Soil + run fourSAIL
## ----------------------------------------------------------------------------

data <- ToolsRTM::dataSpec_PDB
Rsoil1 <- data[, 11]; Rsoil2 <- data[, 12]
rsoil <- 0.5 * Rsoil1 + 0.5 * Rsoil2
wl <- 400:2500

refl <- t(sapply(seq_len(n.samples), function(i) {
  sim <- suppressMessages(foursail(inputLUT = LUT[i, ], rsoil = rsoil, LeafModel = "PROSPECT-PRO"))
  ToolsRTM::Compute_BRF(rdot = sim$rdot, rsot = sim$rsot, tts = LUT$tts[i], data.light = ToolsRTM::dataSpec_PDB)
}))
colnames(refl) <- paste0("R.", wl)

## ----------------------------------------------------------------------------
## 3. Convolve to Sentinel-2A
## ----------------------------------------------------------------------------

refl_X <- as.data.frame(refl)
colnames(refl_X) <- paste0("X", wl)
refl_X <- cbind(id = seq_len(nrow(refl_X)), refl_X)
se2a <- suppressMessages(get.spectra.convolved(rfl = refl_X, sensor = "Sentinel2a", plot.spectra = FALSE))
names(se2a) <- c("id", "B1","B2","B3","B4","B5","B6","B7","B8","B8A","B9","B10","B11","B12")
bands <- c("B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12")

## ----------------------------------------------------------------------------
## 4. Johnson indices per band
## ----------------------------------------------------------------------------

candidate_traits <- c("Cab", "Car", "Anth", "LMA", "Cbrown", "N", "LIDFa", "LAI", "EWT")
X <- LUT[, candidate_traits]
# Drop any trait getLUT() didn't actually vary (e.g. LMA is constant 0 with
# PROSPECT-PRO, which uses CBC/Prot for dry matter instead) -- johnson()'s
# correlation-matrix step fails outright on a zero-variance column.
zero_var <- sapply(X, function(v) stats::sd(v) == 0)
if (any(zero_var)) {
  cat("Dropping constant (zero-variance) trait(s) from the predictor set: ",
      paste(names(X)[zero_var], collapse = ", "), "\n", sep = "")
  X <- X[, !zero_var, drop = FALSE]
}

indices_df <- data.frame(Banda = character(), Parametro = character(), Indice = numeric())
for (banda in bands) {
  y <- se2a[[banda]]
  x <- johnson(X, y)
  indices_df <- rbind(indices_df, data.frame(Banda = banda, Parametro = rownames(x$johnson), Indice = x$johnson$original))
}

indices_df$Banda <- factor(indices_df$Banda, levels = bands)
write.csv(indices_df, file.path(out_dir, "3-Johnson_indices.csv"), row.names = FALSE)

## ----------------------------------------------------------------------------
## 5. Plot
## ----------------------------------------------------------------------------

p <- ggplot(indices_df, aes(x = Banda, y = Indice, fill = Parametro)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Johnson relative-importance indices by Sentinel-2A band",
       x = "Band", y = "Johnson index", fill = "Trait") +
  theme_bw()
print(p)
ggsave(file.path(out_dir, "3-Johnson_byband.png"), plot = p, width = 18, height = 10, dpi = 300, units = "cm")

dev.off()

cat("\nSaved to '", out_dir, "/': 3-Johnson_indices.csv, 3-Johnson_byband.png\n", sep = "")
