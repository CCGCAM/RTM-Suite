# ==============================================================================
# Per-band Sobol sensitivity indices (first-order Si and total Ti) for
# Sentinel-2A bands, using the sensobol package's proper quasi-random Sobol
# design (sobol_matrices() -> evaluate model -> sobol_indices()).
#
# Originally read a pre-simulated 20k-row LUT CSV that isn't bundled with
# this repo -- rewritten to build and simulate its own LUT via ToolsRTM
# instead, so it's fully self-contained and actually runs.
# ==============================================================================

rm(list = ls())

if (!require("sensobol")) { install.packages("sensobol"); require("sensobol") }
if (!require("ggplot2")) { install.packages("ggplot2"); require("ggplot2") }

if (!requireNamespace("ToolsRTM", quietly = TRUE)) {
  remotes::install_gitlab("caminoccg/toolsrtm")
}
library(ToolsRTM)

out_dir <- "outs/Sensibility"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
pdf(file.path(out_dir, "Rplots.pdf"))  # catches any stray plot()/print() call, keeps it out of Scripts/

## ----------------------------------------------------------------------------
## 1. Sobol design: quasi-random sample in [0,1], one column per parameter
## ----------------------------------------------------------------------------

params <- c("Cab", "Car", "Anth", "LMA", "Cbrown", "N", "LIDFa", "LAI")
N <- 100  # base sample size; total fourSAIL evaluations = N * (length(params) + 2) = 1000

set.seed(1234)
mat01 <- sobol_matrices(N = N, params = params)  # every column in [0,1]

# Rescale each [0,1] column to that trait's realistic range (same ranges used
# elsewhere in this repo's PROSAIL-based scripts, e.g. Scripts/Sensibility/v0/1-Sobol_indices.R).
ranges <- list(Cab = c(5, 90), Car = c(0.5, 25), Anth = c(0, 10), LMA = c(0.0005, 0.03),
               Cbrown = c(0, 0.5), N = c(1.3, 2.5), LIDFa = c(30, 80), LAI = c(0.5, 7))

# Baseline LUT from getLUT() first, for every column fourSAIL needs that isn't
# a Sobol parameter (LIDFb, TypeLidf, tts, tto, psi, hspot, etc.) -- building
# the LUT purely from the Sobol matrix (as an earlier draft of this script
# did) leaves those columns missing entirely and foursail() errors.
LUT <- as.data.frame(getLUT(inputs = ToolsRTM::inputsPROSAIL, nLUT = nrow(mat01), setseed = 1234))
for (p in params) LUT[[p]] <- ranges[[p]][1] + mat01[, p] * diff(ranges[[p]])

# Fixed defaults for everything fourSAIL/PROSPECT-PRO needs that isn't a Sobol parameter.
LUT$EWT <- 0.02; LUT$Prot <- 0.0015; LUT$CBC <- 0.008
LUT$hspot <- 0.01; LUT$tts <- 30; LUT$tto <- 0; LUT$psi <- 0
LUT$Cs <- 0; LUT$fqe <- 0.01; LUT$Cx <- 0
LUT$cell.d <- 40; LUT$inter.c <- 0.045; LUT$baseline.abs <- 0.0006
LUT$leaf.thick <- 1.6; LUT$albino.abs <- 0; LUT$lign.cell <- 2; LUT$Nitrogen <- 1

n.samples <- nrow(LUT)
cat(sprintf("=== Sobol design: %d parameters, N=%d, %d total simulations ===\n", length(params), N, n.samples))

## ----------------------------------------------------------------------------
## 2. Soil + run fourSAIL for every design row
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
cat(sprintf("Simulated %d x %d reflectance matrix\n", nrow(refl), ncol(refl)))

## ----------------------------------------------------------------------------
## 3. Convolve to Sentinel-2A
## ----------------------------------------------------------------------------

refl_X <- as.data.frame(refl)
colnames(refl_X) <- paste0("X", wl)
refl_X <- cbind(id = seq_len(nrow(refl_X)), refl_X)
se2a <- suppressMessages(get.spectra.convolved(rfl = refl_X, sensor = "Sentinel2a", plot.spectra = FALSE))
names(se2a) <- c("id", "B1","B2","B3","B4","B5","B6","B7","B8","B8A","B9","B10","B11","B12")
bands <- c("B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12")
cat(sprintf("Convolved to Sentinel-2A: %d bands\n", length(bands)))

## ----------------------------------------------------------------------------
## 4. Sobol indices per band
## ----------------------------------------------------------------------------

indices_df <- data.frame(Banda = character(), Parametro = character(), Indice = numeric(), Tipo = character())
for (banda in bands) {
  Y <- se2a[[banda]]
  ind <- sobol_indices(Y = Y, N = N, params = params)
  indices_df <- rbind(indices_df, data.frame(Banda = banda, Parametro = ind$results$parameters,
                                              Indice = ind$results$original, Tipo = ind$results$sensitivity))
}
write.csv(indices_df, file.path(out_dir, "2-Sobol_perband_indices.csv"), row.names = FALSE)

## ----------------------------------------------------------------------------
## 5. Plots: first-order (Si) and total (Ti) indices, stacked by band
## ----------------------------------------------------------------------------

indices_df$Banda <- factor(indices_df$Banda, levels = bands)

p_si <- ggplot(subset(indices_df, Tipo == "Si"), aes(x = Banda, y = Indice, fill = Parametro)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "First-order Sobol indices (Si) by Sentinel-2A band", x = "Band", y = "Sobol index", fill = "Trait") +
  theme_bw()
print(p_si)
ggsave(file.path(out_dir, "2-Sobol_Si_byband.png"), plot = p_si, width = 18, height = 10, dpi = 300, units = "cm")

p_ti <- ggplot(subset(indices_df, Tipo == "Ti"), aes(x = Banda, y = Indice, fill = Parametro)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Total Sobol indices (Ti) by Sentinel-2A band", x = "Band", y = "Sobol index", fill = "Trait") +
  theme_bw()
print(p_ti)
ggsave(file.path(out_dir, "2-Sobol_Ti_byband.png"), plot = p_ti, width = 18, height = 10, dpi = 300, units = "cm")

dev.off()

cat("\nSaved to '", out_dir, "/': 2-Sobol_perband_indices.csv, 2-Sobol_Si_byband.png, 2-Sobol_Ti_byband.png\n", sep = "")
