# Generates the static SCOPEinR figures shown under the "How in R" / "How in
# Python -- SCOPEinR" tutorial row on the RTM-Suite landing page
# (docs/index.html #tutorials). Mirrors generate_tutorial_figures.R's own
# ToolsRTM/PROSAIL figures, but for the SCOPE energy-balance/fluorescence
# pipeline -- reuses How-in-R-SCOPEinR.Rmd's own verified simulation code.
#
# Run with this file's own directory (Apps/RTMs/) as the working directory.

suppressPackageStartupMessages({
  library(ToolsRTM); library(SCOPEinR); library(ggplot2); library(patchwork); library(randomForest)
})

dir.create("../../Tutorials/images", showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path("../../Tutorials/images", f)
set.seed(42)

theme_scope <- theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 12))

## ----------------------------------------------------------------------------
## 1. One simulation: reflectance components + fluorescence spectrum
## ----------------------------------------------------------------------------

table.with.opts <- read.table(system.file("input", "setoptions.csv", package = "SCOPEinR"), header = TRUE, sep = ",")
Table.LUT <- read.table(system.file("input", "LUT_input.csv", package = "SCOPEinR"), header = TRUE, sep = ",")

db.sim <- get.SCOPE(LUT = Table.LUT[1, ], options.SCOPE = table.with.opts,
                     optipar = SCOPEinR::optipar2021.Pro.CX, leaf.model = "fluspect-CX",
                     canopy.model = "fourSAIL", get.outputs = "ALL", get.plots = FALSE)
res <- db.sim[[1]]

wl_optical <- 400:2400
n <- length(wl_optical)
rad <- res$data.rad

refl_df <- data.frame(
  Wavelength = rep(wl_optical, 3),
  Value = c(rad$refl[1:n], rad$rdd[1:n], rad$rso[1:n]),
  Component = rep(c("refl (total)", "rdd (bi-hemispherical)", "rso (bidirectional)"), each = n)
)
p_refl <- ggplot(refl_df, aes(Wavelength, Value, color = Component)) +
  geom_line(linewidth = 0.8) + theme_scope +
  labs(title = "TOC reflectance components", x = "Wavelength (nm)", y = "Reflectance") +
  theme(legend.position = "bottom", legend.title = element_blank())

wlF <- 640:850
fluor_df <- data.frame(Wavelength = wlF, Fluorescence = rad$LoF_)
p_fluor <- ggplot(fluor_df, aes(Wavelength, Fluorescence)) +
  geom_line(color = "#B2182B", linewidth = 1) + theme_scope +
  labs(title = sprintf("TOC fluorescence (F685=%.2f, F740=%.2f)", rad$F685, rad$F740),
       x = "Wavelength (nm)", y = "Fluorescence radiance (mW m-2 nm-1 sr-1)")

p1 <- p_refl + p_fluor
ggsave(out("scope_reflectance_fluorescence.png"), p1, width = 11, height = 4.5, dpi = 110)

## ----------------------------------------------------------------------------
## 2. 100-row LUT: Actot vs Vcmax25, SIF vs Cab
## ----------------------------------------------------------------------------

inputLUT <- read.table(system.file("input", "inputs_SCOPE.csv", package = "SCOPEinR"), header = TRUE, sep = ",")
N_SAMPLES <- 100
Table.LUT.many <- getLUT.SCOPE(inputLUT = inputLUT, nLUT = N_SAMPLES)
db.sim.many <- get.SCOPE(LUT = Table.LUT.many, n.LUT = N_SAMPLES, options.SCOPE = table.with.opts,
                          optipar = SCOPEinR::optipar2021.Pro.CX, leaf.model = "fluspect-CX",
                          canopy.model = "fourSAIL", get.outputs = "ALL", get.plots = FALSE)

Actot_vals <- sapply(db.sim.many, function(r) r$data.fluxes$Actot)
EoutF_vals <- sapply(db.sim.many, function(r) r$data.rad$EoutF)

df_actot <- data.frame(Vcmax25 = Table.LUT.many$Vcmax25, Actot = Actot_vals)
df_sif   <- data.frame(Cab = Table.LUT.many$Cab, EoutF = EoutF_vals)

p_actot <- ggplot(df_actot, aes(Vcmax25, Actot)) +
  geom_point(color = "#2166AC", size = 1.8) + theme_scope +
  labs(title = "Actot vs Vcmax25", x = "Vcmax25 (umol m-2 s-1)", y = "Actot (umol m-2 s-1)")
p_sif <- ggplot(df_sif, aes(Cab, EoutF)) +
  geom_point(color = "#B2182B", size = 1.8) + theme_scope +
  labs(title = "SIF vs Cab", x = "Cab", y = "EoutF (W/m2)")

p2 <- p_actot + p_sif
ggsave(out("scope_photosynthesis_sif.png"), p2, width = 11, height = 4.5, dpi = 110)

## ----------------------------------------------------------------------------
## 3. LAI retrieved from Sentinel-2A bands (SCOPE reflapp)
## ----------------------------------------------------------------------------

band_refl <- t(sapply(db.sim.many, function(r) {
  rfl_i <- r$data.rad$reflapp[1:n]
  bad <- !is.finite(rfl_i)
  if (any(bad)) rfl_i[bad] <- approx(wl_optical[!bad], rfl_i[!bad], xout = wl_optical[bad])$y
  df_i <- data.frame(wave = wl_optical, rfl = rfl_i)
  get.spectral.convolution.srf(df_i, ToolsRTM::srf.sentinel2a)$RFL
}))
colnames(band_refl) <- paste0("B", seq_len(ncol(band_refl)))
ml_data <- data.frame(band_refl, LAI = Table.LUT.many$LAI)

train_idx <- sample(seq_len(N_SAMPLES), size = round(0.7 * N_SAMPLES))
rf <- randomForest(LAI ~ ., data = ml_data[train_idx, ], ntree = 300)
pred <- predict(rf, ml_data[-train_idx, ])
obs  <- ml_data$LAI[-train_idx]

df_inv <- data.frame(Observed = obs, Predicted = pred)
r2 <- 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
p3 <- ggplot(df_inv, aes(Observed, Predicted)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(color = "#2166AC", size = 2) + theme_scope +
  labs(title = sprintf("LAI from Sentinel-2A bands (R2=%.2f)", r2), x = "Observed LAI", y = "Predicted LAI")
ggsave(out("scope_lai_inversion.png"), p3, width = 6.5, height = 5.5, dpi = 110)

cat("Saved: scope_reflectance_fluorescence.png, scope_photosynthesis_sif.png, scope_lai_inversion.png\n")
