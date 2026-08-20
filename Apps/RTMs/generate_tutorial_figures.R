# Generates the static figures shown in the "How in R" / "How in Python"
# tutorial tabs of app.R (www/tutorial_*.png). Re-run this whenever the
# tutorial narrative's example numbers change. Not sourced by the app itself
# -- the app only displays the PNGs this script writes to www/.
#
# Mirrors the RTM tab's own workflow (fourSAIL + PROSPECT-D as the worked
# example -- the same 6 steps apply to any of the other 14 canopy x leaf
# model combinations the RTM tab supports, and to PROSAIL/PROSAIL-WithSatellite,
# which are just this workflow with the model fixed to fourSAIL+PROSPECT-PRO).

root <- "C:/Users/camin001/OneDrive - Wageningen University & Research/Workspace/0-RTM-Suite"
setwd(file.path(root, "Apps/RTMs"))
suppressPackageStartupMessages({
  library(ToolsRTM); library(ggplot2); library(patchwork); library(randomForest)
})
# Reuse convolve_to_sensor()/.sensor_choices without launching the app.
source("app.R", local = (tutorial_env <- new.env()))
convolve_to_sensor <- tutorial_env$convolve_to_sensor
theme_prosail <- tutorial_env$theme_prosail

dir.create("www", showWarnings = FALSE)
set.seed(42)

wavelength <- ToolsRTM::dataSpec_PDB[, 1]
rsoil <- 0.5 * ToolsRTM::dataSpec_PDB[, 11] + 0.5 * ToolsRTM::dataSpec_PDB[, 12]

baseline <- data.frame(
  LAI = 3, hspot = 0.01, LIDFa = -0.35, LIDFb = -0.15, TypeLidf = 1,
  tts = 30, tto = 0, psi = 0,
  N = 1.5, Cab = 40, Car = 8, Anth = 2, Cbrown = 0, EWT = 0.009, LMA = 0.009, alpha = 40,
  Prot = 0.0045, CBC = 0.005
)

run_one <- function(row, rsoil_vec = rsoil) {
  sail <- ToolsRTM::foursail(inputLUT = row, rsoil = rsoil_vec, LeafModel = "PROSPECT-D")
  ToolsRTM::Compute_BRF(rdot = sail$rdot, rsot = sail$rsot, tts = row$tts, data.light = ToolsRTM::dataSpec_PDB)
}

# ---- Step 1: a single simulation ------------------------------------------
t0 <- Sys.time()
refl_one <- run_one(baseline)
df1 <- data.frame(wavelength = wavelength, reflectance = refl_one)
p1 <- ggplot(df1, aes(wavelength, reflectance)) +
  geom_line(color = "#1b7837", linewidth = 1.3) +
  labs(x = "Wavelength (nm)", y = "Reflectance", title = "One simulation: fourSAIL + PROSPECT-D") +
  theme_prosail(legend_position = "none")
ggsave("www/tutorial_1sim.png", p1, width = 8, height = 4.5, dpi = 110)
cat("step 1 done,", round(as.numeric(Sys.time() - t0, units = "secs"), 2), "s\n")

# ---- Step 2: n_samples = 500, many simulations -----------------------------
t0 <- Sys.time()
n_samples <- 500
lut <- data.frame(
  LAI = runif(n_samples, 0.5, 7), hspot = runif(n_samples, 0, 0.3),
  Cab = runif(n_samples, 10, 70), Car = runif(n_samples, 4, 15),
  Anth = runif(n_samples, 0, 3), Cbrown = runif(n_samples, 0, 0.4),
  EWT = runif(n_samples, 0.005, 0.03), LMA = runif(n_samples, 0.003, 0.015),
  N = runif(n_samples, 1.2, 2.2),
  LIDFa = -0.35, LIDFb = -0.15, TypeLidf = 1, tts = 30, tto = 0, psi = 0, alpha = 40,
  Prot = 0.0045, CBC = 0.005
)
spectra <- matrix(NA_real_, nrow = n_samples, ncol = length(wavelength))
for (i in seq_len(n_samples)) spectra[i, ] <- run_one(lut[i, ])
cat("step 2 (", n_samples, "sims) done,", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")

wl_thin <- seq(1, length(wavelength), by = 4)  # thin for plotting only -- keeps the PNG light
df2 <- do.call(rbind, lapply(seq_len(n_samples), function(i) {
  data.frame(wavelength = wavelength[wl_thin], reflectance = spectra[i, wl_thin], LAI = lut$LAI[i], sim = i)
}))
p2 <- ggplot(df2, aes(wavelength, reflectance, group = sim, color = LAI)) +
  geom_line(alpha = 0.35, linewidth = 0.4) +
  scale_color_viridis_c(name = "LAI", direction = -1) +
  labs(x = "Wavelength (nm)", y = "Reflectance",
       title = paste0(n_samples, " simulations (random LUT), colored by LAI")) +
  theme_prosail(legend_position = "right")
ggsave("www/tutorial_many_spectra.png", p2, width = 8, height = 4.5, dpi = 110)

# ---- Step 3: sensitivity analysis (one trait at a time) -------------------
# Same per-trait-independent-color-scale approach as Exercise-1.Rmd's own
# sensitivity chunk (a shared scale_color_viridis_c() across facets would
# make only the largest-range trait, Cab, show visible color variation).
t0 <- Sys.time()
sweep_trait <- function(trait, values) {
  do.call(rbind, lapply(values, function(v) {
    row <- baseline; row[[trait]] <- v
    data.frame(wavelength = wavelength[wl_thin], reflectance = run_one(row)[wl_thin], value = v)
  }))
}
sens <- list(
  Cab    = sweep_trait("Cab", seq(10, 70, length.out = 6)),
  LAI    = sweep_trait("LAI", seq(0.5, 7, length.out = 6)),
  EWT    = sweep_trait("EWT", seq(0.005, 0.03, length.out = 6)),
  Cbrown = sweep_trait("Cbrown", seq(0, 0.5, length.out = 6))
)
plot_trait <- function(name, df) {
  ggplot(df, aes(wavelength, reflectance, color = value, group = value)) +
    geom_line(linewidth = 0.9) +
    scale_color_viridis_c(name = name, direction = -1) +
    scale_x_continuous(breaks = seq(500, 2500, by = 1000)) +
    labs(x = "Wavelength (nm)", y = "Reflectance", title = name) +
    theme_prosail(base_size = 13, legend_position = "right")
}
p3 <- (plot_trait("Cab", sens$Cab) + plot_trait("LAI", sens$LAI)) /
  (plot_trait("EWT", sens$EWT) + plot_trait("Cbrown", sens$Cbrown)) +
  plot_annotation(
    title = "Sensitivity: one trait varied at a time, others held at baseline",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )
ggsave("www/tutorial_sensitivity.png", p3, width = 12, height = 9, dpi = 110)
cat("step 3 done,", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")

# ---- Step 4: sensor convolution (Sentinel-2A) ------------------------------
t0 <- Sys.time()
df_sat <- convolve_to_sensor(wavelength, refl_one, "Sentinel2A.MSI")
p4 <- ggplot() +
  geom_line(data = df1, aes(wavelength, reflectance, color = "Native (1nm)"), linewidth = 1.4) +
  geom_line(data = df_sat, aes(wl, reflectance, color = "Sentinel-2A bands"), linewidth = 1) +
  geom_point(data = df_sat, aes(wl, reflectance, color = "Sentinel-2A bands"), size = 3) +
  scale_color_manual(values = c("Native (1nm)" = "#1b7837", "Sentinel-2A bands" = "#b2182b")) +
  labs(x = "Wavelength (nm)", y = "Reflectance", color = "",
       title = "Convolved to Sentinel-2A's 13 bands") +
  theme_prosail(legend_position = "top")
ggsave("www/tutorial_sensor_convolution.png", p4, width = 8, height = 4.5, dpi = 110)
cat("step 4 done,", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")

# ---- Step 5: ML inversion (RandomForest, Cab from Sentinel-2A bands) ------
t0 <- Sys.time()
band_refl <- t(vapply(seq_len(n_samples), function(i) {
  convolve_to_sensor(wavelength, spectra[i, ], "Sentinel2A.MSI")$reflectance
}, numeric(nrow(df_sat))))
colnames(band_refl) <- paste0("B", df_sat$band)
ml_df <- data.frame(band_refl, Cab = lut$Cab)

set.seed(1)
train_idx <- sample(seq_len(n_samples), size = round(0.7 * n_samples))
rf <- randomForest::randomForest(Cab ~ ., data = ml_df[train_idx, ], ntree = 300)
pred <- predict(rf, ml_df[-train_idx, ])
obs <- ml_df$Cab[-train_idx]
r2 <- cor(pred, obs)^2
rmse <- sqrt(mean((pred - obs)^2))

df5 <- data.frame(observed = obs, predicted = pred)
p5 <- ggplot(df5, aes(observed, predicted)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(color = "#2166ac", size = 2.5, alpha = 0.7) +
  labs(x = "Observed Cab", y = "Predicted Cab",
       title = "RandomForest inversion of Cab from Sentinel-2A bands",
       subtitle = sprintf("Test set (n=%d): R\u00b2 = %.2f, RMSE = %.2f \u00b5g/cm\u00b2", length(obs), r2, rmse)) +
  theme_prosail(base_size = 16, legend_position = "none")
ggsave("www/tutorial_ml_inversion.png", p5, width = 9, height = 6.5, dpi = 110)
cat("step 5 done,", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s -- R2:", round(r2, 3), "RMSE:", round(rmse, 2), "\n")

cat("\nAll tutorial figures written to www/\n")
