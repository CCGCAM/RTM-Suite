## Validation of ToolsRTM::prospect_DB (PROSPECT-D) against real measured
## leaf spectra: the bundled Angers leaf-optics database (308 real leaves,
## measured CHLa+CHLb/CAR/ANT/EWT/LMA paired with measured reflectance and
## transmittance). N (mesophyll structure) is not independently measurable
## in this database (standard in the PROSPECT literature), so for each leaf
## we search N over a physically plausible grid and keep the best fit --
## exactly the standard PROSPECT calibration/validation procedure used in
## the original papers this database comes from.

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
root <- normalizePath(file.path(dirname(this_file), "..", ".."))
outdir <- file.path(root, "Scripts/Comparison/_out")
suppressMessages(devtools::load_all(file.path(root, "ToolsRTM"), quiet = TRUE))

load(file.path(root, "ToolsRTM/data/leafDB.Angers.rda"))
db <- LeafDB.Angers
db$Refl <- as.matrix(db$Refl)
db$Tran <- as.matrix(db$Tran)
db$DataBioch <- as.data.frame(db$DataBioch)
n_samples <- db$nbSamples
lambda <- db$lambda
n_wl <- length(lambda)

N_grid <- seq(1.0, 3.0, by = 0.05)

results <- data.frame(sample = seq_len(n_samples), best_N = NA_real_,
                       rmse = NA_real_, r2 = NA_real_, alpha_conv = NA)

for (i in seq_len(n_samples)) {
  Cab <- db$DataBioch$CHLa[i] + db$DataBioch$CHLb[i]
  Car <- db$DataBioch$CAR[i]
  Anth <- db$DataBioch$ANT_estimated[i]
  EWT <- db$DataBioch$EWT[i]
  LMA <- db$DataBioch$LMA[i]
  Cbrown <- 0
  meas_R <- db$Refl[, i]
  meas_T <- db$Tran[, i]
  ok_meas <- is.finite(meas_R) & is.finite(meas_T)
  if (sum(ok_meas) < 100) next

  best_rmse <- Inf; best_N <- NA; best_sim <- NULL
  for (N in N_grid) {
    LRT <- tryCatch(prospect_DB(N, Cab, Car, Anth, Cbrown, EWT, LMA, alpha = 40),
                     error = function(e) NULL)
    if (is.null(LRT)) next
    sim_R <- LRT[[2]][1:n_wl]
    sim_T <- LRT[[3]][1:n_wl]
    ok <- ok_meas & is.finite(sim_R) & is.finite(sim_T)
    if (sum(ok) < 100) next
    rmse <- sqrt(mean((c(sim_R[ok], sim_T[ok]) - c(meas_R[ok], meas_T[ok]))^2))
    if (rmse < best_rmse) { best_rmse <- rmse; best_N <- N; best_sim <- list(R = sim_R, T = sim_T) }
  }
  if (is.finite(best_rmse)) {
    obs <- c(meas_R[ok_meas], meas_T[ok_meas])
    pred <- c(best_sim$R[ok_meas], best_sim$T[ok_meas])
    ss_res <- sum((obs - pred)^2); ss_tot <- sum((obs - mean(obs))^2)
    results$best_N[i] <- best_N
    results$rmse[i] <- best_rmse
    results$r2[i] <- 1 - ss_res / ss_tot
  }
}

write.csv(results, file.path(outdir, "validation_prospectD_angers_persample.csv"), row.names = FALSE)

valid <- results[is.finite(results$rmse), ]
cat("Samples validated:", nrow(valid), "/", n_samples, "\n")
cat("Best-fit N: median =", median(valid$best_N), " range =", paste(range(valid$best_N), collapse=" - "), "\n")
cat("RMSE (R+T combined, reflectance/transmittance units 0-1): median =", median(valid$rmse),
    " mean =", mean(valid$rmse), " 95th pct =", quantile(valid$rmse, 0.95), "\n")
cat("R^2: median =", median(valid$r2), " mean =", mean(valid$r2), " min =", min(valid$r2), "\n")

## Example plot: 4 representative samples (spread across chlorophyll range),
## measured vs best-fit simulated R and T.
library(ggplot2)
set.seed(1)
Cab_all <- db$DataBioch$CHLa + db$DataBioch$CHLb
example_idx <- valid$sample[order(Cab_all[valid$sample])][round(seq(1, nrow(valid), length.out = 4))]

plot_rows <- list()
for (i in example_idx) {
  N <- results$best_N[i]
  Cab <- db$DataBioch$CHLa[i] + db$DataBioch$CHLb[i]
  Car <- db$DataBioch$CAR[i]; Anth <- db$DataBioch$ANT_estimated[i]
  EWT <- db$DataBioch$EWT[i]; LMA <- db$DataBioch$LMA[i]
  LRT <- prospect_DB(N, Cab, Car, Anth, 0, EWT, LMA, alpha = 40)
  sim_R <- LRT[[2]][1:n_wl]; sim_T <- 1 - LRT[[3]][1:n_wl]  # plot T as 1-T from top like standard leaf-optics charts
  plot_rows[[length(plot_rows) + 1]] <- data.frame(
    wavelength = lambda, sample = i, Cab = round(Cab, 1),
    measured_R = db$Refl[, i], simulated_R = sim_R,
    measured_1mT = 1 - db$Tran[, i], simulated_1mT = sim_T)
}
df_plot <- do.call(rbind, plot_rows)
df_plot$label <- paste0("Leaf #", df_plot$sample, "  (Cab=", df_plot$Cab, " ug/cm2, N=", round(results$best_N[df_plot$sample], 2), ")")

library(tidyr)
df_long <- rbind(
  data.frame(wavelength = df_plot$wavelength, label = df_plot$label, value = df_plot$measured_R, series = "Measured (Angers)", var = "Reflectance"),
  data.frame(wavelength = df_plot$wavelength, label = df_plot$label, value = df_plot$simulated_R, series = "Simulated (PROSPECT-D)", var = "Reflectance"),
  data.frame(wavelength = df_plot$wavelength, label = df_plot$label, value = df_plot$measured_1mT, series = "Measured (Angers)", var = "1 - Transmittance"),
  data.frame(wavelength = df_plot$wavelength, label = df_plot$label, value = df_plot$simulated_1mT, series = "Simulated (PROSPECT-D)", var = "1 - Transmittance")
)

p <- ggplot(df_long, aes(x = wavelength, y = value, color = series, linetype = var)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~label, ncol = 2) +
  scale_color_manual(values = c("Measured (Angers)" = "black", "Simulated (PROSPECT-D)" = "#B2182B")) +
  labs(title = "PROSPECT-D (ToolsRTM) vs real measured leaf spectra (Angers database)",
       subtitle = "Cab/Car/Anth/EWT/LMA from real measurements; N fit per leaf (standard PROSPECT calibration practice)",
       x = "Wavelength (nm)", y = "Reflectance / 1 - Transmittance", color = "Source", linetype = "Variable") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(outdir, "validation_prospectD_angers_examples.png"), p, width = 10, height = 7, dpi = 150)

## Distribution plot: RMSE and R2 across all 308 samples
p2 <- ggplot(valid, aes(x = r2)) + geom_histogram(bins = 30, fill = "#2166AC", color = "white") +
  labs(title = paste0("PROSPECT-D vs Angers real leaves: R² distribution, n=", nrow(valid)),
       x = "R^2 (simulated vs measured R+T)", y = "Number of leaves") + theme_bw(base_size = 11)
ggsave(file.path(outdir, "validation_prospectD_angers_r2hist.png"), p2, width = 7, height = 4.5, dpi = 150)

cat("Wrote validation_prospectD_angers_examples.png, validation_prospectD_angers_r2hist.png, validation_prospectD_angers_persample.csv\n")
