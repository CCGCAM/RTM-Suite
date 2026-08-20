# ==============================================================================
# Compare TRADITIONAL inversion (LUT nearest-neighbour matching by a merit
# function -- RMSE, FGE) against MACHINE-LEARNING inversion (Random Forest),
# on the exact same held-out test set, convolved to real Sentinel-2A bands --
# so all three approaches see literally the same observations, just invert
# them differently:
#   - Traditional (RMSE / FGE): get.inversionOpt() ranks every TRAIN spectrum
#     by how well it matches each TEST spectrum under the chosen merit
#     function and averages the n.opt best matches' trait values. No model
#     fitting at all -- pure LUT search, one call per merit function.
#   - ML (Random Forest): get.inversion() fits a model on TRAIN (band
#     reflectance -> trait), then predicts on the SAME TEST set.
# All three are scored against the SAME held-out Cab values (R², RMSE), so
# the comparison is apples-to-apples.
# ==============================================================================

rm(list = ls())  # avoid leftover objects from a previous run/session leaking in

library(ToolsRTM)
library(ggplot2)

devtools::load_all("../../../ToolsRTM")

out_dir <- "../../../outs/Comparison"
dir.create(out_dir, showWarnings = FALSE)

set.seed(1)

r2_f   <- function(obs, pred) 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
rmse_f <- function(obs, pred) sqrt(mean((obs - pred)^2))

## ----------------------------------------------------------------------------
## 1. Simulate a LUT (PROSPECT-PRO + fourSAIL) and convolve to Sentinel-2A
## ----------------------------------------------------------------------------

n.samples <- 300
inputs <- ToolsRTM::inputsPROSAIL
LUT <- as.data.frame(getLUT(inputs = inputs, nLUT = n.samples, setseed = 1))

wl <- 400:2500
rsoil <- rep(0.15, length(wl))

cat(sprintf("Simulating %d PROSPECT-PRO + fourSAIL spectra ...\n", n.samples))
refl <- t(sapply(seq_len(n.samples), function(i) {
  sim <- suppressMessages(simulate_RTM(inputLUT = LUT[i, ], rsoil = rsoil,
                                        leaf.model = "PROSPECT-PRO", canopy.model = "fourSAIL"))
  sim$rsot
}))

refl_X <- as.data.frame(refl)
colnames(refl_X) <- paste0("X", wl)
refl_X <- cbind(id = seq_len(nrow(refl_X)), refl_X)

cat("Convolving to Sentinel-2A (13 bands) ...\n")
se2a <- suppressMessages(get.spectra.convolved(rfl = refl_X, sensor = "Sentinel2a", plot.spectra = FALSE))
wl_bands <- as.numeric(names(se2a)[-1])  # actual SRF-weighted band-center wavelengths
band_names <- c("B1","B2","B3","B4","B5","B6","B7","B8","B8A","B9","B10","B11","B12")
names(se2a) <- c("id", band_names)

## ----------------------------------------------------------------------------
## 2. Train/test split (70/30) -- SAME split feeds all three methods
## ----------------------------------------------------------------------------

train_idx <- sample(seq_len(n.samples), size = round(0.7 * n.samples))
test_idx  <- setdiff(seq_len(n.samples), train_idx)

LUT_train <- LUT[train_idx, ]
LUT_test  <- LUT[test_idx, ]
se2a_mat  <- as.matrix(se2a[, band_names])
se2a_train_mat <- se2a_mat[train_idx, ]
se2a_test_mat  <- se2a_mat[test_idx, ]

train_df <- cbind(LUT_train, se2a[train_idx, band_names])
test_df  <- cbind(LUT_test,  se2a[test_idx,  band_names])

cat(sprintf("Train: %d spectra, Test (held out): %d spectra\n", length(train_idx), length(test_idx)))

## ----------------------------------------------------------------------------
## 3. Traditional inversion -- LUT nearest-neighbour matching (RMSE, FGE)
## ----------------------------------------------------------------------------

cat("Traditional inversion (merit-RMSE) ...\n")
opt_rmse <- get.inversionOpt(rfl.sensor = se2a_test_mat, rfl.rtm = se2a_train_mat, LUT = LUT_train,
                              wave = wl_bands, method = "merit-RMSE", nOpt = 5)

cat("Traditional inversion (merit-FGE) ...\n")
opt_fge <- get.inversionOpt(rfl.sensor = se2a_test_mat, rfl.rtm = se2a_train_mat, LUT = LUT_train,
                             wave = wl_bands, method = "merit-FGE", nOpt = 5)

## ----------------------------------------------------------------------------
## 4. ML inversion -- Random Forest (get.inversion), same train/test split
## ----------------------------------------------------------------------------

cat("ML inversion (Random Forest) ...\n")
ml_res <- get.inversion(data = train_df, depVar = "Cab", inputs = band_names,
                         algorithm = "RF", n.cores = 1, n.samples = nrow(train_df))
pred_ml_test <- as.numeric(stats::predict(ml_res$model, newdata = test_df[, c("Cab", band_names)]))

## ----------------------------------------------------------------------------
## 5. Compare: R² and RMSE on the held-out test set, for all three methods
## ----------------------------------------------------------------------------

obs_cab <- LUT_test$Cab

results <- data.frame(
  method = c("Traditional (merit-RMSE)", "Traditional (merit-FGE)", "ML (Random Forest)"),
  R2   = c(r2_f(obs_cab, opt_rmse[[2]]$Cab), r2_f(obs_cab, opt_fge[[2]]$Cab), r2_f(obs_cab, pred_ml_test)),
  RMSE = c(rmse_f(obs_cab, opt_rmse[[2]]$Cab), rmse_f(obs_cab, opt_fge[[2]]$Cab), rmse_f(obs_cab, pred_ml_test))
)
cat("\n=== Cab retrieval on held-out Sentinel-2A test set (n = ", length(test_idx), ") ===\n", sep = "")
print(results, row.names = FALSE, digits = 3)
write.csv(results, file.path(out_dir, "inversion_traditional_vs_ML_stats.csv"), row.names = FALSE)

## ----------------------------------------------------------------------------
## 6. Plot: predicted vs observed Cab, all three methods
## ----------------------------------------------------------------------------

plot_df <- rbind(
  data.frame(method = "Traditional (merit-RMSE)", observed = obs_cab, predicted = opt_rmse[[2]]$Cab),
  data.frame(method = "Traditional (merit-FGE)",  observed = obs_cab, predicted = opt_fge[[2]]$Cab),
  data.frame(method = "ML (Random Forest)",       observed = obs_cab, predicted = pred_ml_test)
)

p <- ggplot(plot_df, aes(x = observed, y = predicted)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = method), alpha = 0.7, size = 2) +
  facet_wrap(~method) +
  labs(x = "Observed Cab (ug/cm2)", y = "Predicted Cab (ug/cm2)",
       title = "Traditional (LUT merit-function) vs ML inversion",
       subtitle = sprintf("PROSPECT-PRO + fourSAIL -> Sentinel-2A, n.train=%d, n.test=%d",
                           length(train_idx), length(test_idx))) +
  theme_bw() + theme(legend.position = "none", plot.title = element_text(face = "bold"))

print(p)
plot_path <- file.path(out_dir, "inversion_traditional_vs_ML.png")
ggsave(plot_path, plot = p, width = 24, height = 9, dpi = 300, units = "cm")
cat("\nSaved plot to '", plot_path, "' and stats to 'inversion_traditional_vs_ML_stats.csv'\n", sep = "")
