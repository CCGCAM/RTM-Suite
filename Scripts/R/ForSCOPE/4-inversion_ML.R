# ==============================================================================
# ForSCOPE course pipeline, step 2/3: classic ML trait inversion
# (ToolsRTM::get.inversion(), 11 algorithms) for Cab, LAI, EWT, Vcmax25.
#
# get.inversion(save.model=TRUE, save.path=...) already saves the fitted
# model (.rds), per-algorithm statistics (.csv) and a predicted-vs-observed
# scatter plot (.png) -- used here rather than reimplementing that.
#
# Extra experiment (section at the bottom): does adding SIF687/SIF760 as
# predictors improve Vcmax25 inversion specifically, beyond reflectance +
# indices alone? Fluorescence is a direct byproduct of photosynthetic
# electron transport, so there's real physical reason to expect it helps.
# ==============================================================================

rm(list = ls())

if (!requireNamespace("ToolsRTM", quietly = TRUE)) {
  remotes::install_gitlab("caminoccg/toolsrtm")
}
library(ToolsRTM)

out_dir <- "outs/ForSCOPE"
sim <- readRDS(file.path(out_dir, "SCOPE-1-datasets.rds"))
datasets <- sim$datasets

target_traits <- c("Cab", "LAI", "EWT", "Vcmax25")
algorithms <- c("PLSR", "SVM", "RF", "GB", "Bayesian",
                "AdaBag", "BRNN", "xGB", "RVM", "qLASSO", "Ensemble")

# TRUE = each algorithm's caret tuning runs on its own parallel::makeCluster();
# FALSE = fully sequential, no cluster created (safer in sandboxed environments).
use.parallel <- FALSE
n.cores <- if (use.parallel) max(1, parallel::detectCores() - 2) else 1

dataset <- datasets$native
all_results <- list()

for (target_trait in target_traits) {
  cat(sprintf("\n\n=== ML inversion: %s ===\n", target_trait))

  model_dir <- file.path(out_dir, "models_ML", target_trait)
  dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

  inputs_inv <- setdiff(names(dataset), names(sim$LUT))  # bands + indices only, not other traits
  inputs_inv <- c(target_trait, inputs_inv)

  results <- list()
  for (algo in algorithms) {
    cat(sprintf("[%s] ... ", algo))
    fit <- tryCatch(
      suppressMessages(suppressWarnings(get.inversion(
        data = dataset[, unique(c(target_trait, inputs_inv))], depVar = target_trait,
        inputs = setdiff(inputs_inv, target_trait),
        algorithm = algo, n.cores = n.cores, n.samples = nrow(dataset), seed = 123,
        save.model = TRUE, save.path = model_dir
      ))),
      error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
    )
    if (is.null(fit) || is.null(fit$model)) {
      results[[algo]] <- data.frame(trait = target_trait, algorithm = algo, R2 = NA, RMSE = NA, status = "failed")
      next
    }
    r2 <- fit$statistics$R2[2]
    rmse <- fit$statistics$RMSE[2]
    cat(sprintf("R2 (test) = %.3f, RMSE (test) = %.3f\n", r2, rmse))
    results[[algo]] <- data.frame(trait = target_trait, algorithm = algo, R2 = r2, RMSE = rmse, status = "ok")
  }

  df_results <- do.call(rbind, results)
  rownames(df_results) <- NULL
  write.csv(df_results, file.path(model_dir, sprintf("metrics_ML_%s.csv", target_trait)), row.names = FALSE)
  all_results[[target_trait]] <- df_results
}

df_all <- do.call(rbind, all_results)
rownames(df_all) <- NULL
cat("\n\n=== Summary: ML inversion, all traits x algorithms ===\n")
print(df_all)
write.csv(df_all, file.path(out_dir, "4-inversion_ML_summary.csv"), row.names = FALSE)

## ----------------------------------------------------------------------------
## Summary figure: R2 per trait x algorithm
## ----------------------------------------------------------------------------

library(ggplot2)
summary_plot <- ggplot(subset(df_all, status == "ok"), aes(x = algorithm, y = R2, fill = trait)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(x = "Algorithm", y = expression(R^2~"(test)"), fill = "Trait",
       title = "ML inversion accuracy by trait and algorithm (ForSCOPE)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 11))
print(summary_plot)
ggsave(file.path(out_dir, "4-inversion_ML_summary.png"), plot = summary_plot,
       width = 20, height = 12, dpi = 300, units = "cm")

cat("\nModels saved under '", out_dir, "/models_ML/<trait>/model-<algorithm>.rds'\n", sep = "")
cat("Metrics saved to '", out_dir, "/4-inversion_ML_summary.csv'\n", sep = "")
cat("Summary figure saved to '", out_dir, "/4-inversion_ML_summary.png'\n", sep = "")

## ----------------------------------------------------------------------------
## Extra: does adding SIF687/SIF760 improve Vcmax25 inversion specifically?
## Re-run a subset of algorithms for Vcmax25 with SIF687/SIF760 added to the
## predictor set, and compare R2 against the SIF-free run already done above.
## ----------------------------------------------------------------------------

cat("\n\n=== SIF-vs-no-SIF comparison for Vcmax25 ===\n")

sif_algorithms <- c("PLSR", "SVM", "RF", "GB")  # a fast representative subset, not all 11 again
inputs_no_sif <- setdiff(names(dataset), names(sim$LUT))
inputs_with_sif <- c(inputs_no_sif, "SIF687", "SIF760")

sif_results <- list()
for (algo in sif_algorithms) {
  for (variant in c("no_SIF", "with_SIF")) {
    inputs_v <- if (variant == "no_SIF") inputs_no_sif else inputs_with_sif
    cat(sprintf("[Vcmax25 / %s / %s] ... ", algo, variant))
    fit <- tryCatch(
      suppressMessages(suppressWarnings(get.inversion(
        data = dataset[, unique(c("Vcmax25", inputs_v))], depVar = "Vcmax25",
        inputs = inputs_v, algorithm = algo, n.cores = n.cores, n.samples = nrow(dataset), seed = 123,
        save.model = FALSE
      ))),
      error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
    )
    if (is.null(fit) || is.null(fit$model)) {
      sif_results[[paste(algo, variant)]] <- data.frame(algorithm = algo, variant = variant, R2 = NA, RMSE = NA)
      next
    }
    cat(sprintf("R2 (test) = %.3f\n", fit$statistics$R2[2]))
    sif_results[[paste(algo, variant)]] <- data.frame(algorithm = algo, variant = variant,
                                                        R2 = fit$statistics$R2[2], RMSE = fit$statistics$RMSE[2])
  }
}
df_sif <- do.call(rbind, sif_results)
rownames(df_sif) <- NULL
print(df_sif)
write.csv(df_sif, file.path(out_dir, "4b-SIF_vs_noSIF_Vcmax25.csv"), row.names = FALSE)

sif_plot <- ggplot(df_sif, aes(x = algorithm, y = R2, fill = variant)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(x = "Algorithm", y = expression(R^2~"(test)"), fill = "Predictors",
       title = "Does SIF687/SIF760 improve Vcmax25 inversion?") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
print(sif_plot)
ggsave(file.path(out_dir, "4b-SIF_vs_noSIF_Vcmax25.png"), plot = sif_plot, width = 16, height = 10, dpi = 300, units = "cm")
cat("\nSIF comparison saved to '", out_dir, "/4b-SIF_vs_noSIF_Vcmax25.{csv,png}'\n", sep = "")
