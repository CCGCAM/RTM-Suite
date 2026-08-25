# ==============================================================================
# ForSPART course pipeline, step 2/3: classic ML trait inversion
# (ToolsRTM::get.inversion(), 12 algorithms) for Cab, LAI, EWT.
#
# get.inversion(save.model=TRUE, save.path=...) already saves the fitted
# model (.rds), per-algorithm statistics (.csv) and a predicted-vs-observed
# scatter plot (.png) -- used here rather than reimplementing that.
# ==============================================================================

rm(list = ls())

if (!requireNamespace("ToolsRTM", quietly = TRUE)) {
  remotes::install_gitlab("caminoccg/toolsrtm")
}
library(ToolsRTM)

out_dir <- "outs/ForSPART"
sim <- readRDS(file.path(out_dir, "1-datasets.rds"))
datasets <- sim$datasets

target_traits <- c("Cab", "LAI", "EWT")
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
       title = "ML inversion accuracy by trait and algorithm (ForSPART, SPART)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 11))
print(summary_plot)
ggsave(file.path(out_dir, "4-inversion_ML_summary.png"), plot = summary_plot,
       width = 20, height = 12, dpi = 300, units = "cm")

cat("\nModels saved under '", out_dir, "/models_ML/<trait>/model-<algorithm>.rds'\n", sep = "")
cat("Metrics saved to '", out_dir, "/4-inversion_ML_summary.csv'\n", sep = "")
cat("Summary figure saved to '", out_dir, "/4-inversion_ML_summary.png'\n", sep = "")
