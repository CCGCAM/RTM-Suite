# ==============================================================================
# SCOPEinR pipeline step 2: trait inversion with ToolsRTM::get.inversion()'s
# 12 classic ML algorithms, on the SCOPE-simulated dataset from SCOPE-1-simulate.R.
# Same algorithm/toggle as Scripts/Pipeline/2-inversion_ML.R -- get.inversion()
# is a ToolsRTM function used here on SCOPEinR-generated data.
# ==============================================================================

rm(list = ls())

if (!requireNamespace("ToolsRTM", quietly = TRUE)) {
  remotes::install_gitlab("caminoccg/toolsrtm")
}
library(ToolsRTM)

out_dir <- "outs/rtm_sims/pipeline_scope"
sim <- readRDS(file.path(out_dir, "SCOPE-1-datasets.rds"))
datasets <- sim$datasets

# Any column of inputs_SCOPE.csv works here -- not just leaf pigments/structure.
# E.g. target_trait <- "Vcmax25" inverts the biochemistry parameter (maximum
# carboxylation rate) straight from simulated reflectance, same mechanism as
# Cab/LAI below (verified: mechanically works end-to-end with real data;
# expect it to need more training samples than Cab for a good R2, since
# Vcmax's spectral signal is weaker/more indirect than pigment absorption).
target_trait <- "Cab"    # change to invert a different trait, e.g. "LAI" or "Vcmax25"
algorithms <- c("PLSR", "SVM", "RF", "GB", "NN", "Bayesian",
                "AdaBag", "BRNN", "xGB", "RVM", "qLASSO", "Ensemble")

# get.inversion()'s n.cores previously did nothing (fixed in ToolsRTM/R/get.inversion.R
# this session). TRUE = each algorithm's caret tuning runs on its own
# parallel::makeCluster(); FALSE = fully sequential, no cluster created.
use.parallel <- TRUE
n.cores <- if (use.parallel) max(1, parallel::detectCores() - 2) else 1

dataset <- datasets$native
inputs_inv <- setdiff(names(dataset), names(sim$LUT))  # bands + indices only, not other traits
inputs_inv <- c(target_trait, inputs_inv)  # get.inversion() needs the target column present too

cat(sprintf("=== Inverting '%s' from %d predictors (bands + indices), %d simulations, %d algorithms ===\n",
            target_trait, length(inputs_inv) - 1, nrow(dataset), length(algorithms)))

results <- list()
for (algo in algorithms) {
  cat(sprintf("\n[%s] ... ", algo))
  fit <- tryCatch(
    suppressMessages(suppressWarnings(get.inversion(
      data = dataset[, unique(c(target_trait, inputs_inv))], depVar = target_trait,
      inputs = setdiff(inputs_inv, target_trait),
      algorithm = algo, n.cores = n.cores, n.samples = nrow(dataset), seed = 123
    ))),
    error = function(e) { cat("FAILED:", conditionMessage(e)); NULL }
  )
  if (is.null(fit) || is.null(fit$model)) {
    results[[algo]] <- data.frame(algorithm = algo, R2 = NA, RMSE = NA, status = "failed")
    next
  }
  r2 <- fit$statistics$R2[2]
  rmse <- fit$statistics$RMSE[2]
  cat(sprintf("R2 (test) = %.3f, RMSE (test) = %.3f", r2, rmse))
  results[[algo]] <- data.frame(algorithm = algo, R2 = r2, RMSE = rmse, status = "ok")
}

df_results <- do.call(rbind, results)
rownames(df_results) <- NULL

cat("\n\n=== Summary: ", target_trait, " inversion (SCOPE-simulated data), all algorithms ===\n", sep = "")
print(df_results)

write.csv(df_results, file.path(out_dir, sprintf("SCOPE-2-inversion_ML_%s.csv", target_trait)), row.names = FALSE)
cat("\nSaved to '", out_dir, "/SCOPE-2-inversion_ML_", target_trait, ".csv'\n", sep = "")
