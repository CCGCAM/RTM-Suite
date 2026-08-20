# ==============================================================================
# Pipeline step 2: trait inversion with ToolsRTM::get.inversion()'s 12
# classic ML algorithms (PLSR, SVM, RF, GB, NN, Bayesian, AdaBag, BRNN, xGB,
# RVM, qLASSO, Ensemble), on the native-resolution dataset built in step 1.
#
# Not every algorithm needs an extra package that may not be installed on
# every machine (e.g. Bayesian needs 'bartMachine' + a working rJava/JVM) --
# each is wrapped in tryCatch so one missing dependency doesn't stop the
# others; failures are reported, not hidden.
# ==============================================================================

rm(list = ls())

if (!requireNamespace("ToolsRTM", quietly = TRUE)) {
  remotes::install_gitlab("caminoccg/toolsrtm")
}
library(ToolsRTM)

out_dir <- "outs/rtm_sims/pipeline"
sim <- readRDS(file.path(out_dir, "1-datasets.rds"))
datasets <- sim$datasets

target_trait <- "Cab"    # any LUT column works, e.g. "Car", "LAI", "Anth", "CBC", "EWT"
algorithms <- c("PLSR", "SVM", "RF", "GB", "NN", "Bayesian",
                "AdaBag", "BRNN", "xGB", "RVM", "qLASSO", "Ensemble")

# get.inversion()'s n.cores previously did nothing (always used detectCores()-2
# internally regardless of what was passed -- fixed in ToolsRTM/R/get.inversion.R).
# TRUE = each algorithm's caret tuning runs on its own parallel::makeCluster();
# FALSE = n.cores = 1, no cluster created at all, fully sequential (needed e.g.
# in sandboxed environments where spawning worker Rscript processes is unstable).
use.parallel <- FALSE
n.cores <- if (use.parallel) max(1, parallel::detectCores() - 2) else 1

dataset <- datasets$native  # or datasets$se2a (Sentinel-2A bands) / datasets$prisma (PRISMA bands)
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

cat("\n\n=== Summary: ", target_trait, " inversion, all algorithms ===\n", sep = "")
print(df_results)

write.csv(df_results, file.path(out_dir, sprintf("2-inversion_ML_%s.csv", target_trait)), row.names = FALSE)
cat("\nSaved to '", out_dir, "/2-inversion_ML_", target_trait, ".csv'\n", sep = "")
