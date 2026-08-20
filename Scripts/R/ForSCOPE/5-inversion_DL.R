# ==============================================================================
# ForSCOPE course pipeline, step 3/3: deep-learning trait inversion
# (ToolsRTM::getMLmodel(), 'CNN' and 'Hidden-layers' architectures,
# TensorFlow/Keras) for Cab, LAI, EWT, Vcmax25. Run
# Scripts/Pipeline/0-setup_python_env.R once first if you haven't already.
#
# getMLmodel(save.model=TRUE, path.model=...) already saves the fitted model
# (.hdf5) -- used here rather than reimplementing that.
# ==============================================================================

rm(list = ls())

# If R's HOME resolves to a network drive, force reticulate's virtualenv/cache
# location onto the local disk -- see Scripts/Pipeline/0-setup_python_env.R.
# Run in a fresh R session: reticulate caches this on first use.
local_venv_root <- file.path(Sys.getenv("LOCALAPPDATA"), "r-reticulate-venvs")
Sys.setenv(RETICULATE_VIRTUALENV_ROOT = local_venv_root)
Sys.setenv(WORKON_HOME = local_venv_root)
Sys.setenv(TF_USE_LEGACY_KERAS = "1")
library(reticulate)
py_require(packages = c("tensorflow", "tf-keras"))
library(keras)  # getMLmodel() calls callback_early_stopping() unnamespaced -- needs this attached, not just installed

use.dev.source <- TRUE
if (use.dev.source) {
  devtools::load_all("../../../ToolsRTM/R")
} else {
  library(ToolsRTM)
}
library(ggplot2)

out_dir <- "../../../outs/ForSCOPE"
sim <- readRDS(file.path(out_dir, "SCOPE-1-datasets.rds"))
datasets <- sim$datasets

target_traits <- c("Cab", "LAI", "EWT", "Vcmax25")
architectures <- c("Hidden-layers", "CNN")

dataset <- datasets$native
all_results <- list()

for (target_trait in target_traits) {
  cat(sprintf("\n\n=== Deep-learning inversion: %s ===\n", target_trait))

  model_dir <- file.path(out_dir, "models_DL", target_trait)
  dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

  inputs_inv <- setdiff(names(dataset), names(sim$LUT))
  model_dataset <- dataset[, c(target_trait, inputs_inv)]

  for (arch in architectures) {
    cat(sprintf("[%s] ... ", arch))
    fit <- tryCatch(
      getMLmodel(dataset = model_dataset, depVar = target_trait, model = arch,
                 optimizer = "adam", batch.size = 16, n.epochs = 50,
                 save.model = TRUE, path.model = model_dir,
                 prop.split = c(0.8, 0.2), data.trans = "preProcess",
                 method.preProcess = "Normalize", depVar.trans = FALSE),
      error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
    )
    if (is.null(fit)) {
      all_results[[paste(target_trait, arch)]] <- data.frame(trait = target_trait, architecture = arch, status = "failed")
      next
    }
    cat("trained OK\n")
    plot.val <- fit[['plot.val']]
    if (!is.null(plot.val)) {
      print(plot.val)
      ggsave(file.path(model_dir, sprintf("scatter_%s.png", arch)), plot = plot.val,
             width = 10, height = 10, dpi = 300, units = "cm")
    }
    all_results[[paste(target_trait, arch)]] <- data.frame(trait = target_trait, architecture = arch, status = "ok")
  }
}

df_all <- do.call(rbind, all_results)
rownames(df_all) <- NULL
cat("\n\n=== Summary: deep-learning inversion, all traits x architectures ===\n")
print(df_all)
write.csv(df_all, file.path(out_dir, "5-inversion_DL_summary.csv"), row.names = FALSE)

cat("\nModels saved under '", out_dir, "/models_DL/<trait>/model_", "*.hdf5'\n", sep = "")
cat("Metrics saved to '", out_dir, "/5-inversion_DL_summary.csv'\n", sep = "")
