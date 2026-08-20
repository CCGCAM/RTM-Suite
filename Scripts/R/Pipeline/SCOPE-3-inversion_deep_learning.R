# ==============================================================================
# SCOPEinR pipeline step 3: trait inversion with deep learning, via
# ToolsRTM::getMLmodel(), on the SCOPE-simulated dataset from SCOPE-1-simulate.R.
# Run Scripts/Pipeline/0-setup_python_env.R once first.
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
library(keras)

use.dev.source <- TRUE

if (use.dev.source) {
  devtools::load_all("ToolsRTM/R")
} else {
  library(ToolsRTM)
}

out_dir <- "outs/rtm_sims/pipeline_scope"
sim <- readRDS(file.path(out_dir, "SCOPE-1-datasets.rds"))
datasets <- sim$datasets

target_trait <- "Cab"
dataset <- datasets$native
inputs_inv <- setdiff(names(dataset), names(sim$LUT))
model_dataset <- dataset[, c(target_trait, inputs_inv)]

cat(sprintf("=== Deep-learning inversion of '%s' (SCOPE-simulated data), %d predictors, %d simulations ===\n",
            target_trait, length(inputs_inv), nrow(model_dataset)))

architectures <- c("Hidden-layers", "CNN")

results <- list()
for (arch in architectures) {
  cat(sprintf("\n=== %s ===\n", arch))
  fit <- tryCatch(
    getMLmodel(dataset = model_dataset, depVar = target_trait, model = arch,
               optimizer = "adam", batch.size = 32, n.epochs = 50,
               save.model = TRUE, path.model = file.path(out_dir, "models_DL"),
               prop.split = c(0.8, 0.2), data.trans = "preProcess",
               method.preProcess = "Normalize", depVar.trans = FALSE),
    error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL }
  )
  results[[arch]] <- fit
}

cat("\n=== Summary ===\n")
for (arch in architectures) {
  if (is.null(results[[arch]])) {
    cat(sprintf("%s: failed (see error above)\n", arch))
  } else {
    cat(sprintf("%s: trained OK, model saved under '%s/models_DL'\n", arch, out_dir))
  }
}
