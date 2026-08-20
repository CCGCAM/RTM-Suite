# ==============================================================================
# One-time setup: a reproducible Python environment for TensorFlow/Keras, so
# ToolsRTM::getMLmodel() (deep-learning trait inversion, used by
# Scripts/Pipeline/5-inversion_deep_learning.R) works out of the box on any
# machine -- no manual Python/conda install, nothing touching a system Python.
#
# Uses reticulate's built-in `py_require()` mechanism (reticulate >= 1.41):
# declare which Python packages you need, and reticulate resolves/downloads
# an isolated, reproducible environment for them on demand (via `uv` under
# the hood) the first time Python is actually used (e.g. the first
# `library(keras)` call) -- no conda/Miniconda install needed at all, and
# nothing is shared with any other Python on the machine.
#
# Run this once per machine to confirm it works; every script in this
# pipeline that needs TensorFlow (5-inversion_deep_learning.R) repeats the
# `py_require()` + `Sys.setenv(TF_USE_LEGACY_KERAS = "1")` lines itself, so
# this file is just a standalone check, not a dependency the others load.
#
# Background: the R `keras` package expects the classic Keras 2 API, but
# modern TensorFlow ships Keras 3 by default -- `tf-keras` is the official
# compatibility shim that restores the Keras 2 API TensorFlow still exposes
# it under. TF_USE_LEGACY_KERAS=1 tells `keras`/`tensorflow` (R side) to use
# that shim. See HANDOFF_Claude_Code.md for where this recipe originated
# (there, written against an older reticulate that needed a manual conda env
# -- py_require() supersedes that, but the TF_USE_LEGACY_KERAS reasoning
# still applies either way).
# ==============================================================================

rm(list = ls())

# IMPORTANT: run this in a FRESH R session (restart R first: Session > Restart R
# in RStudio, or Ctrl+Shift+F10). reticulate caches its Python/virtualenv
# location the first time it's touched in a session -- setting the env vars
# below *after* that point has no effect, which is why re-running this script
# in the same session as a previous failed attempt won't pick up the fix.

# If R's HOME resolves to a network drive (common on institutional/managed
# Windows profiles, e.g. a WUR network home under \\WURNET.NL\...), reticulate
# can fail to provision a Python environment there with "Access is denied" --
# force its virtualenv/cache location onto the local disk instead, which is
# always writable and doesn't depend on network/VPN state. Safe no-op if HOME
# is already local. Setting both env vars reticulate checks for this (newer
# versions prefer RETICULATE_VIRTUALENV_ROOT; WORKON_HOME is the older
# virtualenvwrapper-style fallback it also honors) so this works regardless
# of exactly which reticulate version is installed.
local_venv_root <- file.path(Sys.getenv("LOCALAPPDATA"), "r-reticulate-venvs")
Sys.setenv(RETICULATE_VIRTUALENV_ROOT = local_venv_root)
Sys.setenv(WORKON_HOME = local_venv_root)

if (!requireNamespace("reticulate", quietly = TRUE)) install.packages("reticulate")
if (!requireNamespace("keras", quietly = TRUE)) install.packages("keras")
if (!requireNamespace("tensorflow", quietly = TRUE)) install.packages("tensorflow")

library(reticulate)
py_require(packages = c("tensorflow", "tf-keras"))
Sys.setenv(TF_USE_LEGACY_KERAS = "1")

cat("Importing tensorflow (this provisions the environment on first use -- ",
    "downloads ~350MB the very first time, cached after that) ...\n", sep = "")

tf_ok <- tryCatch({
  tf <- reticulate::import("tensorflow")
  cat("\nTensorFlow version:", as.character(tf$`__version__`), "\n")
  TRUE
}, error = function(e) {
  cat("\nTensorFlow import failed:", conditionMessage(e), "\n")
  FALSE
})

if (tf_ok) {
  cat("\n=== Setup OK. Scripts/Pipeline/5-inversion_deep_learning.R already includes\n")
  cat("the same py_require()/TF_USE_LEGACY_KERAS lines and will just work. ===\n")
} else {
  cat("\n=== Setup did not complete cleanly -- see the error above before running\n")
  cat("Scripts/Pipeline/5-inversion_deep_learning.R. ===\n")
}
