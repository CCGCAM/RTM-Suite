
# ============================================================================
# Running SCOPE v2.1 with SCOPEinR + ToolsRTM
# Soil Canopy Observation, Photochemistry and Energy balance model
# Example script for students
# Author: Carlos Camino (adapted for teaching)
# ============================================================================

rm(list = ls())   # Clear workspace

# ----------------------------------------------------------------------------
# 0. Load or install required packages
# ----------------------------------------------------------------------------

pkgs <- c(
  "parallel", "doParallel",
  "ToolsRTM", "SCOPEinR",
  "dplyr", "tidyr", "ggplot2"
)

for (p in pkgs) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

# ------------------------------------------------------------------------------
# Step 1: Install & Load ToolsRTM from GitLab
# ------------------------------------------------------------------------------

if (!requireNamespace("ToolsRTM", quietly = TRUE)) {
  # NOTE: do NOT auto-install via remotes::install_gitlab() here -- on a dev
  # machine that would silently overwrite the local source build with the
  # published GitLab version and hit the network. Install explicitly instead:
  #   devtools::document("ToolsRTM"); devtools::install("ToolsRTM", upgrade = "never")
  stop("ToolsRTM is not installed. Install the local dev build with ",
       "devtools::install('ToolsRTM', upgrade = 'never') before running this script.")
}
library(ToolsRTM)

cat("\n ToolsRTM is ready: ", as.character(packageVersion("ToolsRTM")), "\n", sep = "")
#Check the version is 0.65

# ------------------------------------------------------------------------------
# Step 2: Install & Load SCOPE from GitLab
# ------------------------------------------------------------------------------
if (!requireNamespace("SCOPEinR", quietly = TRUE)) {
  # NOTE: do NOT auto-install via remotes::install_gitlab() here -- see note above.
  stop("SCOPEinR is not installed. Install the local dev build with ",
       "devtools::install('SCOPEinR', upgrade = 'never') before running this script.")
}
library(SCOPEinR)

cat("\n SCOPEinR is ready: ", as.character(packageVersion("SCOPEinR")), "\n", sep = "")

#Check the version is 0.48


# ----------------------------------------------------------------------------
# 1. Load SCOPE model options
#    These options control radiative transfer, energy balance, and model numerics.
# ----------------------------------------------------------------------------

table.with.opts <- read.table(
  system.file("input", "setoptions.csv", package = "SCOPEinR"),
  header = TRUE, sep = ","
)

# ----------------------------------------------------------------------------
# 2. Build the LUT (Look-Up Table) with model input parameters
#    The LUT contains one row per SCOPE simulation.
# ----------------------------------------------------------------------------

N.Samples <- 100
file.LUT  <- "LUT"   # choose "default" or "LUT"

start_time <- Sys.time()

# ----------------------------------------------------------------------------
# 2A. Option A: Use predefined LUT (one simulation)
# ----------------------------------------------------------------------------
if (file.LUT == "default") {

  Table.LUT <- read.table(
    system.file("input", "LUT_input.csv", package = "SCOPEinR"),
    header = TRUE, sep = ","
  )

  db.sim <- SCOPEinR::get.SCOPE(
    LUT           = Table.LUT[1, ],
    options.SCOPE = table.with.opts,
    optipar       = SCOPEinR::optipar2021.Pro.CX,
    leaf.model    = "fluspect-CX",
    canopy.model  = "fourSAIL",
    get.outputs   = "ALL",
    get.plots     = FALSE
  )

  # ----------------------------------------------------------------------------
  # 2B. Option B: Generate a random LUT with N.Samples simulations
  # ----------------------------------------------------------------------------
} else {

  inputLUT <- read.table(
    system.file("input", "inputs_SCOPE.csv", package = "SCOPEinR"),
    header = TRUE, sep = ","
  )

  Table.LUT <- getLUT.SCOPE(
    inputLUT = inputLUT,
    nLUT     = N.Samples
  )

  db.sim <- get.SCOPE(
    LUT           = Table.LUT,
    n.LUT         = N.Samples,
    options.SCOPE = table.with.opts,
    optipar       = SCOPEinR::optipar2021.Pro.CX,
    leaf.model    = "fluspect-CX",
    canopy.model  = "fourSAIL",
    get.outputs   = "ALL",
    get.plots     = FALSE
  )
}

end_time <- Sys.time()

cat("\nExecution time:", end_time - start_time, "\n")

# ----------------------------------------------------------------------------
# 3. Extract and save SCOPE outputs
#    This saves reflectance, fluorescence, radiance, fluxes, and derived variables.
# ----------------------------------------------------------------------------

get.SCOPE.outputs(
  data.sim        = db.sim,
  N.sims          = N.Samples,
  LUT             = Table.LUT,
  get.outputs = 'ALL',
  path.out        = "../../../outs/rtm_sims/scope_runs/",
  get.more.inputs = c("refl", "lidf", "LIDFb", "Ft_Fo", "rdo"),
  get.plots       = TRUE
)
