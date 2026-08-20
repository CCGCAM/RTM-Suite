# ============================================================================
# Running SCOPE v2.1 with SCOPEinR + ToolsRTM
# Soil Canopy Observation, Photochemistry and Energy balance model
# Example script for students
# Author: Carlos Camino (adapted for teaching)
# ============================================================================

# OPTIONAL: clear environment (comment out if you don't like this behaviour)
rm(list = ls())

# ----------------------------------------------------------------------------
# 0. Load / install required packages
# ----------------------------------------------------------------------------

load_packages <- function(pkgs) {
  for (p in pkgs) {
    if (!require(p, character.only = TRUE)) {
      install.packages(p)
      library(p, character.only = TRUE)
    }
  }
}

load_packages(c(
  "parallel", "doParallel",          # parallel processing
  "ToolsRTM", "SCOPEinR",'copula',            # RTM + SCOPE wrapper
  "dplyr", "tidyr", "ggplot2",       # data handling + plots
  "signal", "pracma", "expint",
  "copula", "progress", "data.table" # required by ToolsRTM / SCOPEinR
))


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

#Check the version is 0.47

# ----------------------------------------------------------------------------
# 1. User settings
# ----------------------------------------------------------------------------

# Number of SCOPE simulations to run
n_samples   <- 200

# Number of chunks to split the LUT into (for memory management)
n_chunks    <- 5

# Paths to input tables
path_setoptions <- system.file("input", "setoptions.csv", package = "SCOPEinR")
path_inputLUT   <- system.file("input", "inputs_SCOPE.csv", package = "SCOPEinR")

# Path where SCOPE outputs are stored (created by SCOPEinR) -- project-level
# outputs folder, never inside Scripts/
path_out <- "../../../outs/rtm_sims/scope_runs"

# Seed for reproducibility
set.seed(123)

# ----------------------------------------------------------------------------
# 2. Read SCOPE options and base LUT
# ----------------------------------------------------------------------------

# Table with SCOPE run options (radiative transfer, energy balance, etc.)
scope_options <- read.table(path_setoptions, header = TRUE, sep = ",")

# Base LUT structure (columns + default values)
inputLUT <- read.table(path_inputLUT, header = TRUE, sep = ",")

# Build an initial LUT with n_samples rows using ToolsRTM helper
LUT <- getLUT.SCOPE(inputLUT = inputLUT, nLUT = n_samples)

# ----------------------------------------------------------------------------
# 3. Randomize key plant traits for the LUT
#    (biochemical, structural, and photosynthetic parameters)
# ----------------------------------------------------------------------------

# Biochemical traits
LUT$Cab   <- runif(n_samples, min = 5,  max = 70)   # chlorophyll (µg cm-2)
LUT$Car   <- runif(n_samples, min = 1,  max = 15)   # carotenoids
LUT$Anth  <- runif(n_samples, min = 0,  max = 5)    # anthocyanins
LUT$Cbrown<- runif(n_samples, min = 0,  max = 1)    # brown pigments
LUT$Cx    <- runif(n_samples, min = 0,  max = 1)    # additional absorber
LUT$LMA   <- runif(n_samples, min = 0.003723, max = 0.006104) # leaf mass area

# Photosynthetic capacity
LUT$Vcmax25 <- runif(n_samples, min = 5, max = 90)  # Rubisco capacity at 25°C

# Water content (fixed in this example)
LUT$EWT <- 0.01   # equivalent water thickness (m)

# Viewing / illumination geometry + structure
LUT$tts   <- runif(n_samples, min = 15, max = 75)   # solar zenith
LUT$tto   <- runif(n_samples, min = 0,  max = 30)   # sensor zenith
LUT$hspot <- runif(n_samples, min = 0,  max = 1)    # hot-spot parameter
LUT$LAI   <- runif(n_samples, min = 0.05, max = 6)  # leaf area index

# Correlated leaf angle distribution parameters (LIDFa, LIDFb)
n_seed <- sample(1:n_samples, 1)
LUT_lidf <- ToolsRTM::getCor( n_inputs   = 2, setseed    = n_seed,
                              distribution = "Uniform",
                              nLUT       = n_samples, rho        = 0.20,
                              Varnames   = c("LIDFa", "LIDFb"),
                              MinRange   = c(-0.5, -0.5), MaxRange   = c( 0.2,  0.2))

LUT$LIDFa <- LUT_lidf$LUT$LIDFa
LUT$LIDFb <- LUT_lidf$LUT$LIDFb

# Quick check: dimension + correlation plot
print(dim(LUT))
plot(LUT$LIDFa, LUT$LIDFb,
     xlab = "LIDFa", ylab = "LIDFb",
     main = "Correlated leaf angle parameters")

# ----------------------------------------------------------------------------
# 4. Run SCOPE simulations in parallel (chunked)
# ----------------------------------------------------------------------------

n_rows     <- nrow(LUT)
chunk_size <- n_rows / n_chunks
cat("Total simulations:", n_rows, " |  chunk size:", chunk_size, "\n")

for (i in seq_len(n_chunks)) {

  cat("\nRunning chunk", i, "of", n_chunks, "...\n")

  start_row <- ((i - 1) * chunk_size) + 1
  end_row   <- min(i * chunk_size, n_rows)

  LUT_chunk <- LUT[start_row:end_row, ]

  # Main SCOPEinR call:
  #   - leaf.model    = 'fluspect-CX'
  #   - canopy.model  = 'fourSAIL'
  #   - parallel      = TRUE  (internal parallelisation)
  #   - get.outputs   = 'ALL' (all variables)
  #   - get.csv       = TRUE  (write outputs as CSV)
  # BUG (fixed): path.out was never passed here, so get.SCOPE.parallel()
  # fell back to its own default ('outs/', relative to the current working
  # directory) instead of path_out defined above -- with cwd = Scripts/,
  # every chunk silently wrote into Scripts/outs/ (violating the "never
  # write outputs inside Scripts/" rule), while the merge/plotting steps
  # further down read from path_out and picked up whatever unrelated SCOPE
  # run folder happened to already be there instead of these chunks.
  SCOPEinR::get.SCOPE.parallel(
    LUT          = LUT_chunk,
    options.SCOPE= scope_options,
    optipar      = SCOPEinR::optipar2017.ProspectD,
    leaf.model   = "fluspect-CX",
    canopy.model = "fourSAIL",
    path.out     = path_out,
    parallel     = TRUE,
    get.outputs  = "ALL",
    get.plots    = FALSE,
    get.csv      = TRUE
  )
}

# ----------------------------------------------------------------------------
# 5. Merge all SCOPE CSV outputs into single files
#    (uses the getCSV helper from ToolsRTM)
# ----------------------------------------------------------------------------

# Here we tell getCSV to:
#  - look in 'outs/'
#  - use the last 5 folders (the 5 chunks above)
#  - merge all available variables ("All")
df.scope<-getCSV( path.out    = path_out, n.folders   = n_chunks,
                            files.names = "All")

# Example: access merged fluxes
# fluxes_all <- merged_outputs$data[[which(merged_outputs$names == "fluxes")]]

# ----------------------------------------------------------------------------
# 6. Generate additional SCOPE diagnostic plots
# ----------------------------------------------------------------------------

# Get the last SCOPE output folder (most recent run)
subdirs          <- list.dirs(path_out, recursive = FALSE)
last_subdir      <- subdirs[length(subdirs)]
cat("\nUsing folder for plots:", last_subdir, "\n")

# Plant traits to explore in the plots
traits <- c("Vcmax25", "Anth")

# Available plot types: "reflectance", "fluorescence", "radiance"
SCOPEinR::get.SCOPE.plots(
  path.files  = last_subdir,
  plant.trait = traits,
  get.plots   = "reflectance"
)

SCOPEinR::get.SCOPE.plots(
  path.files  = last_subdir,
  plant.trait = traits,
  get.plots   = "fluorescence"
)

SCOPEinR::get.SCOPE.plots(
  path.files  = last_subdir,
  plant.trait = traits,
  get.plots   = "radiance"
)

cat("\nSCOPE simulations and plots finished.\n")
