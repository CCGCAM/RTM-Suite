# ==============================================================================
# Compare SCOPE simulations across the leaf.model/canopy.model parameters
# that SCOPEinR::get.SCOPE() accepts, using the SAME LUT row for every run.
#
# FINDING (verified below): get.SCOPE()'s `leaf.model`/`canopy.model`
# arguments are accepted but not wired into the simulation -- SCOPE always
# runs its own multi-layer Fluspect leaf optics + its own multi-layer
# canopy engine (getRTMo), which is a different, purpose-built implementation
# from ToolsRTM::foursail()/foursail2()/inform() (different interface, no
# vertical-profile/gap-probability outputs the rest of the SCOPE pipeline
# needs). Swapping those in would mean reimplementing them inside SCOPE's
# framework, not adding an if/else -- decided NOT to attempt that; SCOPE's
# leaf/canopy model is meant to stay unique/fixed.
#
# What WAS fixed: get.SCOPE()/get.SCOPE.ind() now warn() if you pass
# anything other than the actual defaults, instead of silently ignoring the
# request (see get.SCOPE.R). This script still demonstrates the underlying
# behavior (all 5 combinations produce identical reflectance) and now also
# shows the warnings firing.
# ==============================================================================

rm(list = ls())  # avoid leftover objects from a previous run/session leaking in

devtools::load_all("../../../ToolsRTM")
devtools::load_all("../../../SCOPEinR")
library(SCOPEinR)
library(ToolsRTM)
library(ggplot2)
library(dplyr)

out_dir <- "../../../outs/Comparison"  # project-level outputs folder, never inside Scripts/
dir.create(out_dir, showWarnings = FALSE)
pdf(file.path(out_dir, "Rplots.pdf"))  # catches any stray plot()/print() call, keeps it out of Scripts/

opts <- read.table(system.file("input", "setoptions.csv", package = "SCOPEinR"), header = TRUE, sep = ",")
lut  <- read.table(system.file("input", "LUT_input.csv", package = "SCOPEinR"), header = TRUE, sep = ",")

combos <- list(
  list(leaf.model = "fluspect-CX", canopy.model = "fourSAIL"),
  list(leaf.model = "fluspect-B",  canopy.model = "fourSAIL"),
  list(leaf.model = "PROSPECT",    canopy.model = "fourSAIL"),
  list(leaf.model = "fluspect-CX", canopy.model = "INFORM"),
  list(leaf.model = "fluspect-CX", canopy.model = "foursail2")
)

cat("=== Running get.SCOPE() with each leaf.model/canopy.model combination ===\n\n")

results <- lapply(combos, function(cb) {
  cat(sprintf("leaf.model = %-12s canopy.model = %-10s ... ", cb$leaf.model, cb$canopy.model))
  # Not wrapped in suppressWarnings(): we want the leaf.model/canopy.model
  # no-op warning (see script header) to actually show.
  sim <- suppressMessages(get.SCOPE(
    LUT = lut[1, ], options.SCOPE = opts, optipar = SCOPEinR::optipar2021.Pro.CX,
    leaf.model = cb$leaf.model, canopy.model = cb$canopy.model,
    get.outputs = "ALL", get.plots = FALSE
  ))
  cat("done\n")
  data.frame(wavelength = 1:length(sim[[1]]$data.rad$refl), refl = sim[[1]]$data.rad$refl,
             combo = paste(cb$leaf.model, cb$canopy.model, sep = " + "))
})

df <- bind_rows(results)

## ----------------------------------------------------------------------------
## Verify empirically whether the combinations actually differ
## ----------------------------------------------------------------------------

wide <- tidyr::pivot_wider(df, names_from = combo, values_from = refl, id_cols = wavelength)
combo_cols <- setdiff(names(wide), "wavelength")
max_diff <- max(sapply(combo_cols[-1], function(cn) max(abs(wide[[cn]] - wide[[combo_cols[1]]]), na.rm = TRUE)),
                na.rm = TRUE)

cat("\nMax absolute difference between any two combinations:", max_diff, "\n")
if (max_diff < 1e-9) {
  cat("=> CONFIRMED: leaf.model/canopy.model had NO effect on the output in this run.\n")
  cat("   (Known issue -- get.SCOPE() does not currently wire these parameters through\n")
  cat("   to the simulation. All 5 combinations produced bit-identical reflectance.)\n")
} else {
  cat("=> Combinations DID produce different results (max diff above) -- if you're\n")
  cat("   seeing this, the leaf.model/canopy.model dispatch has since been fixed;\n")
  cat("   update the comment above.\n")
}

## ----------------------------------------------------------------------------
## Plot -- if all lines overlap exactly, that IS the finding
## ----------------------------------------------------------------------------

p <- ggplot(df, aes(x = wavelength, y = refl, color = combo)) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  labs(title = "SCOPE reflectance across leaf.model/canopy.model combinations",
       subtitle = if (max_diff < 1e-9) {
         "All lines overlap exactly -- get.SCOPE() currently ignores these parameters (see script header)"
       } else "Differences are real (see console output for max absolute difference)",
       x = "Wavelength index", y = "Reflectance", color = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "top")

ggsave(file.path(out_dir, "compare_SCOPE_leafcanopy.png"), p, width = 9, height = 5.5, dpi = 150)
cat("\nSaved to '", out_dir, "/compare_SCOPE_leafcanopy.png'\n", sep = "")


dev.off()  # close the Rplots.pdf capture device opened near the top
