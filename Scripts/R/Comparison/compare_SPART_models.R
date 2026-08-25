# ==============================================================================
# Verify ToolsRTM::SPART() (Soil-Plant-Atmosphere Radiative Transfer model:
# TOC reflectance -> TOA reflectance/radiance through a sensor's SRF/SMAC
# atmospheric correction coefficients) actually runs and that its LeafModel
# argument really changes the output (unlike get.SCOPE(), SPART's LeafModel
# is passed straight through to foursail(), which is a real dispatcher).
# ==============================================================================

rm(list = ls())  # avoid leftover objects from a previous run/session leaking in

if (!requireNamespace("ToolsRTM", quietly = TRUE)) {
  remotes::install_gitlab("caminoccg/toolsrtm")
}
library(ToolsRTM)
library(ggplot2)
library(dplyr)

out_dir <- "outs/Comparison"  # project-level outputs folder, never inside Scripts/
dir.create(out_dir, showWarnings = FALSE)
pdf(file.path(out_dir, "Rplots.pdf"))  # catches any stray plot()/print() call, keeps it out of Scripts/

## ----------------------------------------------------------------------------
## 1. Build a single-row LUT from ToolsRTM's own reference ranges (inputsSPART)
## ----------------------------------------------------------------------------

LUT <- as.data.frame(getLUT(inputs = ToolsRTM::inputsSPART, nLUT = 1, setseed = 1234))

# inputsSPART only ships PROSPECT-PRO/-D columns -- add the extra columns
# Liberty and Fluspect-B/-Cx need too (same requirement as foursail(); see
# Scripts/compare_RTM_models.R for where these typical values come from).
LUT$Cs <- 0; LUT$fqe <- 0.01; LUT$Cx <- 0
LUT$cell.d <- 40; LUT$inter.c <- 0.045; LUT$baseline.abs <- 0.0006
LUT$leaf.thick <- 1.6; LUT$albino.abs <- 0; LUT$lign.cell <- 2; LUT$Nitrogen <- 1

cat("LUT columns:", paste(names(LUT), collapse = ", "), "\n")

leaf_models <- c("PROSPECT-PRO", "PROSPECT-D", "Liberty", "Fluspect-B", "Fluspect-B-Cx")

## ----------------------------------------------------------------------------
## 2. Run SPART for each leaf model (TerraAqua.MODIS sensor, default irradiance)
## ----------------------------------------------------------------------------

cat("=== Running SPART() across all 5 leaf models (sensor = TerraAqua.MODIS) ===\n")

results <- lapply(leaf_models, function(lm) {
  cat(sprintf("LeafModel = %-14s ... ", lm))
  out <- tryCatch(
    suppressMessages(suppressWarnings(SPART(
      inputLUT = LUT, CanopyModel = "fourSAIL", LeafModel = lm,
      sensor.i = ToolsRTM::TerraAqua.MODIS, get.plots = FALSE
    ))),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(out)) return(NULL)
  cat("OK\n")
  out
})
names(results) <- leaf_models
ok <- !sapply(results, is.null)
cat(sprintf("\n%d/%d leaf models ran successfully: %s\n",
            sum(ok), length(ok), paste(names(results)[ok], collapse = ", ")))

## ----------------------------------------------------------------------------
## 3. Confirm LeafModel actually changes the TOC reflectance (unlike get.SCOPE())
## ----------------------------------------------------------------------------

toc_list <- results[ok]
cat("\nstr() of one result (to see what SPART() actually returns):\n")
str(toc_list[[1]], max.level = 1)

## ----------------------------------------------------------------------------
## 4. Plot TOC BRDF reflectance for every leaf model -- confirms LeafModel
##    genuinely changes the output here (unlike get.SCOPE(), see
##    Scripts/compare_SCOPE_models.R)
## ----------------------------------------------------------------------------

df_toc <- bind_rows(lapply(names(toc_list), function(lm) {
  d <- toc_list[[lm]]$rfl.toc.brdf
  data.frame(wavelength = d[[1]], reflectance = d[[2]], leaf_model = lm)
}))

leaf_colors <- c("PROSPECT-PRO" = "#0072B2", "PROSPECT-D" = "#56B4E9", "Liberty" = "#009E73",
                 "Fluspect-B" = "#E69F00", "Fluspect-B-Cx" = "#D55E00")

p_toc <- ggplot(df_toc, aes(x = wavelength, y = reflectance, color = leaf_model)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = leaf_colors) +
  labs(title = "SPART(): TOC BRDF reflectance across leaf models",
       subtitle = "Same LUT, CanopyModel = fourSAIL, sensor = TerraAqua.MODIS",
       x = "Wavelength (nm)", y = "TOC BRDF reflectance", color = "Leaf model") +
  theme_bw(base_size = 12) + theme(legend.position = "top")

out_dir <- "outs/Comparison"  # project-level outputs folder, never inside Scripts/
ggsave(file.path(out_dir, "compare_SPART_leafmodels.png"), p_toc, width = 8, height = 5, dpi = 150)
cat("\nSaved to '", out_dir, "/compare_SPART_leafmodels.png'\n", sep = "")

dev.off()  # close the Rplots.pdf capture device opened near the top
