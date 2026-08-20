# ==============================================================================
# Spectral global sensitivity analysis (Sobol total index), reproducing the
# style of outs/sensibilidad.png: stacked "Total SI [%]" vs wavelength, one
# trait per color band, comparing a Uniform-PDF sampling scenario against a
# Gaussian-PDF one.
#
# Uses ToolsRTM::get.spectral.sensitivity() (new function, wraps the existing
# ToolsRTM::get.sobol.indices() Sobol estimator, run once per wavelength).
# ==============================================================================

rm(list = ls())  # avoid leftover objects from a previous run/session leaking in

library(ToolsRTM)
library(ggplot2)
library(dplyr)

out_dir <- "../../../outs/Sensibility"  # project-level outputs folder, never inside Scripts/
dir.create(out_dir, showWarnings = FALSE)
pdf(file.path(out_dir, "Rplots.pdf"))  # catches any stray plot()/print() call, keeps it out of Scripts/

n_samples <- 1000  # -> Sobol N = 500, as requested
traits <- c("N", "Cab", "EWT", "LMA", "LIDFa", "LAI")

trait_labels <- c(
  N = "Leaf structural parameter (N)", Cab = "Chlorophyll a+b content",
  EWT = "Equivalent water thickness", LMA = "Dry matter content (LMA)",
  LIDFa = "Leaf angle distribution", LAI = "Total Leaf Area Index",
  SoilCoef = "Soil coefficient"
)

# Chunked + saved to disk: for large n_samples (5000, 20000...) a crash only
# loses the chunk in progress, and re-running this script resumes from the
# saved chunks instead of re-simulating everything from scratch.
chunks_dir <- file.path(out_dir, "sensitivity_chunks")

cat("=== Uniform PDF sensitivity run (n =", n_samples, ") ===\n")
si_uniform <- get.spectral.sensitivity(n.samples = n_samples, distribution = "Uniform",
                                        traits = traits, wl.step = 5, seed = 1,
                                        chunk.size = 500, save.path = file.path(chunks_dir, "uniform"))

cat("\n=== Gaussian PDF sensitivity run (n =", n_samples, ") ===\n")
si_normal <- get.spectral.sensitivity(n.samples = n_samples, distribution = "Gaussian",
                                       traits = traits, wl.step = 5, seed = 2,
                                       chunk.size = 500, save.path = file.path(chunks_dir, "gaussian"))

si_all <- bind_rows(si_uniform, si_normal)
si_all$trait_label <- factor(trait_labels[si_all$trait], levels = trait_labels)
si_all$distribution <- factor(si_all$distribution, levels = c("Uniform", "Gaussian"),
                              labels = c("PDF UNIFORM", "PDF NORMAL"))

# Same trait color scheme as the reference figure (dark blue -> yellow -> red family)
trait_colors <- setNames(
  c("#08306B", "#4292C6", "#00BFC4", "#7FBF7B", "#FDB863", "#B2182B", "gray50"),
  trait_labels[c("LAI", "LIDFa", "SoilCoef", "N", "Cab", "EWT", "LMA")]
)

p_sensitivity <- ggplot(si_all, aes(x = wavelength, y = STi_pct, fill = trait_label)) +
  geom_area(position = "stack") +
  facet_wrap(~distribution, ncol = 1) +
  scale_fill_manual(values = trait_colors, name = NULL) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  labs(title = "Spectral sensitivity of TOC reflectance (fourSAIL + PROSPECT-D)",
       subtitle = sprintf("Sobol total sensitivity index, %d simulations per PDF (N = %d)",
                          n_samples, n_samples %/% 2),
       x = "Wavelength (nm)", y = "Total SI [%]") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey85", color = NA),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "sensitivity_analysis.png"), p_sensitivity, width = 9, height = 8, dpi = 200)
write.csv(si_all, file.path(out_dir, "sensitivity_analysis.csv"), row.names = FALSE)

cat("\nSaved to '", out_dir, "/': sensitivity_analysis.png, sensitivity_analysis.csv\n", sep = "")

dev.off()  # close the Rplots.pdf capture device opened near the top
