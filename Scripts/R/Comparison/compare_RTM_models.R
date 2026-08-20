# ==============================================================================
# Compare canopy models (fourSAIL, foursail2, INFORM) and leaf models
# (PROSPECT-PRO, PROSPECT-D, Liberty, Fluspect-B, Fluspect-B-Cx) side by side,
# using the SAME LUT (structural/geometric parameters identical, only the
# parameters each specific model actually needs differ).
#
# Three things are produced:
#   1. Spectral comparison plots (canopy models; leaf models; full grid).
#   2. A difference/agreement heatmap: deviation of each model combination
#      from the ensemble mean at every wavelength -- shows at a glance where
#      models agree (near white) and where they diverge (saturated color).
#   3. A numeric RMSE table (pairwise + vs. ensemble mean), printed and saved
#      as CSV, so agreement/disagreement isn't just visual.
#
# Output metric: rsot -- the bi-directional reflectance factor (BRDF term),
# i.e. the actual top-of-canopy reflectance seen by a sensor for the chosen
# sun/view geometry.
# ==============================================================================

rm(list = ls())  # avoid leftover objects from a previous run/session leaking in

library(ToolsRTM)
library(ggplot2)
library(dplyr)
library(tidyr)


devtools::load_all("../../../ToolsRTM")
devtools::load_all("../../../SCOPEinR")

out_dir <- "../../../outs/Comparison"  # project-level outputs folder, never inside Scripts/
dir.create(out_dir, showWarnings = FALSE)
pdf(file.path(out_dir, "Rplots.pdf"))  # catches any stray plot()/print() call, keeps it out of Scripts/

# Consistent palette across every plot in this script (colorblind-safe, from
# the R "Okabe-Ito" set) -- extend if more leaf models are added later.
leaf_colors <- c(
  "PROSPECT-PRO"  = "#0072B2",
  "PROSPECT-D"    = "#56B4E9",
  "Liberty"       = "#009E73",
  "Fluspect-B"    = "#E69F00",
  "Fluspect-B-Cx" = "#D55E00"
)
canopy_colors <- c("fourSAIL" = "#0072B2", "foursail2" = "#009E73", "INFORM" = "#D55E00")

theme_rtm <- theme_bw(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90", color = NA),
        plot.title = element_text(face = "bold"))

## ----------------------------------------------------------------------------
## 1. Shared LUT -- same geometry/structure for every model, plus every
##    column any of the 5 leaf models could need (unused columns are simply
##    ignored by whichever model doesn't need them).
## ----------------------------------------------------------------------------

common_lut <- data.frame(
  # Leaf biochemistry (PROSPECT-PRO/-D)
  N = 1.5, Cab = 40, Car = 8, Anth = 1, Cbrown = 0, EWT = 0.01, LMA = 0.009, alpha = 40,
  Prot = 0.002, CBC = 0.007,
  # Fluspect-B/-Cx extras
  Cs = 0, fqe = 0.01, Cx = 0,
  # Liberty
  cell.d = 40, inter.c = 0.045, baseline.abs = 0.0006, leaf.thick = 1.6,
  albino.abs = 0, lign.cell = 2, Nitrogen = 1,
  # Canopy structure / geometry (fourSAIL, foursail2, INFORM)
  LIDFa = -0.35, LIDFb = -0.15, TypeLidf = 1, LAI = 3, hspot = 0.01,
  tts = 30, tto = 0, psi = 0,
  # foursail2-only (two-layer green/brown canopy)
  fraction_brown = 0.1, diss = 0.5, Cv = 1, Zeta = 0,
  # INFORM-only (forest structure)
  LAIu = 0.5, sd = 650, cd = 4.5, h = 20, skyl = 0.1
)

rsoil_2101 <- rep(0.15, 2101)  # flat soil reflectance, 400-2500 nm (PROSPECT/Liberty range)
rsoil_2001 <- rsoil_2101[1:2001]  # 400-2400 nm (Fluspect/INFORM range)
wl_2101 <- 400:2500
wl_2001 <- 400:2400

leaf_models <- c("PROSPECT-PRO", "PROSPECT-D", "Liberty", "Fluspect-B", "Fluspect-B-Cx")
canopy_models <- c("fourSAIL", "foursail2", "INFORM")

## ----------------------------------------------------------------------------
## 2. Comparison A: canopy models, fixed leaf model (PROSPECT-D)
## ----------------------------------------------------------------------------

cat("=== Canopy model comparison (leaf model = PROSPECT-D) ===\n")

sim_foursail <- foursail(inputLUT = common_lut, rsoil = rsoil_2101, LeafModel = "PROSPECT-D")
sim_foursail2 <- foursail2(inputLUT = common_lut, rsoil = rsoil_2101, LeafModel = "PROSPECT-D")
sim_inform <- inform(inputLUT = common_lut, rsoil = rsoil_2101, LeafModel = "PROSPECT-D")

df_canopy <- bind_rows(
  data.frame(wavelength = wl_2101, rsot = sim_foursail$rsot, model = "fourSAIL"),
  data.frame(wavelength = wl_2101, rsot = sim_foursail2$rsot, model = "foursail2"),
  data.frame(wavelength = wl_2101, rsot = sim_inform, model = "INFORM")  # inform() returns r_forest directly
)

p_canopy <- ggplot(df_canopy, aes(x = wavelength, y = rsot, color = model)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = canopy_colors) +
  labs(title = "Canopy model comparison", subtitle = "Leaf model = PROSPECT-D, same LUT for all three",
       x = "Wavelength (nm)", y = "TOC bidirectional reflectance (rsot)", color = "Canopy model") +
  theme_rtm

## ----------------------------------------------------------------------------
## 3. Comparison B: leaf models, fixed canopy model (fourSAIL -- supports all 5)
## ----------------------------------------------------------------------------

cat("\n=== Leaf model comparison (canopy model = fourSAIL) ===\n")

df_leaf <- bind_rows(lapply(leaf_models, function(lm) {
  short_domain <- lm %in% c("Fluspect-B", "Fluspect-B-Cx")
  sim <- foursail(inputLUT = common_lut,
                   rsoil = if (short_domain) rsoil_2001 else rsoil_2101,
                   LeafModel = lm, spectrum.all = !short_domain)
  data.frame(wavelength = if (short_domain) wl_2001 else wl_2101,
             rsot = sim$rsot, leaf_model = lm)
}))

p_leaf <- ggplot(df_leaf, aes(x = wavelength, y = rsot, color = leaf_model)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = leaf_colors) +
  labs(title = "Leaf model comparison", subtitle = "Canopy model = fourSAIL, same LUT for all five",
       x = "Wavelength (nm)", y = "TOC bidirectional reflectance (rsot)", color = "Leaf model") +
  theme_rtm

## ----------------------------------------------------------------------------
## 4. Comparison C: leaf models x canopy models that support them (grid)
## ----------------------------------------------------------------------------

cat("\n=== Full grid: leaf model x canopy model ===\n")

df_grid <- bind_rows(lapply(canopy_models, function(cm) {
  bind_rows(lapply(leaf_models, function(lm) {
    short_domain <- lm %in% c("Fluspect-B", "Fluspect-B-Cx")
    # foursail()/foursail2() expect a pre-truncated rsoil for Fluspect (400-2400nm);
    # inform() always wants the full 400-2500nm rsoil and truncates internally --
    # passing it the short one causes an internal length-recycling mismatch.
    rsoil <- if (short_domain && cm != "INFORM") rsoil_2001 else rsoil_2101
    sim <- tryCatch({
      if (cm == "fourSAIL") {
        foursail(inputLUT = common_lut, rsoil = rsoil, LeafModel = lm, spectrum.all = !short_domain)$rsot
      } else if (cm == "foursail2") {
        foursail2(inputLUT = common_lut, rsoil = rsoil, LeafModel = lm)$rsot
      } else {
        inform(inputLUT = common_lut, rsoil = rsoil, LeafModel = lm)
      }
    }, error = function(e) { message(cm, " + ", lm, " failed: ", conditionMessage(e)); NULL })
    if (is.null(sim)) return(NULL)
    data.frame(wavelength = if (short_domain) wl_2001 else wl_2101,
               rsot = sim, leaf_model = lm, canopy_model = cm)
  }))
}))
df_grid$combo <- paste(df_grid$canopy_model, df_grid$leaf_model, sep = " + ")

p_grid <- ggplot(df_grid, aes(x = wavelength, y = rsot, color = leaf_model)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~canopy_model) +
  scale_color_manual(values = leaf_colors) +
  labs(title = "All canopy models x all leaf models", subtitle = "Same LUT throughout",
       x = "Wavelength (nm)", y = "TOC bidirectional reflectance (rsot)", color = "Leaf model") +
  theme_rtm

print(p_grid)
## ----------------------------------------------------------------------------
## 5. Agreement/difference heatmap: deviation from the ensemble mean at each
##    wavelength, common domain only (400-2400nm, where all combos overlap)
## ----------------------------------------------------------------------------

cat("\n=== Agreement heatmap (deviation from ensemble mean) ===\n")

df_common <- df_grid[df_grid$wavelength <= 2400, ]
df_common <- df_common %>%
  group_by(wavelength) %>%
  mutate(ensemble_mean = mean(rsot), deviation = rsot - ensemble_mean) %>%
  ungroup()

p_heatmap <- ggplot(df_common, aes(x = wavelength, y = combo, fill = deviation)) +
  geom_raster() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                        name = "rsot - ensemble mean") +
  labs(title = "Where models agree and where they don't",
       subtitle = "Deviation from the across-model ensemble mean at each wavelength (white = agreement)",
       x = "Wavelength (nm)", y = NULL) +
  theme_rtm + theme(legend.position = "right", axis.text.y = element_text(size = 8))

print(p_heatmap)
## ----------------------------------------------------------------------------
## 6. Numeric agreement table: RMSE of each combo vs. the ensemble mean, and
##    the full pairwise RMSE matrix between all combos
## ----------------------------------------------------------------------------

rmse_vs_mean <- df_common %>%
  group_by(combo) %>%
  summarise(RMSE_vs_ensemble_mean = sqrt(mean(deviation^2)), .groups = "drop") %>%
  arrange(desc(RMSE_vs_ensemble_mean))

cat("\nRMSE of each model combination vs. the ensemble mean (largest first):\n")
print(rmse_vs_mean, n = Inf)

wide <- df_common %>% select(wavelength, combo, rsot) %>% pivot_wider(names_from = combo, values_from = rsot)
combo_names <- setdiff(names(wide), "wavelength")
pairwise_rmse <- outer(combo_names, combo_names, Vectorize(function(a, b) {
  sqrt(mean((wide[[a]] - wide[[b]])^2))
}))
dimnames(pairwise_rmse) <- list(combo_names, combo_names)

cat("\nPairwise RMSE matrix between all model combinations:\n")
print(round(pairwise_rmse, 4))

write.csv(rmse_vs_mean, file.path(out_dir, "rmse_vs_ensemble_mean.csv"), row.names = FALSE)
write.csv(as.data.frame(pairwise_rmse), file.path(out_dir, "pairwise_rmse_matrix.csv"))

## ----------------------------------------------------------------------------
## 7. Save plots
## ----------------------------------------------------------------------------

ggsave(file.path(out_dir, "compare_canopy_models.png"), p_canopy, width = 8, height = 5, dpi = 150)
ggsave(file.path(out_dir, "compare_leaf_models.png"), p_leaf, width = 8, height = 5, dpi = 150)
ggsave(file.path(out_dir, "compare_grid_canopy_x_leaf.png"), p_grid, width = 10, height = 6, dpi = 150)
ggsave(file.path(out_dir, "agreement_heatmap.png"), p_heatmap, width = 9, height = 5, dpi = 150)

cat("\nSaved to '", out_dir, "/': compare_canopy_models.png, compare_leaf_models.png, ",
    "compare_grid_canopy_x_leaf.png, agreement_heatmap.png, ",
    "rmse_vs_ensemble_mean.csv, pairwise_rmse_matrix.csv\n", sep = "")

dev.off()  # close the Rplots.pdf capture device opened near the top
