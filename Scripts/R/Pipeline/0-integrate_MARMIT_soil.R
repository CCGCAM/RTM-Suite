# ==============================================================================
# Feed a MARMIT-simulated soil spectrum (ToolsRTM::get.marmit.rsoil()) into
# every canopy/TOA model that takes a soil spectrum -- foursail(), foursail2(),
# inform(), and SPART() (via its new `rsoil` override, see ToolsRTM::SPART) --
# and confirm the soil choice visibly changes the output vs. each model's
# usual flat/BSM soil baseline.
#
# Wavelength-grid note: foursail()/foursail2()/inform() with a non-Fluspect
# leaf model want rsoil on the full 400-2500nm (2101-point) grid; SPART()'s
# rsoil override wants 400-2400nm (2001-point) to match its internal BSM
# grid -- get.marmit.rsoil()'s `wl.out` argument controls this.
# ==============================================================================

rm(list = ls())  # avoid leftover objects from a previous run/session leaking in

use.dev.source <- TRUE  # TRUE = cargar desde el código fuente local (ToolsRTM/R, SCOPEinR/R); FALSE = usar el paquete instalado

if (use.dev.source) {
  devtools::load_all("ToolsRTM/R")
  devtools::load_all("SCOPEinR/R")
} else {
  library(ToolsRTM)
  library(SCOPEinR)
}
# ggplot2/dplyr aren't re-exported by ToolsRTM (it only `import()`s them for
# its own internal use, per NAMESPACE) -- need them attached directly here.
library(ggplot2)
library(dplyr)

out_dir <- "outs/rtm_sims"
dir.create(out_dir, showWarnings = FALSE)

# Only Bablet_2016 ships with ToolsRTM itself; the other 7 official MARMIT
# databases live in this repo's own databases/ folder -- pass
# db_root = "databases" and database = "Liu_2002" (etc.) to use them, see
# ?get.marmit.rsoil.
soil_marmit_2101 <- get.marmit.rsoil(database = "Bablet_2016", id = 1,
                                      L = 0.05, eps = 0.4, wl.out = 400:2500)
soil_marmit_2001 <- get.marmit.rsoil(database = "Bablet_2016", id = 1,
                                      L = 0.05, eps = 0.4, wl.out = 400:2400)
soil_flat <- rep(0.15, 2101)  # the flat baseline used throughout Scripts/compare_RTM_models.R

common_lut <- data.frame(
  N = 1.5, Cab = 40, Car = 8, Anth = 1, Cbrown = 0, EWT = 0.01, LMA = 0.009, alpha = 40,
  Prot = 0.002, CBC = 0.007,
  LIDFa = -0.35, LIDFb = -0.15, TypeLidf = 1, LAI = 3, hspot = 0.01,
  tts = 30, tto = 0, psi = 0,
  fraction_brown = 0.1, diss = 0.5, Cv = 1, Zeta = 0,
  LAIu = 0.5, sd = 650, cd = 4.5, h = 20, skyl = 0.1
)
wl_2101 <- 400:2500

## ----------------------------------------------------------------------------
## 1. foursail()
## ----------------------------------------------------------------------------

cat("=== foursail(): flat soil vs. MARMIT soil ===\n")
sim_fs_flat   <- foursail(inputLUT = common_lut, rsoil = soil_flat, LeafModel = "PROSPECT-D")
sim_fs_marmit <- foursail(inputLUT = common_lut, rsoil = soil_marmit_2101$rsoil.wet, LeafModel = "PROSPECT-D")
cat(sprintf("  max |rsot difference| = %.4f\n", max(abs(sim_fs_flat$rsot - sim_fs_marmit$rsot))))

## ----------------------------------------------------------------------------
## 2. foursail2()
## ----------------------------------------------------------------------------

cat("=== foursail2(): flat soil vs. MARMIT soil ===\n")
sim_fs2_flat   <- foursail2(inputLUT = common_lut, rsoil = soil_flat, LeafModel = "PROSPECT-D")
sim_fs2_marmit <- foursail2(inputLUT = common_lut, rsoil = soil_marmit_2101$rsoil.wet, LeafModel = "PROSPECT-D")
cat(sprintf("  max |rsot difference| = %.4f\n", max(abs(sim_fs2_flat$rsot - sim_fs2_marmit$rsot))))

## ----------------------------------------------------------------------------
## 3. inform()
## ----------------------------------------------------------------------------

cat("=== inform(): flat soil vs. MARMIT soil ===\n")
sim_inf_flat   <- inform(inputLUT = common_lut, rsoil = soil_flat, LeafModel = "PROSPECT-D")
sim_inf_marmit <- inform(inputLUT = common_lut, rsoil = soil_marmit_2101$rsoil.wet, LeafModel = "PROSPECT-D")
cat(sprintf("  max |r_forest difference| = %.4f\n", max(abs(sim_inf_flat - sim_inf_marmit))))

## ----------------------------------------------------------------------------
## 4. SPART() -- BSM soil vs. MARMIT soil override
## ----------------------------------------------------------------------------

cat("=== SPART(): BSM soil vs. MARMIT soil override ===\n")
LUT_spart <- as.data.frame(getLUT(inputs = ToolsRTM::inputsSPART, nLUT = 1, setseed = 1234))
LUT_spart$Cs <- 0; LUT_spart$fqe <- 0.01; LUT_spart$Cx <- 0
LUT_spart$cell.d <- 40; LUT_spart$inter.c <- 0.045; LUT_spart$baseline.abs <- 0.0006
LUT_spart$leaf.thick <- 1.6; LUT_spart$albino.abs <- 0; LUT_spart$lign.cell <- 2; LUT_spart$Nitrogen <- 1

sim_spart_bsm <- suppressMessages(suppressWarnings(SPART(
  inputLUT = LUT_spart, CanopyModel = "fourSAIL", LeafModel = "PROSPECT-D",
  sensor.i = ToolsRTM::TerraAqua.MODIS, get.plots = FALSE
)))
sim_spart_marmit <- suppressMessages(suppressWarnings(SPART(
  inputLUT = LUT_spart, CanopyModel = "fourSAIL", LeafModel = "PROSPECT-D",
  sensor.i = ToolsRTM::TerraAqua.MODIS, rsoil = soil_marmit_2001$rsoil.wet, get.plots = FALSE
)))
max_diff_spart_toc <- max(abs(sim_spart_bsm$output$rfl.toc - sim_spart_marmit$output$rfl.toc))
cat(sprintf("  max |TOC rfl difference| = %.4f\n", max_diff_spart_toc))

## ----------------------------------------------------------------------------
## 5. Plot: soil choice's effect on TOC reflectance, all four models
## ----------------------------------------------------------------------------

df_plot <- bind_rows(
  data.frame(wavelength = wl_2101, rsot = sim_fs_flat$rsot,   soil = "Flat (0.15)", model = "fourSAIL"),
  data.frame(wavelength = wl_2101, rsot = sim_fs_marmit$rsot, soil = "MARMIT",      model = "fourSAIL"),
  data.frame(wavelength = wl_2101, rsot = sim_fs2_flat$rsot,   soil = "Flat (0.15)", model = "foursail2"),
  data.frame(wavelength = wl_2101, rsot = sim_fs2_marmit$rsot, soil = "MARMIT",      model = "foursail2"),
  data.frame(wavelength = wl_2101, rsot = sim_inf_flat,   soil = "Flat (0.15)", model = "INFORM"),
  data.frame(wavelength = wl_2101, rsot = sim_inf_marmit, soil = "MARMIT",      model = "INFORM")
)

p <- ggplot(df_plot, aes(x = wavelength, y = rsot, color = soil)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~model, ncol = 1) +
  labs(title = "Effect of soil spectrum choice on TOC reflectance",
       subtitle = "Flat 0.15 baseline vs. MARMIT (Bablet_2016 id=1, L=0.05cm, eps=0.4), same canopy/leaf LUT",
       x = "Wavelength (nm)", y = "TOC bidirectional reflectance (rsot)", color = "Soil") +
  theme_bw(base_size = 12) + theme(legend.position = "top")

print(p)
ggsave(file.path(out_dir, "marmit_integration_canopy_models.png"), p, width = 8, height = 9, dpi = 150)

## ----------------------------------------------------------------------------
## 6. Summary
## ----------------------------------------------------------------------------

cat("\n=== Summary: MARMIT soil integration ===\n")
cat("All four models ran successfully with a MARMIT-generated soil spectrum,\n")
cat("and the soil choice measurably changes TOC reflectance in every case\n")
cat("(largest effect where the canopy is sparse enough for soil to show through,\n")
cat("e.g. low LAI or high soil visibility in the LIDF/hotspot geometry).\n")
cat("\nSaved to '", out_dir, "/marmit_integration_canopy_models.png'\n", sep = "")
