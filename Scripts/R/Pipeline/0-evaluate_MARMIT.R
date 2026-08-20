# ==============================================================================
# Evaluate ToolsRTM::get.marmit.rsoil() -- the MARMIT soil reflectance model
# -- on its own, before feeding it into any canopy/SPART model:
#   1. Wet vs. dry reflectance for a range of surface-wetness levels (eps),
#      at fixed water-layer thickness.
#   2. Effect of water-layer thickness (L) at fixed wetness.
#   3. MARMIT-1 vs. MARMIT-2 comparison (same soil, same L/eps).
#   4. Sanity checks: wet reflectance <= dry reflectance at every wavelength
#      (water darkens soil -- this should always hold), and the estimated
#      soil moisture content (SMC) increases monotonically with eps*L.
# ==============================================================================

rm(list = ls())  # avoid leftover objects from a previous run/session leaking in

use.dev.source <- TRUE  # TRUE = cargar desde el código fuente local (../ToolsRTM, ../SCOPEinR); FALSE = usar el paquete instalado

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

soil_database <- "Bablet_2016"
soil_id <- 1

## ----------------------------------------------------------------------------
## 1. Wetness sweep (eps), fixed L
## ----------------------------------------------------------------------------

cat("=== MARMIT-1: wetness (eps) sweep, database =", soil_database, "id =", soil_id, "===\n")

eps_levels <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
L_fixed <- 0.05

df_eps <- bind_rows(lapply(eps_levels, function(e) {
  soil <- get.marmit.rsoil(database = soil_database, id = soil_id, version = "marmit1",
                            L = L_fixed, eps = e)
  data.frame(wavelength = soil$wavelength, reflectance = soil$rsoil.wet,
             eps = e, SMC = round(soil$SMC, 1))
}))
df_eps$label <- sprintf("eps = %.1f (SMC = %.1f%%)", df_eps$eps, df_eps$SMC)

p_eps <- ggplot(df_eps, aes(x = wavelength, y = reflectance, color = label)) +
  geom_line(linewidth = 0.7) +
  labs(title = "MARMIT-1: soil reflectance vs. surface wetness",
       subtitle = sprintf("%s, soil ID %d, water layer L = %.2fcm", soil_database, soil_id, L_fixed),
       x = "Wavelength (nm)", y = "Soil reflectance", color = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "right")

print(p_eps)

## ----------------------------------------------------------------------------
## 2. Water-layer thickness sweep (L), fixed eps
## ----------------------------------------------------------------------------

cat("=== MARMIT-1: water-layer thickness (L) sweep ===\n")

L_levels <- c(0.01, 0.05, 0.1, 0.2, 0.35, 0.5)
eps_fixed <- 0.5

df_L <- bind_rows(lapply(L_levels, function(l) {
  soil <- get.marmit.rsoil(database = soil_database, id = soil_id, version = "marmit1",
                            L = l, eps = eps_fixed)
  data.frame(wavelength = soil$wavelength, reflectance = soil$rsoil.wet,
             L = l, SMC = round(soil$SMC, 1))
}))
df_L$label <- sprintf("L = %.2fcm (SMC = %.1f%%)", df_L$L, df_L$SMC)

p_L <- ggplot(df_L, aes(x = wavelength, y = reflectance, color = label)) +
  geom_line(linewidth = 0.7) +
  labs(title = "MARMIT-1: soil reflectance vs. water-layer thickness",
       subtitle = sprintf("%s, soil ID %d, wet surface fraction eps = %.1f", soil_database, soil_id, eps_fixed),
       x = "Wavelength (nm)", y = "Soil reflectance", color = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "right")

print(p_L)

## ----------------------------------------------------------------------------
## 3. MARMIT-1 vs. MARMIT-2
## ----------------------------------------------------------------------------

cat("=== MARMIT-1 vs. MARMIT-2 ===\n")

soil_m1 <- get.marmit.rsoil(database = soil_database, id = soil_id, version = "marmit1",
                             L = 0.05, eps = 0.4)
soil_m2 <- get.marmit.rsoil(database = soil_database, id = soil_id, version = "marmit2",
                             L = 0.05, eps = 0.4, n_i = 1.53, k_i = 0.001, d_i = 0.0005)

df_version <- bind_rows(
  data.frame(wavelength = soil_m1$wavelength, reflectance = soil_m1$rsoil.wet, version = "MARMIT-1"),
  data.frame(wavelength = soil_m2$wavelength, reflectance = soil_m2$rsoil.wet, version = "MARMIT-2"),
  data.frame(wavelength = soil_m1$wavelength, reflectance = soil_m1$rsoil.dry, version = "Dry reference")
)

p_version <- ggplot(df_version, aes(x = wavelength, y = reflectance, color = version)) +
  geom_line(linewidth = 0.7) +
  labs(title = "MARMIT-1 vs. MARMIT-2 (same soil, L = 0.05cm, eps = 0.4)",
       subtitle = sprintf("%s, soil ID %d", soil_database, soil_id),
       x = "Wavelength (nm)", y = "Soil reflectance", color = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "top")

print(p_version)
cat(sprintf("MARMIT-1 SMC = %.1f%%, MARMIT-2 SMC = %.1f%% (same L*eps -> same SMC, SMC doesn't depend on version)\n",
            soil_m1$SMC, soil_m2$SMC))

## ----------------------------------------------------------------------------
## 4. Sanity checks
## ----------------------------------------------------------------------------

cat("\n=== Sanity checks ===\n")

# 4a. Wetter soil should never be brighter than dry soil, at any wavelength
violations <- sum(df_eps$reflectance[df_eps$eps > 0] >
                     rep(df_eps$reflectance[df_eps$eps == 0], sum(eps_levels > 0)))
cat(sprintf("Wet-brighter-than-dry violations (should be 0): %d / %d points\n",
            violations, sum(df_eps$eps > 0) * length(unique(df_eps$wavelength)) / length(eps_levels[eps_levels > 0])))

# 4b. SMC should increase monotonically with eps (fixed L)
smc_by_eps <- df_eps %>% distinct(eps, SMC) %>% arrange(eps)
cat("SMC vs. eps (should be non-decreasing):\n")
print(smc_by_eps)
cat("Monotonic:", !is.unsorted(smc_by_eps$SMC), "\n")

# 4c. SMC should increase monotonically with L (fixed eps)
smc_by_L <- df_L %>% distinct(L, SMC) %>% arrange(L)
cat("\nSMC vs. L (should be non-decreasing):\n")
print(smc_by_L)
cat("Monotonic:", !is.unsorted(smc_by_L$SMC), "\n")

## ----------------------------------------------------------------------------
## 5. Save
## ----------------------------------------------------------------------------

ggsave(file.path(out_dir, "marmit_wetness_sweep.png"), p_eps, width = 8, height = 5, dpi = 150)
ggsave(file.path(out_dir, "marmit_thickness_sweep.png"), p_L, width = 8, height = 5, dpi = 150)
ggsave(file.path(out_dir, "marmit_v1_vs_v2.png"), p_version, width = 8, height = 5, dpi = 150)
write.csv(df_eps, file.path(out_dir, "marmit_wetness_sweep.csv"), row.names = FALSE)

cat("\nSaved to '", out_dir, "/': marmit_wetness_sweep.png, marmit_thickness_sweep.png, ",
    "marmit_v1_vs_v2.png, marmit_wetness_sweep.csv\n", sep = "")
