# ----------------------------------------------------------------------------
# Inspect SCOPE output folder: list CSV files and their dimensions
# ----------------------------------------------------------------------------
# OPTIONAL: clear environment (comment out if you don't like this behaviour)
rm(list = ls())

library(dplyr)
library(readr)   # faster than base::read.csv but you can use read.csv if you prefer
library(ggplot2)
library(tidyr)
library(SCOPEinR)
library(ToolsRTM)


# Path to ONE SCOPE run (change this to your run folder)

# Get the last SCOPE output folder (most recent run) -- project-level outputs
# folder, never inside Scripts/
path_out = '../../outs/rtm_sims/scope_runs/'
subdirs          <- list.dirs(path_out, recursive = FALSE)
last_subdir      <- subdirs[length(subdirs)]
cat("\nUsing folder for plots:", last_subdir, "\n")
path_run <- last_subdir   # <- adjust!

# List all CSV files in the main folder (not Parameters / Plots)
csv_files <- list.files(path_run, pattern = "\\.csv$", full.names = TRUE)
print(csv_files)

# ----------------------------------------------------------------------------
# 0. Load LUT from Parameters/inputLUT.csv and plot trait variability
# ----------------------------------------------------------------------------


# Path to LUT
lut_file <- file.path(path_run, "Parameters", "inputLUT.csv")

lut <- read.csv(lut_file)
head(lut)
# Choose traits to inspect (students can modify this vector)
traits_to_plot <- c("Cab", "Car", "Anth", "LAI", 'EWT',"Vcmax25" )

# Keep only traits that actually exist in the LUT
traits_to_plot <- traits_to_plot[traits_to_plot %in% names(lut)]

if (length(traits_to_plot) == 0) {
  stop("None of the selected traits are present in inputLUT.csv")
}

# Long format for ggplot
lut_long <- lut |>
  select(all_of(traits_to_plot)) |>
  pivot_longer(cols = everything(),
               names_to = "Trait",
               values_to = "Value")

# Basic boxplot per trait
ggplot(lut_long, aes(x = Trait, y = Value)) +
  geom_boxplot() +
  theme_bw() +
  ylab("Value") +
  ggtitle("Variability of SCOPE input traits in LUT")


# ----------------------------------------------------------------------------
# Density and histogram plots per trait (helps visualise distributions)
# ----------------------------------------------------------------------------

ggplot(lut_long, aes(x = Value)) +
  geom_density(fill = "grey80") +
  facet_wrap(~ Trait, scales = "free") +
  theme_bw() +
  xlab("Value") +
  ggtitle("Distribution of LUT input traits")



# Histogram (barplot) per trait
ggplot(lut_long, aes(x = Value)) +
  geom_histogram(fill = "grey70", color = "black", bins = 20) +
  facet_wrap(~ Trait, scales = "free") +
  theme_bw() +
  xlab("Value") +
  ylab("Frequency") +
  ggtitle("Histogram of LUT input traits (SCOPE simulations)")



# ----------------------------------------------------------------------------
# 1. Reflectance components: rdd, rdo, rsd, rso, refl, reflapp
# ----------------------------------------------------------------------------

# rdd  = diffuse-in  / diffuse-out reflectance
#        Light arriving as diffuse radiation (sky light)
#        and leaving the canopy also as diffuse radiation.
#        Describes multiple scattering inside the canopy.

# rdo  = diffuse-in  / directional-out reflectance
#        Diffuse incoming radiation (sky light)
#        leaving the canopy in a specific viewing direction.
#        Important for BRDF effects under overcast conditions.

# rsd  = direct-in (sun beam) / diffuse-out reflectance
#        Direct sunlight entering the canopy
#        and emerging as diffuse scattered light.
#        Shows sun–canopy interactions and internal scattering.

# rso  = direct-in (sun beam) / directional-out reflectance
#        Direct sun illumination leaving the canopy
#        in the exact viewing direction of the sensor.
#        Contains hot-spot effects and strong geometry interactions.

# refl = total canopy reflectance
#        The sum of all four components (rdd + rdo + rsd + rso).
#        This is the “true” physical reflectance of the canopy.

# reflapp = apparent reflectance
#          Reflectance derived from outgoing radiance at sensor level.
#          This accounts for illumination geometry and directional effects.
#          This is what a real sensor (e.g., Sentinel-2, PRISMA) observes.
# ----------------------------------------------------------------------------


# Load wavelength definitions from SCOPEinR
bands <- SCOPEinR::define.bands()
wlS   <- bands$wlS[1:2001]   # 400–2500 nm (shortwave)
wlF   <- bands$wlF           # 640–850 nm (fluorescence)


# ----------------------------------------------------------------------------
# 1.1 Read reflectance CSVs
# ----------------------------------------------------------------------------

refl    <- read.csv(file.path(path_run, "refl.csv"))
reflapp <- read.csv(file.path(path_run, "reflapp.csv"))
rdd     <- read.csv(file.path(path_run, "rdd.csv"))
rdo     <- read.csv(file.path(path_run, "rdo.csv"))
rsd     <- read.csv(file.path(path_run, "rsd.csv"))
rso     <- read.csv(file.path(path_run, "rso.csv"))

# Each row = one simulation, each column = wavelength (X1, X2, ...)

# ----------------------------------------------------------------------------
# 1.2 Plot mean reflectance spectrum (all simulations)
# ----------------------------------------------------------------------------

# Helper to convert one matrix-like data.frame to long format with wavelengths
to_long_spectrum <- function(df, wl, var_name) {
  df |>
    mutate(nsim = dplyr::row_number()) |>
    pivot_longer(cols = starts_with("X"),
                 names_to  = "band",
                 values_to = "value") |>
    mutate(Wavelength = rep(wl, times = length(unique(nsim))),
           Variable   = var_name)
}

refl_long <- to_long_spectrum(refl, wlS, "refl")

refl_mean <- refl_long |>
  group_by(Wavelength, Variable) |>
  summarize(mean_refl = mean(value, na.rm = TRUE),
            .groups = "drop")

ggplot(refl_mean, aes(x = Wavelength, y = mean_refl)) +
  geom_line() +
  theme_bw() +
  xlab("Wavelength (nm)") +
  ylab("Reflectance") +
  ggtitle("Mean canopy reflectance (refl) across all simulations") +
  xlim(400, 2500)

# ----------------------------------------------------------------------------
# 1.3 Compare reflectance components for one simulation (e.g. nsim = 1)
# ----------------------------------------------------------------------------

nsim_sel <- 1

extract_sim <- function(df, wl, var_name, nsim = 1) {
  # select row nsim and transpose to a vector
  v <- as.numeric(df[nsim, ])
  tibble(
    Wavelength = wl,
    Value      = v,
    Variable   = var_name
  )
}

df_rdd <- extract_sim(rdd, wlS, "rdd", nsim_sel)
df_rdo <- extract_sim(rdo, wlS, "rdo", nsim_sel)
df_rsd <- extract_sim(rsd, wlS, "rsd", nsim_sel)
df_rso <- extract_sim(rso, wlS, "rso", nsim_sel)
df_refl<- extract_sim(refl, wlS, "refl", nsim_sel)
df_reflApp<- extract_sim(reflapp, wlS, "reflapp", nsim_sel)
refl_components <- bind_rows(df_rdd, df_rdo, df_rsd, df_rso, df_refl, df_reflApp)

ggplot(refl_components, aes(x = Wavelength, y = Value, color = Variable)) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  xlab("Wavelength (nm)") +
  ylab("Reflectance") +
  ggtitle(paste("Reflectance components (simulation", nsim_sel, ")")) +
  xlim(400, 2500)

# ----------------------------------------------------------------------------
# 2. Radiance and fluorescence outputs
# ----------------------------------------------------------------------------
# Lo_spectrum.csv
#   Top-of-canopy radiance WITHOUT fluorescence.
#   Only reflected sunlight is included (no SIF signal).
#   Useful as the “reference” radiance.

# Lo_spectrum_includingF.csv
#   Top-of-canopy radiance WITH fluorescence included.
#   This is what a real spectrometer would measure (radiance + SIF).

# Eout_spectrum.csv
#   Hemispherically integrated upwelling shortwave radiation.
#   Total energy leaving the canopy across all directions (W m⁻² µm⁻¹).
# fluorescence.csv
#   Total fluorescence radiance at the top of the canopy (LoF).
#   Sum of all canopy fluorescence contributions.

# fluorescence_sunlit.csv
#   Fluorescence emitted by sunlit leaves (LoF_sunlit).
#   Typically the strongest contribution to canopy SIF.

# fluorescence_shaded.csv
#   Fluorescence emitted by shaded leaves (LoF_shaded).
#   Important in dense canopies and under low sun angles.

# fluorescence_scattered.csv
#   Scattered fluorescence inside the canopy (LoF_scattered).
#   Portion of SIF that is re-emitted after multiple scattering events.

# fluorescence_soil.csv
#   Soil fluorescence contribution (LoF_soil).
#   Usually small, but relevant at low LAI.
# ----------------------------------------------------------------------------

# fluorescence_All_Leaves.csv
#   Leaf-level emitted fluorescence from all leaves (Fem_leaves).
#   Pure emission from leaves before canopy interactions.

# fluorescence_hemis.csv
#   Hemispherically integrated fluorescence (EoutF).
#   Total fluorescence flux leaving the canopy into the hemisphere
#     (W m⁻² µm⁻¹), independent of viewing direction.

# fluorescence_ReabsCorr.csv
#   Hemispheric fluorescence corrected for re-absorption (EoutF_rc).
#   Removes the effect of leaves re-absorbing part of their own
#     fluorescence.
#   Closer to the “true emitted fluorescence” escaping the canopy.
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# 2.1 Read radiance and fluorescence CSVs
# ----------------------------------------------------------------------------

Lo_noF  <- read.csv(file.path(path_run, "Lo_spectrum.csv"))
Lo_withF<- read.csv(file.path(path_run, "Lo_spectrum_includingF.csv"))

F_tot        <- read.csv(file.path(path_run, "fluorescence.csv"))
F_sunlit     <- read.csv(file.path(path_run, "fluorescence_sunlit.csv"))
F_shaded     <- read.csv(file.path(path_run, "fluorescence_shaded.csv"))
F_scattered  <- read.csv(file.path(path_run, "fluorescence_scattered.csv"))
F_soil       <- read.csv(file.path(path_run, "fluorescence_soil.csv"))
F_all_leaves <- read.csv(file.path(path_run, "fluorescence_All_Leaves.csv"))
F_hemis      <- read.csv(file.path(path_run, "fluorescence_hemis.csv"))
F_reabs      <- read.csv(file.path(path_run, "fluorescence_ReabsCorr.csv"))

# ----------------------------------------------------------------------------
# 2.2 Compare Lo with and without fluorescence for one simulation
# ----------------------------------------------------------------------------

Lo_noF_sim   <- extract_sim(Lo_noF,  wlS, "Lo_noF",   nsim_sel)
Lo_withF_sim <- extract_sim(Lo_withF, wlS, "Lo_withF", nsim_sel)

Lo_compare <- bind_rows(Lo_noF_sim, Lo_withF_sim)

ggplot(Lo_compare, aes(x = Wavelength, y = Value, color = Variable)) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  xlab("Wavelength (nm)") +
  ylab("Radiance (W m-2 sr-1 µm-1)") +
  ggtitle(paste("Top-of-canopy radiance with / without fluorescence (sim", nsim_sel, ")")) +
  xlim(640, 850)

# ----------------------------------------------------------------------------
# 2.3 Fluorescence components for one simulation in the F window (640–850 nm)
# ----------------------------------------------------------------------------


F_tot_sim       <- extract_sim(F_tot,       wlF, "Total_F",     nsim_sel)
F_sunlit_sim    <- extract_sim(F_sunlit,    wlF, "Sunlit_F",    nsim_sel)
F_shaded_sim    <- extract_sim(F_shaded,    wlF, "Shaded_F",    nsim_sel)
F_scattered_sim <- extract_sim(F_scattered, wlF, "Scattered_F", nsim_sel)
F_soil_sim      <- extract_sim(F_soil,      wlF, "Soil_F",      nsim_sel)

F_components <- bind_rows(
  F_tot_sim, F_sunlit_sim, F_shaded_sim, F_scattered_sim, F_soil_sim
)

ggplot(F_components, aes(x = Wavelength, y = Value, color = Variable)) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  xlab("Wavelength (nm)") +
  ylab("Fluorescence radiance (W m-2 sr-1 µm-1)") +
  ggtitle(paste("Fluorescence components (simulation", nsim_sel, ")")) +
  xlim(640, 850)

# ----------------------------------------------------------------------------
# 2.4 Hemispheric fluorescence (EoutF and re-absorption corrected)
# ----------------------------------------------------------------------------

F_hemis_sim <- extract_sim(F_hemis, wlF, "EoutF",    nsim_sel)
F_reabs_sim <- extract_sim(F_reabs, wlF, "EoutFrc",  nsim_sel)

F_hemis_compare <- bind_rows(F_hemis_sim, F_reabs_sim)

ggplot(F_hemis_compare, aes(x = Wavelength, y = Value, color = Variable)) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  xlab("Wavelength (nm)") +
  ylab("Fluorescence (W m-2 µm-1)") +
  ggtitle(paste("Hemispheric fluorescence (sim", nsim_sel, ")")) +
  xlim(640, 850)


# ----------------------------------------------------------------------------
# 3. Irradiance inputs (what lights the canopy)
# ----------------------------------------------------------------------------
# Esun.csv
#   Spectral direct solar irradiance at the top of the canopy (Esun_).
#   → Collimated beam from the sun (W m⁻² µm⁻¹).
#   → Dominant under clear-sky conditions and strong sun.

# Esky.csv
#   Spectral diffuse sky irradiance at the top of the canopy (Esky_).
#   → Light scattered in the atmosphere, coming from all directions.
#   → Important under cloudy or hazy conditions.

# Eout_spectrum.csv
#   Spectral upwelling hemispheric radiation from the canopy (Eout_).
#   → Total shortwave energy leaving the canopy (W m⁻² µm⁻¹).
#   → Result of Esun + Esky interacting with soil + vegetation.


bands <- SCOPEinR::define.bands()
wlS   <- bands$wlS[1:2001]

Esun <- read.csv(file.path(path_run, "Esun.csv"))
Esky <- read.csv(file.path(path_run, "Esky.csv"))

nsim_sel <- 1  # escogemos una simulación para visualizar

Esun_sim <- as.numeric(Esun[nsim_sel, ])
Esky_sim <- as.numeric(Esky[nsim_sel, ])
Etot_sim <- Esun_sim + Esky_sim

irr_df <- tibble(
  Wavelength = wlS,
  Esun       = Esun_sim,
  Esky       = Esky_sim,
  Etot       = Etot_sim
) |>
  pivot_longer(cols = c(Esun, Esky, Etot),
               names_to = "Component",
               values_to = "Irradiance")

ggplot(irr_df, aes(x = Wavelength, y = Irradiance, color = Component)) +
  geom_line(linewidth = 0.5) +
  theme_bw() +
  xlab("Wavelength (nm)") +
  ylab("Irradiance (W m-2 µm-1)") +
  ggtitle(paste("Incoming irradiance (simulation", nsim_sel, ")")) +
  xlim(400, 2500)

# ----------------------------------------------------------------------------
# 4. Energy and carbon fluxes (fluxes.csv)
# ----------------------------------------------------------------------------
# This file contains the surface energy balance and CO₂ assimilation
# simulated by SCOPE.
#
# Each ROW   = one simulation in your LUT
# Each COLUMN = one energy or carbon flux from canopy or soil
#
# ----------------------------- CANOPY FLUXES --------------------------------
# Rnctot  = Net radiation absorbed by the canopy (incoming – reflected – emitted)
#           Main driver of leaf heating and photosynthesis.          [W m⁻²]
#
# lEctot  = Latent heat flux from canopy transpiration
#           Energy used to evaporate water from leaves.              [W m⁻²]
#
# Hctot   = Sensible heat flux from canopy
#           Heat transferred from warm leaves to the air.            [W m⁻²]
#
# Actot   = Net CO₂ assimilation by the canopy
#           Canopy-level photosynthetic carbon uptake.         [µmol m⁻² s⁻¹]
#
# Tcave   = Average canopy temperature
#           Mean leaf temperature across all canopy layers.          [°C]
#
# ----------------------------- SOIL FLUXES ----------------------------------
# Rnstot  = Net radiation absorbed by the soil                           [W m⁻²]
#
# lEstot  = Latent heat flux from soil evaporation                       [W m⁻²]
#
# Hstot   = Sensible heat flux from the soil                             [W m⁻²]
#
# Gtot    = Ground heat flux (energy conducted into deeper soil)         [W m⁻²]
#
# Tsave   = Soil surface temperature                                     [°C]
#
# ----------------------------- TOTAL FLUXES ---------------------------------
# Rntot   = Total net radiation (canopy + soil)                           [W m⁻²]
#
# lEtot   = Total latent heat flux (transpiration + soil evaporation)     [W m⁻²]
#
# Htot    = Total sensible heat flux (canopy + soil)                      [W m⁻²]
# ----------------------------------------------------------------------------



# ----------------------------------------------------------------------------
# 4.1 Load fluxes and vegetation outputs
# ----------------------------------------------------------------------------

fluxes <- read.csv(file.path(path_run, "fluxes.csv"))
veg    <- read.csv(file.path(path_run, "vegatation.csv"))  # nombre con typo del script

# Remove n.sim column if present
if ("n.sim" %in% names(fluxes)) fluxes <- fluxes[ , !(names(fluxes) %in% "n.sim")]
if ("n.sim" %in% names(veg))    veg    <- veg[ , !(names(veg)    %in% "n.sim")]

# ----------------------------------------------------------------------------
# 4.2 Boxplots for flux variables (distribución entre simulaciones)
# ----------------------------------------------------------------------------

flux_long <- fluxes |>
  pivot_longer(cols = everything(),
               names_to = "Flux",
               values_to = "Value")

ggplot(flux_long, aes(x = Flux, y = Value)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Value") +
  ggtitle("Distribution of SCOPE energy and carbon fluxes")


# ----------------------------------------------------------------------------
# 5. Physiological and thermal vegetation variables (vegatation.csv)
# ----------------------------------------------------------------------------
# This file contains canopy-level physiological indicators simulated by SCOPE.
# Each row = one simulation; each column = a canopy property related to
# photosynthesis, energy dissipation, or thermal state.
#
# A
#   = Canopy-averaged photosynthesis rate
#   Net CO₂ assimilation of the entire canopy.           [µmol m⁻² s⁻¹]
#
# Ja
#   = Electron transport rate
#   Represents the rate of energy flow through the photosynthetic
#     electron transport chain (linked to light reactions). [µmol m⁻² s⁻¹]
#
# ENPQ
#   = Non-photochemical quenching energy flux
#   Heat dissipated by the plant as a photoprotective mechanism
#     when excess light cannot be used for photosynthesis.   [W m⁻²]
#
# PNPQ
#   = Non-photochemical quenching per photon
#   Efficiency of photoprotective heat dissipation relative
#     to absorbed photons.                                   [µmol photons⁻¹]
#
# LST
#   = Land Surface Temperature (canopy)
#   Effective temperature of the canopy top, influenced by
#     transpiration, radiation, and stomatal behavior.       [K]
#
# emis
#   = Thermal emissivity of the canopy
#   Capacity of the canopy to emit thermal radiation;
#     values close to 1 indicate nearly perfect emitters.    [unitless]
# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# 5.1 Boxplots for vegetation variables (A, Ja, ENPQ, PNPQ, LST, emis)
# ----------------------------------------------------------------------------

veg_long <- veg |>
  pivot_longer(cols = everything(),
               names_to = "Variable",
               values_to = "Value")

ggplot(veg_long, aes(x = Variable, y = Value)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Value") +
  ggtitle("Distribution of SCOPE vegetation outputs")




# ----------------------------------------------------------------------------
# 6. Relating LUT traits (Cab, LAI, Vcmax25, etc.) with SCOPE fluxes and vegetation outputs
# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# 6.1. Load LUT, fluxes, and vegetation outputs
# ----------------------------------------------------------------------------

lut_file     <- file.path(path_run, "Parameters", "inputLUT.csv")
fluxes_file  <- file.path(path_run, "fluxes.csv")
veg_file     <- file.path(path_run, "vegatation.csv")   # (name saved by SCOPE)

lut    <- read.csv(lut_file)
fluxes <- read.csv(fluxes_file)
veg    <- read.csv(veg_file)

# Remove n.sim column if present
if ("n.sim" %in% names(fluxes)) fluxes <- fluxes[, !(names(fluxes) %in% "n.sim")]
if ("n.sim" %in% names(veg))    veg    <- veg[,    !(names(veg)    %in% "n.sim")]

# ----------------------------------------------------------------------------
# 6.2. Create a combined dataframe (one row = one simulation)
# ----------------------------------------------------------------------------

n_sim <- nrow(lut)

df_all <- lut |>
  mutate(sim = 1:n_sim) |>
  bind_cols(
    fluxes |> mutate(sim = 1:n_sim),
    veg    |> mutate(sim = 1:n_sim)
  )

# df_all now contains:
# - Inputs: Cab, Car, Anth, EWT, LAI, Vcmax25, LIDFa, LIDFb, etc.
# - Fluxes: Actot, Rntot, lEtot, Htot, etc.
# - Vegetation: A, Ja, ENPQ, PNPQ, LST, emis


# ----------------------------------------------------------------------------
# 6.3. Load LUT, fluxes, and vegetation outputs
# ----------------------------------------------------------------------------

lut_file     <- file.path(path_run, "Parameters", "inputLUT.csv")
fluxes_file  <- file.path(path_run, "fluxes.csv")
veg_file     <- file.path(path_run, "vegatation.csv")   # (name saved by SCOPE)

lut    <- read.csv(lut_file)
fluxes <- read.csv(fluxes_file)
veg    <- read.csv(veg_file)

# Remove n.sim column if present
if ("n.sim" %in% names(fluxes)) fluxes <- fluxes[, !(names(fluxes) %in% "n.sim")]
if ("n.sim" %in% names(veg))    veg    <- veg[,    !(names(veg)    %in% "n.sim")]

# ----------------------------------------------------------------------------
# 6.4. Create a combined dataframe (one row = one simulation)
# ----------------------------------------------------------------------------

n_sim <- nrow(lut)

df_all <- lut |>
  mutate(sim = 1:n_sim) |>
  bind_cols(
    fluxes |> mutate(sim = 1:n_sim),
    veg    |> mutate(sim = 1:n_sim)
  )

# df_all now contains:
# - Inputs: Cab, Car, Anth, EWT, LAI, Vcmax25, LIDFa, LIDFb, etc.
# - Fluxes: Actot, Rntot, lEtot, Htot, etc.
# - Vegetation: A, Ja, ENPQ, PNPQ, LST, emis

# ----------------------------------------------------------------------------
# 6.5. Scatterplot example: LAI vs Actot (CO₂ assimilation)
# ----------------------------------------------------------------------------
# LAI    = Leaf Area Index (m² leaf / m² ground)
# Actot  = Total canopy net CO₂ assimilation (µmol m⁻² s⁻¹)
# ----------------------------------------------------------------------------

ggplot(df_all, aes(x = LAI, y = Actot)) +
  geom_point(alpha = 0.6) +
  theme_bw() +
  xlab("LAI (m² m⁻²)") +
  ylab("Actot (µmol m⁻² s⁻¹)") +
  ggtitle("Relationship between LAI and canopy CO₂ assimilation (Actot)")

# ----------------------------------------------------------------------------
# 6.6. Vcmax25 vs A (photosynthetic capacity vs canopy photosynthesis)
# ----------------------------------------------------------------------------

# ----------------------------------------------------------------------------
# Vcmax25 = Maximum Rubisco carboxylation rate at 25°C (µmol m⁻² s⁻¹)
# A       = Average canopy photosynthesis (µmol m⁻² s⁻¹)
# Color indicates LAI to show structural effects.
# ----------------------------------------------------------------------------

ggplot(df_all, aes(x = Vcmax25, y = A, color = LAI)) +
  geom_point(alpha = 0.7) +
  scale_color_viridis_c(name = "LAI") +
  theme_bw() +
  xlab("Vcmax25 (µmol m⁻² s⁻¹)") +
  ylab("A (µmol m⁻² s⁻¹)") +
  ggtitle("Vcmax25 vs canopy photosynthesis (A)\ncolored by LAI")

# ----------------------------------------------------------------------------
# Vcmax25 = Maximum Rubisco carboxylation capacity at 25°C (µmol m⁻² s⁻¹)
#           → This is an INPUT trait controlling biochemical photosynthesis.
#
# Ja       = Electron transport rate (µmol m⁻² s⁻¹)
#           → This is an OUTPUT from SCOPE, reflecting light-driven electron flow.
#
# Plot goal:
#   Examine how biochemical capacity (Vcmax25) influences electron transport (Ja),
#   while color-coding by LAI to show how canopy structure modifies the response.
# ----------------------------------------------------------------------------

ggplot(df_all, aes(x = Vcmax25, y = Ja, color = LAI)) +
  geom_point(alpha = 0.7) +
  scale_color_viridis_c(name = "LAI") +
  theme_bw() +
  xlab("Vcmax25 (µmol m⁻² s⁻¹)\nBiochemical capacity for CO₂ fixation") +
  ylab("Ja (µmol m⁻² s⁻¹)\nElectron transport rate") +
  ggtitle("Relationship between Vcmax25 and Electron Transport (Ja)\nwith canopy structure (LAI) as color gradient")



# ----------------------------------------------------------------------------
# 6.7. LAI vs LST (how canopy structure affects surface temperature)
# ----------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# LST = Land Surface Temperature (K)
# ----------------------------------------------------------------------------

ggplot(df_all, aes(x = LAI, y = LST)) +
  geom_point(alpha = 0.6) +
  theme_bw() +
  xlab("LAI (m² m⁻²)") +
  ylab("LST (K)") +
  ggtitle("Relationship between LAI and canopy temperature (LST)")


 # ----------------------------------------------------------------------------
 # 6.8. Cab vs reflectance at 550 nm
 # ----------------------------------------------------------------------------
 # Biological relevance:
 #   • Chlorophyll strongly absorbs red and blue light but reflects green
 #   • Reflectance at 550 nm increases when Cab decreases
 # ----------------------------------------------------------------------------

 refl <- read.csv(file.path(path_run, "refl.csv"))
 refl_550 <- refl[, which(abs(SCOPEinR::define.bands()$wlS - 550) < 0.3)]

 df_all$refl_550 <- refl_550

 ggplot(df_all, aes(x = Cab, y = refl_550)) +
   geom_point(alpha = 0.6) +
   theme_bw() +
   xlab("Cab (µg cm⁻²)\nChlorophyll content") +
   ylab("Reflectance at 550 nm") +
   ggtitle("Chlorophyll decreases green reflectance (550 nm)")

 # ----------------------------------------------------------------------------
 # 6.9. Cab vs F740 /761
 # ----------------------------------------------------------------------------
 # Biological relevance:
 #   • F740 peak increases with chlorophyll concentration
 #   • A key variable for SIF-based photosynthesis estimation
 # ----------------------------------------------------------------------------

 sif <- read.csv(file.path(path_run, "fluorescence_scalar.csv"))
 df_all$F740 <- sif$F740
 df_all$F761 <- sif$F761

 ggplot(df_all, aes(x = Cab, y = F740)) +
   geom_point(alpha = 0.6, color = "purple") +
   theme_bw() +
   xlab("Cab (µg cm⁻²)") +
   ylab("F740 (W m⁻² µm⁻¹ sr⁻¹)") +
   ggtitle("Red/Far-Red fluorescence (F740) increases with chlorophyll")


 # ----------------------------------------------------------------------------
 # Plot Cab vs F761
 # ----------------------------------------------------------------------------
 ggplot(df_all, aes(x = Cab, y = F761)) +
   geom_point(alpha = 0.6, color = "red") +
   theme_bw() +
   xlab("Cab (µg cm⁻²)\nLeaf chlorophyll content") +
   ylab("F761 (W m⁻² µm⁻¹ sr⁻¹)\nFluorescence emission at 761 nm") +
   ggtitle("Increase in F761 fluorescence with increasing chlorophyll")


 # ----------------------------------------------------------------------------
 # 7. Absorbed PAR (aPAR.csv)
 # ----------------------------------------------------------------------------

 # ----------------------------------------------------------------------------
 # SCOPE separates absorbed PAR into components associated with different
 # canopy elements. All values are in µmol m⁻² s⁻¹.
 #
 # aPARTot    = Total absorbed PAR by the canopy + soil
 #              → Primary driver of photosynthesis and fluorescence.
 #
 # aPARsun    = Absorbed PAR by sunlit leaves
 #              → Highest light load; strong driver of electron transport (Ja).
 #
 # aPARshaded = Absorbed PAR by shaded leaves
 #              → Lower irradiance; controls understory/photosynthesis stability.
 #
 # aPARsoil   = Absorbed PAR by the soil surface
 #              → Important when LAI is low (open canopies).
 # ----------------------------------------------------------------------------

 apar <- read.csv(file.path(path_run, "aPAR.csv"))

 df_all$aPARTot    <- apar$aPARTot
 df_all$aPARsun    <- apar$aPARsun
 df_all$aPARshaded <- apar$aPARshaded
 df_all$aPARsoil   <- apar$aPARsoil


 # ----------------------------------------------------------------------------
 # 7.1. Relationship 1: Total aPAR vs Photosynthesis (A)
 # ----------------------------------------------------------------------------
 # Biological meaning:
 #   • Photosynthesis responds directly to absorbed PAR, not incoming sunlight.
 #   • At low aPAR: relationship is linear (light-limited).
 #   • At high aPAR: saturates due to biochemical constraints (Vcmax, Jmax).
 # ----------------------------------------------------------------------------

 ggplot(df_all, aes(x = aPARTot, y = A)) +
   geom_point(alpha = 0.6, color = "blue") +
   theme_bw() +
   xlab("aPARTot (µmol m⁻² s⁻¹)\nAbsorbed PAR") +
   ylab("A (µmol m⁻² s⁻¹)\nNet CO₂ assimilation") +
   ggtitle("Light-response relationship: Photosynthesis (A) vs Absorbed PAR")


 # ----------------------------------------------------------------------------
 # 7.2. Relationship 2: Total aPAR vs Fluorescence (F740, F761)
 # ----------------------------------------------------------------------------
 # Why this matters:
 #   • Fluorescence is proportional to excitation energy (absorbed photons).
 #   • F740 and F761 are the main far-red SIF peaks used by FLEX/TROPOMI.
 #   • Higher aPAR → more excitation → more fluorescence emission.
 # ----------------------------------------------------------------------------

 ggplot(df_all, aes(x = aPARTot, y = F740)) +
   geom_point(alpha = 0.6, color = "purple") +
   theme_bw() +
   xlab("aPARTot (µmol m⁻² s⁻¹)") +
   ylab("F740 (W m⁻² µm⁻¹ sr⁻¹)") +
   ggtitle("Fluorescence peak F740 increases with absorbed PAR")

 ggplot(df_all, aes(x = aPARTot, y = F761)) +
   geom_point(alpha = 0.6, color = "red") +
   theme_bw() +
   xlab("aPARTot (µmol m⁻² s⁻¹)") +
   ylab("F761 (W m⁻² µm⁻¹ sr⁻¹)") +
   ggtitle("Fluorescence peak F761 increases with absorbed PAR")


 # ----------------------------------------------------------------------------
 # 7.3. Relationship 3: LAI vs aPARTot
 # ----------------------------------------------------------------------------
 # Interpretation:
 #   • LAI determines how much light is intercepted by the canopy.
 #   • aPARTot increases with LAI until light saturation/self-shading occurs.
 # ----------------------------------------------------------------------------

 ggplot(df_all, aes(x = LAI, y = aPARTot)) +
   geom_point(alpha = 0.6, color = "darkgreen") +
   theme_bw() +
   xlab("LAI (m² m⁻²)") +
   ylab("aPARTot (µmol m⁻² s⁻¹)") +
   ggtitle("Absorbed PAR increases with canopy density (LAI)")


 # ----------------------------------------------------------------------------
 # 7.4. Relationship 4: Net Radiation (Rnctot) vs aPARTot
 # ----------------------------------------------------------------------------
 # Why this matters:
 #   • Net radiation provides the energy budget at the canopy surface.
 #   • A strong correlation is expected: more net radiation → more absorbed PAR.
 #   • Non-linearities appear depending on leaf angle distribution and LAI.
 # ----------------------------------------------------------------------------

 ggplot(df_all, aes(x = Rnctot, y = aPARTot)) +
   geom_point(alpha = 0.6, color = "orange") +
   theme_bw() +
   xlab("Rnctot (W m⁻²)\nNet radiation absorbed by canopy") +
   ylab("aPARTot (µmol m⁻² s⁻¹)") +
   ggtitle("Link between canopy net radiation and absorbed PAR")



 # ----------------------------------------------------------------------------
 # 8. Aerodynamic and stomatal resistances (resistance.csv)
 # ----------------------------------------------------------------------------

 # ----------------------------------------------------------------------------
 # This file contains resistances that control how heat and water are exchanged
 # between the canopy, soil and the atmosphere.
 #
 # Each row = one SCOPE simulation (matching n.sim in other files).
 #
 # VARIABLES
 # ----------------------------------------------------------------------------
 # raws      = Aerodynamic resistance above the canopy                  [s m⁻¹]
 #             → Describes how efficiently heat and water vapour are
 #               transported from the canopy to the reference height
 #               (affected by wind speed and canopy roughness).
 #
 # raforSoil = Aerodynamic resistance for the soil surface               [s m⁻¹]
 #             → Controls how easily heat and vapour are transported
 #               from the soil surface to the air above.
 #
 # rss       = Soil surface resistance                                   [s m⁻¹]
 #             → Represents how difficult it is for water to evaporate
 #               from the soil (strongly linked to soil moisture).
 #             → High rss = dry soil (low evaporation),
 #               Low rss = wet soil (high evaporation).
 #
 # ustar     = Friction velocity                                         [m s⁻¹]
 #             → A turbulence metric; higher u* means stronger mixing
 #               between the surface (canopy/soil) and the atmosphere.
 # ----------------------------------------------------------------------------

 res <- read.csv(file.path(path_run, "resistance.csv"))

 # Merge resistance variables into df_all (assumes same simulation order)
 df_all$raws      <- res$raws
 df_all$raforSoil <- res$raforSoil
 df_all$rss       <- res$rss
 df_all$ustar     <- res$ustar


 # ----------------------------------------------------------------------------
 # 8.1. Plot 1: Aerodynamic resistance above canopy (raws) vs sensible heat (Hctot)
 # ----------------------------------------------------------------------------
 # Interpretation:
 #   • Lower raws → more efficient turbulent transfer → higher Hctot.
 #   • Higher raws → canopy more "decoupled" from the air → lower Hctot.
 # ----------------------------------------------------------------------------

 ggplot(df_all, aes(x = raws, y = Hctot)) +
   geom_point(alpha = 0.6, color = "orange") +
   theme_bw() +
   xlab("Aerodynamic resistance above canopy (raws, s m⁻¹)") +
   ylab("Canopy sensible heat flux Hctot (W m⁻²)") +
   ggtitle("Lower aerodynamic resistance (raws) increases canopy sensible heat flux")


 # ----------------------------------------------------------------------------
 # 8.2. Plot 2: Soil surface resistance (rss) vs soil latent heat (lEstot)
 # ----------------------------------------------------------------------------
 # Interpretation:
 #   • High rss = dry soil → evaporation is suppressed → low lEstot.
 #   • Low rss = wet soil → evaporation allowed → high lEstot.
 # ----------------------------------------------------------------------------

 ggplot(df_all, aes(x = rss, y = lEstot)) +
   geom_point(alpha = 0.6, color = "dodgerblue4") +
   theme_bw() +
   xlab("Soil surface resistance rss (s m⁻¹)") +
   ylab("Soil latent heat flux lEstot (W m⁻²)") +
   ggtitle("Dry soil (high rss) reduces soil evaporation (lEstot)")


 # ----------------------------------------------------------------------------
 # 8.3. Plot 3: Aerodynamic resistance for soil (raforSoil) vs soil sensible heat (Hstot)
 # ----------------------------------------------------------------------------
 # Interpretation:
 #   • Low raforSoil → efficient transfer of heat from soil to air → larger Hstot.
 #   • High raforSoil → soil more isolated → smaller sensible heat flux.
 # ----------------------------------------------------------------------------

 ggplot(df_all, aes(x = raforSoil, y = Hstot)) +
   geom_point(alpha = 0.6, color = "brown4") +
   theme_bw() +
   xlab("Aerodynamic resistance for soil (raforSoil, s m⁻¹)") +
   ylab("Soil sensible heat flux Hstot (W m⁻²)") +
   ggtitle("Soil sensible heat flux decreases with increasing aerodynamic resistance")


 # ----------------------------------------------------------------------------
 # 8.4. Plot 4: Friction velocity (ustar) vs total sensible heat (Htot)
 # ----------------------------------------------------------------------------
 # Interpretation:
 #   • ustar measures turbulence intensity.
 #   • Higher u* usually enhances turbulent exchange → can increase Htot.
 # ----------------------------------------------------------------------------

 ggplot(df_all, aes(x = ustar, y = Htot)) +
   geom_point(alpha = 0.6, color = "darkred") +
   theme_bw() +
   xlab("Friction velocity u* (m s⁻¹)") +
   ylab("Total sensible heat flux Htot (W m⁻²)") +
   ggtitle("Stronger turbulence (higher u*) tends to enhance sensible heat exchange")


 # ----------------------------------------------------------------------------
 # 9. General-purpose function to create any trait–output scatterplot
 # ----------------------------------------------------------------------------

 # ----------------------------------------------------------------------------
 # Usage examples:
 #   plot_trait_vs_output(df_all, "Cab", "Actot")
 #   plot_trait_vs_output(df_all, "Anth", "A", color_by = "LAI")
 #   plot_trait_vs_output(df_all, "Vcmax25", "LST", color_by = "LAI")
 # ----------------------------------------------------------------------------

 plot_trait_vs_output <- function(df, xvar, yvar, color_by = NULL) {

   if (!xvar %in% names(df)) stop(paste("Variable", xvar, "is not in the dataframe"))
   if (!yvar %in% names(df)) stop(paste("Variable", yvar, "is not in the dataframe"))

   p <- ggplot(df, aes_string(x = xvar, y = yvar))

   if (!is.null(color_by)) {
     if (!color_by %in% names(df)) stop(paste("Variable", color_by, "is not in the dataframe"))
     p <- p + aes_string(color = color_by) +
       scale_color_viridis_c(name = color_by)
   }

   p + geom_point(alpha = 0.7) +
     theme_bw() +
     ggtitle(paste(xvar, "vs", yvar)) +
     xlab(xvar) +
     ylab(yvar)
 }

# Example calls:
plot_trait_vs_output(df_all, "Cab", "Actot")
plot_trait_vs_output(df_all, "Anth", "A", color_by = "LAI")
plot_trait_vs_output(df_all, "Vcmax25", "LST", color_by = "LAI")#
plot_trait_vs_output(df_all, "Vcmax25", "LST", color_by = "LAI")#
plot_trait_vs_output(df_all, "Rnctot", "aPARTot", color_by = "LAI")#









