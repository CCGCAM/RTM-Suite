
# PROSAIL / PROSAIL-WithSatellite / Model Explorer online reflectance simulator -- GRS32306 Advanced Earth Observation
#
# Five tabs:
#  - "PROSAIL": native-resolution (1nm) canopy reflectance simulator (as before).
#  - "PROSAIL-WithSatellite": the same simulation, PLUS a satellite-band overlay -- pick a sensor
#    (Sentinel-2A/B, Landsat 5/7/8, Terra/Aqua MODIS, Sentinel-3A/B, PRISMA), see the native
#    curve and the sensor-band-convolved values (points + connecting line) together, with each
#    band's FWHM shown as a shaded grey polygon.
#  - "Model Explorer" (formerly "RTM"): a general-purpose explorer -- pick ANY canopy model
#    (fourSAIL, fourSAIL2, INFORM) x ANY leaf model (PROSPECT-D, PROSPECT-PRO, Fluspect-B,
#    Fluspect-B-Cx, Liberty) and the input panel switches to exactly the parameters that
#    combination needs. Three sub-tabs: Simulation (the canopy reflectance plot), SIF
#    (leaf-level fluorescence, Fluspect models only), Sensitivity (sweep one trait at a time).
#  - "How in R" / "How in Python": a static tutorial mirroring this app's own workflow in code.

# 1. load libraries (fast path: try require() first, only install what's actually
#    missing -- avoids the slow installed.packages() full-library scan on every
#    startup that the original version always paid, even when nothing needed installing) ----

required_packages <- c("shiny", "shinythemes", "ggplot2", "reshape2", "patchwork")
for (pkg in required_packages) {
  if (!suppressWarnings(requireNamespace(pkg, quietly = TRUE))) {
    install.packages(pkg)
  }
}
suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(ggplot2)
  library(reshape2)
  library(patchwork)
  library(ToolsRTM)
})

# 2. Satellite sensor helpers ----------------------------------------------------------------

# All sensors get.coef.SMAC()-style objects support, plus PRISMA (handled
# separately below, via ToolsRTM::srf.prisma/fwhm.prisma -- PRISMA isn't part
# of the 9-sensor SMAC atmospheric-correction catalogue, it's a plain SRF table).
.sensor_choices <- c(
  "Sentinel-2A" = "Sentinel2A.MSI", "Sentinel-2B" = "Sentinel2B.MSI",
  "Landsat 5 TM" = "LANDSAT5.TM", "Landsat 7 ETM+" = "LANDSAT7.ETM", "Landsat 8 OLI" = "LANDSAT8.OLI",
  "Terra/Aqua MODIS" = "TerraAqua.MODIS",
  "Sentinel-3A OLCI" = "Sentinel3A.OLCI", "Sentinel-3B OLCI" = "Sentinel3B.OLCI",
  "PRISMA (hyperspectral)" = "PRISMA"
)

# For the "RTM" tab: which canopy model calls which ToolsRTM function, and
# which leaf models each one accepts (all 3 canopy models support the same
# 5 leaf models -- see toolsrtm.canopy._leaf_optics-equivalent dispatch in
# foursail()/foursail2()/inform(), all via a `LeafModel` argument).
.canopy_model_choices <- c("fourSAIL" = "foursail", "fourSAIL2" = "foursail2", "INFORM" = "inform")
.leaf_model_choices <- c(
  "PROSPECT-D" = "PROSPECT-D", "PROSPECT-PRO" = "PROSPECT-PRO",
  "Fluspect-B" = "Fluspect-B", "Fluspect-B-Cx" = "Fluspect-B-Cx", "Liberty" = "Liberty"
)

#' Leaf-level simulated SIF (sun-induced fluorescence) for the Fluspect-B /
#' Fluspect-B-Cx leaf models, plus the "apparent reflectance" a sensor would
#' see once that fluorescence is added to the reflected signal.
#'
#' getFluspect.B()/getFluspect.Cx() return excitation-emission matrices
#' (MbI/MbII, backward/PSI+PSII -- the upward, sensor-facing direction for a
#' leaf illuminated and viewed from above), not a single emission spectrum:
#' rows are emission wavelengths (wlF, 640-850nm), columns are excitation
#' wavelengths (wlE, 400-750nm). Weighting each excitation column by the
#' incident irradiance at that wavelength and summing across excitation
#' wavelengths (a matrix-vector product) gives the emitted SIF spectrum
#' under natural illumination -- the standard leaf-level fluorescence
#' integral (the same first step SCOPE/RTMf's own canopy-coupled version
#' builds on, before scaling through canopy structure/LAI, which this app's
#' canopy models (fourSAIL/fourSAIL2/INFORM) don't include -- so this is
#' LEAF-level SIF, not the full canopy-transported quantity RTMf would give).
#'
#' "Apparent reflectance" follows directly from at-sensor radiance = leaf
#' reflected radiance + emitted SIF radiance, divided by incident
#' irradiance at the (emission) wavelength: apparent = true_refl + SIF/E.
#'
#' Verified against known chlorophyll fluorescence physiology: this produces
#' the expected double-peaked emission (red peak ~685nm from PSII, far-red
#' peak ~740nm from PSI+PSII combined) and a small (few-percent) apparent-
#' reflectance bump concentrated in the red-edge/near-infrared, matching the
#' textbook SIF emission spectrum shape.
compute_leaf_sif <- function(leaf_row, leaf_model) {
  fluspect_fn <- if (leaf_model == "Fluspect-B") ToolsRTM::getFluspect.B else ToolsRTM::getFluspect.Cx
  version <- if (leaf_model == "Fluspect-B") "D" else "Cx"
  out <- fluspect_fn(inputsLeaf = leaf_row, inputsOptipar = ToolsRTM::optipar, version = version)

  wlE <- seq(400, 750, 1)
  wlF <- seq(640, 850, 1)
  d <- ToolsRTM::dataSpec_PDB
  incident <- d$direct_light + d$diffuse_light  # W m-2 nm-1, 400-2500nm grid

  E_wle <- incident[match(wlE, d$wavelength)]
  # Fluspect-B splits backward fluorescence by photosystem (MbI/MbII);
  # Fluspect-B-Cx returns a single combined matrix (Mb) instead -- same
  # (wlF x wlE) shape either way, just under a different field name.
  Mb_total <- if (leaf_model == "Fluspect-B") out$MbI + out$MbII else out$Mb
  sif <- as.numeric(Mb_total %*% E_wle)

  refl_wlf <- out$refl[match(wlF, out$lambda)]
  E_wlf <- incident[match(wlF, d$wavelength)]
  apparent_refl <- refl_wlf + sif / E_wlf

  data.frame(wavelength = wlF, reflectance = refl_wlf, apparent_reflectance = apparent_refl, sif = sif)
}

#' Dispatch a fully-built LUT row to the chosen canopy model and return its
#' reflectance spectrum. Shared by the RTM tab's main plot (rtm_result) and
#' its sensitivity sweep (rtm_sensitivity_result), so both use the exact
#' same, already-tested dispatch logic instead of two copies drifting apart.
run_canopy_model <- function(row, canopy_model, leaf_model, rsoil) {
  if (canopy_model == "foursail") {
    sail <- ToolsRTM::foursail(inputLUT = row, rsoil = rsoil, LeafModel = leaf_model)
    wavelength <- ToolsRTM::dataSpec_PDB[, 1]
    reflectance <- ToolsRTM::Compute_BRF(rdot = sail$rdot, rsot = sail$rsot, tts = row$tts, data.light = ToolsRTM::dataSpec_PDB)
  } else if (canopy_model == "foursail2") {
    sail <- ToolsRTM::foursail2(inputLUT = row, rsoil = rsoil, LeafModel = leaf_model)
    wavelength <- ToolsRTM::dataSpec_PDB[, 1]
    reflectance <- ToolsRTM::Compute_BRF(rdot = sail$rdot, rsot = sail$rsot, tts = row$tts, data.light = ToolsRTM::dataSpec_PDB)
  } else {
    out <- ToolsRTM::inform(inputLUT = row, rsoil = rsoil, LeafModel = leaf_model)
    # inform() returns TOC reflectance directly (already the BRF-equivalent
    # quantity), on the same 400-2500nm grid as fourSAIL/fourSAIL2 -- no
    # separate Compute_BRF() step needed.
    wavelength <- ToolsRTM::dataSpec_PDB[, 1]
    reflectance <- as.numeric(out)
  }
  data.frame(wavelength = wavelength[seq_along(reflectance)], reflectance = reflectance)
}

# Which traits can be swept in the RTM tab's Sensitivity sub-tab, per leaf
# model -- mirrors rtm_leaf_ui()'s own slider definitions (same input IDs,
# min/max) so the sweep range always matches what's actually adjustable.
# "LAI" (canopy-level) is available for every leaf model.
.rtm_sweep_traits <- list(
  "PROSPECT-D" = list(
    LAI = c(0.5, 8), N = c(1, 3), Cab = c(0, 80), Car = c(0, 20),
    Anth = c(0, 7), Cbrown = c(0, 1), EWT = c(0.0001, 0.05), LMA = c(0.0001, 0.03)
  ),
  "PROSPECT-PRO" = list(
    LAI = c(0.5, 8), N = c(1, 3), Cab = c(0, 80), Car = c(0, 20),
    Anth = c(0, 7), Cbrown = c(0, 1), EWT = c(0.0001, 0.05),
    Prot = c(0.0001, 0.03), CBC = c(0.0001, 0.03)
  ),
  "Fluspect-B" = list(
    LAI = c(0.5, 8), N = c(1, 3), Cab = c(0, 80), Car = c(0, 20), Anth = c(0, 7),
    EWT = c(0.0001, 0.05), LMA = c(0.0001, 0.03), Cs = c(0, 1), fqe = c(0, 0.05)
  ),
  "Fluspect-B-Cx" = list(
    LAI = c(0.5, 8), N = c(1, 3), Cab = c(0, 80), Car = c(0, 20), Anth = c(0, 7),
    EWT = c(0.0001, 0.05), LMA = c(0.0001, 0.03), Cs = c(0, 1), fqe = c(0, 0.05), Cx = c(0, 1)
  ),
  # cell.d/inter.c narrowed from the full slider range (20-100 / 0.01-0.1) --
  # liberty()'s internal fixed-point solver for R (Compute_BRF-equivalent
  # inter-layer reflectance) can fail to converge near those extremes,
  # producing NaN reflectance. Confirmed via direct sweep (Python port, but
  # the same iterative solver is a direct port of ToolsRTM::liberty, so the
  # same instability applies to R). baseline.abs dropped entirely (not just
  # narrowed) -- its own safe range is so small (0.0004-0.001) that a
  # one-trait-at-a-time sweep's color legend ends up all but unreadable
  # (labels like "4e-05"), and it stays fixed at its one_sim_code_r/python()
  # default instead (0.0006, itself safely inside the narrow stable range).
  # Even the remaining narrower bounds don't fully eliminate the instability
  # under COMBINED random sampling (about 4% of rows in a 500-sample LUT,
  # many_sims_code_r/python()'s own use of this table) since several
  # parameters can still compound near their own edges -- hence the
  # is.finite()/np.isfinite() filter in that generated code, needed only for
  # Liberty in practice but applied generically.
  "Liberty" = list(
    LAI = c(0.5, 8), `cell.d` = c(30, 100), `inter.c` = c(0.01, 0.09),
    `leaf.thick` = c(1, 10), `albino.abs` = c(0, 4),
    Cab = c(0, 80), EWT = c(0.0001, 0.05), `lign.cell` = c(1, 8), Nitrogen = c(0.3, 2)
  )
)

# Canopy-level traits sweepable in the Sensitivity sub-tab, regardless of
# leaf model: geometry/hotspot (every canopy model) plus each canopy
# model's own extra structural parameters (empty for plain fourSAIL).
# Merged into .rtm_sweep_traits[[leaf_model]] at lookup time -- see
# rtm_sweep_trait_choices() below, the single place both the Sensitivity
# sub-tab's own dropdown and its sweep logic read from.
# psi (relative azimuth) deliberately excluded -- it's mathematically
# degenerate whenever tto == 0 (nadir view, this tab's own default): with no
# "direction" to a straight-down view, azimuth cannot affect reflectance at
# all (confirmed via direct repro: rdot/rsot bit-for-bit identical across
# the full 0-180 psi range at tto=0; off-nadir, e.g. tto=30, it has a real,
# substantial effect, especially near the hotspot at psi=0). Since the
# sweep always starts from whatever tto currently is, a psi sweep looks
# like nothing happened by default, which is confusing rather than useful --
# dropped rather than kept with a caveat.
.rtm_geometry_traits <- list(
  hspot = c(0, 0.5), tts = c(0, 80), tto = c(0, 80), psoil = c(0, 1),
  # Special case (see rtm_sensitivity_result()'s own handling): the "Type
  # LIDF" dropdown picks one of 6 preset (LIDFa, LIDFb) shape-pairs
  # (TypeLidf=1 mode); sweeping LIDFa directly instead means switching to
  # the alternative ellipsoidal/Campbell parameterization (TypeLidf=2),
  # where LIDFa IS the mean leaf inclination angle in degrees (0=planophile/
  # horizontal, 90=erectophile/vertical) and LIDFb is unused (fixed at 0)
  # -- confirmed numerically stable across the full 0-90 range.
  LIDFa = c(0, 90)
)
.rtm_canopy_extra_traits <- list(
  "foursail" = list(),
  "foursail2" = list(fraction_brown = c(0, 1), diss = c(0, 1), Cv = c(0, 1), Zeta = c(0.1, 5)),
  # LAIu starts at 0.1, not 0 -- foursail.inform()'s hotspot integral divides
  # by (ks * lai) (sumint <- (1 - tss) / (ks * lai)); at lai == 0 exactly
  # both numerator and denominator are 0, giving R's 0/0 = NaN instead of
  # the correct removable-singularity limit (confirmed via direct repro:
  # LAIu=0 -> NaN reflectance, any LAIu>0 -> fine). Same reason the main LAI
  # trait above starts at 0.5, not 0.
  "inform" = list(LAIu = c(0.1, 3), cd = c(1, 15), h = c(2, 40), sd = c(50, 2000), skyl = c(0, 1))
)

#' All traits sweepable in the Sensitivity sub-tab for a given canopy x leaf
#' combination -- that leaf model's own params (.rtm_sweep_traits),
#' PLUS geometry (.rtm_geometry_traits), PLUS that canopy model's own extra
#' structural params (.rtm_canopy_extra_traits). Same registry backs both
#' the sub-tab's trait dropdown and the sweep computation itself, so they
#' can never disagree on what's available or what range to sweep it over.
rtm_sweep_trait_choices <- function(canopy_model, leaf_model) {
  c(.rtm_sweep_traits[[leaf_model]], .rtm_geometry_traits, .rtm_canopy_extra_traits[[canopy_model]])
}

#' Same universe of traits as rtm_sweep_trait_choices(), grouped into
#' <optgroup> categories for the Sensitivity dropdown -- LAI/hspot/LIDFa/
#' canopy-extras all live in .rtm_sweep_traits/.rtm_geometry_traits mixed in
#' with leaf params (many_sims_code_r/python() need that flat structure for
#' LUT sampling, so it isn't reorganized at the source), this just
#' re-buckets the same names by what kind of thing they are for display.
rtm_sweep_trait_groups <- function(canopy_model, leaf_model) {
  leaf_names <- setdiff(names(.rtm_sweep_traits[[leaf_model]]), "LAI")
  canopy_names <- c("LAI", "hspot", "LIDFa", names(.rtm_canopy_extra_traits[[canopy_model]]))
  list(
    "Leaf biochemistry" = leaf_names,
    "Canopy structure" = canopy_names,
    "Viewing & illumination angles" = c("tts", "tto"),
    "Soil" = "psoil"
  )
}

#' Generate R/Python code reproducing the CURRENT sensitivity sweep exactly
#' -- LIVE slider values (not the tutorial's fixed defaults), so a student
#' can copy/paste and immediately reproduce what they're looking at, then
#' change one number and re-run to explore further themselves. `row` is the
#' live rtm_row() (or an equivalent one-row data.frame); `psoil` is the
#' current soil-brightness slider value (not a row column, see rsoil_for()
#' in rtm_sensitivity_result() for why psoil needs separate handling).
sensitivity_code_r <- function(canopy, leaf, trait, row, psoil, rng) {
  fn <- .canopy_fn_name[[canopy]]
  # LIDFa/LIDFb/TypeLidf are always fixed (never a sweep target) but aren't
  # in .leaf_model_params/.canopy_model_params/geometry at all -- included
  # unconditionally here (from the CURRENT TypeLIDF dropdown choice, not a
  # hardcoded default like the tutorial's simpler one_sim_code_r() uses) so
  # every leaf-angle-distribution-dispatching call (foursail()/foursail2()/
  # inform() all read these) actually has them.
  # LIDFa sweep forces the ellipsoidal (Campbell) leaf-angle
  # parameterization -- TypeLidf<-2, LIDFb<-0 -- for every swept row (must
  # match rtm_sensitivity_result()'s row_for() exactly), so LIDFb/TypeLidf
  # are excluded from the fixed block here rather than emitted at their
  # (irrelevant, possibly TypeLidf==1-preset) current slider values.
  excl <- if (trait == "LIDFa") c(trait, "LIDFb", "TypeLidf") else trait
  fixed_names <- setdiff(c("LAI", "hspot", "tts", "tto", "psi", "LIDFa", "LIDFb", "TypeLidf",
                            names(.leaf_model_params[[leaf]]),
                            names(.canopy_model_params[[canopy]])), excl)
  fixed_vals <- vapply(fixed_names, function(n) unname(row[[n]]), numeric(1))
  fixed_args <- .fmt_params_r(setNames(fixed_vals, fixed_names))

  rsoil_lines <- if (trait == "psoil") {
    "# psoil itself is swept -- rsoil is rebuilt per value below, not fixed"
  } else {
    sprintf("rsoil <- %s * dataSpec_PDB[, 11] + %s * dataSpec_PDB[, 12]  # current soil brightness", psoil, 1 - psoil)
  }
  loop_row_line <- if (trait == "psoil") {
    "  rsoil_v <- v * dataSpec_PDB[, 11] + (1 - v) * dataSpec_PDB[, 12]"
  } else if (trait == "LIDFa") {
    "  row_v <- row; row_v$TypeLidf <- 2; row_v$LIDFb <- 0; row_v$LIDFa <- v  # ellipsoidal (Campbell) mode"
  } else {
    sprintf("  row_v <- row; row_v$%s <- v", trait)
  }
  # tts always read from row_v (has every column, including the swept one if
  # trait == "tts"), never the outer fixed `row` (missing that one column
  # whenever it's the trait being swept) -- for the psoil branch, `row_v`
  # isn't defined (only rsoil_v varies), so `row` itself is used there
  # instead, correctly, since psoil sweeps never remove `tts` from `row`.
  tts_ref <- if (trait == "psoil") "row$tts" else "row_v$tts"
  call_line <- if (canopy == "inform") {
    sprintf('  %s(inputLUT = %s, rsoil = %s, LeafModel = "%s")',
            fn, if (trait == "psoil") "row" else "row_v", if (trait == "psoil") "rsoil_v" else "rsoil", leaf)
  } else {
    sprintf(
      '  sail_v <- %s(inputLUT = %s, rsoil = %s, LeafModel = "%s")\n  Compute_BRF(rdot = sail_v$rdot, rsot = sail_v$rsot,\n              tts = %s, data.light = dataSpec_PDB)',
      fn, if (trait == "psoil") "row" else "row_v", if (trait == "psoil") "rsoil_v" else "rsoil", leaf, tts_ref
    )
  }

  paste0(
    "row <- data.frame(\n  ", fixed_args, "\n)\n",
    rsoil_lines, "\n\n",
    sprintf("values <- seq(%s, %s, length.out = 6)  # sweeping %s\n", rng[1], rng[2], trait),
    "spectra <- t(sapply(values, function(v) {\n",
    loop_row_line, "\n",
    call_line, "\n",
    "}))"
  )
}

sensitivity_code_python <- function(canopy, leaf, trait, row, psoil, rng) {
  fn <- .canopy_fn_name[[canopy]]
  excl <- if (trait == "LIDFa") c(trait, "LIDFb", "TypeLidf") else trait
  fixed_names <- setdiff(c("LAI", "hspot", "tts", "tto", "psi", "LIDFa", "LIDFb", "TypeLidf",
                            names(.leaf_model_params[[leaf]]),
                            names(.canopy_model_params[[canopy]])), excl)
  fixed_vals <- vapply(fixed_names, function(n) unname(row[[n]]), numeric(1))
  fixed_args <- .fmt_params_py(setNames(fixed_vals, fixed_names))

  import_line <- sprintf("from toolsrtm import %s, compute_brf", fn)
  rsoil_lines <- if (trait == "psoil") {
    "# psoil itself is swept -- rsoil is rebuilt per value below, not fixed"
  } else {
    sprintf("rsoil = %s * dry_soil + %s * wet_soil  # current soil brightness (see toolsrtm._data.data_spec_pdb)",
            psoil, 1 - psoil)
  }
  loop_body <- if (trait == "psoil") {
    "    rsoil_v = v * dry_soil + (1 - v) * wet_soil"
  } else if (trait == "LIDFa") {
    '    row_v = dict(row); row_v["TypeLidf"] = 2; row_v["LIDFb"] = 0; row_v["LIDFa"] = v  # ellipsoidal (Campbell) mode'
  } else {
    sprintf('    row_v = dict(row); row_v["%s"] = v', trait)
  }
  tts_ref <- if (trait == "psoil") 'row["tts"]' else 'row_v["tts"]'
  call_line <- if (canopy == "inform") {
    sprintf('    spectra[i] = np.asarray(%s(%s, %s, leaf_model="%s"))',
            fn, if (trait == "psoil") "row" else "row_v", if (trait == "psoil") "rsoil_v" else "rsoil", leaf)
  } else {
    sprintf(
      '    sail_v = %s(%s, %s, leaf_model="%s")\n    spectra[i] = compute_brf(sail_v.rdot, sail_v.rsot, %s)',
      fn, if (trait == "psoil") "row" else "row_v", if (trait == "psoil") "rsoil_v" else "rsoil", leaf, tts_ref
    )
  }

  paste0(
    "import numpy as np\n", import_line, "\n\n",
    "row = {\n    ", fixed_args, ",\n}\n",
    rsoil_lines, "\n\n",
    sprintf("values = np.linspace(%s, %s, 6)  # sweeping %s\n", rng[1], rng[2], trait),
    "spectra = np.empty((6, ", if (leaf %in% c("Fluspect-B", "Fluspect-B-Cx")) "2001" else "2101", "))\n",
    "for i, v in enumerate(values):\n",
    loop_body, "\n",
    call_line
  )
}

#' Convolve a native-resolution reflectance spectrum onto a sensor's bands.
#'
#' Uses each sensor's own spectral response function (SRF) table
#' (`wl_srf_smac`/`p_srf_smac`, the same data `ToolsRTM::get.spectral.convolution`
#' uses) rather than the coarser band-center interpolation `SPART()` itself
#' uses for reflectance -- this gives a real SRF-weighted band value, and lets
#' us compute an empirical FWHM per band directly from the SRF weight profile
#' (`band_width` is NaN for several bundled sensors, e.g. Landsat 8 and MODIS,
#' so it can't be used directly).
#'
#' **A real indexing bug found while building this**: for any sensor whose
#' `id_smac_in_all` is not the identity permutation (MODIS, Landsat, OLCI --
#' true for Sentinel-2 only by coincidence, since its mapping happens to be
#' `1:13`), `wl_srf_smac`/`p_srf_smac`'s COLUMN order follows the sensor's own
#' natural band numbering (parsed from `band_id_smac`, e.g. "Band 17  (NIR)"),
#' NOT `wl_smac`'s array position. Naively pairing column `b` with
#' `wl_smac[b]` (as `ToolsRTM::get.spectral.convolution` does for its `wave`
#' column, and as an early draft of this function also did) puts each band's
#' SRF weights against the WRONG center wavelength for those sensors --
#' confirmed empirically: MODIS "band 8" pairs with `wl_smac[8]` = 416nm, but
#' its actual SRF weight profile peaks at ~905nm (MODIS's real band 17).
#' Corrected here by re-deriving each column's true center wavelength via
#' `id_smac_in_all` + the natural band numbers parsed from `band_id_smac`.
convolve_smac_sensor <- function(wl_native, refl_native, sensor_obj) {
  wl_srf <- sensor_obj$wl_srf_smac
  p_srf <- sensor_obj$p_srf_smac
  nbands <- ncol(wl_srf)

  id_map <- as.numeric(sensor_obj$id_smac_in_all)
  band_labels <- sensor_obj$band_id_smac
  has_labels <- length(band_labels) == nbands && all(grepl("Band[[:space:]]*[0-9]+", band_labels))
  natural_band_no <- if (has_labels) {
    as.numeric(regmatches(band_labels, regexpr("[0-9]+", band_labels)))
  } else {
    # No descriptive labels (Sentinel-2): id_smac_in_all is a direct 1:nbands
    # identity permutation into wl_smac already.
    seq_len(nbands)
  }
  true_center <- as.numeric(sensor_obj$wl_smac)[match(id_map, natural_band_no)]

  out <- data.frame(band = seq_len(nbands), wl = true_center, fwhm = NA_real_, reflectance = NA_real_)
  for (b in seq_len(nbands)) {
    wl_b <- wl_srf[, b]; p_b <- p_srf[, b]
    valid <- !is.na(wl_b) & !is.na(p_b) & p_b > 0
    if (!any(valid)) next
    wl_v <- wl_b[valid]; p_v <- p_b[valid]

    half <- max(p_v) / 2
    above <- wl_v[p_v >= half]
    out$fwhm[b] <- if (length(above) >= 2) max(above) - min(above) else NA

    idx <- match(round(wl_v), round(wl_native))
    ok <- !is.na(idx)
    if (sum(ok) == 0) next
    out$reflectance[b] <- sum(p_v[ok] * refl_native[idx[ok]]) / sum(p_v[ok])
  }
  out <- out[!is.na(out$wl), ]
  out[order(out$wl), ]
}

#' Convolve onto PRISMA/Sentinel-2A/Sentinel-2B via `ToolsRTM::get.spectral.convolution.srf()`
#' -- the package's own plain-per-band-SRF-table convolution function,
#' shared by exactly these 3 sensors (PRISMA has no SMAC atmospheric-
#' correction coefficients at all; Sentinel-2A/B additionally ship a
#' second, plain, publisher-original SRF table alongside their SMAC bundle
#' -- see `srf.prisma`/`srf.sentinel2a`/`srf.sentinel2b`'s own docs). This
#' used to be a private, app-only helper (`convolve_prisma()`) reimplementing
#' this exact SRF-weighting logic; promoted into ToolsRTM itself so the app
#' depends on one well-tested, documented package function instead of its
#' own copy (PRISMA's convolved reflectance is byte-identical to the old
#' in-app version, cross-checked before switching over).
convolve_srf_sensor <- function(wl_native, refl_native, srf_table, fwhm_table = NULL) {
  df <- data.frame(wave = wl_native, rfl = refl_native)
  out <- ToolsRTM::get.spectral.convolution.srf(df, srf_table, fwhm = fwhm_table)
  data.frame(band = seq_along(out$wl), wl = out$wl, fwhm = out$fwhm, reflectance = out$RFL)
}

#' The 6 remaining bundled sensors (Landsat 4/5/7/8, Sentinel-3A/B OLCI,
#' Terra/Aqua MODIS) only ship SMAC-bundled SRF data, no plain per-band
#' table -- still routed through this app's own convolve_smac_sensor()
#' rather than `ToolsRTM::get.spectral.convolution.rfl()` (also fixed this
#' round, see its own docstring) because that package function doesn't
#' compute FWHM, and Sentinel-3/MODIS's SRF matrix needs the UNCOLLAPSED,
#' per-detector-resolved profile for FWHM (get.coef.SMAC() averages that
#' away via colMeans() for exactly those 3 sensors) -- a real data-shape
#' complication left for future work rather than rushed here.
convolve_to_sensor <- function(wl_native, refl_native, sensor_key) {
  if (sensor_key == "PRISMA") {
    return(convolve_srf_sensor(wl_native, refl_native, ToolsRTM::srf.prisma, ToolsRTM::fwhm.prisma))
  }
  if (sensor_key == "Sentinel2A.MSI") {
    return(convolve_srf_sensor(wl_native, refl_native, ToolsRTM::srf.sentinel2a))
  }
  if (sensor_key == "Sentinel2B.MSI") {
    return(convolve_srf_sensor(wl_native, refl_native, ToolsRTM::srf.sentinel2b))
  }
  sensor_obj <- getExportedValue("ToolsRTM", sensor_key)
  convolve_smac_sensor(wl_native, refl_native, sensor_obj)
}

# Reference tables (Reference tab) -----------------------------------------------------------

#' Every leaf/canopy model parameter this app exposes, in one place, cross-
#' checked against the actual slider definitions in leaf_canopy_sidebar()/
#' rtm_leaf_ui()/rtm_canopy_extra_ui() so the reference table can never
#' silently drift from what the sliders themselves allow.
.rtm_reference_leaf_params <- rbind(
  data.frame(Model = "PROSPECT-D", Parameter = c("N", "Cab", "Car", "Anth", "Cbrown", "EWT", "LMA"),
             Description = c("Leaf structure parameter (mesophyll layers)", "Chlorophyll a+b content",
                              "Carotenoid content", "Anthocyanin content", "Brown pigment content",
                              "Equivalent water thickness", "Leaf mass per area (dry matter)"),
             Range = c("1 - 3", "0 - 80", "0 - 20", "0 - 7", "0 - 1", "0.0001 - 0.05", "0.0001 - 0.03"),
             Units = c("-", "ug/cm2", "ug/cm2", "ug/cm2", "-", "g/cm2", "g/cm2")),
  data.frame(Model = "PROSPECT-PRO", Parameter = c("N", "Cab", "Car", "Anth", "Cbrown", "EWT", "Prot *", "CBC *"),
             Description = c("Leaf structure parameter (mesophyll layers)", "Chlorophyll a+b content",
                              "Carotenoid content", "Anthocyanin content", "Brown pigment content",
                              "Equivalent water thickness",
                              "Leaf protein content -- see LMA = Prot + CBC table below.",
                              "Carbon-based constituents (cellulose+lignin+starch) -- see LMA = Prot + CBC table below."),
             Range = c("1 - 3", "0 - 80", "0 - 20", "0 - 7", "0 - 1", "0.0001 - 0.05", "0.0001 - 0.03", "0.0001 - 0.03"),
             Units = c("-", "ug/cm2", "ug/cm2", "ug/cm2", "-", "g/cm2", "g/cm2", "g/cm2")),
  data.frame(Model = "Fluspect-B", Parameter = c("N", "Cab", "Car", "Anth", "EWT", "LMA", "Cs", "fqe"),
             Description = c("Leaf structure parameter (mesophyll layers)", "Chlorophyll a+b content",
                              "Carotenoid content", "Anthocyanin content", "Equivalent water thickness",
                              "Leaf mass per area (dry matter)", "Senescent/brown pigment content",
                              "Fluorescence quantum efficiency"),
             Range = c("1 - 3", "0 - 80", "0 - 20", "0 - 7", "0.0001 - 0.05", "0.0001 - 0.03", "0 - 1", "0 - 0.05"),
             Units = c("-", "ug/cm2", "ug/cm2", "ug/cm2", "g/cm2", "g/cm2", "-", "-")),
  data.frame(Model = "Fluspect-B-Cx", Parameter = c("N", "Cab", "Car", "Anth", "EWT", "LMA", "Cs", "fqe", "Cx"),
             Description = c("Leaf structure parameter (mesophyll layers)", "Chlorophyll a+b content",
                              "Carotenoid content", "Anthocyanin content", "Equivalent water thickness",
                              "Leaf mass per area (dry matter)", "Senescent/brown pigment content",
                              "Fluorescence quantum efficiency", "Violaxanthin-zeaxanthin transition state"),
             Range = c("1 - 3", "0 - 80", "0 - 20", "0 - 7", "0.0001 - 0.05", "0.0001 - 0.03", "0 - 1", "0 - 0.05", "0 - 1"),
             Units = c("-", "ug/cm2", "ug/cm2", "ug/cm2", "g/cm2", "g/cm2", "-", "-", "-")),
  data.frame(Model = "Liberty", Parameter = c("cell.d", "inter.c", "baseline.abs", "leaf.thick",
                                               "albino.abs", "Cab", "EWT", "lign.cell", "Nitrogen"),
             Description = c("Average leaf cell diameter", "Intercellular air space fraction",
                              "Wavelength-independent baseline absorption",
                              "Leaf thickness (Benford layer-stacking parameter)",
                              "Visible-region absorption due to lignin (albino-leaf proxy)",
                              "Chlorophyll content", "Equivalent water thickness",
                              "Combined lignin+cellulose content", "Nitrogen content"),
             Range = c("30 - 100 (20 - 100 slider)", "0.01 - 0.09 (0.01 - 0.1 slider)", "0.0001 - 0.001",
                       "1 - 10", "0 - 4", "0 - 80", "0.0001 - 0.05", "1 - 8", "0.3 - 2"),
             Units = c("um", "-", "-", "-", "-", "ug/cm2", "g/cm2", "ug/cm2", "g/m2"))
)

#' * = PROSPECT-PRO's Prot/CBC vs LMA relationship, spelled out separately
#' instead of crammed into the leaf-params table's Description cell (same
#' "put the long explanation in its own small table, point to it with an
#' asterisk" pattern as .rtm_reference_lidf_presets below).
.rtm_reference_lma_prot_cbc <- data.frame(
  Quantity = c("LMA (PROSPECT-D)", "Prot + CBC (PROSPECT-PRO)"),
  Meaning = c("Leaf mass per area (total dry matter) as ONE lumped number.",
              "The SAME dry matter split into its two chemical components: Prot (protein) and CBC (cellulose+lignin+starch, i.e. carbon-based constituents)."),
  Relationship = c("LMA = Prot + CBC", "Prot + CBC = LMA"),
  Note = c("Used by PROSPECT-D and Fluspect (no Prot/CBC split available).",
           "This app always drives PROSPECT-PRO with Prot/CBC, not LMA -- the model itself falls back to LMA only when Prot=CBC=0, which never happens here.")
)

.rtm_reference_canopy_params <- rbind(
  data.frame(Model = "All (fourSAIL/fourSAIL2/INFORM)",
             Parameter = c("LAI", "hspot", "TypeLidf", "LIDFa *", "LIDFb *", "psoil", "tts", "tto", "psi"),
             Description = c(
               "Leaf area index", "Hotspot parameter (leaf-size / canopy-height ratio)",
               "Which leaf-angle-distribution parameterization LIDFa/LIDFb are read as: 1 = one of 6 preset canonical shapes, 2 = ellipsoidal (Campbell) mode. See the leaf-angle-distribution table below.",
               "Leaf-angle-distribution parameter #1 -- meaning depends on TypeLidf. See table below.",
               "Leaf-angle-distribution parameter #2 -- meaning depends on TypeLidf. See table below.",
               "Soil brightness (0 = wet/dark soil, 1 = dry/bright soil)",
               "Solar zenith angle", "Observer (view) zenith angle", "Relative azimuth angle"),
             Range = c("0.001 - 10", "0 - 1", "1 (preset) or 2 (ellipsoidal)", "see table below *", "see table below *",
                       "0 - 1", "0 - 90", "0 - 90", "0 - 180"),
             Units = c("m2/m2", "-", "-", "-", "-", "-", "deg", "deg", "deg")),
  data.frame(Model = "fourSAIL2 (extra)", Parameter = c("fraction_brown", "diss", "Cv", "Zeta"),
             Description = c("Fraction of brown (senescent) leaf area", "Green/brown layer dissociation factor",
                              "Vertical crown cover fraction", "Tree shape factor (crown diameter / height)"),
             Range = c("0 - 1", "0 - 1", "0 - 1", "0.1 - 5"),
             Units = c("-", "-", "-", "-")),
  data.frame(Model = "INFORM (extra)", Parameter = c("LAIu", "cd", "h", "sd", "skyl"),
             Description = c("Understorey leaf area index", "Tree crown diameter", "Tree height",
                              "Stem density", "Diffuse-radiation fraction"),
             Range = c("0.1 - 3 (0 - 3 slider)", "1 - 15", "2 - 40", "50 - 2000", "0 - 1"),
             Units = c("m2/m2", "m", "m", "ha-1", "-"))
)

#' * = LIDFa/LIDFb's meaning depends on TypeLidf -- this table spells out
#' both modes so the "*" note in .rtm_reference_canopy_params has somewhere
#' concrete to point to instead of just re-explaining TypeLidf inline. The
#' 6 TypeLidf=1 preset codes are the fixed (LIDFa, LIDFb) pairs the app's
#' own TypeLIDF dropdown maps each named shape to (see rtm_row(), the
#' `switch(input$rtm_TypeLIDF, ...)` calls) -- they are NOT physical
#' angles, just canonical shape-selector codes from the SAIL leaf-angle
#' distribution literature (Verhoef 1998 / Bunnik 1978).
.rtm_reference_lidf_presets <- rbind(
  data.frame(TypeLidf = "1", Name = "Planophile", LIDFa = "1", LIDFb = "0",
             Description = "Mostly horizontal leaves (e.g. prostrate/rosette canopies)"),
  data.frame(TypeLidf = "1", Name = "Erectophile", LIDFa = "-1", LIDFb = "0",
             Description = "Mostly vertical/erect leaves (e.g. many grasses, cereals)"),
  data.frame(TypeLidf = "1", Name = "Plagiophile", LIDFa = "0", LIDFb = "-1",
             Description = "Mostly oblique leaves, peaked around ~45 deg"),
  data.frame(TypeLidf = "1", Name = "Extremophile", LIDFa = "0", LIDFb = "1",
             Description = "Bimodal mix of mostly-horizontal + mostly-vertical leaves, few oblique"),
  data.frame(TypeLidf = "1", Name = "Uniform", LIDFa = "0", LIDFb = "0",
             Description = "All inclination angles equally likely"),
  data.frame(TypeLidf = "1", Name = "Spherical", LIDFa = "-0.35", LIDFb = "-0.15",
             Description = "Leaf normals distributed as if uniformly random on a sphere -- no dominant leaf angle"),
  data.frame(TypeLidf = "2", Name = "Ellipsoidal / Campbell", LIDFa = "0 - 90 (mean leaf angle, deg)", LIDFb = "0 (unused)",
             Description = "LIDFa is used directly as the canopy's mean leaf inclination angle in degrees (0 = horizontal, 90 = vertical), continuously variable rather than one of 6 fixed shapes. This app's Sensitivity tab's LIDFa sweep always uses this mode, since it's the one where LIDFa varies continuously with a direct physical meaning.")
)

#' Actual leaf-angle-distribution SHAPES for the table above -- not just
#' codes in a table, the real frequency-vs-inclination-angle curve each
#' preset produces, computed via the SAME dladgen()/campbell() functions
#' foursail()/foursail2()/inform() call internally (never re-derived here).
#' litab is identical (bin midpoints 5,15,...,89 deg) between dladgen() and
#' campbell(), confirmed by reading both functions' source, so the two
#' halves stack into one long data frame safely.
.rtm_lidf_preset_codes <- list(
  Planophile = c(1, 0), Erectophile = c(-1, 0), Plagiophile = c(0, -1),
  Extremophile = c(0, 1), Uniform = c(0, 0), Spherical = c(-0.35, -0.15)
)
.rtm_lidf_distribution <- rbind(
  do.call(rbind, lapply(names(.rtm_lidf_preset_codes), function(nm) {
    ab <- .rtm_lidf_preset_codes[[nm]]
    ld <- ToolsRTM::dladgen(ab[1], ab[2])
    data.frame(Angle = ld$litab, Frequency = ld$lidf, Shape = nm)
  })),
  do.call(rbind, lapply(c(20, 45, 70), function(ala) {
    ld <- ToolsRTM::campbell(ala)
    data.frame(Angle = ld$litab, Frequency = ld$lidf, Shape = sprintf("Ellipsoidal (LIDFa=%d deg)", ala))
  }))
)
.rtm_lidf_distribution$Shape <- factor(.rtm_lidf_distribution$Shape,
  levels = c(names(.rtm_lidf_preset_codes), "Ellipsoidal (LIDFa=20 deg)", "Ellipsoidal (LIDFa=45 deg)", "Ellipsoidal (LIDFa=70 deg)"))
# dladgen()/campbell()'s 13 angle bins are NOT evenly spaced -- 10 deg wide
# up to 80 deg, then 2 deg wide near grazing incidence (81/83/85/87/89) --
# a fixed geom_col(width=) would visually overlap those narrow tail bins
# (and trips ggplot2's position_stack() overlap warning); this per-row
# width instead matches each bin's real extent.
.rtm_lidf_distribution$BinWidth <- rep(c(rep(10, 8), rep(2, 5)), 9)

#' Band center wavelength + FWHM for every band of every bundled sensor --
#' reuses convolve_to_sensor() itself (the SAME, already-tested SRF-weighted
#' logic the PROSAIL-WithSatellite tab's own plot uses), called once at
#' app-load time with a flat dummy reflectance spectrum (only the wl/fwhm
#' columns are used here, not the resulting -- meaningless for a flat input
#' -- convolved reflectance).
.rtm_reference_sensor_bands <- do.call(rbind, lapply(names(.sensor_choices), function(sensor_name) {
  sensor_key <- .sensor_choices[[sensor_name]]
  wl_native <- ToolsRTM::dataSpec_PDB[, 1]
  dummy_refl <- rep(0.2, length(wl_native))
  out <- convolve_to_sensor(wl_native, dummy_refl, sensor_key)
  data.frame(
    SensorKey = sensor_key,
    Sensor = sensor_name,
    Band = if ("band" %in% colnames(out)) paste0("B", out$band) else seq_len(nrow(out)),
    `Center wavelength (nm)` = round(out$wl, 1),
    `FWHM (nm)` = if ("fwhm" %in% colnames(out)) round(out$fwhm, 1) else NA,
    check.names = FALSE
  )
}))

# Shared plot theme, sized for large/high-DPI screens -- plotOutput containers
# are already vh-based (grow with the browser window), but ggplot2's default
# point-based text stays a fixed pixel size regardless of canvas size, so on a
# big screen the (bigger) plot area makes the (same-size) axis text/labels
# look proportionally tiny. base_size=20 (vs. the original 13-14) plus larger
# explicit axis-text/legend sizes fixes that; render_plot_res (used as
# renderPlot's own `res=` argument) sharpens text on high-DPI displays.
render_plot_res <- 110
theme_prosail <- function(base_size = 20, legend_position = "top") {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size + 4),
      plot.subtitle = element_text(hjust = 0.5, size = base_size - 3, face = "italic"),
      axis.title = element_text(face = "bold", size = base_size + 2),
      axis.text = element_text(size = base_size - 1),
      legend.position = legend_position,
      legend.text = element_text(face = "bold", size = base_size - 1),
      legend.title = element_text(size = base_size),
      legend.key.width = unit(1.6, "cm"),
      panel.grid.minor = element_blank()
    )
}

# "How in R" / "How in Python" tutorial tabs -- static walkthrough (code +
# pre-generated figures, run once via generate_tutorial_figures.R and saved
# to www/) of the same 6 steps for BOTH languages: load the package, run one
# simulation, run a LUT of many simulations, sweep one trait at a time
# (sensitivity), convolve to a specific sensor, and invert a trait with ML.
# Worked example uses fourSAIL + PROSPECT-D (the RTM tab's default) -- see
# each tab's closing note for what's different for PROSAIL/PROSAIL-WithSatellite.
# Monotonically-increasing counter so every code_block() call gets its own
# DOM id -- needed because renderUI() re-invokes code_block() reactively
# (a new id each time is fine, always exactly one <pre> with that id exists
# at once) and because multiple code blocks appear on the same page at once
# (Sensitivity's R/Python tabs, each tutorial step, etc. -- colliding ids
# would make the copy button copy the WRONG block's text).
.code_block_id <- local({
  i <- 0
  function() { i <<- i + 1; paste0("codeblock_", i) }
})

code_block <- function(code) {
  id <- .code_block_id()
  tags$div(
    style = "position: relative;",
    tags$button(
      "Copy", id = paste0(id, "_btn"), class = "btn btn-default btn-xs",
      onclick = sprintf(
        "navigator.clipboard.writeText(document.getElementById('%s').innerText); var b=document.getElementById('%s_btn'); b.innerText='Copied!'; setTimeout(function(){b.innerText='Copy';}, 1500);",
        id, id
      ),
      style = "position:absolute; top:8px; right:8px; z-index:2;"
    ),
    tags$pre(
      id = id,
      style = paste(
        "background:#f7f7f9; border:1px solid #ddd; border-radius:6px; margin:0;",
        "padding:12px 16px; padding-top:36px; font-size:15px; line-height:1.45; overflow-x:auto;"
      ),
      code
    )
  )
}

# Per-model "one simulation" code, generated (not hand-typed x15) from the
# SAME parameter sets the RTM tab's own rtm_row()/rtm_result() reactives use
# -- one leaf-model parameter block x one canopy-model call, in both
# languages, so picking a combination in the dropdowns below shows exactly
# the code that combination's own RTM-tab simulation runs. Values match
# rtm_row()'s per-leaf-model defaults (app.R's own reactive) exactly.
.leaf_model_params <- list(
  "PROSPECT-D"    = c(N = 1.5, Cab = 40, Car = 8, Anth = 2, Cbrown = 0, EWT = 0.009, LMA = 0.009, alpha = 40),
  "PROSPECT-PRO"  = c(N = 1.5, Cab = 40, Car = 8, Anth = 2, Cbrown = 0, EWT = 0.009, LMA = 0, alpha = 40, Prot = 0.0045, CBC = 0.005),
  # Cx is read unconditionally by getFluspect.B() even for plain Fluspect-B
  # (version='D'), so it needs a harmless default here too; alpha is only
  # needed when this row also feeds inform()'s always-run understory step
  # (foursail.inform()), not by Fluspect itself, but is included on every
  # leaf block so the SAME row works for all 3 canopy models.
  "Fluspect-B"    = c(N = 1.8, Cab = 40, Car = 10, Anth = 0, EWT = 0.015, LMA = 0.01, Cs = 0.1, fqe = 0.01, Cx = 0, alpha = 40),
  "Fluspect-B-Cx" = c(N = 1.8, Cab = 40, Car = 10, Anth = 0, EWT = 0.015, LMA = 0.01, Cs = 0.1, fqe = 0.01, Cx = 0.5,
                       Prot = 0.0045, CBC = 0.005, alpha = 40),
  "Liberty"       = c(cell.d = 40, inter.c = 0.045, baseline.abs = 0.0006, leaf.thick = 1.6,
                       albino.abs = 2, Cab = 40, EWT = 0.009, lign.cell = 4, Nitrogen = 1, alpha = 40)
)
.canopy_model_params <- list(
  "foursail"  = c(),
  "foursail2" = c(fraction_brown = 0.2, diss = 0.5, Cv = 1, Zeta = 1),
  "inform"    = c(LAIu = 0.1, cd = 4.5, h = 20, sd = 650, skyl = 0.1)
)
.canopy_fn_name <- c("foursail" = "foursail", "foursail2" = "foursail2", "inform" = "inform")

#' Specific absorption coefficient curves for the Absorption coefficients
#' sub-tab -- genuinely different data per leaf model family, not the same
#' table relabeled: PROSPECT-D uses dataSpec_PDB's SC_* columns,
#' PROSPECT-PRO uses dataSpec_PRO's spA_* columns (a different
#' parameterization -- protein/carbon-based split instead of one combined
#' dry-matter term), Fluspect-B/-B-Cx use ToolsRTM::optipar's K* columns
#' (shared table; -Cx additionally exposes Kp/Kcbc, matching that it's the
#' one Fluspect variant that also accepts Prot/CBC inputs). Liberty has no
#' equivalent at all -- its absorption model is parameterized by cell
#' geometry/intercellular air space, not a per-unit-content specific
#' absorption spectrum, so it returns NULL (handled explicitly by the
#' server, not silently plotting nonsense).
#' `cx`: the CURRENT Cx slider value (only meaningful for Fluspect models,
#' ignored otherwise) -- used to compute the live-interpolated Kca curve for
#' panel3 (below), rather than showing optipar's own static bundled `Kca`
#' column mislabeled as "current".
.abs_coef_source <- function(leaf_model, cx = 0) {
  if (leaf_model == "PROSPECT-D") {
    d <- ToolsRTM::dataSpec_PDB
    list(wave = d$wavelength,
         panel1 = list(`Chlorophyll content` = d$SC_chl, `Carotenoids content` = d$SC_car,
                        `Anthocyanin content` = d$SC_Anth, `Brown pigments` = d$SC_Brwon),
         panel2 = list(`Water content` = d$SC_Cw, `Dry matter` = d$SC_Cm))
  } else if (leaf_model == "PROSPECT-PRO") {
    d <- ToolsRTM::dataSpec_PRO
    list(wave = d$wave,
         panel1 = list(`Chlorophyll content` = d$spA_Cab, `Carotenoids content` = d$spA_Car,
                        `Anthocyanin content` = d$spA_Atn, `Brown pigments` = d$spA_Br),
         # spA_Cm (dry matter, LMA-parameterized) included alongside
         # Protein/nonPro -- PROSPECT-PRO's own Kall formula uses EITHER
         # LMA*Kdm OR Prot*Kp + CBC*Kcbc (foursail.R zeroes out whichever
         # wasn't supplied), so both pathways' coefficients are genuinely
         # part of this model, not redundant.
         panel2 = list(`Water content` = d$spA_Cw, `Dry matter (LMA)` = d$spA_Cm,
                        `Protein content` = d$spA_Protein, `Carbon-based compounds` = d$spA_nonPro))
  } else if (leaf_model %in% c("Fluspect-B", "Fluspect-B-Cx")) {
    d <- ToolsRTM::optipar
    panel2 <- list(`Water content` = d$Kw, `Dry matter` = d$Kdm)
    if (leaf_model == "Fluspect-B-Cx") {
      panel2[["Protein content"]] <- d$Kp
      panel2[["Carbon-based compounds"]] <- d$Kcbc
    }
    # Both Fluspect variants (not just -Cx) actually interpolate carotenoid
    # absorption between a violaxanthin state (KcaV, Cx=0) and zeaxanthin
    # state (KcaZ, Cx=1) -- getFluspect.B itself does this
    # (Kca <- (1-Cx)*KcaV + Cx*KcaZ) unless the legacy Cx=-999 sentinel is
    # used, which this app never does. Shown as its own panel since it's
    # genuinely new information (the *mechanism* behind the single
    # "Carotenoids content" line in panel 1), not available at all for
    # PROSPECT (no analogous xanthophyll-cycle state in that model).
    current_kca <- (1 - cx) * d$KcaV + cx * d$KcaZ
    list(wave = d$wl,
         panel1 = list(`Chlorophyll content` = d$Kab, `Carotenoids content` = d$Kca,
                        `Anthocyanin content` = d$Kant, `Senescence (Cs)` = d$Ks),
         panel2 = panel2,
         panel3 = list(`Violaxanthin state (Cx=0)` = d$KcaV, `Zeaxanthin state (Cx=1)` = d$KcaZ,
                        `Current Kca (interpolated)` = current_kca))
  } else {
    NULL
  }
}

#' Format a named numeric vector as R `data.frame(...)` / Python `dict(...)`
#' constructor arguments, one per line, matching this tutorial's own style.
.fmt_params_r <- function(params) {
  nm <- names(params)
  nm_quoted <- ifelse(grepl("\\.", nm), paste0("`", nm, "`"), nm)
  paste(paste0(nm_quoted, " = ", params), collapse = ", ")
}
.fmt_params_py <- function(params) {
  # Dict-LITERAL syntax ("key": value), not dict(key=value) kwargs -- Liberty's
  # parameter names (e.g. "cell.d") contain a dot, which Python's kwarg syntax
  # can't express at all (dict(cell.d=40) is a syntax error: it parses as
  # attribute assignment, not a keyword argument) but a plain string dict key
  # handles for every leaf model uniformly.
  nm <- names(params)
  paste(paste0('"', nm, '": ', params), collapse = ", ")
}

one_sim_code_r <- function(canopy, leaf) {
  leaf_p <- .leaf_model_params[[leaf]]
  canopy_p <- .canopy_model_params[[canopy]]
  fn <- .canopy_fn_name[[canopy]]
  row_args <- paste(c(
    "LAI = 3, hspot = 0.01, LIDFa = -0.35, LIDFb = -0.15, TypeLidf = 1",
    "tts = 30, tto = 0, psi = 0",
    .fmt_params_r(leaf_p),
    if (length(canopy_p) > 0) .fmt_params_r(canopy_p)
  ), collapse = ",\n  ")
  call_line <- sprintf('sail <- %s(inputLUT = row, rsoil = rsoil, LeafModel = "%s")', fn, leaf)
  brf_lines <- if (canopy == "inform") {
    "# inform() returns TOC reflectance directly -- no separate Compute_BRF() step\nreflectance <- as.numeric(sail)"
  } else {
    'reflectance <- Compute_BRF(rdot = sail$rdot, rsot = sail$rsot,\n                            tts = row$tts, data.light = dataSpec_PDB)'
  }
  paste0(
    "row <- data.frame(\n  ", row_args, "\n)\n",
    "rsoil <- 0.5 * dataSpec_PDB[, 11] + 0.5 * dataSpec_PDB[, 12]  # dry/wet soil blend\n\n",
    call_line, "\n", brf_lines
  )
}

one_sim_code_python <- function(canopy, leaf) {
  leaf_p <- .leaf_model_params[[leaf]]
  canopy_p <- .canopy_model_params[[canopy]]
  fn <- .canopy_fn_name[[canopy]]
  row_args <- paste(c(
    '"LAI": 3, "hspot": 0.01, "LIDFa": -0.35, "LIDFb": -0.15, "TypeLidf": 1',
    '"tts": 30, "tto": 0, "psi": 0',
    .fmt_params_py(leaf_p),
    if (length(canopy_p) > 0) .fmt_params_py(canopy_p)
  ), collapse = ",\n    ")
  import_line <- if (canopy == "foursail") {
    "from toolsrtm import foursail, compute_brf"
  } else {
    sprintf("from toolsrtm import %s, compute_brf", fn)
  }
  call_line <- sprintf('sail = %s(row, rsoil, leaf_model="%s")', fn, leaf)
  brf_lines <- if (canopy == "inform") {
    "# inform() returns TOC reflectance directly -- no separate compute_brf() step\nreflectance = np.asarray(sail)"
  } else {
    'reflectance = compute_brf(sail.rdot, sail.rsot, row["tts"])'
  }
  paste0(
    "import numpy as np\n", import_line, "\n\n",
    # dict-literal syntax ({"key": value}), not dict(key=value) kwargs --
    # Liberty's parameter names (e.g. "cell.d") aren't valid Python
    # identifiers, so they can't be written as kwargs at all.
    "row = {\n    ", row_args, ",\n}\n",
    "wavelength = np.arange(400, 2501)\n",
    "rsoil = np.full(2101, 0.15)  # flat soil (see toolsrtm.marmit/BSM for a modeled soil)\n\n",
    call_line, "\n", brf_lines
  )
}

#' Per-model "n_samples = 500" code (step 3) -- same idea as one_sim_code_r/
#' _python, but sampling a LUT instead of one fixed row. Samples every trait
#' this leaf model actually supports over the SAME ranges the RTM tab's own
#' Sensitivity sub-tab uses (.rtm_sweep_traits, including LAI); any column a
#' leaf model needs but doesn't expose as a slider (e.g. Fluspect-B's Cx,
#' Fluspect-B-Cx's Prot/CBC, every model's alpha for inform()'s always-run
#' understorey step) stays fixed at its one_sim default instead of being
#' sampled.
many_sims_code_r <- function(canopy, leaf) {
  traits <- .rtm_sweep_traits[[leaf]]
  leaf_p <- .leaf_model_params[[leaf]]
  canopy_p <- .canopy_model_params[[canopy]]
  fn <- .canopy_fn_name[[canopy]]

  sample_args <- vapply(names(traits), function(t) {
    rng <- traits[[t]]
    nm <- if (grepl("\\.", t)) paste0("`", t, "`") else t
    sprintf("%s = runif(n_samples, %s, %s)", nm, rng[1], rng[2])
  }, character(1))
  fixed_extra <- leaf_p[setdiff(names(leaf_p), names(traits))]

  lut_args <- paste(c(
    paste(sample_args, collapse = ", "),
    "hspot = 0.01, LIDFa = -0.35, LIDFb = -0.15, TypeLidf = 1, tts = 30, tto = 0, psi = 0",
    if (length(fixed_extra) > 0) .fmt_params_r(fixed_extra),
    if (length(canopy_p) > 0) .fmt_params_r(canopy_p)
  ), collapse = ",\n  ")

  brf_line <- if (canopy == "inform") {
    sprintf('%s(inputLUT = lut[i, ], rsoil = rsoil, LeafModel = "%s")', fn, leaf)
  } else {
    sprintf(
      'sail_i <- %s(inputLUT = lut[i, ], rsoil = rsoil, LeafModel = "%s")\n  Compute_BRF(rdot = sail_i$rdot, rsot = sail_i$rsot,\n              tts = lut$tts[i], data.light = dataSpec_PDB)',
      fn, leaf
    )
  }

  paste0(
    "n_samples <- 500\nset.seed(42)\n",
    "lut <- data.frame(\n  ", lut_args, "\n)\n\n",
    "spectra <- t(sapply(seq_len(n_samples), function(i) {\n  ", brf_line, "\n}))\n\n",
    "# A few random parameter draws can be numerically unstable for some leaf\n",
    "# models (e.g. Liberty's inter-layer solver, near the edges of its own\n",
    "# range) and come back as NaN -- drop those rows rather than propagating them.\n",
    "ok <- apply(spectra, 1, function(r) all(is.finite(r)))\n",
    "spectra <- spectra[ok, , drop = FALSE]\n",
    "lut <- lut[ok, ]\n",
    "# spectra: one valid simulated spectrum per row (length depends on the leaf model's own domain --\n",
    "# 2101 for PROSPECT-D/-PRO/Liberty, 2001 for Fluspect-B/-B-Cx)"
  )
}

many_sims_code_python <- function(canopy, leaf) {
  traits <- .rtm_sweep_traits[[leaf]]
  leaf_p <- .leaf_model_params[[leaf]]
  canopy_p <- .canopy_model_params[[canopy]]
  fn <- .canopy_fn_name[[canopy]]

  sample_args <- vapply(names(traits), function(t) {
    rng <- traits[[t]]
    sprintf('"%s": rng.uniform(%s, %s, n_samples)', t, rng[1], rng[2])
  }, character(1))
  fixed_extra <- leaf_p[setdiff(names(leaf_p), names(traits))]

  lut_args <- paste(c(
    paste(sample_args, collapse = ", "),
    '"hspot": 0.01, "LIDFa": -0.35, "LIDFb": -0.15, "TypeLidf": 1, "tts": 30, "tto": 0, "psi": 0',
    if (length(fixed_extra) > 0) .fmt_params_py(fixed_extra),
    if (length(canopy_p) > 0) .fmt_params_py(canopy_p)
  ), collapse = ",\n    ")

  import_line <- sprintf("from toolsrtm import %s, compute_brf", fn)
  brf_line <- if (canopy == "inform") {
    sprintf('spectra[i] = np.asarray(%s(row_i, rsoil, leaf_model="%s"))', fn, leaf)
  } else {
    sprintf(
      'sail_i = %s(row_i, rsoil, leaf_model="%s")\n    spectra[i] = compute_brf(sail_i.rdot, sail_i.rsot, row_i["tts"])',
      fn, leaf
    )
  }

  paste0(
    "import numpy as np\n", import_line, "\n\n",
    "n_samples = 500\nrng = np.random.default_rng(42)\n",
    "lut = {\n    ", lut_args, ",\n}\n\n",
    "spectra = np.empty((n_samples, ", if (leaf %in% c("Fluspect-B", "Fluspect-B-Cx")) "2001" else "2101", "))\n",
    "for i in range(n_samples):\n",
    "    row_i = {k: (v[i] if hasattr(v, \"__len__\") else v) for k, v in lut.items()}\n    ",
    brf_line, "\n\n",
    "# A few random parameter draws can be numerically unstable for some leaf\n",
    "# models (e.g. Liberty's inter-layer solver, near the edges of its own\n",
    "# range) and come back as NaN -- drop those rows rather than propagating them.\n",
    "ok = np.isfinite(spectra).all(axis=1)\n",
    "spectra = spectra[ok]\n",
    "lut = {k: (np.asarray(v)[ok] if hasattr(v, \"__len__\") else v) for k, v in lut.items()}\n",
    "# spectra: one valid simulated spectrum per row"
  )
}

.tutorial_r_code <- list(
  setup = 'library(ToolsRTM)',
  sensitivity = 'sweep_trait <- function(trait, values) {
  sapply(values, function(v) {
    row2 <- row; row2[[trait]] <- v
    sail_i <- foursail(inputLUT = row2, rsoil = rsoil, LeafModel = "PROSPECT-D")
    Compute_BRF(rdot = sail_i$rdot, rsot = sail_i$rsot,
                tts = row2$tts, data.light = dataSpec_PDB)
  })
}
cab_curves <- sweep_trait("Cab", seq(10, 70, length.out = 6))     # 2101 x 6
ewt_curves <- sweep_trait("EWT", seq(0.005, 0.03, length.out = 6))
# repeat for LAI, Cbrown, ... -- plot each trait in ITS OWN panel with its
# own color scale (a scale shared across traits hides every trait except
# the one with the largest absolute range -- see Exercise-1.Rmd\'s own fix
# for this exact bug, "sensitivity-plot" chunk)',
  convolution = 'band_values <- get.spectral.convolution.srf(
  df  = data.frame(wave = dataSpec_PDB[, 1], rfl = reflectance),
  srf = ToolsRTM::srf.sentinel2a)
# band_values$wl / $fwhm / $RFL -- one row per Sentinel-2A band, sorted by
# wavelength. For the other 6 bundled sensors (Landsat/OLCI/MODIS), swap in
# get.spectral.convolution.rfl(df, sensor.i = ToolsRTM::LANDSAT8.OLI) instead
# -- see this step\'s own explanation above for which function each sensor needs.',
  convolution_own_sensor = 'df <- data.frame(wave = dataSpec_PDB[, 1], rfl = reflectance)

# Your own 15-band, 3-camera synchronized system -- center + FWHM (nm)
# for every band, e.g. from a calibration sheet:
own_centers <- c(444, 475, 502, 531, 550, 560, 570, 650, 668, 678, 705, 717, 740, 754, 842)
own_fwhm    <- c(28,  32,  18,  14,  12,  27,  14,  16,  14,  14,  10,  12,  18,  10,  57)

band_values <- get.spectral.convolution.gaussian(df, centers = own_centers, fwhm = own_fwhm)
# No measured SRF for your camera -- each band is approximated as a
# Gaussian response with that center/FWHM.

# A pushbroom hyperspectral camera (e.g. Headwall) where only band centers
# are known (copied from its ENVI header\'s "wavelength = {...}" block, no
# "fwhm = {...}" block) -- FWHM is estimated from band spacing instead:
headwall_centers <- c(398.0166, 400.247, 402.4774, 404.7078, 406.9382)  # (full header has 272)
headwall_bands <- get.spectral.convolution.gaussian(df, centers = headwall_centers)

# A bundled sensor with only nominal characteristics (no measured SRF):
enmap_bands <- get.spectral.convolution.gaussian(df, sensor.i = "EnMAP")
modis_bands <- get.spectral.convolution.gaussian(df, sensor.i = "MODIS")',
  ml_inversion = 'library(randomForest)
band_refl <- t(sapply(seq_len(n_samples), function(i) {
  df_i <- data.frame(wave = dataSpec_PDB[, 1], rfl = spectra[i, ])
  get.spectral.convolution.srf(df_i, ToolsRTM::srf.sentinel2a)$rfl
}))
colnames(band_refl) <- paste0("B", seq_len(ncol(band_refl)))
ml_data <- data.frame(band_refl, Cab = lut$Cab)

train_idx <- sample(seq_len(n_samples), size = round(0.7 * n_samples))
rf <- randomForest(Cab ~ ., data = ml_data[train_idx, ], ntree = 300)
pred <- predict(rf, ml_data[-train_idx, ])
obs  <- ml_data$Cab[-train_idx]
cat("R2:", cor(pred, obs)^2, " RMSE:", sqrt(mean((pred - obs)^2)))'
)

.tutorial_python_code <- list(
  setup = 'import numpy as np
from toolsrtm import foursail, compute_brf',
  sensitivity = 'def sweep_trait(trait, values):
    curves = np.empty((len(values), len(wavelength)))
    for j, v in enumerate(values):
        row_j = dict(row); row_j[trait] = v
        sail_j = foursail(row_j, rsoil, leaf_model="PROSPECT-D")
        curves[j] = compute_brf(sail_j.rdot, sail_j.rsot, row_j["tts"])
    return curves

cab_curves = sweep_trait("Cab", np.linspace(10, 70, 6))
ewt_curves = sweep_trait("EWT", np.linspace(0.005, 0.03, 6))
# repeat for LAI, Cbrown, ... -- same per-trait-own-scale reasoning as the R
# version (see this tutorial\'s R tab for why a shared scale is wrong here).',
  convolution = 'from toolsrtm.srf import srf_sentinel2a, spectral_convolution_srf

sensor = srf_sentinel2a()
band_values = spectral_convolution_srf(wavelength, reflectance, sensor)
# band_values.wl / .fwhm / .rfl / .band_names -- one entry per Sentinel-2A
# band, sorted by wavelength. For the other 6 bundled sensors (Landsat/OLCI/
# MODIS), use toolsrtm.smac.spectral_convolution() with a SmacSensor instead
# -- see this step\'s own explanation above for which function each sensor needs.',
  convolution_own_sensor = 'from toolsrtm.srf import spectral_convolution_gaussian

# Your own 15-band, 3-camera synchronized system -- center + FWHM (nm)
# for every band, e.g. from a calibration sheet:
own_centers = [444, 475, 502, 531, 550, 560, 570, 650, 668, 678, 705, 717, 740, 754, 842]
own_fwhm    = [28,  32,  18,  14,  12,  27,  14,  16,  14,  14,  10,  12,  18,  10,  57]

band_values = spectral_convolution_gaussian(wavelength, reflectance, centers=own_centers, fwhm=own_fwhm)
# No measured SRF for your camera -- each band is approximated as a
# Gaussian response with that center/FWHM.

# A pushbroom hyperspectral camera (e.g. Headwall) where only band centers
# are known (copied from its ENVI header\'s "wavelength = {...}" block, no
# "fwhm = {...}" block) -- FWHM is estimated from band spacing instead:
headwall_centers = [398.0166, 400.247, 402.4774, 404.7078, 406.9382]  # (full header has 272)
headwall_bands = spectral_convolution_gaussian(wavelength, reflectance, centers=headwall_centers)

# A bundled sensor with only nominal characteristics (no measured SRF):
enmap_bands = spectral_convolution_gaussian(wavelength, reflectance, sensor="EnMAP")
modis_bands = spectral_convolution_gaussian(wavelength, reflectance, sensor="MODIS")',
  ml_inversion = 'from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, root_mean_squared_error

band_refl = np.array([spectral_convolution_srf(wavelength, spectra[i], sensor).rfl
                       for i in range(n_samples)])
X_train, X_test, y_train, y_test = train_test_split(
    band_refl, lut["Cab"], test_size=0.3, random_state=1)

rf = RandomForestRegressor(n_estimators=300, random_state=1)
rf.fit(X_train, y_train)
pred = rf.predict(X_test)
print("R2:", r2_score(y_test, pred),
      "RMSE:", root_mean_squared_error(y_test, pred))'
)

.tutorial_sif_note_r <- 'sif_df <- compute_leaf_sif(row, "Fluspect-B")   # or "Fluspect-B-Cx"
# sif_df$reflectance:           true leaf reflectance (unaffected by fqe)
# sif_df$apparent_reflectance:  reflectance + SIF / incident irradiance --
#                                what a sensor actually measures, since
#                                at-sensor radiance = reflected + emitted
# sif_df$sif:                   simulated SIF emission (W m-2 nm-1), 640-850nm
# compute_leaf_sif() is defined at the top of this app\'s own app.R -- see
# its docstring for the full physics (leaf-level, not canopy-transported
# SIF), or the Model Explorer tab\'s own SIF sub-tab for a live, interactive
# version (same computation, plotted for you).'

.tutorial_sif_note_python <- 'from toolsrtm import fluspect_b  # or fluspect_cx for Fluspect-B-Cx

wlE, wlF = np.arange(400, 751), np.arange(640, 851)
incident = np.full(len(wlE), 1.5)  # simplified flat placeholder irradiance (W m-2 nm-1)

lrt = fluspect_b(Cab=40, Car=10, EWT=0.015, LMA=0.01, Cs=0.1, N=1.8, fqe=0.01, Cx=0)
sif = (lrt.MbI + lrt.MbII) @ incident                       # simulated SIF emission, 640-850nm
refl_wlf = lrt.refl[(wlF[0] - 400):(wlF[-1] - 400 + 1)]      # true reflectance, same window
apparent_reflectance = refl_wlf + sif / incident[0]          # what a sensor actually measures
# toolsrtm doesn\'t (yet) have a public helper for this -- shown inline here;
# see this app\'s own R app.R, compute_leaf_sif(), for the version with a
# real (non-flat) incident irradiance spectrum, and the Model Explorer
# tab\'s SIF sub-tab for a live, interactive version.'

#' Builds the shared narrative for one "How in <language>" tab: 6 numbered
#' steps (setup, one simulation, n=500 simulations, sensitivity, sensor
#' convolution, ML inversion), each with a code block and -- where relevant
#' -- the corresponding pre-generated figure from www/ (same figures for
#' both languages: the two implementations are verified numerically
#' equivalent, see python/README.md's parity test suite).
tutorial_tab_content <- function(lang, code, id_prefix) {
  step_title <- function(text) h4(strong(text))
  tagList(
    p(style = "font-size:16px;",
      strong(paste0("Reproducing this app's workflow in ", lang, ".")),
      " These 6 steps mirror the ", strong("Model Explorer"), " tab's own pipeline.",
      " Steps 4-6 use fourSAIL + PROSPECT-D as one fixed worked example; steps",
      " 2 and 3 below both adapt live to whichever canopy/leaf model you",
      " pick -- exactly like the ", strong("Model Explorer"), " tab's own",
      " dropdowns."),

    step_title("1. Load the package"),
    p("Just the package import -- nothing computed yet. Every later step",
      " assumes this ran first."),
    code_block(code$setup),

    step_title("2. Run one simulation"),
    p("Build a single-row parameter set (one plant/canopy's worth of trait",
      " values), run the canopy model you pick below, and combine its",
      " direct/diffuse-light-weighted outputs into ONE reflectance spectrum",
      " -- the most basic thing this whole app does, just without a plot."),
    fluidRow(
      column(width = 6, selectInput(paste0(id_prefix, "_canopy_model"), "Canopy model",
                                     choices = .canopy_model_choices, selected = "foursail")),
      column(width = 6, selectInput(paste0(id_prefix, "_leaf_model"), "Leaf model",
                                     choices = .leaf_model_choices, selected = "PROSPECT-D"))
    ),
    uiOutput(paste0(id_prefix, "_one_sim_code")),
    div(style = "max-width: 700px;", img(src = "tutorial_1sim.png", style = "width:100%;")),
    p(em("The code above updates live for whichever combination you pick --",
         " it's a reference for how to call that specific model, on its own.",
         " The figure and steps 4-6 below always continue with the fourSAIL +",
         " PROSPECT-D worked example, regardless of what's selected here (use",
         " the interactive ", strong("Model Explorer"), " tab to actually plot any of the",
         " other 14 combinations, and re-run step 2's own code for one of",
         " them in a fresh session if you want to carry it through the rest",
         " of the pipeline yourself).")),
    uiOutput(paste0(id_prefix, "_fluspect_sif_note")),

    p(strong("A note on BRDF and TOC correction: "), "the ", code("Compute_BRF()"),
      "/", code("compute_brf()"), " call above blends a canopy's",
      " hemispherical-directional (", code("rdot"), ") and bi-directional",
      " (", code("rsot"), ") reflectance factors using the sun/view geometry",
      " (", code("tts"), "/", code("tto"), "/", code("psi"), ") -- this IS the BRDF correction:",
      " it's needed whenever illumination or viewing geometry isn't fixed",
      " nadir-sun-normalized (off-nadir sensors, multi-angle acquisitions,",
      " comparing images from different overpass times/dates), so that the",
      " resulting top-of-canopy (TOC) reflectance is comparable across",
      " geometries rather than confounded by them. For a flat, purely",
      " Lambertian assumption (no angular effects) you'd only need",
      " ", code("rsot"), " itself; using the full BRDF blend is what this",
      " app -- and the real satellites it's teaching you to interpret --",
      " actually do."),

    step_title("3. Many simulations (n_samples = 500)"),
    p("Now do step 2 five hundred times, with the trait values sampled",
      " randomly (not fixed at one baseline) -- this builds a LUT",
      " (look-up table): one row of traits, one simulated spectrum, times",
      " 500. Everything from here on (sensitivity, ML inversion) works from",
      " this LUT, not from single simulations one at a time. Uses the",
      " ", strong("same dropdowns as step 2"), " above, and the same sampling",
      " ranges as the interactive ", strong("Model Explorer"), " tab's own",
      " Sensitivity sub-tab."),
    uiOutput(paste0(id_prefix, "_many_sims_code")),
    div(style = "max-width: 700px;", img(src = "tutorial_many_spectra.png", style = "width:100%;")),
    p(em("Figure above is still the fourSAIL + PROSPECT-D worked example --",
         " steps 4-6 below continue from it regardless of the dropdowns.")),

    step_title("4. Sensitivity: vary one trait at a time"),
    p("Different question from step 3: instead of random combinations,",
      " hold EVERY parameter at one fixed baseline except one, and sweep",
      " just that one trait across its realistic range -- this isolates",
      " which spectral regions respond to which trait (Cab: red/red-edge;",
      " LAI: NIR plateau; EWT: SWIR; Cbrown: visible/red-edge), matching",
      " Exercise-1's own Part 2. The interactive equivalent is the",
      " ", strong("Model Explorer"), " tab's own Sensitivity sub-tab, which",
      " lets you pick the trait (any leaf parameter, LAI, angles, or",
      " canopy-model-specific parameters like fourSAIL2's ", code("fraction_brown"),
      " or INFORM's ", code("LAIu"), ") instead of it being fixed in code."),
    code_block(code$sensitivity),
    div(style = "max-width: 950px;", img(src = "tutorial_sensitivity.png", style = "width:100%;")),

    step_title("5. Convolve to a specific sensor"),
    p("Resample the native 1nm spectrum from step 2 onto a real sensor's",
      " bands (13 for Sentinel-2A here) using its spectral response",
      " function (SRF) -- this is what a satellite actually observes, coarser",
      " and fewer bands than the simulation itself (see the",
      " ", strong("PROSAIL-WithSatellite"), " tab for an interactive version",
      " with FWHM shown too, and any of the other 8 bundled sensors)."),
    p(strong("Three convolution functions, one for each kind of sensor data you might have:")),
    tags$ol(
      tags$li(code("ToolsRTM::get.spectral.convolution.srf()"), " (used below) -- a real, MEASURED",
              " per-nm spectral response function (SRF), no atmospheric-correction coefficients",
              " needed. Covers PRISMA and Sentinel-2A/Sentinel-2B (", code("srf.prisma"), "/",
              code("srf.sentinel2a"), "/", code("srf.sentinel2b"), ")."),
      tags$li(code("ToolsRTM::get.spectral.convolution.rfl()"), " -- a measured SRF bundled",
              " together WITH SMAC atmospheric-correction coefficients. Covers the other 6 bundled",
              " sensors: Landsat 4/5/7/8, Sentinel-3A/B OLCI, Terra/Aqua MODIS."),
      tags$li(code("ToolsRTM::get.spectral.convolution.gaussian()"), " -- no measured SRF at all,",
              " only NOMINAL band characteristics (center wavelength + FWHM, or center + published",
              " band edges), approximated as a Gaussian response. This is what you need for",
              " ", strong("EnMAP"), " (", code('sensor.i = "EnMAP"'), "), for ", strong("ALI, Hyperion,",
              " MODIS (19-band nominal set), Quickbird, RapidEye, WorldView-2"), " and older/nominal",
              " Landsat characteristics (", code('sensor.i = "MODIS"'), ", etc. -- see", " ", code("unique(ToolsRTM::sensor.characteristics$Sensor)"),
              " for the full list), and for ", strong("your own sensor or camera"), " -- pass your own",
              " ", code("centers"), " (and, if you know it, ", code("fwhm"), ") instead of",
              " ", code("sensor.i"), ". If you only know band centers (e.g. copied from an ENVI",
              " header's ", code("wavelength = {...}"), " block, no ", code("fwhm = {...}"), " block",
              " given -- common for pushbroom cameras like Headwall), FWHM is estimated from the",
              " spacing between neighboring bands.")
    ),
    code_block(code$convolution),
    div(style = "max-width: 700px;", img(src = "tutorial_sensor_convolution.png", style = "width:100%;")),
    p(strong("Worked example -- your own sensor: "), "a synchronized 3-camera, 15-band rig with a",
      " known center + FWHM for every band (e.g. from a calibration sheet):"),
    code_block(code$convolution_own_sensor),

    step_title("6. Trait inversion with machine learning"),
    p("The reverse problem from every step above: given ONLY the sensor-band",
      " reflectance (what a satellite actually measures), retrieve the trait",
      " (Cab here) that produced it, without knowing the ground truth. Train",
      " a model on the step-3 LUT's (sensor-band reflectance -> trait) pairs,",
      " then evaluate it on held-out rows the model never saw during",
      " training -- matching Exercise-1's own Part 3, and the whole reason",
      " step 3 built a LUT of many simulations instead of just one."),
    code_block(code$ml_inversion),
    div(style = "max-width: 800px;", img(src = "tutorial_ml_inversion.png", style = "width:100%;")),

    step_title("What changes on the other tabs"),
    tags$ul(
      tags$li(strong("PROSAIL-WithSatellite"), ": steps 1 + 5 together,",
              " interactively, with the model fixed to fourSAIL + PROSPECT-PRO",
              " -- pick any bundled sensor from a dropdown instead of",
              " hardcoding Sentinel-2A."),
      tags$li(strong("Model Explorer"), ": step 1, generalized -- pick any of the 3",
              " canopy models (fourSAIL/fourSAIL2/INFORM) and any of the 5",
              " leaf models interactively; the LeafModel argument in step 1's",
              " code is exactly the dropdown's current value; plus the",
              " interactive SIF, Sensitivity, and Absorption coefficients",
              " sub-tabs described above.")
    ),
    br()
  )
}

# 3. Define UI ---------------------------------------------------------------------------------

leaf_canopy_sidebar <- function() {
  tagList(
    fluidRow(
      h4("Leaf parameters"),
      column(width = 4,
             sliderInput("Cab", "Chlorophyll content", min = 0, max = 80, value = 40),
             sliderInput("Cbrown", "Brown pigments", min = 0, max = 1, step = 0.05, value = 0.0)),
      column(width = 4,
             sliderInput("Car", "Carotenoid content", min = 0, max = 20, value = 8),
             sliderInput("EWT", "Equivalent water thickness", min = 0.0001, max = 0.05, value = 0.009)),
      column(width = 4,
             sliderInput("Anth", "Anthocyanin Content", min = 0, max = 7, step = 0.2, value = 2),
             sliderInput("CBC", "Carbon-based constituent", min = 0.0001, max = 0.03, value = 0.005)),
      column(width = 4, sliderInput("N", "Structure parameter N", min = 1, max = 3, step = 0.2, value = 1.5)),
      column(width = 4, sliderInput("Prot", "Leaf protein content", min = 0.0001, max = 0.03, value = 0.0045))
    ),
    fluidRow(
      h4("Canopy parameters"),
      column(width = 4, sliderInput("LAI", "Leaf Area Index", min = 0.001, max = 10, step = 0.1, value = 4)),
      column(width = 4, sliderInput("hspot", "Hotspot parameter", min = 0, max = 1, value = 0.01)),
      column(width = 4,
             selectInput("TypeLIDF", "Type LIDF",
                         choices = list(Planophile = "plano", Erectophile = "erecto", Plagiophile = "plagio",
                                        Extremophile = "extremo", Uniform = "uniform", Spherical = "sph"),
                         selected = "sph"))
    ),
    fluidRow(
      h4("Soil parameters"),
      column(width = 4, sliderInput("psoil", "Soil brightness", min = 0, max = 1, value = 0.5)),
      column(width = 4, checkboxInput("include_soil_reflectance", "Include soil reflectance", value = FALSE))
    ),
    fluidRow(
      h4("Angle parameters"),
      column(width = 4, sliderInput("tts", "Solar zenith angle", min = 0, max = 90, value = 30)),
      column(width = 4, sliderInput("tto", "Observer zenith angle", min = 0, max = 90, value = 10)),
      column(width = 4, sliderInput("psi", "Relative azimuth angle", min = 0, max = 180, value = 0))
    )
  )
}

# Narrowly scoped on purpose -- an earlier, broader CSS pass (global font
# scaling) made things worse (clipped plot titles, wrapped sidebar labels;
# see AEO-Course/README.md) and was reverted. This one only pins the top
# navbar and makes the sidebar's own scroll independent from the main
# content's, via position:sticky -- nothing about font/text sizing at all.
# Bootstrap 3's default navbar (shinytheme flatly is BS3-based) is 50px
# tall; the sidebar's own sticky offset matches that so it starts right
# below the pinned navbar instead of under it.
.sticky_sidebar_css <- "
.navbar { position: sticky; top: 0; z-index: 1030; margin-bottom: 0; }
.well { position: sticky; top: 50px; max-height: calc(100vh - 70px); overflow-y: auto; }
"

ui <- navbarPage(
  "Online reflectance simulator", theme = shinytheme("flatly"), collapsible = TRUE, fluid = TRUE,
  id = "main_nav",
  header = tags$head(tags$style(HTML(.sticky_sidebar_css))),

  tabPanel(
    title = "Model Explorer",
    sidebarLayout(
      sidebarPanel(
        width = 5,
        fluidRow(
          h4("Model choice"),
          column(width = 6,
                 selectInput("rtm_canopy_model", "Canopy model", choices = .canopy_model_choices,
                             selected = "foursail")),
          column(width = 6,
                 selectInput("rtm_leaf_model", "Leaf model", choices = .leaf_model_choices,
                             selected = "PROSPECT-D"))
        ),
        h4("Leaf parameters"),
        uiOutput("rtm_leaf_ui"),
        h4("Canopy parameters"),
        fluidRow(
          column(width = 4, sliderInput("rtm_LAI", "Leaf Area Index", min = 0.001, max = 10, step = 0.1, value = 3)),
          column(width = 4, sliderInput("rtm_hspot", "Hotspot parameter", min = 0, max = 1, value = 0.01)),
          column(width = 4,
                 selectInput("rtm_lidf_mode", "Leaf angle distribution (TypeLidf)",
                             choices = list("1: named shape (Planophile, Erectophile, ...)" = "1",
                                            "2: set mean leaf angle directly (degrees)" = "2"),
                             selected = "1"),
                 tags$small(
                   "1 = pick one of 6 canonical leaf-orientation shapes. ",
                   "2 = choose the average leaf tilt yourself, in degrees (0 deg = flat, 90 deg = upright). See ",
                   tags$strong("Reference > Canopy models"), " for what each option actually looks like.",
                   style = "color:#666;"))
        ),
        # TypeLidf=1: pick one of the 6 canonical shapes (LIDFa/LIDFb are
        # preset codes, see .rtm_reference_lidf_presets in the Reference
        # tab). TypeLidf=2: LIDFa becomes the mean leaf angle directly (deg)
        # and LIDFb is unused -- same ellipsoidal/Campbell mode the
        # Sensitivity tab's own LIDFa sweep always uses.
        conditionalPanel(
          condition = "input.rtm_lidf_mode == '1'",
          fluidRow(column(width = 6,
            selectInput("rtm_TypeLIDF", "Type LIDF",
                        choices = list(Planophile = "plano", Erectophile = "erecto", Plagiophile = "plagio",
                                       Extremophile = "extremo", Uniform = "uniform", Spherical = "sph"),
                        selected = "sph")))
        ),
        conditionalPanel(
          condition = "input.rtm_lidf_mode == '2'",
          fluidRow(
            column(width = 6, sliderInput("rtm_LIDFa_ellipsoidal", "LIDFa: mean leaf angle (deg)",
                                           min = 0, max = 90, value = 45)),
            column(width = 6, tags$p(tags$em("LIDFb = 0 in this mode (ellipsoidal/Campbell -- see Reference tab)."),
                                      style = "margin-top: 30px;"))
          )
        ),
        uiOutput("rtm_canopy_extra_ui"),
        h4("Soil & angles"),
        fluidRow(
          column(width = 3, sliderInput("rtm_psoil", "Soil brightness", min = 0, max = 1, value = 0.5)),
          column(width = 3, sliderInput("rtm_tts", "Solar zenith", min = 0, max = 90, value = 30)),
          column(width = 3, sliderInput("rtm_tto", "Observer zenith", min = 0, max = 90, value = 0)),
          column(width = 3, sliderInput("rtm_psi", "Relative azimuth", min = 0, max = 180, value = 0))
        ),
        fluidRow(column(width = 12,
          checkboxInput("rtm_include_soil", "Include soil reflectance", value = FALSE)
        )),
        fluidRow(column(
          width = 12,
          div(style = "margin-top: 10px;", downloadButton("downloadRtmData", "Download RTM spectrum"))
        ))
      ),
      mainPanel(
        width = 7,
        tabsetPanel(
          id = "rtm_tabs",
          tabPanel("Simulation",
                   h4("Simulated canopy reflectance"),
                   div(style = "height: 60vh; min-height: 340px; overflow-y: auto;",
                       plotOutput(outputId = "rtm_plot", height = "100%"))),
          tabPanel("SIF (fluorescence)",
                   conditionalPanel(
                     condition = "input.rtm_leaf_model == 'Fluspect-B' || input.rtm_leaf_model == 'Fluspect-B-Cx'",
                     h4("Leaf-level fluorescence (SIF)"),
                     div(style = "height: 60vh; min-height: 340px; overflow-y: auto;",
                         plotOutput(outputId = "rtm_sif_plot", height = "100%"))
                   ),
                   conditionalPanel(
                     condition = "input.rtm_leaf_model != 'Fluspect-B' && input.rtm_leaf_model != 'Fluspect-B-Cx'",
                     div(style = "padding-top: 40px; text-align: center; color: #888; font-size: 16px;",
                         "Select Fluspect-B or Fluspect-B-Cx as the leaf model (left) to see simulated SIF.")
                   )),
          tabPanel("Sensitivity",
                   h4("Vary one trait, keep everything else at its current value"),
                   uiOutput("rtm_sens_trait_ui"),
                   # Taller than most plot containers -- Fluspect leaf models
                   # show 3 stacked panels here (reflectance, apparent
                   # reflectance, SIF), not just 1.
                   div(style = "height: 85vh; min-height: 460px; overflow-y: auto;",
                       plotOutput(outputId = "rtm_sensitivity_plot", height = "100%")),
                   h4("Fixed values used for this sweep"),
                   p("Everything below except the swept trait itself (shown as its own",
                     " range, not a single fixed value) was held at exactly this value."),
                   tableOutput("rtm_sens_fixed_table"),
                   h4("Reproduce this sweep in code"),
                   p("Copy/paste and run -- uses the exact values above, not placeholder defaults."),
                   tabsetPanel(
                     tabPanel("R", br(), uiOutput("rtm_sens_code_r")),
                     tabPanel("Python", br(), uiOutput("rtm_sens_code_py"))
                   )
                   ),
          tabPanel("Absorption coefficients",
                   p("Specific absorption coefficients for whichever leaf model is picked",
                     " on the left -- genuinely different data per model (PROSPECT-D, PROSPECT-PRO,",
                     " Fluspect-B/-B-Cx each have their own), not the same table relabeled."),
                   div(style = "height: 40vh; min-height: 280px;",
                       plotOutput(outputId = "abs_plot1", height = "100%")),
                   div(style = "height: 40vh; min-height: 280px;",
                       plotOutput(outputId = "abs_plot2", height = "100%")),
                   conditionalPanel(
                     condition = "input.rtm_leaf_model == 'Fluspect-B' || input.rtm_leaf_model == 'Fluspect-B-Cx'",
                     h4("Carotenoid absorption states (Fluspect only)"),
                     p("Both Fluspect variants interpolate carotenoid absorption between a",
                       " violaxanthin state (Cx=0) and a zeaxanthin state (Cx=1) -- the",
                       " \"Carotenoids content\" line above is this interpolation at your",
                       " CURRENT Cx value, shown here explicitly."),
                     div(style = "height: 40vh; min-height: 280px;",
                         plotOutput(outputId = "abs_plot3", height = "100%"))
                   ))
        )
      )
    )
  ),

  tabPanel(
    title = "PROSAIL-WithSatellite",
    sidebarLayout(
      sidebarPanel(
        width = 5,
        fluidRow(
          h4("Satellite sensor"),
          column(width = 8,
                 selectInput("sat_sensor", "Compare to satellite", choices = .sensor_choices,
                             selected = "Sentinel2A.MSI")),
          column(width = 4,
                 checkboxInput("show_fwhm", "Show FWHM", value = TRUE))
        ),
        leaf_canopy_sidebar(),
        fluidRow(column(
          width = 12,
          div(style = "margin-top: 10px;", downloadButton("downloadSatData", "Download sensor-band data"))
        ))
      ),
      mainPanel(
        width = 7,
        h4("Native (1nm) vs. satellite-band reflectance"),
        div(style = "height: 60vh; min-height: 320px; overflow-y: auto;",
            plotOutput(outputId = "sat_plot", height = "100%")),
        h4("Sensor band values"),
        div(style = "height: 38vh; min-height: 220px; overflow-y: auto;",
            tableOutput(outputId = "sat_table"))
      )
    )
  ),

  tabPanel(
    title = "How in R",
    fluidRow(column(width = 10, offset = 1,
      tutorial_tab_content("R", .tutorial_r_code, id_prefix = "tut_r")
    ))
  ),

  tabPanel(
    title = "How in Python",
    fluidRow(column(width = 10, offset = 1,
      tutorial_tab_content("Python", .tutorial_python_code, id_prefix = "tut_py")
    ))
  ),

  tabPanel(
    title = "Install",
    fluidRow(column(width = 10, offset = 1,
      p(style = "font-size:16px;",
        strong("Installing ToolsRTM/SCOPEinR (R) and toolsrtm/scopeinpython (Python)."),
        " Two kinds of install for each language -- from the packages' own published repo",
        " (stable, versioned releases), or straight from this monorepo's own source (always",
        " the latest code, e.g. the fixes made throughout this course) -- plus where each",
        " package's manual lives. Same commands whether you're running this app locally or",
        " on a server -- nothing app-specific here."),

      h4(strong("R: ToolsRTM and SCOPEinR")),
      tags$p(strong("Option 1 -- from GitLab (published releases):")),
      code_block('if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_gitlab("caminoccg/toolsrtm")
remotes::install_gitlab("caminoccg/scopeinr")

library(ToolsRTM)
library(SCOPEinR)'),
      tags$p(strong("Option 2 -- from a local clone of this monorepo (0-RTM-Suite), always the",
                     " latest source (e.g. any fixes made while working through this course):")),
      code_block('git clone <0-RTM-Suite GITLAB URL>
cd 0-RTM-Suite

# from R:
devtools::install("ToolsRTM")
devtools::install("SCOPEinR")
# (devtools::document() first if you edited any .R file\'s roxygen comments)

library(ToolsRTM)
library(SCOPEinR)'),

      h4(strong("Python: toolsrtm and scopeinpython")),
      tags$p(strong("Option 1 -- from GitHub (published releases, in progress):"), " each Python",
             " package is being split into its own GitHub repo, mirroring the R packages' GitLab",
             " layout above -- ", code("ToolsRTMinPython"), " and ", code("SCOPEinRinPython"),
             " (exact org/URLs to be confirmed; the install command will be",
             " ", code("pip install git+<repo URL>"), " once published, no cloning needed)."),
      tags$p(strong("Option 2 -- from a local clone of this monorepo (works today), editable so",
                     " local changes take effect immediately:")),
      code_block('git clone <0-RTM-Suite GITLAB/GITHUB URL>
cd 0-RTM-Suite

pip install -e python/toolsrtm
pip install -e python/scopeinpython   # depends on toolsrtm, install this one second

# to also run the test suites:
pip install -e "python/toolsrtm[test]"
pip install -e "python/scopeinpython[test]"'),
      p("Verify either language's install by reproducing one line from this app's own",
        " ", strong("How in R"), "/", strong("How in Python"), " tabs -- if step 1's ",
        code("library(ToolsRTM)"), "/", code("from toolsrtm import foursail, compute_brf"),
        " runs without error, the install worked."),

      h4(strong("Where the manuals are")),
      tags$ul(
        tags$li(strong("ToolsRTM (R): "), code("docs/toolsrtm/index.html"), " -- pkgdown site,",
                " every exported function's full help page, browsable offline (open the file",
                " directly in a browser, no server needed)."),
        tags$li(strong("SCOPEinR (R): "), code("docs/scopeinr/index.html"), " -- same, for SCOPEinR."),
        tags$li(strong("toolsrtm + scopeinpython (Python): "), code("docs/python/index.html"),
                " -- one combined Sphinx/Read-the-Docs-style site for both Python packages,",
                " including the install page and the numerical-verification-against-R writeup",
                " (", code("python/README.md"), " has the same content in plain Markdown, if",
                " you'd rather read it on GitHub/GitLab directly without opening the built site)."),
        tags$li("All three are paths inside this same ", code("0-RTM-Suite"), " repo -- once the",
                " GitLab/GitHub repos above are live, these would typically also be published as",
                " GitLab/GitHub Pages sites, not yet set up.")
      )
    ))
  ),

  tabPanel(
    title = "Reference",
    fluidRow(column(width = 10, offset = 1,
      p(style = "font-size:16px;",
        strong("Full parameter and sensor-band reference."),
        " Every leaf/canopy model parameter this app's sliders expose, with its",
        " description, range, and units, plus every bundled satellite sensor's",
        " own bands (center wavelength and FWHM) -- the same numbers the",
        " sliders and ", strong("PROSAIL-WithSatellite"), " tab themselves use,",
        " not a separately-maintained copy."),
      tabsetPanel(
        tabPanel("Leaf models",
                 br(),
                 selectInput("reference_leaf_model", "Leaf model", choices = .leaf_model_choices,
                             selected = "PROSPECT-D", width = "300px"),
                 tableOutput("reference_leaf_table"),
                 conditionalPanel(
                   condition = "input.reference_leaf_model == 'PROSPECT-PRO'",
                   tags$p(tags$strong("* LMA vs Prot + CBC: "),
                          "PROSPECT-D/Fluspect use one lumped dry-matter number (LMA); PROSPECT-PRO splits it into its two chemical components instead."),
                   tableOutput("reference_lma_prot_cbc")
                 )),
        tabPanel("Canopy models",
                 br(),
                 selectInput("reference_canopy_model", "Canopy model", choices = .canopy_model_choices,
                             selected = "foursail", width = "300px"),
                 tableOutput("reference_canopy_table"),
                 tags$p(tags$strong("* Leaf-angle distribution (LIDFa / LIDFb): "),
                        "these two parameters change meaning depending on TypeLidf -- see the preset/mode table below."),
                 tableOutput("reference_lidf_presets"),
                 tags$p(tags$strong("Leaf inclination angle distributions: "),
                        "how much leaf area falls in each angle bin (0 deg = horizontal, 90 deg = vertical) for the 6 TypeLidf=1 presets, plus 3 example mean angles for TypeLidf=2 (ellipsoidal/Campbell) -- this is what changing LIDFa/LIDFb actually does to the canopy's leaves."),
                 plotOutput("reference_lidf_plot", height = "600px")),
        tabPanel("Satellite sensors",
                 br(),
                 selectInput("reference_sensor", "Sensor", choices = .sensor_choices,
                             selected = "Sentinel2A.MSI", width = "300px"),
                 p("FWHM values come from each sensor's own spectral response function",
                   " (empirical half-max crossing), not a nominal spec sheet number."),
                 tags$p(strong("Which ToolsRTM function convolves each sensor: "), "two functions,",
                   " split by what data each sensor actually has bundled:"),
                 tags$ul(
                   tags$li(strong("Landsat 4/5/7/8, Sentinel-3A/B OLCI, Terra/Aqua MODIS"), ": ",
                           code("ToolsRTM::get.spectral.convolution.rfl()"), " -- uses each sensor's SMAC",
                           " atmospheric-correction-bundle SRF (band-ordering bug for these 6 fixed",
                           " at the source this round: the raw ", code("wl_smac"), " array order does",
                           " NOT match the SRF table's own natural band order for any of them)."),
                   tags$li(strong("PRISMA, Sentinel-2A, Sentinel-2B"), ": ",
                           code("ToolsRTM::get.spectral.convolution.srf()"), " -- these 3 also have a",
                           " plain, publisher-original per-band SRF table with no SMAC coefficients",
                           " attached at all (", code("srf.prisma"), "/", code("srf.sentinel2a"), "/",
                           code("srf.sentinel2b"), "). PRISMA has ", em("only"), " this (no SMAC bundle",
                           " exists for it); Sentinel-2A/B have both, and this app always uses the plain",
                           " table for them.")
                 ),
                 tags$p(strong("Sentinel-2A vs Sentinel-2B: "), "not two copies of the same sensor --",
                        " twin satellites with genuinely different measured SRFs (manufacturing",
                        " tolerances). Band centers agree to within ~1nm for most bands but differ by",
                        " up to ~17nm in B12 (SWIR2); B7/B9/B10/B11 differ by ~2-3nm. Pick whichever",
                        " satellite actually acquired the real imagery you're comparing against; for a",
                        " generic \"Sentinel-2\" exercise either is fine except when B7/B9/B10/B11/B12",
                        " specifically matter."),
                 tableOutput("reference_sensor_table"))
      )
    ))
  )
)

# 4. Define server logic ------------------------------------------------------------------------

server <- function(input, output, session) {

  parameters_model <- reactive({
    c(
      N = as.numeric(input$N), Cab = as.numeric(input$Cab), Car = as.numeric(input$Car),
      Anth = as.numeric(input$Anth), Cbrown = as.numeric(input$Cbrown), EWT = as.numeric(input$EWT),
      LMA = 0, alpha = 40,
      Prot = as.numeric(input$Prot), CBC = as.numeric(input$CBC),
      LAI = as.numeric(input$LAI), TypeLidf = 1, hspot = as.numeric(input$hspot),
      tts = as.numeric(input$tts), tto = as.numeric(input$tto), psi = as.numeric(input$psi),
      psoil = as.numeric(input$psoil)
    )
  })

  lut_data.sim <- reactive({
    lut.to_sim <- data.frame(t(as.data.frame(parameters_model())))
    row.names(lut.to_sim) <- NULL

    data <- ToolsRTM::dataSpec_PDB
    Rsoil.dry <- data[, 11]
    Rsoil.wet <- data[, 12]
    psoil <- lut.to_sim[1, "psoil"]
    rsoil <- psoil * Rsoil.dry + (1 - psoil) * Rsoil.wet

    list(lut.to_sim, rsoil)
  })

  reflectance_data <- reactive({
    df.LUT <- lut_data.sim()[[1]]
    row.names(df.LUT) <- NULL
    rsoil <- lut_data.sim()[[2]]

    df.LUT$LIDFa <- switch(input$TypeLIDF, plano = 1, erecto = -1, plagio = 0, extremo = 0, uniform = 0, sph = -0.35)
    df.LUT$LIDFb <- switch(input$TypeLIDF, plano = 0, erecto = 0, plagio = -1, extremo = 1, uniform = 0, sph = -0.15)

    sim_leaf_values <- ToolsRTM::prospect_PRO(
      N = df.LUT$N, Cab = df.LUT$Cab, Car = df.LUT$Car, Anth = df.LUT$Anth, Cbrown = df.LUT$Cbrown,
      EWT = df.LUT$EWT, LMA = df.LUT$LMA, alpha = 40, Prot = df.LUT$Prot, CBC = df.LUT$CBC
    )
    wavelength <- sim_leaf_values[[1]]

    reflectance_values <- ToolsRTM::foursail(inputLUT = df.LUT[1, ], rsoil = rsoil, LeafModel = "PROSPECT-PRO")
    rdot <- reflectance_values[[1]]
    rsot <- reflectance_values[[2]]

    reflectance_values <- ToolsRTM::Compute_BRF(rdot = rdot, rsot = rsot, tts = df.LUT[1, "tts"], data.light = ToolsRTM::dataSpec_PDB)
    showNotification("Forward simulation done successfully.", type = "message", duration = 3)

    data.frame(wavelength = wavelength, reflectance = reflectance_values)
  })

  # Sensor-band convolution, shared by the PROSAIL-WithSatellite plot/table/download ---------------------
  sat_data <- reactive({
    df.plot <- reflectance_data()
    out <- convolve_to_sensor(df.plot$wavelength, df.plot$reflectance, input$sat_sensor)
    showNotification(
      paste0("Convolved to ", names(.sensor_choices)[.sensor_choices == input$sat_sensor], " bands."),
      type = "message", duration = 3
    )
    out
  })

  # PROSAIL-WithSatellite tab plot -----------------------------------------------------------------------
  output$sat_plot <- renderPlot({
    df.plot <- reflectance_data()
    df.sat <- sat_data()
    validate(need(nrow(df.sat) > 0, "No valid bands found for this sensor."))

    p <- ggplot()

    if (isTRUE(input$show_fwhm) && all(c("wl", "fwhm") %in% colnames(df.sat))) {
      df.rect <- df.sat[!is.na(df.sat$fwhm), ]
      if (nrow(df.rect) > 0) {
        p <- p + geom_rect(
          data = df.rect,
          aes(xmin = wl - fwhm / 2, xmax = wl + fwhm / 2, ymin = -Inf, ymax = Inf),
          fill = "grey60", alpha = 0.25, inherit.aes = FALSE
        )
      }
    }

    # Native curve drawn LAST (on top) and noticeably thicker than the
    # satellite points/line, so it stays clearly visible as the reference
    # curve rather than getting visually lost behind the band markers.
    p <- p +
      geom_line(data = df.sat, aes(x = wl, y = reflectance, color = "Satellite bands"), linewidth = 1.1) +
      geom_point(data = df.sat, aes(x = wl, y = reflectance, color = "Satellite bands"), size = 3.5) +
      geom_line(data = df.plot, aes(x = wavelength, y = reflectance, color = "Native (1nm)"), linewidth = 1.6)

    color_values <- c("Native (1nm)" = "#1b7837", "Satellite bands" = "#b2182b")
    y_max <- max(0.5, max(df.plot$reflectance) * 1.05)

    if (input$include_soil_reflectance) {
      rsoil <- lut_data.sim()[[2]]
      df.soil <- data.frame(wavelength = df.plot$wavelength, reflectance = rsoil)
      p <- p + geom_line(data = df.soil, aes(x = wavelength, y = reflectance, color = "Soil"), linewidth = 1.2)
      color_values <- c(color_values, "Soil" = "#8c510a")
      y_max <- max(y_max, max(rsoil) * 1.05)
      showNotification("Soil reflectance spectrum added.", type = "message", duration = 3)
    }

    p +
      scale_color_manual(values = color_values) +
      labs(x = "Wavelength (nm)", y = "Reflectance", color = "",
           title = paste("Native vs.", names(.sensor_choices)[.sensor_choices == input$sat_sensor],
                          "band reflectance"),
           subtitle = if (isTRUE(input$show_fwhm)) "Grey polygons: each band's full width at half maximum (FWHM)" else NULL) +
      coord_cartesian(ylim = c(0, y_max)) +
      theme_prosail()
  }, res = render_plot_res)

  output$sat_table <- renderTable({
    df.sat <- sat_data()
    df.sat$band_label <- paste0("B", df.sat$band)
    df.sat[, c("band_label", "wl", "fwhm", "reflectance")] |>
      setNames(c("Band", "Center wavelength (nm)", "FWHM (nm)", "Reflectance"))
  }, digits = c(0, 0, 1, 1, 4))

  # RTM tab: dynamic leaf-parameter panel, switches with input$rtm_leaf_model ------------------
  output$rtm_leaf_ui <- renderUI({
    lm <- input$rtm_leaf_model
    if (lm == "PROSPECT-D") {
      tagList(fluidRow(
        column(4, sliderInput("rtm_N", "N", min = 1, max = 3, step = 0.2, value = 1.5)),
        column(4, sliderInput("rtm_Cab", "Cab", min = 0, max = 80, value = 40)),
        column(4, sliderInput("rtm_Car", "Car", min = 0, max = 20, value = 8)),
        column(4, sliderInput("rtm_Anth", "Anth", min = 0, max = 7, step = 0.2, value = 2)),
        column(4, sliderInput("rtm_Cbrown", "Cbrown", min = 0, max = 1, step = 0.05, value = 0)),
        column(4, sliderInput("rtm_EWT", "EWT", min = 0.0001, max = 0.05, value = 0.009)),
        column(4, sliderInput("rtm_LMA", "LMA", min = 0.0001, max = 0.03, value = 0.009))
      ))
    } else if (lm == "PROSPECT-PRO") {
      tagList(fluidRow(
        column(4, sliderInput("rtm_N", "N", min = 1, max = 3, step = 0.2, value = 1.5)),
        column(4, sliderInput("rtm_Cab", "Cab", min = 0, max = 80, value = 40)),
        column(4, sliderInput("rtm_Car", "Car", min = 0, max = 20, value = 8)),
        column(4, sliderInput("rtm_Anth", "Anth", min = 0, max = 7, step = 0.2, value = 2)),
        column(4, sliderInput("rtm_Cbrown", "Cbrown", min = 0, max = 1, step = 0.05, value = 0)),
        column(4, sliderInput("rtm_EWT", "EWT", min = 0.0001, max = 0.05, value = 0.009)),
        column(4, sliderInput("rtm_Prot", "Protein content", min = 0.0001, max = 0.03, value = 0.0045)),
        column(4, sliderInput("rtm_CBC", "Carbon-based constituent", min = 0.0001, max = 0.03, value = 0.005))
      ))
    } else if (lm %in% c("Fluspect-B", "Fluspect-B-Cx")) {
      tagList(fluidRow(
        column(4, sliderInput("rtm_N", "N", min = 1, max = 3, step = 0.2, value = 1.5)),
        column(4, sliderInput("rtm_Cab", "Cab", min = 0, max = 80, value = 40)),
        column(4, sliderInput("rtm_Car", "Car", min = 0, max = 20, value = 8)),
        column(4, sliderInput("rtm_Anth", "Anth", min = 0, max = 7, step = 0.2, value = 2)),
        column(4, sliderInput("rtm_EWT", "EWT", min = 0.0001, max = 0.05, value = 0.009)),
        column(4, sliderInput("rtm_LMA", "LMA", min = 0.0001, max = 0.03, value = 0.009)),
        column(4, sliderInput("rtm_Cs", "Senescence factor Cs", min = 0, max = 1, step = 0.05, value = 0)),
        column(4, sliderInput("rtm_fqe", "Fluorescence quantum yield fqe", min = 0, max = 0.05, step = 0.005, value = 0.01)),
        if (lm == "Fluspect-B-Cx") {
          column(4, sliderInput("rtm_Cx", "Zeaxanthin fraction Cx", min = 0, max = 1, step = 0.05, value = 0))
        }
      ))
    } else if (lm == "Liberty") {
      tagList(fluidRow(
        column(4, sliderInput("rtm_cell.d", "Cell diameter", min = 20, max = 100, value = 40)),
        column(4, sliderInput("rtm_inter.c", "Intercellular air space", min = 0.01, max = 0.1, step = 0.005, value = 0.045)),
        column(4, sliderInput("rtm_baseline.abs", "Baseline absorption", min = 0.0001, max = 0.001, step = 0.0001, value = 0.0006)),
        column(4, sliderInput("rtm_leaf.thick", "Leaf thickness", min = 1, max = 10, step = 0.1, value = 1.6)),
        column(4, sliderInput("rtm_albino.abs", "Albino absorption", min = 0, max = 4, step = 0.1, value = 2)),
        column(4, sliderInput("rtm_Cab", "Cab", min = 0, max = 80, value = 40)),
        column(4, sliderInput("rtm_EWT", "EWT", min = 0.0001, max = 0.05, value = 0.009)),
        column(4, sliderInput("rtm_lign.cell", "Lignin + cellulose", min = 1, max = 8, step = 0.5, value = 4)),
        column(4, sliderInput("rtm_Nitrogen", "Nitrogen content", min = 0.3, max = 2, step = 0.1, value = 1))
      ))
    }
  })

  # RTM tab: extra canopy-parameter panel, only for fourSAIL2/INFORM -----------------------------
  output$rtm_canopy_extra_ui <- renderUI({
    cm <- input$rtm_canopy_model
    if (cm == "foursail2") {
      tagList(fluidRow(
        column(3, sliderInput("rtm_fraction_brown", "Fraction brown", min = 0, max = 1, step = 0.05, value = 0.2)),
        column(3, sliderInput("rtm_diss", "Layer dissociation", min = 0, max = 1, step = 0.05, value = 0.5)),
        column(3, sliderInput("rtm_Cv", "Crown cover Cv", min = 0, max = 1, step = 0.05, value = 1)),
        column(3, sliderInput("rtm_Zeta", "Tree shape Zeta", min = 0.1, max = 5, step = 0.1, value = 1))
      ))
    } else if (cm == "inform") {
      tagList(fluidRow(
        column(4, sliderInput("rtm_LAIu", "Understorey LAI", min = 0, max = 3, step = 0.1, value = 0.1)),
        column(4, sliderInput("rtm_cd", "Crown diameter (m)", min = 1, max = 15, step = 0.5, value = 4.5)),
        column(4, sliderInput("rtm_h", "Tree height (m)", min = 2, max = 40, value = 20)),
        column(4, sliderInput("rtm_sd", "Stem density (ha-1)", min = 50, max = 2000, step = 50, value = 650)),
        column(4, sliderInput("rtm_skyl", "Diffuse fraction skyl", min = 0, max = 1, step = 0.05, value = 0.1))
      ))
    } else {
      NULL
    }
  })

  # RTM tab: build the LUT row + dispatch to foursail/foursail2/inform ---------------------------
  rtm_row <- reactive({
    lm <- input$rtm_leaf_model
    req(input$rtm_Cab, input$rtm_EWT)  # wait until rtm_leaf_ui has actually rendered

    # TypeLidf=1 (preset shape, default) vs TypeLidf=2 (ellipsoidal/Campbell
    # -- LIDFa becomes the mean leaf angle directly, LIDFb unused/0) --
    # rtm_lidf_mode is the sidebar's "Leaf angle distribution (TypeLidf)"
    # selector; conditionalPanel only toggles CSS visibility of the two
    # underlying controls, both stay registered as inputs regardless of
    # which is shown.
    if (identical(input$rtm_lidf_mode, "2")) {
      TypeLidf <- 2
      LIDFa <- input$rtm_LIDFa_ellipsoidal
      LIDFb <- 0
    } else {
      TypeLidf <- 1
      LIDFa <- switch(input$rtm_TypeLIDF, plano = 1, erecto = -1, plagio = 0, extremo = 0, uniform = 0, sph = -0.35)
      LIDFb <- switch(input$rtm_TypeLIDF, plano = 0, erecto = 0, plagio = -1, extremo = 1, uniform = 0, sph = -0.15)
    }

    row <- data.frame(
      LAI = input$rtm_LAI, hspot = input$rtm_hspot, LIDFa = LIDFa, LIDFb = LIDFb, TypeLidf = TypeLidf,
      tts = input$rtm_tts, tto = input$rtm_tto, psi = input$rtm_psi,
      # Union of every leaf model's columns, defaulted to 0/harmless placeholders --
      # each dispatch below only ever reads the subset its LeafModel actually uses.
      N = 1.5, Cab = 40, Car = 8, Anth = 0, Cbrown = 0, EWT = 0.01, LMA = 0, alpha = 40,
      Prot = 0, CBC = 0, Cs = 0, fqe = 0.01, Cx = 0,
      cell.d = 40, inter.c = 0.045, baseline.abs = 0.0006, leaf.thick = 1.6, albino.abs = 2,
      lign.cell = 4, Nitrogen = 1,
      fraction_brown = 0.2, diss = 0.5, Cv = 1, Zeta = 1,
      LAIu = 0.1, cd = 4.5, h = 20, sd = 650, skyl = 0.1
    )

    if (lm %in% c("PROSPECT-D", "PROSPECT-PRO")) {
      row$N <- input$rtm_N; row$Cab <- input$rtm_Cab; row$Car <- input$rtm_Car
      row$Anth <- input$rtm_Anth; row$Cbrown <- input$rtm_Cbrown; row$EWT <- input$rtm_EWT
      if (lm == "PROSPECT-D") {
        row$LMA <- input$rtm_LMA
      } else {
        row$LMA <- 0; row$Prot <- input$rtm_Prot; row$CBC <- input$rtm_CBC
      }
    } else if (lm %in% c("Fluspect-B", "Fluspect-B-Cx")) {
      row$N <- input$rtm_N; row$Cab <- input$rtm_Cab; row$Car <- input$rtm_Car; row$Anth <- input$rtm_Anth
      row$EWT <- input$rtm_EWT; row$LMA <- input$rtm_LMA; row$Cs <- input$rtm_Cs; row$fqe <- input$rtm_fqe
      row$Cx <- if (lm == "Fluspect-B-Cx") input$rtm_Cx else 0
    } else if (lm == "Liberty") {
      row$cell.d <- input$`rtm_cell.d`; row$inter.c <- input$`rtm_inter.c`
      row$baseline.abs <- input$`rtm_baseline.abs`; row$leaf.thick <- input$`rtm_leaf.thick`
      row$albino.abs <- input$`rtm_albino.abs`; row$Cab <- input$rtm_Cab; row$EWT <- input$rtm_EWT
      row$lign.cell <- input$`rtm_lign.cell`; row$Nitrogen <- input$rtm_Nitrogen
    }

    if (input$rtm_canopy_model == "foursail2") {
      req(input$rtm_fraction_brown)
      row$fraction_brown <- input$rtm_fraction_brown; row$diss <- input$rtm_diss
      row$Cv <- input$rtm_Cv; row$Zeta <- input$rtm_Zeta
    } else if (input$rtm_canopy_model == "inform") {
      req(input$rtm_LAIu)
      row$LAIu <- input$rtm_LAIu; row$cd <- input$rtm_cd; row$h <- input$rtm_h
      row$sd <- input$rtm_sd; row$skyl <- input$rtm_skyl
    }

    row
  })

  rtm_rsoil <- reactive({
    input$rtm_psoil * ToolsRTM::dataSpec_PDB[, 11] + (1 - input$rtm_psoil) * ToolsRTM::dataSpec_PDB[, 12]
  })

  rtm_result <- reactive({
    lm <- input$rtm_leaf_model
    cm <- input$rtm_canopy_model
    out <- run_canopy_model(rtm_row(), cm, lm, rtm_rsoil())
    showNotification(paste0(names(.canopy_model_choices)[.canopy_model_choices == cm], " + ", lm, " simulation done."),
                      type = "message", duration = 3)
    out
  })

  output$rtm_plot <- renderPlot({
    df.plot <- rtm_result()
    cm_label <- names(.canopy_model_choices)[.canopy_model_choices == input$rtm_canopy_model]
    lm_label <- input$rtm_leaf_model

    p <- ggplot(df.plot, aes(x = wavelength, y = reflectance)) +
      geom_line(aes(color = "Canopy"), linewidth = 1.4) +
      labs(x = "Wavelength (nm)", y = "Reflectance", color = "",
           title = paste0(cm_label, " + ", lm_label))
    y_max <- max(0.5, max(df.plot$reflectance, na.rm = TRUE) * 1.05)

    # NB (fix): checkbox existed in the UI but was never read anywhere in
    # the server -- "Include soil reflectance" silently did nothing.
    if (isTRUE(input$rtm_include_soil)) {
      rsoil <- rtm_rsoil()
      df.soil <- data.frame(wavelength = df.plot$wavelength, reflectance = rsoil[seq_along(df.plot$wavelength)])
      p <- p + geom_line(data = df.soil, aes(x = wavelength, y = reflectance, color = "Soil"), linewidth = 1.2)
      y_max <- max(y_max, max(df.soil$reflectance) * 1.05)
      showNotification("Soil reflectance spectrum added.", type = "message", duration = 3)
    }

    p + scale_color_manual(values = c("Canopy" = "#2166ac", "Soil" = "#8c510a")) +
      coord_cartesian(ylim = c(0, y_max)) +
      theme_prosail(legend_position = if (isTRUE(input$rtm_include_soil)) "top" else "none")
  }, res = render_plot_res)

  # Leaf-level SIF panel -- only shown (conditionalPanel, UI side) for the two
  # Fluspect leaf models, since MbI/MbII (the fluorescence matrices this
  # needs) only exist for those. Two stacked panels: reflectance vs.
  # apparent (SIF-added) reflectance on top, the raw SIF emission spectrum
  # below -- see compute_leaf_sif()'s docstring for the physics.
  output$rtm_sif_plot <- renderPlot({
    req(input$rtm_leaf_model %in% c("Fluspect-B", "Fluspect-B-Cx"))
    sif_df <- compute_leaf_sif(rtm_row(), input$rtm_leaf_model)
    showNotification("SIF simulation done.", type = "message", duration = 3)
    # Smaller base_size than the default theme_prosail() -- these two panels
    # are stacked (patchwork), so each only gets half the container height.

    p_refl <- ggplot(sif_df, aes(x = wavelength)) +
      geom_line(aes(y = reflectance, color = "Reflectance"), linewidth = 1.2) +
      geom_line(aes(y = apparent_reflectance, color = "Apparent reflectance (+SIF)"), linewidth = 1.2) +
      scale_color_manual(values = c("Reflectance" = "#1b7837", "Apparent reflectance (+SIF)" = "#b2182b")) +
      labs(x = "Wavelength (nm)", y = "Reflectance", color = "",
           title = "Leaf reflectance vs. apparent (SIF-added) reflectance") +
      theme_prosail(base_size = 13)

    p_sif <- ggplot(sif_df, aes(x = wavelength, y = sif)) +
      geom_line(color = "#b2182b", linewidth = 1.2) +
      labs(x = "Wavelength (nm)", y = "SIF (W m-2 nm-1)", title = "Simulated leaf-level SIF emission") +
      theme_prosail(base_size = 13, legend_position = "none")

    p_refl / p_sif
  }, res = render_plot_res)

  # Sensitivity sub-tab: pick ONE trait, sweep it across its slider's own
  # range (6 values) while holding every other parameter at its CURRENT
  # slider value (rtm_row()), re-running the same canopy dispatch each time
  # -- same one-trait-at-a-time approach as the tutorial's own step 4 and
  # Exercise-1.Rmd's Part 2, but interactive and against whatever
  # canopy/leaf model + parameter values are currently selected.
  output$rtm_sens_trait_ui <- renderUI({
    traits <- rtm_sweep_trait_choices(input$rtm_canopy_model, input$rtm_leaf_model)
    req(traits)
    groups <- rtm_sweep_trait_groups(input$rtm_canopy_model, input$rtm_leaf_model)
    # selectInput() renders a named list of vectors as <optgroup> blocks;
    # each vector must itself be named choice-label -> choice-value.
    grouped_choices <- lapply(groups, function(nms) setNames(nms, nms))
    choices <- names(traits)
    tagList(
      fluidRow(
        column(width = 7, selectInput("rtm_sens_trait", "Parameter to vary", choices = grouped_choices,
                                       selected = if ("Cab" %in% choices) "Cab" else choices[1])),
        column(width = 5, selectInput("rtm_sens_palette", "Color palette",
                                       choices = c("Viridis" = "viridis", "Plasma" = "plasma",
                                                    "Magma" = "magma", "Cividis" = "cividis",
                                                    "Inferno" = "inferno"),
                                       selected = "viridis"))
      )
    )
  })

  rtm_sensitivity_result <- reactive({
    trait <- input$rtm_sens_trait
    lm <- input$rtm_leaf_model
    cm <- input$rtm_canopy_model
    traits <- rtm_sweep_trait_choices(cm, lm)
    req(trait, traits, trait %in% names(traits))

    rng <- traits[[trait]]
    values <- seq(rng[1], rng[2], length.out = 6)
    base_row <- rtm_row()
    rsoil <- rtm_rsoil()

    # psoil isn't a row column at all -- rtm_rsoil() derives the actual
    # soil spectrum from input$rtm_psoil separately, so sweeping it means
    # varying the rsoil ARGUMENT passed to run_canopy_model(), not row$psoil
    # (which doesn't exist).
    rsoil_for <- function(v) {
      if (trait == "psoil") v * ToolsRTM::dataSpec_PDB[, 11] + (1 - v) * ToolsRTM::dataSpec_PDB[, 12] else rsoil
    }
    # LIDFa is also special: switches TypeLidf to 2 (ellipsoidal) and pins
    # LIDFb to 0 for the sweep, overriding whatever the "Type LIDF" dropdown
    # currently has (its own TypeLidf=1 preset-shape system) -- see
    # .rtm_geometry_traits's own comment for why.
    row_for <- function(v) {
      row_v <- base_row
      if (trait == "LIDFa") {
        row_v$TypeLidf <- 2; row_v$LIDFb <- 0; row_v$LIDFa <- v
      } else if (trait != "psoil") {
        row_v[[trait]] <- v
      }
      row_v
    }
    refl_list <- lapply(values, function(v) {
      df <- run_canopy_model(row_for(v), cm, lm, rsoil_for(v))
      df$value <- v
      df
    })
    refl_df <- do.call(rbind, refl_list)

    # When sweeping LAI itself, also carry the (LAI-independent) soil
    # reflectance baseline, so low-LAI curves visibly sit close to it
    # (background shows through) and high-LAI curves pull away (canopy
    # dominates). Not shown for the psoil sweep itself -- the swept
    # reflectance curves there ARE the soil-brightness comparison, so a
    # separate dashed reference line was redundant clutter, not useful
    # information (user-caught after trying it).
    soil_df <- NULL
    if (trait == "LAI") {
      wl <- unique(refl_df$wavelength)
      soil_df <- data.frame(wavelength = wl, reflectance = rsoil[seq_along(wl)], group = "soil")
    }

    # SIF/apparent-reflectance only vary with LEAF-level traits
    # (compute_leaf_sif() reads N/Cab/Car/.../fqe/Cx from the row, nothing
    # canopy- or soil-level like LAI/hspot/tts/tto/psoil/fraction_brown/
    # LAIu/etc.) -- sweeping any of those would show a flat, unchanging SIF
    # panel for the same reason psi looked flat in the reflectance panel
    # (varying something the output doesn't actually depend on). Skipped
    # entirely rather than shown-but-misleading.
    sif_df <- NULL
    if (lm %in% c("Fluspect-B", "Fluspect-B-Cx") && trait %in% names(.leaf_model_params[[lm]])) {
      sif_list <- lapply(values, function(v) {
        row_v <- base_row
        row_v[[trait]] <- v
        df <- compute_leaf_sif(row_v, lm)
        df$value <- v
        df
      })
      sif_df <- do.call(rbind, sif_list)
    }

    showNotification(paste0("Sensitivity sweep for ", trait, " done."), type = "message", duration = 3)
    list(refl_df = refl_df, soil_df = soil_df, sif_df = sif_df, trait = trait)
  })

  output$rtm_sensitivity_plot <- renderPlot({
    res <- rtm_sensitivity_result()
    trait_label <- res$trait
    palette <- input$rtm_sens_palette
    if (is.null(palette)) palette <- "viridis"
    # direction = -1: low trait values get the palette's "light/yellow" end,
    # high values get its "dark/blue" end (e.g. high Cab -> blue) -- the
    # opposite of viridis's own default (dark = low, yellow = high).
    has_fluor <- !is.null(res$sif_df)

    p_refl <- ggplot(res$refl_df, aes(x = wavelength, y = reflectance, color = value, group = value)) +
      geom_line(linewidth = 1) +
      scale_color_viridis_c(name = trait_label, option = palette, direction = -1) +
      labs(x = "Wavelength (nm)", y = "Reflectance",
           title = paste0("Canopy reflectance sensitivity to ", trait_label)) +
      theme_prosail(base_size = if (has_fluor) 13 else 20)

    if (!is.null(res$soil_df)) {
      # group: "soil" for LAI's single fixed reference line, "min psoil"/
      # "max psoil" for psoil's two sweep-endpoint lines -- always present,
      # so 2 lines never get zigzag-connected into one.
      p_refl <- p_refl +
        geom_line(data = res$soil_df, aes(x = wavelength, y = reflectance, group = group),
                  inherit.aes = FALSE, color = "#8c510a", linewidth = 1.1, linetype = "dashed")
    }

    if (!has_fluor) {
      p_refl
    } else {
      # fqe (and only fqe) leaves true canopy/leaf reflectance completely
      # unchanged -- it only scales the fluorescence matrices, so sweeping it
      # would otherwise look like "nothing happened" in the reflectance
      # panel above. Always showing apparent (SIF-added) reflectance and raw
      # SIF as their own panels for Fluspect models -- not just when
      # sweeping fqe specifically -- guarantees the sweep's actual effect is
      # visible somewhere, for every trait, not only the ones that happen to
      # move true reflectance.
      p_apparent <- ggplot(res$sif_df, aes(x = wavelength, y = apparent_reflectance, color = value, group = value)) +
        geom_line(linewidth = 1) +
        scale_color_viridis_c(name = trait_label, option = palette, direction = -1) +
        labs(x = "Wavelength (nm)", y = "Apparent reflectance",
             title = paste0("Leaf apparent (SIF-added) reflectance sensitivity to ", trait_label)) +
        theme_prosail(base_size = 13)

      p_sif <- ggplot(res$sif_df, aes(x = wavelength, y = sif, color = value, group = value)) +
        geom_line(linewidth = 1) +
        scale_color_viridis_c(name = trait_label, option = palette, direction = -1) +
        labs(x = "Wavelength (nm)", y = "SIF (W m-2 nm-1)",
             title = paste0("Leaf-level SIF sensitivity to ", trait_label)) +
        theme_prosail(base_size = 13)
      p_refl / p_apparent / p_sif
    }
  }, res = render_plot_res)

  # Fixed-values table + live-value R/Python code for the current sweep --
  # same fixed_names logic as sensitivity_code_r()/_python() themselves (see
  # their own docstring), so the table and the code can never show different
  # numbers for the same sweep.
  output$rtm_sens_fixed_table <- renderTable({
    res <- rtm_sensitivity_result()
    trait <- res$trait
    lm <- input$rtm_leaf_model
    cm <- input$rtm_canopy_model
    row <- rtm_row()
    traits <- rtm_sweep_trait_choices(cm, lm)
    rng <- traits[[trait]]

    # LIDFa sweep forces the ellipsoidal (Campbell) leaf-angle
    # parameterization -- TypeLidf<-2, LIDFb<-0 -- for every swept row (see
    # rtm_sensitivity_result()'s row_for()), overriding whatever TypeLidf/
    # LIDFb the sliders currently hold. Showing those as "fixed at their
    # slider value" would be wrong, so they're excluded here and replaced
    # with a note of the values actually used during the sweep.
    excl <- if (trait == "LIDFa") c(trait, "LIDFb", "TypeLidf") else trait
    fixed_names <- setdiff(c("LAI", "hspot", "tts", "tto", "psi", "LIDFa", "LIDFb", "TypeLidf",
                              names(.leaf_model_params[[lm]]), names(.canopy_model_params[[cm]])), excl)
    fixed_vals <- vapply(fixed_names, function(n) unname(row[[n]]), numeric(1))
    df <- data.frame(Parameter = fixed_names, Value = fixed_vals)
    if (trait != "psoil") {
      df <- rbind(df, data.frame(Parameter = "psoil (soil brightness)", Value = input$rtm_psoil))
    }
    if (trait == "LIDFa") {
      df <- rbind(df, data.frame(
        Parameter = c("TypeLidf (forced for this sweep)", "LIDFb (forced for this sweep)"),
        Value = c("2 (ellipsoidal / Campbell)", "0")
      ))
    }
    rbind(
      data.frame(Parameter = paste0(trait, " (SWEPT, not fixed)"),
                 Value = sprintf("%s values from %s to %s", 6, rng[1], rng[2])),
      df
    )
  })

  output$rtm_sens_code_r <- renderUI({
    res <- rtm_sensitivity_result()
    trait <- res$trait
    lm <- input$rtm_leaf_model
    cm <- input$rtm_canopy_model
    traits <- rtm_sweep_trait_choices(cm, lm)
    code_block(sensitivity_code_r(cm, lm, trait, rtm_row(), input$rtm_psoil, traits[[trait]]))
  })
  output$rtm_sens_code_py <- renderUI({
    res <- rtm_sensitivity_result()
    trait <- res$trait
    lm <- input$rtm_leaf_model
    cm <- input$rtm_canopy_model
    traits <- rtm_sweep_trait_choices(cm, lm)
    code_block(sensitivity_code_python(cm, lm, trait, rtm_row(), input$rtm_psoil, traits[[trait]]))
  })

  output$downloadRtmData <- downloadHandler(
    filename = function() paste0("rtm_", input$rtm_canopy_model, "_", input$rtm_leaf_model, "_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(rtm_result(), file, row.names = FALSE)
      showNotification("Downloaded RTM data.", type = "message", duration = 3)
    }
  )

  # Absorption coefficients: genuinely different data per leaf model
  # (.abs_coef_source()'s own docstring has the full detail), not the same
  # PROSPECT-PRO table shown regardless of what's selected. Liberty has no
  # equivalent data at all -- shown as an explicit message, not a plot.
  .abs_coef_colors <- c(
    "Chlorophyll content" = "darkolivegreen", "Carotenoids content" = "darkorange3",
    "Anthocyanin content" = "firebrick3", "Brown pigments" = "gold4", "Senescence (Cs)" = "gold4",
    "Water content" = "dodgerblue4", "Dry matter" = "khaki3", "Dry matter (LMA)" = "khaki3",
    "Protein content" = "orchid4", "Carbon-based compounds" = "khaki3",
    "Violaxanthin state (Cx=0)" = "darkolivegreen3", "Zeaxanthin state (Cx=1)" = "darkorange3",
    "Current Kca (interpolated)" = "black"
  )
  abs_coef_panel <- function(wave, panel, wave_range = NULL) {
    # Built directly as long-format (one geom_line() per trait, i.e. this IS
    # already "one line per trait" -- as.data.frame(panel) + melt() used to
    # sit in between, but as.data.frame() silently mangles list names with
    # spaces ("Chlorophyll content" -> "Chlorophyll.content"), which then
    # failed to match names(panel) in the factor()/scale_color_manual()
    # calls right after -- every trait's label came back NA, breaking the
    # legend and color mapping (not literally "one merged line", but visibly
    # broken all the same). Skipping the wide/mangle-prone intermediate
    # step entirely avoids the whole class of bug.
    df_long <- do.call(rbind, lapply(names(panel), function(nm) {
      data.frame(wave = wave, SpecificAbsorption = panel[[nm]], Trait = nm)
    }))
    df_long$Trait <- factor(df_long$Trait, levels = names(panel))
    if (!is.null(wave_range)) df_long <- df_long[df_long$wave >= wave_range[1] & df_long$wave <= wave_range[2], ]
    ggplot(df_long, aes(x = wave, y = SpecificAbsorption, color = Trait)) +
      geom_line(linewidth = 1) +
      labs(x = "Wavelength (nm)", y = "Specific Absorption") +
      theme_prosail(legend_position = "right") +
      theme(legend.title = element_blank()) +
      scale_color_manual(values = .abs_coef_colors[names(panel)])
  }

  # Live Cx (only meaningful/settable for Fluspect-B-Cx; Fluspect-B itself
  # always simulates with Cx fixed at 0, per rtm_row()'s own dispatch).
  rtm_current_cx <- reactive({
    if (isTRUE(input$rtm_leaf_model == "Fluspect-B-Cx") && !is.null(input$rtm_Cx)) input$rtm_Cx else 0
  })

  output$abs_plot1 <- renderPlot({
    src <- .abs_coef_source(input$rtm_leaf_model, rtm_current_cx())
    validate(need(!is.null(src),
                   "Liberty has no equivalent specific-absorption-coefficient data -- its absorption model is parameterized by cell geometry and intercellular air space, not a per-unit-content spectrum like PROSPECT/Fluspect."))
    showNotification("Absorption spectrum done successfully.", type = "message", duration = 3)
    abs_coef_panel(src$wave, src$panel1, wave_range = c(400, 900))
  }, res = render_plot_res)

  output$abs_plot2 <- renderPlot({
    src <- .abs_coef_source(input$rtm_leaf_model, rtm_current_cx())
    validate(need(!is.null(src), ""))
    abs_coef_panel(src$wave, src$panel2)
  }, res = render_plot_res)

  output$abs_plot3 <- renderPlot({
    src <- .abs_coef_source(input$rtm_leaf_model, rtm_current_cx())
    req(!is.null(src), !is.null(src$panel3))
    abs_coef_panel(src$wave, src$panel3)
  }, res = render_plot_res)

  output$downloadSatData <- downloadHandler(
    filename = function() paste0("prosail_", input$sat_sensor, "_bands_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(sat_data(), file, row.names = FALSE)
      showNotification("Downloaded sensor-band data.", type = "message", duration = 3)
    }
  )

  # "How in R"/"How in Python" tutorial tabs: step 2's code adapts to
  # whichever canopy/leaf model combination is picked in that tab's own
  # dropdowns (one_sim_code_r()/one_sim_code_python(), verified against all
  # 15 combinations in both languages -- see generate_tutorial_figures.R's
  # sibling verification scripts).
  output$tut_r_one_sim_code <- renderUI({
    code_block(one_sim_code_r(input$tut_r_canopy_model, input$tut_r_leaf_model))
  })
  output$tut_py_one_sim_code <- renderUI({
    code_block(one_sim_code_python(input$tut_py_canopy_model, input$tut_py_leaf_model))
  })
  # Step 3 -- same dropdowns as step 2, same idea (many_sims_code_r()/
  # _python(), sampling ranges shared with the Model Explorer tab's own
  # Sensitivity sub-tab via .rtm_sweep_traits).
  output$tut_r_many_sims_code <- renderUI({
    code_block(many_sims_code_r(input$tut_r_canopy_model, input$tut_r_leaf_model))
  })
  output$tut_py_many_sims_code <- renderUI({
    code_block(many_sims_code_python(input$tut_py_canopy_model, input$tut_py_leaf_model))
  })

  # Fluspect-only SIF/apparent-reflectance note, shown right after step 2
  # when that tab's leaf-model dropdown is Fluspect-B/-B-Cx -- empty
  # (NULL) otherwise.
  tutorial_fluspect_note <- function(leaf_model, note_text) {
    if (isTRUE(leaf_model %in% c("Fluspect-B", "Fluspect-B-Cx"))) {
      tagList(
        p(strong("Using Fluspect -- SIF and apparent reflectance available too: "),
          "Fluspect models (unlike PROSPECT/Liberty) also simulate chlorophyll",
          " fluorescence (SIF). This isn't part of the reflectance computed above,",
          " but is available from the same leaf model:"),
        code_block(note_text)
      )
    } else {
      NULL
    }
  }
  output$tut_r_fluspect_sif_note <- renderUI({
    tutorial_fluspect_note(input$tut_r_leaf_model, .tutorial_sif_note_r)
  })
  output$tut_py_fluspect_sif_note <- renderUI({
    tutorial_fluspect_note(input$tut_py_leaf_model, .tutorial_sif_note_python)
  })

  # Reference tab -- all three tables are static (built once at app-load
  # time from the same constants/helpers the sliders and sensor plots
  # themselves use), so no reactivity needed here at all.
  output$reference_leaf_table <- renderTable({
    req(input$reference_leaf_model)
    subset(.rtm_reference_leaf_params, Model == input$reference_leaf_model)
  })

  output$reference_lma_prot_cbc <- renderTable({
    .rtm_reference_lma_prot_cbc
  })

  # Canopy models: the shared "All" rows apply no matter which model is
  # picked, plus that model's own "(extra)" rows if it has any (fourSAIL
  # itself doesn't -- .canopy_extra_labels[["foursail"]] matches nothing,
  # so only the shared rows show, which is correct).
  .canopy_extra_labels <- c("foursail" = "fourSAIL (extra)", "foursail2" = "fourSAIL2 (extra)",
                             "inform" = "INFORM (extra)")
  output$reference_canopy_table <- renderTable({
    req(input$reference_canopy_model)
    extra_label <- .canopy_extra_labels[[input$reference_canopy_model]]
    tbl <- subset(.rtm_reference_canopy_params,
                   Model == "All (fourSAIL/fourSAIL2/INFORM)" | Model == extra_label)
    # Model column originally distinguishes "shared by every canopy model"
    # rows from model-specific "(extra)" rows -- once the table is already
    # filtered down to the ONE selected model, that distinction is no
    # longer useful information, just a confusing "All (fourSAIL/...)"
    # label; overwrite it with the actually-selected model's own name.
    model_name <- names(.canopy_model_choices)[.canopy_model_choices == input$reference_canopy_model]
    tbl$Model <- model_name
    tbl
  })

  output$reference_lidf_presets <- renderTable({
    .rtm_reference_lidf_presets
  })

  output$reference_lidf_plot <- renderPlot({
    ggplot(.rtm_lidf_distribution, aes(x = Angle, y = Frequency, width = BinWidth)) +
      geom_col(fill = "#2166ac", position = "identity") +
      facet_wrap(~Shape, ncol = 3) +
      labs(x = "Leaf inclination angle (deg)", y = "Fraction of leaf area") +
      theme_prosail(legend_position = "none")
  }, res = render_plot_res)

  output$reference_sensor_table <- renderTable({
    req(input$reference_sensor)
    subset(.rtm_reference_sensor_bands, SensorKey == input$reference_sensor, select = -SensorKey)
  })
}

# 5. Create the Shiny app object -----

shinyApp(ui = ui, server = server)
