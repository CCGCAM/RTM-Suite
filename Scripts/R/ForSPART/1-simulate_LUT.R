# ==============================================================================
# ForSPART course pipeline, step 1/3: simulate a LUT with ToolsRTM::SPART()
# (Soil-Plant-Atmosphere Radiative Transfer -- couples fourSAIL canopy +
# BSM/MARMIT soil + SMAC atmosphere to give TOC *and* TOA reflectance at a
# real sensor's own bands), compute indices, save diagnostic figures.
#
# Unlike simulate_RTM() (fourSAIL/foursail2/INFORM), SPART() outputs directly
# at a chosen sensor's bands -- there's no separate "native then convolve"
# step, so "convolution" here means running SPART() once per sensor instead.
# Always 100 simulations. Feeds 2-inversion_ML.R and 3-inversion_DL.R next.
# ==============================================================================

rm(list = ls())

if (!requireNamespace("ToolsRTM", quietly = TRUE)) {
  remotes::install_gitlab("caminoccg/toolsrtm")
}
library(ToolsRTM)
library(ggplot2)

## ----------------------------------------------------------------------------
## User-configurable settings
## ----------------------------------------------------------------------------

n.samples  <- 100                 # always 100 for this course pipeline
leaf.model <- "PROSPECT-PRO"
soil.source <- "marmit"           # "bsm" (SPART's own built-in soil model) or "marmit"
seed <- 1
sensors <- list(Sentinel2A = ToolsRTM::Sentinel2A.MSI, MODIS = ToolsRTM::TerraAqua.MODIS)

correlate.traits <- list(Car = list(with = "Cab", scale = 1/4, r = 0.8))

out_dir <- "outs/ForSPART"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
pdf(file.path(out_dir, "Rplots.pdf"))  # captures any stray plot() call, see ForPROSAIL/3-simulate_LUT.R for why

cat(sprintf("=== Simulating %d SPART spectra (%s + fourSAIL, soil = %s) ===\n",
            n.samples, leaf.model, soil.source))

## ----------------------------------------------------------------------------
## 1. Build the LUT
## ----------------------------------------------------------------------------

LUT <- as.data.frame(getLUT(inputs = ToolsRTM::inputsSPART, nLUT = n.samples, setseed = seed))

for (trait in names(correlate.traits)) {
  spec <- correlate.traits[[trait]]
  scale_ <- if (is.null(spec$scale)) 1 else spec$scale
  LUT[[trait]] <- correlatedValue(x = LUT[[spec$with]] * scale_, r = spec$r)
  cat(sprintf("Correlation: %s ~ %s (target r=%.2f, realized r=%.2f)\n",
              trait, spec$with, spec$r, cor(LUT[[trait]], LUT[[spec$with]])))
}

# inputsSPART only ships PROSPECT-PRO/-D columns -- add what other leaf
# models need too (same requirement as foursail()/SPART() itself).
LUT$Cs <- 0; LUT$fqe <- 0.01; LUT$Cx <- 0
LUT$cell.d <- 40; LUT$inter.c <- 0.045; LUT$baseline.abs <- 0.0006
LUT$leaf.thick <- 1.6; LUT$albino.abs <- 0; LUT$lign.cell <- 2; LUT$Nitrogen <- 1

key_traits <- intersect(c("Cab", "Car", "Anth", "LAI", "EWT", "LMA", "LIDFa"), names(LUT))

hist_data <- stack(LUT[, key_traits])
colnames(hist_data) <- c("value", "trait")
hist_plot <- ggplot(hist_data, aes(x = value)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  facet_wrap(~trait, scales = "free") +
  theme_bw() +
  labs(x = NULL, y = "count", title = sprintf("Sampled trait distributions (n=%d)", n.samples)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
print(hist_plot)
ggsave(file.path(out_dir, "0-trait_histograms.png"), plot = hist_plot, width = 18, height = 12, dpi = 300, units = "cm")

corr_mat <- cor(LUT[, key_traits])
png(file.path(out_dir, "0-trait_correlations.png"), width = 1400, height = 1400, res = 200)
corrplot::corrplot(corr_mat, method = "color", type = "upper",
                    addCoef.col = "black", tl.col = "black", tl.srt = 45, number.cex = 0.8,
                    title = "Sampled trait correlations", mar = c(0, 0, 2, 0))
dev.off()
corrplot::corrplot(corr_mat, method = "color", type = "upper",
                    addCoef.col = "black", tl.col = "black", tl.srt = 45, number.cex = 0.8,
                    title = "Sampled trait correlations", mar = c(0, 0, 2, 0))

## ----------------------------------------------------------------------------
## 2. Soil spectrum (SPART wants 400-2400nm, 1nm step, length 2001)
## ----------------------------------------------------------------------------

wl_soil <- 400:2400
if (soil.source == "marmit") {
  soil <- get.marmit.rsoil(database = "Bablet_2016", id = 1, L = 0.05, eps = 0.3, wl.out = wl_soil)
  rsoil <- soil$rsoil.wet
  cat(sprintf("Soil: MARMIT (Bablet_2016 id=1, SMC = %.1f%%)\n", soil$SMC))
} else {
  rsoil <- NULL  # SPART's own built-in BSM soil model
  cat("Soil: SPART's built-in BSM model\n")
}

## ----------------------------------------------------------------------------
## 3. Run SPART for every LUT row, once per sensor
## ----------------------------------------------------------------------------

datasets <- list()
for (sensor_name in names(sensors)) {
  cat(sprintf("Running SPART for sensor = %s ... ", sensor_name))
  sensor_obj <- sensors[[sensor_name]]

  sim_list <- lapply(seq_len(n.samples), function(i) {
    tryCatch(
      suppressMessages(suppressWarnings(SPART(
        inputLUT = LUT[i, ], CanopyModel = "fourSAIL", LeafModel = leaf.model,
        sensor.i = sensor_obj, rsoil = rsoil, get.plots = FALSE
      ))),
      error = function(e) NULL
    )
  })
  ok <- !sapply(sim_list, is.null)
  cat(sprintf("%d/%d rows OK\n", sum(ok), n.samples))

  toc <- t(sapply(sim_list[ok], function(s) s$output$rfl.toc.BRDF))
  toa <- t(sapply(sim_list[ok], function(s) s$output$rfl.toa))
  band_wl <- sim_list[[which(ok)[1]]]$output$wave
  colnames(toc) <- paste0("R.", round(band_wl))
  colnames(toa) <- paste0("R.", round(band_wl))

  indices_toc <- tryCatch(
    suppressMessages(getIndices(as.data.frame(toc), pattern.rfl = "R.", spectral.domain = "VNIR-SWIR")),
    error = function(e) { cat("  indices failed:", conditionMessage(e), "\n"); as.data.frame(toc) }
  )
  indices_toc <- indices_toc[, colSums(is.na(indices_toc)) == 0, drop = FALSE]

  dataset <- cbind(LUT[ok, ], indices_toc)
  dataset <- dataset[, !duplicated(names(dataset))]

  datasets[[sensor_name]] <- list(toc = toc, toa = toa, wl = band_wl, dataset = dataset)
}

## ----------------------------------------------------------------------------
## 4. Diagnostic plot: TOC vs TOA reflectance (Sentinel-2A, by Cab bin)
## ----------------------------------------------------------------------------

ref <- datasets[["Sentinel2A"]]
LUT_ok <- LUT[seq_len(nrow(ref$toc)), ]  # rows are already filtered to the OK subset for this sensor above

plot.compare <- rbind(
  data.frame(wave = rep(ref$wl, each = nrow(ref$toc)), RFL = as.vector(ref$toc), Type = "TOC (canopy)"),
  data.frame(wave = rep(ref$wl, each = nrow(ref$toa)), RFL = as.vector(ref$toa), Type = "TOA (satellite)")
)
toc_toa_plot <- ggplot(plot.compare, aes(x = wave, y = RFL, color = Type)) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
  stat_summary(fun = mean, geom = "point", size = 2) +
  theme_bw() +
  labs(x = "Wavelength (nm)", y = "Reflectance",
       title = sprintf("SPART: mean TOC vs TOA reflectance, Sentinel-2A bands (n=%d)", nrow(ref$toc))) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
print(toc_toa_plot)
ggsave(file.path(out_dir, "1-TOC_vs_TOA_reflectance.png"), plot = toc_toa_plot, width = 14, height = 9, dpi = 300, units = "cm")

## ----------------------------------------------------------------------------
## 5. Diagnostic plot: example TOC spectra at Cab's 10th/50th/90th percentile
## ----------------------------------------------------------------------------

pctl_probs <- c(0.1, 0.5, 0.9)
pctl_vals <- quantile(LUT_ok$Cab, probs = pctl_probs)
pctl_rows <- sapply(pctl_vals, function(v) which.min(abs(LUT_ok$Cab - v)))

pctl_long <- do.call(rbind, lapply(seq_along(pctl_rows), function(k) {
  data.frame(wave = ref$wl, RFL = as.numeric(ref$toc[pctl_rows[k], ]),
             Percentile = sprintf("P%d (Cab=%.1f)", pctl_probs[k] * 100, LUT_ok$Cab[pctl_rows[k]]))
}))
pctl_plot <- ggplot(pctl_long, aes(x = wave, y = RFL, color = Percentile)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  theme_bw() +
  labs(x = "Wavelength (nm)", y = "Reflectance (TOC, Sentinel-2A bands)",
       title = "Example SPART TOC spectra at Cab's 10th/50th/90th percentile") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
print(pctl_plot)
ggsave(file.path(out_dir, "2-example_spectra_by_percentile.png"), plot = pctl_plot, width = 14, height = 9, dpi = 300, units = "cm")

## ----------------------------------------------------------------------------
## 6. Save
## ----------------------------------------------------------------------------

saveRDS(list(LUT = LUT, datasets = list(native = datasets$Sentinel2A$dataset, se2a = datasets$Sentinel2A$dataset,
                                         modis = datasets$MODIS$dataset),
             leaf.model = leaf.model, soil.source = soil.source, n.samples = n.samples, seed = seed),
        file.path(out_dir, "1-datasets.rds"))

dev.off()

cat("\nSaved to '", out_dir, "/1-datasets.rds' -- read back with readRDS() in 2-inversion_ML.R / ",
    "3-inversion_DL.R ($datasets$native / $se2a / $modis are ready for inversion).\n", sep = "")
