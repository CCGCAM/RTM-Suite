# ==============================================================================
# ForPROSAIL course pipeline, step 1/3: simulate a LUT (fourSAIL canopy model),
# convolve to Sentinel-2A/2B/PRISMA, compute indices, save a full set of
# diagnostic figures (trait histograms, trait correlations, spectra by trait
# percentile, sensor-convolution comparison).
#
# This is independent of 1-GetSimulationsLUTs.R/1-ReduceLUTs.R/2-Model*.R above
# (those calibrate PROSAIL against Carlos' real field spectrometer data and
# need hsdar-replacement machinery + private field data -- see this folder's
# data/README.md). This script needs neither: it's the self-contained
# simulate -> convolve -> index course pipeline, always 100 simulations,
# feeding 4-inversion_ML.R and 5-inversion_DL.R next.
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

n.samples    <- 100              # always 100 for this course pipeline
leaf.model   <- "PROSPECT-D"
canopy.model <- "fourSAIL"
soil.source  <- "marmit"         # "flat" or "marmit" (ToolsRTM::get.marmit.rsoil())
seed         <- 1
sensors      <- c("Sentinel2a", "Sentinel2b", "PRISMA")

# Make Car co-vary with Cab (real leaf pigments aren't independent) via
# ToolsRTM::correlatedValue() -- same mechanism as Scripts/Pipeline/1-simulate_LUT.R.
correlate.traits <- list(Car = list(with = "Cab", scale = 1/4, r = 0.8))

# Everything lands directly under outs/ForPROSAIL/ -- no extra subfolder, so
# it's obvious where to look.
out_dir <- "outs/ForPROSAIL"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Any plot() / print(ggplot) call that doesn't have its own explicit device
# open (e.g. a stray interactive-viewing call) would otherwise fall back to a
# default 'Rplots.pdf' in the *working* directory (Scripts/ForPROSAIL/) --
# opening this device up front means that fallback lands in outs/ instead.
pdf(file.path(out_dir, "Rplots.pdf"))

cat(sprintf("=== Simulating %d spectra (%s + %s, soil = %s) ===\n",
            n.samples, leaf.model, canopy.model, soil.source))

## ----------------------------------------------------------------------------
## 1. Build the LUT (trait values)
## ----------------------------------------------------------------------------

inputs <- ToolsRTM::inputsPROSAIL
LUT <- as.data.frame(getLUT(inputs = inputs, nLUT = n.samples, setseed = seed))

for (trait in names(correlate.traits)) {
  spec <- correlate.traits[[trait]]
  scale_ <- if (is.null(spec$scale)) 1 else spec$scale
  LUT[[trait]] <- correlatedValue(x = LUT[[spec$with]] * scale_, r = spec$r)
  cat(sprintf("Correlation: %s ~ %s (scale=%.3g, target r=%.2f, realized r=%.2f)\n",
              trait, spec$with, scale_, spec$r, cor(LUT[[trait]], LUT[[spec$with]])))
}

# inputsPROSAIL only ships PROSPECT-PRO/-D columns -- add what other leaf
# models / fourSAIL itself need too (same requirement as foursail() itself).
LUT$Cs <- 0; LUT$fqe <- 0.01; LUT$Cx <- 0
LUT$cell.d <- 40; LUT$inter.c <- 0.045; LUT$baseline.abs <- 0.0006
LUT$leaf.thick <- 1.6; LUT$albino.abs <- 0; LUT$lign.cell <- 2; LUT$Nitrogen <- 1

## ----------------------------------------------------------------------------
## 2. Diagnostic figures: trait distributions and correlations
## ----------------------------------------------------------------------------

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
corrplot::corrplot(corr_mat, method = "color", type = "upper",  # also draw into the Rplots.pdf capture device
                    addCoef.col = "black", tl.col = "black", tl.srt = 45, number.cex = 0.8,
                    title = "Sampled trait correlations", mar = c(0, 0, 2, 0))

## ----------------------------------------------------------------------------
## 3. Soil spectrum
## ----------------------------------------------------------------------------

short_domain <- leaf.model %in% c("Fluspect-B", "Fluspect-B-Cx")
wl <- if (short_domain) 400:2400 else 400:2500

if (soil.source == "marmit") {
  soil <- get.marmit.rsoil(database = "Bablet_2016", id = 1, L = 0.05, eps = 0.3, wl.out = wl)
  rsoil <- soil$rsoil.wet
  cat(sprintf("Soil: MARMIT (Bablet_2016 id=1, SMC = %.1f%%)\n", soil$SMC))
} else {
  rsoil <- rep(0.15, length(wl))
  cat("Soil: flat 0.15\n")
}

## ----------------------------------------------------------------------------
## 4. Run fourSAIL for every LUT row
## ----------------------------------------------------------------------------

refl <- t(sapply(seq_len(n.samples), function(i) {
  sim <- suppressMessages(simulate_RTM(inputLUT = LUT[i, ], rsoil = rsoil,
                                        leaf.model = leaf.model, canopy.model = canopy.model))
  sim$rsot
}))
colnames(refl) <- paste0("R.", wl)
refl <- as.data.frame(refl)

cat(sprintf("Simulated: %d x %d reflectance matrix\n", nrow(refl), ncol(refl)))

## ----------------------------------------------------------------------------
## 5. Convolve to sensors
## ----------------------------------------------------------------------------

refl_X <- refl
colnames(refl_X) <- paste0("X", wl)
refl_X <- cbind(id = seq_len(nrow(refl_X)), refl_X)

convolved <- list()
for (s in sensors) {
  cat(sprintf("Convolving to %s ... ", s))
  convolved[[s]] <- suppressMessages(get.spectra.convolved(rfl = refl_X, sensor = s, plot.spectra = FALSE))
  cat(sprintf("%d bands\n", ncol(convolved[[s]]) - 1))
}

## ----------------------------------------------------------------------------
## 6. Vegetation indices: native + per sensor
## ----------------------------------------------------------------------------

cat("Computing indices: native ... ")
indices_native <- suppressMessages(getIndices(refl, pattern.rfl = "R.", spectral.domain = "VNIR"))
indices_native <- indices_native[, colSums(is.na(indices_native)) == 0, drop = FALSE]
cat(sprintf("%d columns (reflectance + indices)\n", ncol(indices_native)))

se2a_bands <- convolved$Sentinel2a
names(se2a_bands) <- c("id", "B1","B2","B3","B4","B5","B6","B7","B8","B8A","B9","B10","B11","B12")
cat("Computing indices: Sentinel-2A ... ")
indices_se2a <- suppressMessages(getIndicesSE2(se2a_bands[, -1], sensor = "Sentinel-2a", df.data = NULL, fast.process = TRUE))
cat(sprintf("%d indices\n", ncol(indices_se2a)))

prisma_bands <- convolved$PRISMA
wl_prisma <- as.numeric(names(prisma_bands)[-1])
names(prisma_bands) <- c("id", paste0("R.", round(wl_prisma)))
cat("Computing indices: PRISMA ... ")
indices_prisma <- tryCatch(
  suppressMessages(getIndices(prisma_bands[, -1], pattern.rfl = "R.", spectral.domain = "VNIR-SWIR")),
  error = function(e) { cat("failed:", conditionMessage(e), " "); NULL }
)
if (!is.null(indices_prisma)) {
  indices_prisma <- indices_prisma[, colSums(is.na(indices_prisma)) == 0, drop = FALSE]
  cat(sprintf("%d columns\n", ncol(indices_prisma)))
}

dataset_native <- cbind(LUT, indices_native)
dataset_native <- dataset_native[, !duplicated(names(dataset_native))]

dataset_se2a <- cbind(LUT, se2a_bands[, -1], indices_se2a)
dataset_se2a <- dataset_se2a[, !duplicated(names(dataset_se2a))]

dataset_prisma <- NULL
if (!is.null(indices_prisma)) {
  dataset_prisma <- cbind(LUT, prisma_bands[, -1], indices_prisma)
  dataset_prisma <- dataset_prisma[, !duplicated(names(dataset_prisma))]
}

## ----------------------------------------------------------------------------
## 7. Diagnostic plot: reflectance by Cab bin (mean spectra)
## ----------------------------------------------------------------------------

refl.cols.plot <- colnames(refl)[seq(1, ncol(refl), by = max(1, ncol(refl) %/% 400))]
plot.rfl <- cbind(LUT[, "Cab", drop = FALSE], stack(refl[, refl.cols.plot]))
colnames(plot.rfl) <- c("Cab", "RFL", "wave")
plot.rfl$wave <- as.numeric(sub("R.", "", plot.rfl$wave, fixed = TRUE))
plot.rfl$GroupID <- cut(plot.rfl$Cab, breaks = c(0, 20, 30, 40, 60, 80, Inf),
                         labels = c('0-20', '20-30', '30-40', '40-60', '60-80', '>80'))
plot.rfl.mean <- aggregate(RFL ~ GroupID + wave, data = plot.rfl, FUN = "mean")

spectral_plot <- ggplot(plot.rfl.mean, aes(x = wave, y = RFL, group = GroupID)) +
  theme_bw() +
  geom_line(aes(color = GroupID), linetype = "dashed", linewidth = 0.7) +
  labs(color = 'Cab (ug/cm2)', x = 'Wavelength (nm)', y = 'Reflectance') +
  ggtitle(sprintf("%s + %s reflectance by Cab (n=%d)", canopy.model, leaf.model, n.samples)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

print(spectral_plot)
ggsave(file.path(out_dir, "1-reflectance_spectra_byCab.png"), plot = spectral_plot,
       width = 14, height = 9, dpi = 300, units = "cm")

## ----------------------------------------------------------------------------
## 8. Diagnostic plot: example spectra at Cab's 10th/50th/90th percentile
## ----------------------------------------------------------------------------
## Individual simulated spectra (not bin-averaged), at the LUT row whose Cab
## is closest to each percentile -- shows what a single low/typical/high-Cab
## canopy actually looks like, complementing the binned-mean plot above.

pctl_probs <- c(0.1, 0.5, 0.9)
pctl_vals <- quantile(LUT$Cab, probs = pctl_probs)
pctl_rows <- sapply(pctl_vals, function(v) which.min(abs(LUT$Cab - v)))

pctl_long <- do.call(rbind, lapply(seq_along(pctl_rows), function(k) {
  data.frame(wave = wl, RFL = as.numeric(refl[pctl_rows[k], ]),
             Percentile = sprintf("P%d (Cab=%.1f)", pctl_probs[k] * 100, LUT$Cab[pctl_rows[k]]))
}))

pctl_plot <- ggplot(pctl_long, aes(x = wave, y = RFL, color = Percentile)) +
  geom_line(linewidth = 0.8) +
  theme_bw() +
  labs(x = "Wavelength (nm)", y = "Reflectance",
       title = "Example spectra at Cab's 10th/50th/90th percentile") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
print(pctl_plot)
ggsave(file.path(out_dir, "2-example_spectra_by_percentile.png"), plot = pctl_plot,
       width = 14, height = 9, dpi = 300, units = "cm")

## ----------------------------------------------------------------------------
## 9. Diagnostic plot: sensor convolution comparison (median-Cab spectrum)
## ----------------------------------------------------------------------------

i_med <- pctl_rows[2]  # the 50th-percentile row picked above
se2a_wl <- c(443,490,560,665,705,740,783,842,865,945,1375,1610,2190)  # Sentinel-2A central wavelengths, same B1..B12 order as se2a_bands
compare_long <- rbind(
  data.frame(wave = wl, RFL = as.numeric(refl[i_med, ]), Sensor = "Native"),
  data.frame(wave = se2a_wl, RFL = as.numeric(se2a_bands[i_med, -1]), Sensor = "Sentinel-2A"),
  data.frame(wave = wl_prisma, RFL = as.numeric(prisma_bands[i_med, -1]), Sensor = "PRISMA")
)

convolution_plot <- ggplot(compare_long, aes(x = wave, y = RFL, color = Sensor)) +
  geom_line(data = subset(compare_long, Sensor == "Native"), linewidth = 0.4, alpha = 0.5) +
  geom_point(data = subset(compare_long, Sensor != "Native"), size = 2) +
  geom_line(data = subset(compare_long, Sensor != "Native"), linewidth = 0.6) +
  theme_bw() +
  labs(x = "Wavelength (nm)", y = "Reflectance",
       title = sprintf("Sensor convolution comparison (Cab=%.1f spectrum)", LUT$Cab[i_med])) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
print(convolution_plot)
ggsave(file.path(out_dir, "3-sensor_convolution_comparison.png"), plot = convolution_plot,
       width = 14, height = 9, dpi = 300, units = "cm")

## ----------------------------------------------------------------------------
## 10. Save
## ----------------------------------------------------------------------------

saveRDS(list(LUT = LUT, refl = refl, wl = wl, convolved = convolved,
             datasets = list(native = dataset_native, se2a = dataset_se2a, prisma = dataset_prisma),
             leaf.model = leaf.model, canopy.model = canopy.model,
             soil.source = soil.source, n.samples = n.samples, seed = seed),
        file.path(out_dir, "1-datasets.rds"))

dev.off()  # close the Rplots.pdf capture device opened at the top

cat("\nSaved to '", out_dir, "/1-datasets.rds' -- read back with readRDS() in 4-inversion_ML.R / ",
    "5-inversion_DL.R ($datasets$native / $se2a / $prisma are ready for inversion).\n", sep = "")
