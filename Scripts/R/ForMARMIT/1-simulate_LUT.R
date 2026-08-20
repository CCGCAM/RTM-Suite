# ==============================================================================
# ForMARMIT course pipeline, step 1/3: simulate soil reflectance spectra with
# ToolsRTM::get.marmit.rsoil() across random surface-water combinations,
# convolve to Sentinel-2A/PRISMA, compute indices, save diagnostic figures.
#
# MARMIT is soil-only (no vegetation) -- there's no Cab/LAI/EWT here. The
# target trait is SMC (gravimetric soil moisture content, %), MARMIT's own
# physical output, driven by two knobs: L (surface water film thickness, cm)
# and eps (fraction of the surface that's wet). Both are randomized here;
# soil `id` (which of Bablet_2016's 17 reference soils) is held fixed so the
# LUT varies wetness on one consistent soil, not soil type + wetness at once.
# Always 100 simulations. Feeds 2-inversion_ML.R and 3-inversion_DL.R next.
# ==============================================================================

rm(list = ls())

use.dev.source <- TRUE
if (use.dev.source) {
  devtools::load_all("../../../ToolsRTM/R")
} else {
  library(ToolsRTM)
}
library(ggplot2)

## ----------------------------------------------------------------------------
## User-configurable settings
## ----------------------------------------------------------------------------

n.samples <- 100          # always 100 for this course pipeline
database  <- "Bablet_2016" # soil database -- see db.root below
soil.id   <- 1            # which reference soil within `database` (Bablet_2016: 1-17)
version   <- "marmit1"    # "marmit1" or "marmit2"
seed      <- 1
sensors   <- c("Sentinel2a", "Sentinel2b", "PRISMA")

# Only Bablet_2016 ships with ToolsRTM itself (keeps the package install
# small). The other 7 official MARMIT databases (Dupiau_2020, Humper_2015,
# Lesaignoux_2008, Liu_2002, Lobell_2002, Marcq_2012, Philpot_2014 -- see
# https://pss-gitlab.math.univ-paris-diderot.fr/marmit/marmit) live in this
# monorepo's own databases/ folder (repo root, ~200MB total). Set db.root
# below and change `database` above to use any of them -- no download/copy
# needed, just point at the folder:
#   database <- "Liu_2002"
#   db.root  <- "../../../databases"
db.root <- NULL          # NULL = only Bablet_2016 (bundled); or "../../../databases"

out_dir <- "../../../outs/ForMARMIT"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
pdf(file.path(out_dir, "Rplots.pdf"))  # captures any stray plot() call, see ForPROSAIL/3-simulate_LUT.R for why

cat(sprintf("=== Simulating %d MARMIT soil spectra (soil id=%d, %s) ===\n", n.samples, soil.id, version))

## ----------------------------------------------------------------------------
## 1. Sample L (water film thickness) and eps (wetness fraction)
## ----------------------------------------------------------------------------

set.seed(seed)
L   <- runif(n.samples, min = 0.001, max = 0.15)   # cm
eps <- runif(n.samples, min = 0, max = 1)           # fraction wet

LUT <- data.frame(ID = seq_len(n.samples), L = L, eps = eps)

key_traits <- c("L", "eps")
hist_data <- stack(LUT[, key_traits])
colnames(hist_data) <- c("value", "trait")
hist_plot <- ggplot(hist_data, aes(x = value)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  facet_wrap(~trait, scales = "free") +
  theme_bw() +
  labs(x = NULL, y = "count", title = sprintf("Sampled MARMIT input distributions (n=%d)", n.samples)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
print(hist_plot)
ggsave(file.path(out_dir, "0-input_histograms.png"), plot = hist_plot, width = 14, height = 8, dpi = 300, units = "cm")

## ----------------------------------------------------------------------------
## 2. Run MARMIT for every (L, eps) combination
## ----------------------------------------------------------------------------

wl <- 400:2400
marmit_runs <- lapply(seq_len(n.samples), function(i) {
  get.marmit.rsoil(database = database, id = soil.id, version = version,
                    L = L[i], eps = eps[i], wl.out = wl, db_root = db.root)
})

refl <- t(sapply(marmit_runs, function(r) r$rsoil.wet))
colnames(refl) <- paste0("R.", wl)
refl <- as.data.frame(refl)
LUT$SMC <- sapply(marmit_runs, function(r) r$SMC)

cat(sprintf("Simulated: %d x %d reflectance matrix. SMC range: %.1f-%.1f%%\n",
            nrow(refl), ncol(refl), min(LUT$SMC), max(LUT$SMC)))

corr_mat <- cor(LUT[, c("L", "eps", "SMC")])
png(file.path(out_dir, "0-input_correlations.png"), width = 1200, height = 1200, res = 200)
corrplot::corrplot(corr_mat, method = "color", type = "upper",
                    addCoef.col = "black", tl.col = "black", tl.srt = 45, number.cex = 0.9,
                    title = "L / eps / SMC correlations", mar = c(0, 0, 2, 0))
dev.off()
corrplot::corrplot(corr_mat, method = "color", type = "upper",
                    addCoef.col = "black", tl.col = "black", tl.srt = 45, number.cex = 0.9,
                    title = "L / eps / SMC correlations", mar = c(0, 0, 2, 0))

## ----------------------------------------------------------------------------
## 3. Convolve to sensors
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
## 4. Indices (native + per sensor) -- note: most vegetation indices (NDVI
## etc.) aren't physically meaningful for bare soil, but they're still just
## spectral-ratio math, so they're computed as-is and used as ML predictors
## alongside the raw bands (some, like brightness/soil-moisture-sensitive SWIR
## ratios, do carry real signal for this use case).
## ----------------------------------------------------------------------------

cat("Computing indices: native ... ")
indices_native <- suppressMessages(getIndices(refl, pattern.rfl = "R.", spectral.domain = "VNIR-SWIR"))
indices_native <- indices_native[, colSums(is.na(indices_native)) == 0, drop = FALSE]
cat(sprintf("%d columns (reflectance + indices)\n", ncol(indices_native)))

se2a_bands <- convolved$Sentinel2a
names(se2a_bands) <- c("id", "B1","B2","B3","B4","B5","B6","B7","B8","B8A","B9","B10","B11","B12")

dataset_native <- cbind(LUT, indices_native)
dataset_native <- dataset_native[, !duplicated(names(dataset_native))]

dataset_se2a <- cbind(LUT, se2a_bands[, -1])
dataset_se2a <- dataset_se2a[, !duplicated(names(dataset_se2a))]

## ----------------------------------------------------------------------------
## 5. Diagnostic plot: reflectance by SMC bin (mean spectra)
## ----------------------------------------------------------------------------

refl.cols.plot <- colnames(refl)[seq(1, ncol(refl), by = max(1, ncol(refl) %/% 400))]
plot.rfl <- cbind(LUT[, "SMC", drop = FALSE], stack(refl[, refl.cols.plot]))
colnames(plot.rfl) <- c("SMC", "RFL", "wave")
plot.rfl$wave <- as.numeric(sub("R.", "", plot.rfl$wave, fixed = TRUE))
plot.rfl$GroupID <- cut(plot.rfl$SMC, breaks = quantile(LUT$SMC, probs = seq(0, 1, 0.2)),
                         include.lowest = TRUE, labels = c('Q1 (driest)', 'Q2', 'Q3', 'Q4', 'Q5 (wettest)'))
plot.rfl.mean <- aggregate(RFL ~ GroupID + wave, data = plot.rfl, FUN = "mean")

spectral_plot <- ggplot(plot.rfl.mean, aes(x = wave, y = RFL, group = GroupID)) +
  theme_bw() +
  geom_line(aes(color = GroupID), linetype = "dashed", linewidth = 0.7) +
  labs(color = 'SMC quintile', x = 'Wavelength (nm)', y = 'Reflectance') +
  ggtitle(sprintf("MARMIT soil reflectance by SMC quintile (n=%d)", n.samples)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
print(spectral_plot)
ggsave(file.path(out_dir, "1-reflectance_spectra_bySMC.png"), plot = spectral_plot,
       width = 14, height = 9, dpi = 300, units = "cm")

## ----------------------------------------------------------------------------
## 6. Diagnostic plot: example spectra at SMC's 10th/50th/90th percentile
## ----------------------------------------------------------------------------

pctl_probs <- c(0.1, 0.5, 0.9)
pctl_vals <- quantile(LUT$SMC, probs = pctl_probs)
pctl_rows <- sapply(pctl_vals, function(v) which.min(abs(LUT$SMC - v)))

pctl_long <- do.call(rbind, lapply(seq_along(pctl_rows), function(k) {
  data.frame(wave = wl, RFL = as.numeric(refl[pctl_rows[k], ]),
             Percentile = sprintf("P%d (SMC=%.1f%%)", pctl_probs[k] * 100, LUT$SMC[pctl_rows[k]]))
}))
pctl_plot <- ggplot(pctl_long, aes(x = wave, y = RFL, color = Percentile)) +
  geom_line(linewidth = 0.8) +
  theme_bw() +
  labs(x = "Wavelength (nm)", y = "Reflectance",
       title = "Example soil spectra at SMC's 10th/50th/90th percentile") +
  theme(plot.title = element_text(hjust = 0.5, size = 12))
print(pctl_plot)
ggsave(file.path(out_dir, "2-example_spectra_by_percentile.png"), plot = pctl_plot,
       width = 14, height = 9, dpi = 300, units = "cm")

## ----------------------------------------------------------------------------
## 7. Save
## ----------------------------------------------------------------------------

saveRDS(list(LUT = LUT, refl = refl, wl = wl, convolved = convolved,
             datasets = list(native = dataset_native, se2a = dataset_se2a),
             soil.id = soil.id, version = version, n.samples = n.samples, seed = seed),
        file.path(out_dir, "1-datasets.rds"))

dev.off()

cat("\nSaved to '", out_dir, "/1-datasets.rds' -- read back with readRDS() in 2-inversion_ML.R / ",
    "3-inversion_DL.R ($datasets$native / $se2a are ready for inversion of SMC).\n", sep = "")
