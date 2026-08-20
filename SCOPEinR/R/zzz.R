#' @importFrom ggplot2 aes geom_line ggplot ggsave guide_legend guides labs scale_color_gradientn scale_color_manual theme theme_bw xlab xlim ylab ylim
#' @importFrom dplyr group_by summarize ungroup %>%
#' @importFrom tidyr pivot_longer starts_with gather
#' @importFrom foreach %dopar%
#' @importFrom stats aggregate approx integrate median profile runif spectrum
#' @importFrom utils head install.packages read.csv read.table setTxtProgressBar tail txtProgressBar write.csv write.table
#' @importFrom grDevices terrain.colors
NULL

# Non-standard-evaluation column names / loop-local symbols referenced inside
# ggplot2/dplyr pipelines and data.frame constructions throughout the package
# (not actual undefined functions/variables) - registered here to silence the
# "no visible binding for global variable" R CMD check NOTE.
utils::globalVariables(c(
  "..inputs", "Average_Spectrum", "Date", "Emin_1", "Emin_5", "EoutF_", "EoutFrc_",
  "Eout_", "Eoutte_", "Eplu_1", "Eplu_5", "Esky", "Esky_", "Esun_", "F.all.Leaves",
  "Gs", "Jmax", "LoF_", "Lot_", "Lotot_", "Lototf_", "Mb", "Mf", "Pn1d", "Pn1d_CabP",
  "Spectrum", "T1", "T15", "WN", "Wavelength", "chunk", "constant", "data.timeseries",
  "db.sim", "emisi.soil", "group.trait", "k", "kCarrel", "kChlrel", "layer",
  "leafbio", "nr.", "phiI", "phiII", "phiSum", "pvertical", "rdd", "rdo", "refl",
  "reflapp", "rfl", "rfl.complete", "rfl.soil", "rsd", "rso", "sigmaF", "sind", "soil",
  "spectral", "value", "variable", "wave", "wave.f", "wavelength"
))

# Package-level call counter for get.computeA(), mirroring the original SCOPE
# model's global call-count diagnostic. Stored in a mutable environment
# (rather than a plain top-level binding) because the package namespace is
# locked once installed, which would make a `<<-` reassignment fail.
.scopeinr_state <- new.env(parent = emptyenv())
.scopeinr_state$fcount <- 0
