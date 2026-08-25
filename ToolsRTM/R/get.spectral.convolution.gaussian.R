#' Convolve reflectance onto a sensor using only NOMINAL band characteristics
#' (center wavelength + FWHM, or center + band edges) -- for sensors with
#' no real measured per-nm SRF curve available at all, approximated as a
#' Gaussian response (optionally truncated to a published band edge range).
#'
#' `get.spectral.convolution.rfl()` needs SMAC atmospheric-correction
#' coefficients bundled; `get.spectral.convolution.srf()` needs a real,
#' measured, per-nm SRF table. Neither exists for most sensors -- often all
#' that is published (or all a student has, e.g. from an instrument's own
#' ENVI header, or their own camera calibration sheet) is a list of band
#' center wavelengths and full widths at half maximum (FWHM). This function
#' covers that case, from three possible inputs:
#'
#' 1. `sensor.i = "EnMAP"` -- `ToolsRTM::EnMap.characteristics`'s 242
#'    channels (center + FWHM already given).
#' 2. `sensor.i` = one of `unique(ToolsRTM::sensor.characteristics$Sensor)`
#'    ("ALI", "Hyperion", "Landsat4", "Landsat5", "Landsat7", "Landsat8",
#'    "MODIS", "Quickbird", "RapidEye", "Sentinel2a", "Sentinel2b",
#'    "WorldView2-4", "WorldView2-8") -- these ship published band edges
#'    (`lb`/`ub`), not FWHM directly; FWHM is derived as `ub - lb` and the
#'    Gaussian response is additionally hard-truncated to `[lb, ub]` (same
#'    convention `get.srf.from_fwhm()` already uses for these sensors).
#'    (Sentinel-2A/B and PRISMA have a REAL measured SRF table bundled
#'    instead -- use `get.spectral.convolution.srf()` for those, it's more
#'    accurate than this Gaussian approximation.)
#' 3. Your OWN sensor/camera: supply `centers` yourself, in nm (e.g. copied
#'    straight out of an ENVI header's `wavelength = {...}` block). `fwhm`
#'    is optional -- if omitted, each band's width is approximated from its
#'    distance to its neighboring bands (the standard assumption for a
#'    CONTIGUOUS pushbroom imaging spectrometer, e.g. a Headwall camera,
#'    where the true per-band SRF calibration isn't available). If you do
#'    know each band's real FWHM (e.g. from a camera datasheet, or an ENVI
#'    header that includes an `fwhm = {...}` block), pass it explicitly for
#'    a more accurate result.
#'
#' @param df A data frame with a `wave` column (wavelength, nm) and an
#'   `rfl` column (reflectance) -- same convention as
#'   `get.spectral.convolution.rfl()`/`get.spectral.convolution.srf()`. Used
#'   for a SINGLE spectrum; for many spectra at once (e.g. an entire LUT),
#'   use `rfl`/`wave` instead (see below) -- much faster, vectorized as one
#'   matrix multiplication instead of one call per row.
#' @param sensor.i Optional: a bundled sensor name (see above). If given,
#'   `centers`/`fwhm` are looked up automatically and any values you also
#'   pass for them are ignored.
#' @param centers Optional (required if `sensor.i` is not given): your own
#'   sensor's band center wavelengths, nm.
#' @param fwhm Optional: your own sensor's per-band FWHM, nm, same length
#'   and order as `centers`. Derived from band spacing if omitted.
#' @param get.plots logical, plot the convolved spectrum? Default `FALSE`.
#'   Ignored when `rfl`/`wave` (bulk mode) is used instead of `df`.
#' @param rfl Optional: a matrix or data.frame of MANY reflectance spectra
#'   at once, one row per spectrum, columns matching `wave` (e.g. a LUT's
#'   simulated reflectance matrix). When supplied (together with `wave`),
#'   `df` is ignored and the whole matrix is convolved in one vectorized
#'   pass -- this is the fast path for large tables.
#' @param wave Required together with `rfl`: the wavelength grid (nm) `rfl`'s
#'   columns are on.
#'
#' @return Single-spectrum mode (`df`): a data frame with one row per band:
#'   `band` (index), `wl` (center wavelength, nm), `fwhm` (nm), `RFL`
#'   (convolved reflectance) -- sorted by `wl`. Bulk mode (`rfl`/`wave`): a
#'   numeric matrix, one row per input spectrum (same row order as `rfl`),
#'   one column per band (named by center wavelength, sorted by `wl`); the
#'   per-band `wl`/`fwhm` vectors are attached as `attr(out, "wl")` /
#'   `attr(out, "fwhm")`.
#' @export
#' @examples
#' df <- data.frame(wave = ToolsRTM::dataSpec_PDB[, 1],
#'                   rfl = 0.05 + 0.3 * ToolsRTM::dataSpec_PDB[, 1] / 2500)
#'
#' # A bundled sensor with only nominal characteristics (no measured SRF):
#' enmap_bands <- get.spectral.convolution.gaussian(df, sensor.i = "EnMAP")
#' modis_bands <- get.spectral.convolution.gaussian(df, sensor.i = "MODIS")
#'
#' # Your OWN sensor -- e.g. a 3-camera 15-band synchronized rig, with
#' # known band center + FWHM for every band:
#' own_centers <- c(444, 475, 502, 531, 550, 560, 570, 650, 668, 678, 705, 717, 740, 754, 842)
#' own_fwhm    <- c(28,  32,  18,  14,  12,  27,  14,  16,  14,  14,  10,  12,  18,  10,  57)
#' own_bands <- get.spectral.convolution.gaussian(df, centers = own_centers, fwhm = own_fwhm)
#'
#' # A pushbroom imaging spectrometer's own band list from its ENVI header
#' # (e.g. a Headwall camera) -- only wavelength centers known, no FWHM:
#' headwall_centers <- c(398.02, 400.25, 402.48, 404.71, 406.94)  # (truncated example)
#' headwall_bands <- get.spectral.convolution.gaussian(df, centers = headwall_centers)
#'
#' # Bulk mode: an entire LUT's reflectance matrix at once (fast, vectorized)
#' wave <- 400:2500
#' rfl_matrix <- matrix(0.05 + 0.3 * wave / 2500, nrow = 50, ncol = length(wave), byrow = TRUE)
#' bulk_bands <- get.spectral.convolution.gaussian(centers = own_centers, fwhm = own_fwhm,
#'                                                  rfl = rfl_matrix, wave = wave)
get.spectral.convolution.gaussian <- function(df = NULL, sensor.i = NULL, centers = NULL, fwhm = NULL,
                                               get.plots = FALSE, rfl = NULL, wave = NULL) {
  lb <- ub <- NULL

  if (!is.null(sensor.i)) {
    if (identical(sensor.i, "EnMAP")) {
      centers <- ToolsRTM::EnMap.characteristics$center
      fwhm <- ToolsRTM::EnMap.characteristics$fwhm
    } else {
      all_sensors <- unique(ToolsRTM::sensor.characteristics$Sensor)
      rows <- ToolsRTM::sensor.characteristics[ToolsRTM::sensor.characteristics$Sensor == sensor.i, ]
      if (nrow(rows) == 0) {
        stop('Unknown sensor.i "', sensor.i, '". Bundled options: "EnMAP", ',
             paste(all_sensors, collapse = ", "),
             '. For any other sensor (including your own camera), supply `centers` ',
             '(and optionally `fwhm`) directly instead of `sensor.i`.')
      }
      centers <- rows$average
      fwhm <- rows$ub - rows$lb
      lb <- rows$lb
      ub <- rows$ub
    }
  }
  if (is.null(centers)) {
    stop("Supply either `sensor.i` (a bundled sensor name) or your own `centers`.")
  }
  # A single scalar fwhm (one uniform width for every band, e.g. a regular
  # output grid) is recycled to match `centers` -- convenient for the common
  # "resample onto a uniform-FWHM grid" case, e.g. `wl <- seq(468, 842, 2)`.
  if (!is.null(fwhm) && length(fwhm) == 1 && length(centers) > 1) {
    fwhm <- rep(fwhm, length(centers))
  }

  ord <- order(centers)
  centers <- centers[ord]
  if (!is.null(fwhm)) fwhm <- fwhm[ord]
  if (!is.null(lb)) { lb <- lb[ord]; ub <- ub[ord] }

  if (is.null(fwhm)) {
    d <- diff(centers)
    fwhm <- c(d[1], (d[-length(d)] + d[-1]) / 2, d[length(d)])
  }
  stopifnot(length(centers) == length(fwhm))

  nbands <- length(centers)

  ## ---- bulk mode: many spectra at once, one vectorized matrix multiply ----
  if (!is.null(rfl)) {
    if (is.null(wave)) {
      stop("get.spectral.convolution.gaussian(): 'wave' is required together with 'rfl' (bulk mode).")
    }
    rfl <- as.matrix(rfl)
    wl_native <- wave

    W <- matrix(0, nrow = length(wl_native), ncol = nbands)
    for (i in seq_len(nbands)) {
      sigma <- fwhm[i] / (2 * sqrt(2 * log(2)))
      w <- exp(-0.5 * ((wl_native - centers[i]) / sigma)^2)
      if (!is.null(lb)) w[wl_native < lb[i] | wl_native > ub[i]] <- 0
      W[, i] <- w
    }
    col_sums <- colSums(W)
    col_sums[col_sums == 0] <- NA  # band has no overlap with wl_native -> output column is NA, not div-by-zero
    W <- sweep(W, 2, col_sums, "/")

    out <- rfl %*% W
    colnames(out) <- paste0("wl_", round(centers, 1))
    attr(out, "wl") <- centers
    attr(out, "fwhm") <- fwhm
    return(out)
  }

  ## ---- single-spectrum mode (unchanged behavior) ----
  if (is.null(df)) {
    stop("get.spectral.convolution.gaussian(): supply either 'df' (one spectrum) or 'rfl'+'wave' (many spectra at once).")
  }
  wl_native <- df$wave
  refl_native <- df$rfl

  out <- data.frame(band = seq_len(nbands), wl = centers, fwhm = fwhm, RFL = NA_real_)
  for (i in seq_len(nbands)) {
    sigma <- fwhm[i] / (2 * sqrt(2 * log(2)))
    w <- exp(-0.5 * ((wl_native - centers[i]) / sigma)^2)
    if (!is.null(lb)) w[wl_native < lb[i] | wl_native > ub[i]] <- 0
    if (sum(w) == 0) next
    out$RFL[i] <- sum(w * refl_native) / sum(w)
  }

  if (isTRUE(get.plots)) {
    plot.conv <- ggplot2::ggplot(data = out, ggplot2::aes(x = wl, y = RFL)) +
      ggplot2::labs(y = "Reflectance", x = "Wavelength (nm)") +
      ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::theme_bw()
    print(plot.conv)
  }
  out
}
