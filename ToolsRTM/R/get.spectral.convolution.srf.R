#' Convolve reflectance onto a sensor using a plain per-band SRF table
#'
#' `get.spectral.convolution.rfl()` needs a sensor object with SMAC
#' atmospheric-correction coefficients bundled (see its own docs for the 9
#' supported sensors) -- fine for Landsat/OLCI/MODIS, but PRISMA has no SMAC
#' coefficients at all, and Sentinel-2A/B additionally ship a second, plain,
#' publisher-original spectral response function table alongside their SMAC
#' bundle. This function convolves onto any of those "just a per-band SRF
#' table, nothing else" sensors instead -- one shared implementation for
#' `ToolsRTM::srf.prisma`, `ToolsRTM::srf.sentinel2a`, and
#' `ToolsRTM::srf.sentinel2b` (generalizes what used to be two separate,
#' app-only helper functions, `convolve_prisma()`/`convolve_smac_sensor()`,
#' in the AEO-Course PROSAIL Shiny app's own app.R).
#'
#' @param df A data frame with a `wave` column (wavelength, nm) and an `rfl`
#'   column (reflectance) -- same convention as `get.spectral.convolution.rfl()`.
#' @param srf A plain per-band SRF table: its FIRST column is wavelength
#'   (nm, any column name), every other column is one sensor band's SRF
#'   weight at that wavelength (e.g. `ToolsRTM::srf.prisma`,
#'   `ToolsRTM::srf.sentinel2a`, `ToolsRTM::srf.sentinel2b`).
#' @param fwhm Optional data frame with a `fwhm` column, one row per band in
#'   `srf` (same order), e.g. `ToolsRTM::fwhm.prisma` -- when supplied, its
#'   (precisely interpolated) FWHM values are used instead of the coarser
#'   half-max-crossing estimate this function would otherwise derive
#'   directly from `srf`'s own sampled weight profile. Sentinel-2A/B have no
#'   equivalent bundled FWHM table, so leave this `NULL` for them.
#' @param get.plots logical, plot the convolved spectrum? Default `FALSE`.
#'
#' @return A data frame with one row per SRF column: `band` (index), `wl`
#'   (SRF-weighted mean center wavelength, nm), `fwhm` (full width at half
#'   maximum, nm), `RFL` (convolved reflectance).
#' @export
#' @examples
#' df <- data.frame(wave = ToolsRTM::dataSpec_PDB[, 1],
#'                   rfl = 0.05 + 0.3 * ToolsRTM::dataSpec_PDB[, 1] / 2500)
#' prisma_bands <- get.spectral.convolution.srf(df, ToolsRTM::srf.prisma,
#'                                               fwhm = ToolsRTM::fwhm.prisma)
#' s2a_bands <- get.spectral.convolution.srf(df, ToolsRTM::srf.sentinel2a)
get.spectral.convolution.srf <- function(df, srf, fwhm = NULL, get.plots = FALSE) {
  wl_native <- df$wave
  refl_native <- df$rfl

  wl_srf_all <- srf[[1]]
  band_cols <- colnames(srf)[-1]
  nbands <- length(band_cols)

  use_bundled_fwhm <- !is.null(fwhm) && nrow(fwhm) == nbands

  out <- data.frame(band = seq_len(nbands), wl = NA_real_, fwhm = NA_real_, RFL = NA_real_)
  for (i in seq_len(nbands)) {
    p_v_all <- srf[[band_cols[i]]]
    valid <- !is.na(p_v_all) & p_v_all > 0
    if (!any(valid)) next
    wl_v <- wl_srf_all[valid]; p_v <- p_v_all[valid]

    out$wl[i] <- sum(wl_v * p_v) / sum(p_v)

    if (use_bundled_fwhm) {
      out$fwhm[i] <- fwhm$fwhm[i]
    } else {
      half <- max(p_v) / 2
      above <- wl_v[p_v >= half]
      out$fwhm[i] <- if (length(above) >= 2) max(above) - min(above) else NA_real_
    }

    idx <- match(round(wl_v), round(wl_native))
    ok <- !is.na(idx)
    if (sum(ok) == 0) next
    out$RFL[i] <- sum(p_v[ok] * refl_native[idx[ok]]) / sum(p_v[ok])
  }
  out <- out[!is.na(out$wl), ]
  out <- out[order(out$wl), ]

  if (isTRUE(get.plots)) {
    plot.conv <- ggplot2::ggplot(data = out, ggplot2::aes(x = wl, y = RFL)) +
      ggplot2::labs(y = "Reflectance", x = "Wavelength (nm)") +
      ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::theme_bw()
    print(plot.conv)
  }
  out
}
