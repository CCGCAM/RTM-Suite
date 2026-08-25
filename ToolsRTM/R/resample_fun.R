

#' Generate a normalized Gaussian spectral response function
#'
#' @param center numeric. Center wavelength (nm) of the spectral band.
#' @param wl numeric vector. Wavelengths (nm) at which to evaluate the response function.
#' @param fwhm numeric. Full width at half maximum (nm) of the spectral band, controlling the width of the Gaussian.
#'
#' @return A numeric vector of relative spectral response values, normalized to the range 0-1, one per element of \code{wl}.
#' @export
#'
resample_fun<-function(center, wl, fwhm)
{
  a <- dnorm(wl, mean = center, sd = fwhm/2)
  a <- (a-min(a))/(max(a) - min(a))
  return(a)
}
