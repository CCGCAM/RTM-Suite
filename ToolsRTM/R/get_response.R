
#' spectral response function
#'
#' @param center a vector with the center wavelength for your sensor
#' @param wl  a vector with the wavelength from the simulations
#' @param fwhm a  vector with the fwhm from each band
#'
#' @return spectral response
#' @export
#'
#' 
get_response.R<-function(center, wl, fwhm)
{
  a <- dnorm(wl, mean = center, sd = fwhm/2)
  a <- (a-min(a))/(max(a) - min(a))
  return(a)
}
