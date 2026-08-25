#' Planck function
#'
#' \code{get.Planck} estimates the spectral radiance emitted by a blackbody
#' (or greybody, if \code{em} is supplied) at given wavelength(s) and
#' temperature(s), following Planck's radiation law.
#'
#' @param wl numeric vector. Wavelength(s) of interest, in nanometers (nm).
#' @param Tb numeric vector. Temperature(s) of the emitting object, in Kelvin (K).
#' @param em numeric vector, optional. Emissivity of the object at each wavelength/temperature (dimensionless, 0-1). If missing, emissivity is set to 1 (ideal blackbody) for every element of \code{Tb}.
#'
#' @return numeric vector. Spectral radiance \code{Lb} emitted at each wavelength/temperature pair (W m^-2 sr^-1 um^-1).
#' @export
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' Lb <- get.Planck(wl, Tb)
#' }
get.Planck<-function(wl,Tb,em = NULL){

  c1 = 1.191066e-22
  c2 = 14388.33

  if (missing(em) == T){
    em = rep(1,length(Tb))
  }


  #em_wl <- outer(em, (wl * 1e-9)^(-5), FUN = "*")

  #second_<- outer(Tb, (wl * 1e-3), FUN = "*")

  #Lb = em_wl / (exp(c2 /second_ )-1);

  Lb <- em * c1 * (wl * 1e-9)^(-5) / (exp(c2 / (wl * 1e-3 * Tb)) - 1)






  return(Lb)

}
