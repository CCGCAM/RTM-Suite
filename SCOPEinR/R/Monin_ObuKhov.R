
#' @title get.Monin.Obukhov function
#'
#' @description
#' \code{get.Monin.Obukhov} computes the Monin-Obukhov length, a stability
#' parameter used to characterize atmospheric stability and to correct the
#' aerodynamic resistances above the canopy for non-neutral conditions.
#' @param data.meteo list or data.frame. Meteorological forcing data; must contain \code{ustar} (friction velocity, m/s) and \code{Ta} (air temperature, deg C).
#' @param H numeric vector. Sensible heat flux from the canopy/soil to the atmosphere (W m^-2).
#'
#' @return numeric vector. Monin-Obukhov length \code{L} (m). Values that are \code{NA} (e.g. when \code{H} is zero) are set to \code{-1e6}, representing near-neutral stability.
#' @export
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' L <- get.Monin.Obukhov(data.meteo, H)
#' }

get.Monin.Obukhov<-function(data.meteo,H){


  rhoa =  subset( SCOPEinR::constants,constant == 'rhoa')[[2]]
  cp =  subset( SCOPEinR::constants,constant == 'cp')[[2]]
  kappa =  subset( SCOPEinR::constants,constant == 'kappa')[[2]]
  g =  subset( SCOPEinR::constants,constant == 'g')[[2]]

  L = -rhoa * cp * data.meteo[['ustar']]^3 * (data.meteo[['Ta']]  + 273.15) / (kappa * g * H)
  L[is.na(L)] <- -1e6
  return(L)

}
