#'  Stefan-Boltzmann equation
#'
#' \code{get.Stefan_Boltzmann} computes the total radiant exitance of a
#' blackbody at a given temperature, using the Stefan-Boltzmann law
#' (\eqn{H = \sigma T^4}).
#'
#' @param T_C numeric vector. Temperature(s) of the emitting surface, in degrees Celsius.
#'
#' @return numeric vector. Blackbody radiant exitance \code{H} (W m^-2) at each temperature in \code{T_C}.
#' @export
#'
#' @author Carlos Camino
#'
#' @examples
#' \dontrun{
#' H <- get.Stefan_Boltzmann(T_C = 20)
#' }
get.Stefan_Boltzmann<-function(T_C){

  if (!require("SCOPEinR")) {
    message("The 'SCOPEinR' package is not installed. Please install it using install.packages('ggplot2')")
  }
  Kelvin_temp <- subset(SCOPEinR::constants,constant == 'C2K')[[2]]
  sigmaSB <- subset(SCOPEinR::constants,constant == 'sigmaSB')[[2]];
  H <-sigmaSB * (T_C + Kelvin_temp)^4


  return(H)
}




