#' soil respiration
#'
#' @param Ts numeric. Soil temperature.
#'
#' @return R: soil respiration rate, in umol m-2 s-1. This R port always returns 0 (soil respiration is not simulated) - kept for interface compatibility with the original SCOPE model.
#' @export
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' soil_respiration(20)
soil_respiration<-function(Ts){
  R = 0 + 0 * Ts;   # in umol m-2 s-1
  return(R)
}
