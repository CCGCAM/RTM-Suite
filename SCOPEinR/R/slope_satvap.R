#' calculates the saturated vapour pressure at temperature T (degrees C)
#' and the derivative of es to temperature s (kPa/C)
#'
#' @param Temp numeric. Air/leaf temperature, in degrees Celsius.
#'
#' @return The derivative of saturated vapour pressure with respect to temperature (slope), in kPa/degree C.
#' @export
#' @details Date: 2003.
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' slope_satvap(20)
slope_satvap<-function(Temp){

  slope_satvap<-list()

  ## constants
  a = 7.5;
  b = 237.3;         #degrees C
  log10 = 2.3026;

  # the output is in mbar or hPa. The approximation formula that is used is:
  # es(T) = es(0)*10^(aT/(b+T));
  # where es(0) = 6.107 mb, a = 7.5 and b = 237.3 degrees C
  # and s(T) = es(T)*ln(10)*a*b/(b+T)^2

  ## calculations
  es = 6.107 * 10^(7.5 * Temp / (b + Temp));
  slope_satvap$es = es
  s = es * log10 * a * b / (b + Temp)^2;
  slope_satvap$s = s
  return(slope_satvap)
}
