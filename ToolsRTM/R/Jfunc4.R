
#' J4 function for treating (near) conservative scattering
#' edit 2021 11 28: change sampling of angles to match with Jfunc4
#' @param m numeric. Extinction coefficient for direct (solar or observer) flux
#' @param t numeric. Leaf Area Index
#' @return Jout numeric.
#' @export
Jfunc4 <- function(m,t){
  
  del <- m * t
  out <- 0 * del
  out[del> 1e-3] <- (1 - exp(-del))/(m*(1 + exp(-del)))
  out[del <= 1e-3] <- 0.5 * t *(1.-del * del/ 12.)
  return(out)
}
