
#' Title
#'

#*#' J1 function with avoidance of singularity problem
#' edit 2017 12 28: change sampling of angles to match with  Jfunc1.m
#' @param k_para numeric. Extinction coefficient for direct (solar or observer) flux
#' @param l_para numeric.
#' @param t numeric. Leaf Area Index
#' @return Jout numeric.
#' @export

#' 
#' 
jfunc1_ <- function(k_para, l_para, t) {
  kl <- k_para - l_para
  Del <- kl * t
  minlt <- exp(-l_para * t)
  minkt <- exp(-k_para * t)
  Jout <- ifelse(abs(Del) > 1e-3,
                 (minlt - minkt) / kl,
                 0.5 * t * (minkt + minlt) * (1.0 - Del * Del / 12))
  return(list(Jout = Jout, minkt = minkt))
}