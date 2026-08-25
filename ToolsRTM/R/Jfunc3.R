#' J3 function with avoidance of singularity problem
#' edit 2017 12 28: change sampling of angles to match with  Jfunc3.m

#' @param k numeric. Extinction coefficient for direct (solar or observer) flux
#' @param l numeric.
#' @param t numeric. Leaf Area Index
#' @return Jout numeric.
#' @export
Jfunc3 <- function(k,l,t){
out <- (1 - exp(-(k + l) * t))/(k + l)
}
