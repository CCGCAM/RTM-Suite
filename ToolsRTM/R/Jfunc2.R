
#' J2 function with avoidance of singularity problem
#' edit 2017 12 28: change sampling of angles to match with  Jfunc2.m

#' @param k numeric. Extinction coefficient for direct (solar or observer) flux
#' @param l numeric.
#' @param t numeric. Leaf Area Index
#' @return Jout numeric.
#' @export
Jfunc2 <- function(k,l,t){
#	J2 function
Jout <- (1 - exp(-(k + l) * t))/(k + l)
}

