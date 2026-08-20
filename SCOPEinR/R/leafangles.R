






#'Subroutine FluorSail_dladgen (Version 2.3)
#'
#' \code{leafangles} computes the leaf inclination distribution function
#' (LIDF) frequencies for the 13 standard leaf angle classes, from the two
#' Verhoef ellipsoidal LIDF shape parameters, by numerically inverting the
#' cumulative distribution via \code{\link{dcum}}.
#'
#' @param a numeric. LIDFa parameter, controls the average leaf slope (bimodality parameter of the ellipsoidal distribution).
#' @param b numeric. LIDFb parameter, controls the distribution's bimodality/spread.
#'
#' @return A list with: \code{lidf} (numeric vector of length 13, the fraction of leaf area in each leaf angle class) and \code{litab} (numeric vector of length 13, the leaf angle class centers, degrees: \code{5,15,...,89}).
#' @export
#'
#' @author 	Joris Timmermans, Christiaan van der Tol  (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @references For more information look to page 128 of "theory of radiative transfer models applied in optical remote sensing of
#' vegetation canopies"
#'
#' @examples
#' \dontrun{
#' LeafDistribution <- leafangles(a = -0.35, b = -0.15)
#' }
leafangles<-function(a,b){

  # Subroutine FluorSail_dladgen taken from PROSAIL
  # FluorSail for Matlab
  # FluorSail is created by Wout Verhoef,

  litab <- c(5,15,25,35,45,55,65,75,81,83,85,87,89)
  freq=c() ##added
  for (i1 in c(1:8)){
    t   <- i1 * 10
    freq[i1] <- ToolsRTM::dcum(a,b,t)
  }

  for (i2 in c(9:12)){
    t   <- 80 + (i2 - 8) * 2
    freq[i2] <- ToolsRTM::dcum(a,b,t)
  }

  freq[13] <- 1
  for (i   in seq(13,2,-1)){
    freq[i] <- freq[i] - freq[i-1]
  }
LeafDistribution <- list("lidf" = freq,"litab" =litab) ##added
return(LeafDistribution) ##added
}



#' dcum function
#'
#' \code{dcum} computes the cumulative leaf inclination distribution value at
#' a given angle \code{theta}, following the ellipsoidal LIDF model of
#' Verhoef; for \code{a >= 1} an analytical solution is used, otherwise the
#' value is found by fixed-point iteration.
#'
#' @param a numeric. LIDFa parameter, controls the average leaf slope.
#' @param b numeric. LIDFb parameter, controls the distribution's bimodality.
#' @param theta numeric. Leaf inclination angle (degrees).
#' @return numeric. Cumulative LIDF value \code{f} at angle \code{theta}, dimensionless, in \[0,1\].
#' @export
#'
#' @author 	Joris Timmermans, Christiaan van der Tol  (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @details edit 2017-12-28: change sampling of angles to match with dcum.m
#'
#' @examples
#' \dontrun{
#' f <- dcum(a = -0.35, b = -0.15, theta = 40)
#' }
dcum <- function(a,b,theta) {

  ### this function is the same as ToolsRTM::dcum
  ## SubRoutines
  rd <- pi/180
  if (a >= 1){
    f <- 1 - cos(rd * theta)
  }
  else {
    eps <- 1e-8
    delx <- 1
    x <- 2 * rd * theta
    p <- x
    while (delx >=  eps){
      y <- a * sin(x) + 0.5 * b * sin(2.0 * x)
      dx <- 0.5*(y - x + p)
      x <- x + dx
      delx <- abs(dx)
    }
    f <- (2.0 * y + p) / pi
  }
  return(f) ##added
}


