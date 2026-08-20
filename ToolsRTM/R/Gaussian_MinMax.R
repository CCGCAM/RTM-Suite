
#' Gaussion distribution
#'
#' @param n number of cases in numeric format
#' @param m  mean value in numeric format
#' @param s  standard deviation in numeric format
#' @param lwr  min value in numeric format
#' @param upr max value in numeric format
#' @param nnorm  number of cases * a number (numeric format). This value is taken by the function to
#' take values in the selected range (min-max) upt to have the complete the number of cases
#'
#' @return a gaussian dist
#' @export
#'
#' 
gauss_byMin_Max <- function(n, m, s, lwr, upr, nnorm) {
  #set.seed(42)
  samp <- rnorm(nnorm, m, s)
  samp <- samp[samp >= lwr & samp <= upr]
  if (length(samp) >= n) {
    return(sample(samp, n))
  }  
  stop(simpleError("Not enough values to sample from. Try increasing nnorm."))
}