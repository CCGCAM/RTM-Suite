#' Simpson-like trapezoidal integration
#'
#' \code{Sint} numerically integrates the vector \code{y} over \code{x} using
#' the trapezoidal rule, i.e. the sum of the areas of the trapezoids formed by
#' consecutive point pairs.
#'
#' @param y numeric vector. Values of the integrand, same length as \code{x}.
#' @param x numeric vector. Values of the integration variable (e.g. wavelength), same length as \code{y}; must be a monotonically increasing series.
#'
#' @return numeric value. The integral of \code{y} with respect to \code{x}.
#' @export
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @references WV Jan. 2013, for SCOPE 1.40
#'
#' @examples
#' \dontrun{
#' int <- Sint(y, x)
#' }
Sint <- function(y, x) {

  if (all(diff(x) > 0)) {
   # print("x is monotonically increasing")
  } else {
  #  print("x is not monotonically increasing")
  }

  nx <- length(x)

  if (length(x) == 1) {
    x <- t(x)
  }

  if (length(y) == length(x)) {
    y <- t(y)
  }

  step <- x[2:nx] - x[1:nx - 1]

  mean <- 0.5 * (y[1:nx - 1] + y[2:nx])

  ## in matlab retunr as unique value
  int <- sum(mean * step)

  return(int)
}

