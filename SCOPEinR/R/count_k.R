#' count_k
#'
#' \code{count_k} increments a multi-digit "odometer" counter, used to
#' enumerate all combinations of a set of discretized variables (e.g. when
#' building a full-factorial input LUT). Starting from digit \code{id}, every
#' digit that is already at its maximum value is reset to 1 and carried over
#' to the next digit (cyclically, wrapping from \code{nvars} back to 1); the
#' first digit found that is not at its maximum is incremented by one.
#'
#' @param nvars integer. Total number of digits (variables) in the counter.
#' @param v integer vector of length \code{nvars}. Current value of each digit.
#' @param vmax integer vector of length \code{nvars}. Maximum value allowed for each digit.
#' @param id integer. Index of the digit at which to start incrementing.
#'
#' @return integer vector of length \code{nvars}. Updated digit vector \code{vnew} after incrementing.
#' @export
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' vnew <- count_k(nvars, v, vmax, id)
#' }
count_k <- function(nvars, v, vmax, id) {
  i <- id

  # starting at id, set digits which are at its maximum equal to 1
  # first digit that is not at its maximum is incremented
  while (v[i] == vmax[i]) {
    v[i] <- 1
    i <- (i %% nvars) + 1
  }

  v[i] <- (v[i] %% vmax[i]) + 1
  vnew <- v

  return(vnew)
}
