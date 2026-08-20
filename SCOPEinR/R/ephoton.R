#' get.ephoton
#'
#' \code{get.ephoton} calculates the energy content (J) of 1 photon of wavelength lambda (m)
#'
#' @param lambda numeric vector. Wavelength(s) of the photon(s) (m).
#' @param constants data.frame. Physical constants table (as used throughout SCOPEinR), must contain rows for \code{h} (Planck's constant, J s) and \code{c} (speed of light, m s^-1) in columns \code{constant}/\code{value}.
#'
#' @return numeric vector. Energy content \code{E} of one photon at each wavelength \code{lambda} (J).
#' @export
#'
#' @author 	 Wout Verhoef, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' E <- get.ephoton(lambda, constants)
#' }
get.ephoton <- function(lambda, constants) {
  h <- subset(constants , constant == 'h')['value']
  #    Planck's constant [J s]
  h <-rep(h$value,length(lambda))
  c <- subset(constants , constant == 'c')['value'] #speed of light [m s-1]
  c <- rep(c$value,length(lambda))

  E <- h[1] * c[1] / lambda


  # energy of 1 photon [J]
  return(E)
}
