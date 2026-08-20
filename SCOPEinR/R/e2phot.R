#' get.e2phot calculates the number of moles of photons
#'
#' \code{get.e2phot} calculates the number of moles of photons
#' corresponding to E Joules of energy of wavelength lambda (m)
#'
#' @param lambda numeric vector. Wavelength(s) of the photons (m).
#' @param E numeric vector. Radiant energy at each wavelength (J), same length as \code{lambda}.
#' @param constants data.frame. Physical constants table (as used throughout SCOPEinR), must contain rows for \code{h} (Planck's constant, J s), \code{c} (speed of light, m s^-1) and \code{A} (Avogadro's number, mol^-1) in columns \code{constant}/\code{value}.
#'
#' @return numeric vector. Number of moles of photons (mol) corresponding to \code{E} Joules of energy at each wavelength \code{lambda}.
#' @export
#'
#' @author 	 Wout Verhoef, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' molphotons <- get.e2phot(lambda, E, constants)
#' }
get.e2phot <- function(lambda, E, constants) {
  e <- get.ephoton(lambda, constants)
  photons <- E / e
  A <- subset(constants , constant == 'A')['value']

  molphotons <- photons / A$value
  return(molphotons)
}

