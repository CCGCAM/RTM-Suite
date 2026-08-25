#' Brightness-Shape-Moisture soil model
#'
#' \code{getBSM} Calculates the spectrum of fluorescent radiance in the observer's direction and also
#' the TOC spectral hemispherical upward Fs flux.
#' @param soilpar list. Soil parameters: \code{BSMBrightness} (overall soil brightness), \code{BSMlat} (spectral shape "latitude", typical range 20-40 deg), \code{BSMlon} (spectral shape "longitude", typical range 45-65 deg).
#' @param spec list. Spectral inputs: \code{GSV} (Global Soil Vectors basis spectra, matrix nwl x 3), \code{Kw} (water absorption coefficient spectrum, nwl x 1), \code{nw} (refraction index of water spectrum, nwl x 1).
#' @param emp list. Empirical wetting parameters: \code{SMp} (soil moisture volume percentage, range 5-55), \code{SMC} (soil moisture capacity, recommended 0.25), \code{film} (effective optical thickness of a single water film, recommended 0.015).
#'
#' @return
#' 	Wet soil reflectance spectra across 400 nm to 2400 nm
#' @export
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' save.soil <- getBSM(soilpar, spec, emp)
#' }
getBSM <- function(soilpar, spec, emp) {

  # Spectral parameters

  # wl  = spec.wl          # wavelengths
  GSV <- spec[['GSV']]         # Global Soil Vectors spectra (nwl * 3)
  kw  <- spec[['Kw']]          # water absorption spectrum
  nw  <- spec[['nw']]          # water refraction index spectrum

  # Soil parameters

  B   <- soilpar[['BSMBrightness']]       # soil brightness
  lat <- soilpar[['BSMlat']]     # spectral shape latitude (range = 20 - 40 deg)
  lon <- soilpar[['BSMlon']]     # spectral shape longitude (range = 45 - 65 deg)
  SMp <- emp[['SMp']] #* 1E2   # soil moisture volume percentage (5 - 55)

  # Empirical parameters

  SMC  <- emp[['SMC']]  # soil moisture capacity parameter
  film <- emp[['film']]        # single water film optical thickness

  f1 <- B * sin((pi / 180) * lat)
  f2 <- B * cos((pi / 180) * lat) * sin((pi / 180) * lon)
  f3 <- B * cos((pi / 180) * lat) * cos((pi / 180) * lon)

  rdry <- f1 * GSV[, 1] + f2 * GSV[, 2] + f3 * GSV[, 3]

  # Soil moisture effect

  rwet <- soilwat(rdry, nw, kw, SMp, SMC, film)


  return(rwet)
}

#' @title soilwat function
#' \code{soilwat} In this model it is assumed that the water film area is built up
# according to a Poisson process. The fractional areas are as follows:
#' @param rdry dry soil reflectance (NW,1)
#' @param nw refraction index of water  in (NW,1)
#' @param kw absorption coefficient of water   in (NW,1)
#' @param SMp soil moisture volume percentage  (1,NS)
#' @param SMC soil moisture capacity (recommended 0.25) (1,1)
#' @param film effective optical thickness (deleff) of single water film (1,1) (recommended 0.015)
#'
#' @return soil spectrum
#' @export
#' @description

#' Date: Version 1.0, September 2012 Wout Verhoef
#' Update: Version 1.1, Jan 2020, Peiqi Yang

#' @author 	Wout Verhoef, Peiqi Yang, Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' soilwat(rdry, nw, kw, SMp, SMC, film)
#' }

soilwat <- function(rdry, nw, kw, SMp, SMC, film) {

  #

  # P(0)   = dry soil area
  # P(1)   = single water film area
  # P(2)   = double water film area
  # ...
  # et cetera

  # The fractional areas are given by P(k) = mu^k * exp(-mu) / k!

  #For water films of multiple thickness only the transmission loss due to water absorption is modified, since surface reflectance effects
  # are not influenced by the thickness of the film

  k <- c(0:6)
  nk <- length(k)

  mu <- (SMp - 5) / SMC # Mu-parameter of Poisson distribution


  if (mu <= 0) {
    rwet <- rdry
  } else {
    # Lekner & Dorf (1988) modified soil background reflectance
    #for soil refraction index = 2.0; uses the tav-function of PROSPECT

    rbac <- 1 - (1 - rdry) * (rdry * ToolsRTM::calctav(90, 2.0 / nw) / ToolsRTM::calctav(90, 2.0) + 1 - rdry)

    #  total reflectance at bottom of water film surface
    p <- 1 - ToolsRTM::calctav(90,nw) / nw^2# rho21, water to air, diffuse

    #reflectance of water film top surface, use 40 degrees incidence angle, like in PROSPECT
    Rw <- 1 - ToolsRTM::calctav(40, nw) #rho12, air to water, direct

    # fraction of areas
    # P(0) = dry soil area fmul(1)
    # P(1) = single water film area fmul(2)
    # P(2) = double water film area fmul(3)
    # without loop


    # calculate probability matrix
    fmul <- stats::dpois(round(mu), k)

    # calculate two-way transmittance
    tw <- exp(-2*outer(kw, film*k, "*")) # Two-way transmittance, exp(-2kw*k Delta)

    # calculate wet reflectance for each k
    Rwet_k <- Rw + (1-Rw) * (1-p) * tw * rbac / (1 - p * tw * rbac)
    rwet <- rdry * fmul[1] + Rwet_k[,2:nk] %*% fmul[2:nk]


  }
  return(rwet)
}


