#' get.Extraterrestrial_radiance
#'
#' @param Ea0 Solar constant spectral irradiance 
#' @param DOY day of year  integer. 
#' @param tts solar zenith angle
#'
#' @return solar extraterrestrial spectrum 
#' @export 
#'
#' @examples
#' # Example usage:
#' # Assuming Ea0, DOY, and tts are your input values
#' Ea <- get.Extraterrestrial_radiance(Ea0 = 1361, DOY = 180, tts = 45)
#' print(Ea)
#' # 296.1273
get.Extraterrestrial_radiance <- function(Ea0, DOY, tts) {
  # Extraterrestrial radiation
  # Peiqi Yang (p.yang@utwente.nl)
  # peiqiyangweb@gmail.com
  # University of Twente.
  # 17-Nov-2019
  # input:    Ea0, solar constant spectral irradiance 
  #           DOY, day of year  integer. 
  # output:   Ea, solar extraterrestrial spectrum
  
  # NOTE:
  # Extraterrestrial radiation is the intensity (power) of the sun at the top of the Earths atmosphere.
  # It is usually expressed in irradiance units (Watts per square meter) on a plane normal to the sun. 
  # It varies throughout the year because of the Earths elliptical orbit
  # Which results in the Earths-Sun distance varying during the year in a predictable way. 
  # This effect can be represented empirically with the following equations:
  
  b <- 2 * pi * DOY / 365  # radians 
  correct_factor <- 1.00011 + 0.034221 * cos(b) + 0.00128 * sin(b) +
    0.000719 * cos(2 * b) + 0.000077 * sin(2 * b)
  
  La <- Ea0 * correct_factor * cos(tts * pi / 180) / pi
  
  return(La)
}

