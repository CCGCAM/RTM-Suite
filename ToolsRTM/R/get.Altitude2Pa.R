#' Get atmospheric pressure (in hPa) as a function of altitude
#'
#' This function calculates the atmospheric pressure at a given altitude using the barometric formula.
#'
#' @param altitude Altitude in meters.
#' @param Pa0 Air pressure at sea level in hPa (default is 1013.25 hPa).
#' @param Ta0 Air temperature at sea level in Kelvin (default is 288.15 K).
#'
#' @return The atmospheric pressure at the given altitude in hPa.
#' @export
#'
#' @examples
#' # Atmospheric pressure at 1000 meters altitude
#' get.Altitude2Pa(altitude = 1000)
#'
#' # Atmospheric pressure at 5000 meters altitude with custom sea-level pressure
#' get.Altitude2Pa(altitude = 5000, Pa0 = 1010)
#' 
get.Altitude2Pa <- function(altitude, Pa0 = 1013.25, Ta0 = 288.15) {
  
  ## Atmospheric pressure (in hpa) as a function of altitude (in meters)
  # Author: Peiqi Yang (p.yang@utwente.nl)
  #         28-Aug-2019

  if (is.null(Pa0)){
    Pa0 = 1013.25
    
  }
  
  if (is.null(Ta0)){
    Ta0 = 288.15
    
  }
  
  # important constants:
    # normal temperature and pressure at sea level  Pa0 = 101325 (Pa)
    # temperature at sea level                      Ta0 = 15 degrees C or 288.15 degrees absolute
    # changes of temperature per 1km                dT  = 6.5 degrees per 1km
    # Earth-surface gravitational acceleration      g   = 9.80665 m/s2;
    # Molar mass of dry air                         M   = 0.02896968 kg/mol
    # Temperature lapse rate, = g/cp for dry air    L   = ~ 0.00976 K/m
    # Constant-pressure specific heat               Cp  = 1004.68506 J/(kg·K)
  # Universal gas constant                        R0  = 8.314462618 J/(mol·K)
  # Note: The gas constant is equivalent to the Boltzmann constant,
  #       but expressed in units of energy per temperature increment per mole,
  #       i.e. the pressure–volume product,
  #       rather than energy per temperature increment per particle.
  
  # [Barometric formula] (USED in here)
  # https://en.wikipedia.org/wiki/Barometric_formula
  # https://en.wikipedia.org/wiki/Atmospheric_pressure
  
  # Alternative: Standard Atmosphere model (USED in ORIGINAL SMAC)
  # 1. Assumes that temperature is 15 degrees C at sea level (288.15 degrees absolute?
  # 2. Drops 6.5 degrees per 1000 meters of altitude, up to 11000 meters.
  # 3. It assumes that at sea level air pressure is 101325 Newtons per square meter,
  # 4. It assumes air density is 1.225 kilograms per cubic meter.
  
  # Then
  # Pa/Pa0 = (Ta/Ta0)^2.558
  # Ta = Ta0 - 6.5/1000 * altitude;
  
  g <- 9.80665
  M <- 0.02896968
  R0 <- 8.314462618
  # Air pressure (Pa) in hPa
  Pa <- Pa0 * exp(-(g * altitude * M / (Ta0 * R0)))
  
  return(Pa)
}