#' Get the coefficients for the SMAC model
#'
#' This function retrieves the atmospheric correction coefficients for different satellite sensors
#' including Landsat, Sentinel-2, Sentinel-3, and MODIS (Aqua and Terra).
#'
#' @param sensor A character string specifying the sensor name. Options include "Landsat", "Sentinel-2", 
#'               "Sentinel-3", "MODIS_Aqua", and "MODIS_Terra".
#'
#' @return A list of coefficients corresponding to the specified sensor.
#' @export
#' @examples
#' \dontrun{
#' # Get coefficients for Landsat
#' coef_landsat <- get.coef.SMAC("Landsat")
#'
#' # Get coefficients for Sentinel-2
#' coef_sentinel2 <- get.coef.SMAC("Sentinel-2")
#' }
get.coef.SMAC <- function(sensor){

  Sensor.name <- sensor[['mission']]
  coef <-sensor[['SMAC_coef']]
  center.wvl <- sensor[['center_wvl']]
  wl.smac <- sensor[['wl_smac']]
  wl.srf.smac <-sensor[['wl_srf_smac']]
  # Assuming df is your data frame
  wl.srf.smac[is.na(wl.srf.smac)] <- NA
  p.srf.smac <- sensor$p_srf_smac
  Coeffs <- coef
  
  return(list(coefs.SMAC = Coeffs,wl.smac =wl.smac, 
              Sensor.name = Sensor.name,
              center.wvl = center.wvl, 
              wl.smac = wl.smac,
              wl.srf.smac = wl.srf.smac,
              p.srf.smac = p.srf.smac))

   

}
