#' Get spectral bandset for interpolation
#'
#' This function retrieves the spectral bandset for various sensors
#' to be used in interpolation or analysis.
#'
#' @param sensor A character string with the name of the sensor; available options are:
#'               "ALI", "Hyperion", "Landsat4", "Landsat5", "Landsat7", "Landsat8",
#'               "MODIS", "Quickbird", "RapidEye", "Sentinel2a", "Sentinel2b", 
#'               "WorldView2-4", "WorldView2-8".
#'
#' @return A data frame containing the spectral bands for the specified sensor.
#' @export
#'
#' @examples
#' df.sentinel2a <- get.spectral(sensor='Sentinel2a')
get.spectral <- function(sensor='Sentinel2a') {

  data<-ToolsRTM::sensor.characteristics
  sensors<-unique(data$Sensor)
  if (sensor %in% sensor) {
    df.sensor <- subset(data, Sensor == sensor)
    return(df.sensor)
  } else {
    stop(paste('Please use the avalaible sensors : ',sensors,sep = ''))
  }

}
