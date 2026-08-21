#' get.spectral.convolution
#'
#' This function performs spectral convolution on the input irradiance (or radiance) dataframe using the spectral response function of the given sensor.
#'
#' @param df.irradiance A dataframe containing high-resolution irradiance (or radiance) data with wavelengths and corresponding values.
#' @param sensor.i A dataframe or list containing the spectral response function (wavelengths and sensitivity) of the sensor.
#' @param get.plots A boolean indicating whether to generate a plot showing the convolution results. Default is TRUE.
#'
#' @return A dataframe of convolved values for each sensor band, optionally with a plot if `get.plots` is TRUE.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage with sample irradiance data and sensor response function:
#' # Example irradiance data
#' df.irradiance <- data.frame(wavelength = seq(400, 700, 1), irradiance = rnorm(301))
#' # Example sensor response
#' sensor.i <- data.frame(wavelength = seq(400, 700, 1), response = runif(301))
#' get.spectral.convolution(df.irradiance, sensor.i, get.plots = TRUE)
#' 
#' }
get.spectral.convolution <- function(df.irradiance, sensor.i, get.plots=T) {
  
  # Spectral convolution for a given spectral response function
  # input:
  # get.spectral.convolution:    irradiance or radiance in high resolution, to be convoluted
  # sensor.i:       sensor characterisrtics

  sensors.properties = get.coef.SMAC(sensor = sensor.i)
  wlSensor =  sensors.properties[['wl.smac']]
  coefs.SMAC  = sensors.properties[['coefs.SMAC']]
  Sensor.name = sensors.properties[['Sensor.name']] ## this is mission name
  
  bands_df <- data.frame(sensors.properties$wl.srf.smac)
  weights_df <- data.frame(sensors.properties$p.srf.smac)
  
  if ( Sensor.name == 'Sentinel3B' || Sensor.name =='Sentinel3A' || Sensor.name == 'TerraAqua'){
    bands_df <- colMeans(bands_df,na.rm=T)
    bands_df <- round(bands_df, digits = 0)
    weights_df <- colMeans(weights_df,na.rm=T)
  }

  # Function to select wavelengths for each band including weights
  selectWavelengths <- function(band, weights) {
    merged_df <- merge(data.frame(wave = band, weight = weights), df.irradiance, by = "wave", all.x = TRUE)
    return(merged_df)
  }
 
  # Function to calculate convolution for each band
  calculateConvolution <- function(selected_band) {
    if (sum(selected_band$weight) == 0) {
      conv_result <- NA
    } else {
      conv_result <- sum(selected_band$weight * selected_band$irrad, na.rm = TRUE) / sum(selected_band$weight)
    }
    return(conv_result)
  }
  
  # Apply the function to each band with respective weights
  selected_wavelengths <- Map(selectWavelengths, bands_df, weights_df)
  # Apply the convolution function to each selected band
  convolution_results <- lapply(selected_wavelengths, calculateConvolution)
  
  df.conv = data.frame(wave=wlSensor,  Eo=as.numeric(convolution_results))
 
  
  if (get.plots ==  T){
    
    plot.conv <- ggplot2::ggplot(data = df.conv, aes(x = wave, y = Eo)) +
      labs(y= " Extraterrestrial irradiance", x = "") +
      geom_line() + theme_bw()
    
    print(plot.conv)
    
  }
  return(df.conv)
}





