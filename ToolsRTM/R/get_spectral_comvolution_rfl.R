#' Get Spectral Convolution for Reflectance
#'
#' This function performs spectral convolution on reflectance data using a specified sensor's
#' spectral response function. The convolution is essential for simulating how well the sensor
#' captures the true reflectance values based on its spectral characteristics.
#'
#' @param df A data frame containing the high-resolution reflectance data to be convoluted.
#' Must have a `wave` column (wavelength, nm) and an `rfl` column (reflectance) -- e.g.
#' `data.frame(wave = dataSpec_PDB[, 1], rfl = my_reflectance)`.
#' @param sensor.i The sensor to convolve onto -- one of the package's bundled sensor objects,
#' e.g. `ToolsRTM::Sentinel2A.MSI`, `ToolsRTM::LANDSAT8.OLI`, `ToolsRTM::TerraAqua.MODIS`
#' (also accepts the bare sensor name as a character string, e.g. `"LANDSAT8.OLI"`, which is
#' resolved to the matching bundled object automatically). Supported sensors:
#' "LANDSAT4.TM", "LANDSAT5.TM", "LANDSAT7.ETM", "LANDSAT8.OLI", "Sentinel2A.MSI",
#' "Sentinel2B.MSI", "Sentinel3A.OLCI", "Sentinel3B.OLCI", and "TerraAqua.MODIS".
#' @param get.plots A boolean value indicating whether to generate plots of the spectral convolution.
#'                  Default is TRUE.
#'
#' @return A data frame containing the convoluted reflectance values for the specified sensor,
#' with columns `wave` (band center wavelength, nm) and `RFL` (convolved reflectance).
#' @export
#' @examples
#' df <- data.frame(wave = ToolsRTM::dataSpec_PDB[, 1],
#'                   rfl = 0.05 + 0.3 * ToolsRTM::dataSpec_PDB[, 1] / 2500) # toy reflectance
#' convoluted_results <- get.spectral.convolution.rfl(df, sensor.i = "LANDSAT8.OLI", get.plots = FALSE)

get.spectral.convolution.rfl <- function(df, sensor.i, get.plots=T) {

  # Spectral convolution for a given spectral response function
  # input:
  # get.spectral.convolution:    irradiance or radiance in high resolution, to be convoluted
  # sensor.i:       sensor characterisrtics

  # sensor.i is documented as (and, in every other caller in this package --
  # get.smac(), and ToolsRTM.app's own Shiny app -- is always passed as) the
  # actual bundled sensor OBJECT, e.g. ToolsRTM::LANDSAT8.OLI, not a string.
  # A bare sensor-name string (as this function's own @examples used to show)
  # used to crash here with "subscript out of bounds", because
  # get.coef.SMAC() immediately does sensor[['mission']] on it -- resolve it
  # to the real object first so both calling conventions work.
  if (is.character(sensor.i)) {
    sensor.i <- get(sensor.i, envir = asNamespace("ToolsRTM"))
  }

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
    merged_df <- merge(data.frame(wave = band, weight = weights), df, by = "wave", all.x = TRUE)
    return(merged_df)
  }

  # Function to calculate convolution for each band
  calculateConvolution <- function(selected_band) {
    if (sum(selected_band$weight) == 0) {
      conv_result <- NA
    } else {
      conv_result <- sum(selected_band$weight * selected_band$rfl, na.rm = TRUE) / sum(selected_band$weight)
    }
    return(conv_result)
  }

  # Apply the function to each band with respective weights
  selected_wavelengths <- Map(selectWavelengths, bands_df, weights_df)
  # Apply the convolution function to each selected band
  convolution_results <- lapply(selected_wavelengths, calculateConvolution)

  # wlSensor (= sensor.i$wl_smac, in RAW ARRAY ORDER) is only the correct
  # band-center label for convolution_results[[b]] when the sensor's
  # id_smac_in_all mapping happens to be the identity permutation (true for
  # Sentinel-2A/B and Landsat 8 in the currently-bundled sensors, but NOT
  # for Landsat 4/5/7, Sentinel-3A/B OLCI, or Terra/Aqua MODIS -- confirmed
  # by inspecting id_smac_in_all for every bundled sensor). For those,
  # convolution_results[[b]] (computed in wl.srf.smac/p.srf.smac's own
  # natural per-band column order) was silently paired with the WRONG
  # wl.smac entry. Fixed the same way convolve_smac_sensor() (the
  # AEO-Course PROSAIL Shiny app's own, already-corrected version of this
  # exact convolution) does: re-derive each column's true center wavelength
  # from id_smac_in_all + the natural band number parsed out of
  # band_id_smac, instead of assuming column order == wl.smac order.
  id_map <- as.numeric(sensor.i$id_smac_in_all)
  band_labels <- sensor.i$band_id_smac
  has_labels <- length(band_labels) == length(id_map) && all(grepl("Band[[:space:]]*[0-9]+", band_labels))
  if (has_labels) {
    natural_band_no <- as.numeric(regmatches(band_labels, regexpr("[0-9]+", band_labels)))
    wave_out <- as.numeric(wlSensor)[match(id_map, natural_band_no)]
  } else {
    # No descriptive labels (Sentinel-2): id_smac_in_all is already a
    # direct 1:nbands identity permutation into wl_smac.
    wave_out <- as.numeric(wlSensor)
  }

  df.conv = data.frame(wave=wave_out,  RFL=as.numeric(convolution_results))


  if (get.plots ==  T){

    plot.conv <- ggplot2::ggplot(data = df.conv, aes(x = wave, y = RFL)) +
      labs(y= " Extraterrestrial irradiance", x = "") +
      geom_line() + theme_bw()

    print(plot.conv)

  }
  return(df.conv)
}
