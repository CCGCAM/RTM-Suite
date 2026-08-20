#' Generate Spectral Response Function from Full Width at Half Maximum (FWHM)
#'
#' This function generates a spectral response function based on a Gaussian distribution, which is defined by the center wavelength and standard deviation (sd). A normalization based on the maximum (max) value is applied. The sensors considered in this function include:
#' - **Sentinel2a** and **Sentinel2b**: For these sensors, the function uses the spectral response data stored in `ToolsRTM::str.sentinel2a` and `ToolsRTM::str.sentinel2b`, respectively. These datasets provide the necessary characteristics to accurately model their spectral responses.
#' - **PRISMA**: The function utilizes the original resampling function for PRISMA, with the spectral response characteristics stored in `ToolsRTM::srf.prisma`. This allows for precise calculations based on the unique properties of the PRISMA sensor.
#' - Other sensors (e.g., ALI, Hyperion, Landsat4-8, MODIS, Quickbird, RapidEye, WorldView2-4, WorldView2-8) use the Full Width at Half Maximum (FWHM) function for their spectral response calculations, based on their respective characteristics.
#'
#' @param sensor A string specifying the sensor name. Options include:
#'   "ALI", "Hyperion", "Landsat4", "Landsat5", "Landsat7", "Landsat8",
#'   "MODIS", "Quickbird", "RapidEye", "Sentinel2a", "Sentinel2b",
#'   "WorldView2-4", "WorldView2-8".
#' @param path.out A string specifying the output path for the CSV files.
#' @param save Logical indicating whether to save SRF file (default is TRUE).
#' @param get.plot Logical indicating whether to plot the spectral response curves (default is FALSE).
#'
#' @return A data frame containing the spectral response functions.
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate spectral response function for Sentinel-2b
#' srf.matrix <- get.srf.from_fwhm(sensor = "Sentinel2b", save=T, path.out = "Tables", get.plot = TRUE)
#'
#' }
get.srf.from_fwhm <- function( sensor='', save=T, path.out, get.plot = FALSE) {

  # Check if the path.out argument is provided only when save is TRUE
  if (is.null(sensor) | missing(sensor)) {
    stop("The 'sensor' argument is missing. Please specify the output directory.")
  }
  # Check if the path.out argument is provided only when save is TRUE
  if (save && missing(path.out)) {
    stop("The 'path.out' argument is missing. Please specify the output directory.")
  }
  df.data <- ToolsRTM::sensor.characteristics

  # Check if the specified sensor is valid
  if (!sensor %in% unique(df.data$Sensor)) {
    stop("Invalid sensor name. Please select from the following options: ",
         paste(unique(df.data$Sensor), collapse = ", "))
  }

  print('Spectral resampling function will be processed ...')

  df.data$fwhm <- df.data$ub - df.data$lb
  df.data$wavelength <- (df.data$ub + df.data$lb) / 2

  fwhm_col <- as.numeric(grep("fwhm", names(df.data)))
  wl_col <- as.numeric(grep("wavelength", names(df.data)))

  # Subset data for the specified sensor
  df.data_sub <- df.data[df.data$Sensor == sensor, ]

  nm_range <- seq(300, 2600, by = 1)  # Cover full wavelength range (300 nm to 2600 nm)
  srf <- data.frame(wavelength = nm_range)  # Create base dataframe with wavelengths
  # Loop through each channel in df.data_sub

  if (sensor == 'Sentinel2a'){

    srf <- ToolsRTM::srf.sentinel2a

  } else if (sensor == 'Sentinel2b'){

    srf <- ToolsRTM::srf.sentinel2b


  } else if (sensor == 'PRISMA'){

    srf <- ToolsRTM::srf.prisma
    # Rename the columns based on the channels in df.data_sub
    for (i in 2:ncol(srf)) {
      # Create a new column name in the format 'PRISMA_SR_AV_<channel>'
      current_name <- colnames(srf)[i]
      new_col_name <- paste0("PRISMA_SR_AV_", gsub("X", "", current_name))
      # Assign the new column name
      colnames(srf)[i] <- new_col_name
    }


  } else {

    for (i in 1:nrow(df.data_sub)) {
      band_center <- df.data_sub[i, 'wavelength']
      fwhm_value <- df.data_sub[i, 'fwhm']
      lb <- df.data_sub[i, 'lb']  # Lower bound of the band
      ub <- df.data_sub[i, 'ub']  # Upper bound of the band


      # Calculate standard deviation from FWHM
      sd <- fwhm_value / (2 * sqrt(2 * log(2)))

      # Generate Gaussian curve over nm_range
      gaussian_response <- dnorm(nm_range, mean = band_center, sd = sd)

      # Normalize Gaussian response to have peak value of 1 at the center wavelength
      gaussian_response <- gaussian_response / max(gaussian_response)
      # Normalize Gaussian response to ensure sum equals 1
      #gaussian_response <- gaussian_response / sum(gaussian_response)

      # Set values outside the lb-ub range to zero
      gaussian_response[nm_range < lb | nm_range > ub] <- 0

      # Assign column name based on band wavelength
      col_name <- paste0("X", round(band_center, 1))

      # Add Gaussian response to srf
      srf[[col_name]] <- gaussian_response
    }
    # Default case if sensor is unknown
    for (i in 2:ncol(srf)) {
      current_name <- colnames(srf)[i]
      new_col_name <- paste0(sensor, "_SR_AV_", gsub("X", "", current_name))
      colnames(srf)[i] <- new_col_name
    }




  }



  # Write the output to a CSV file only if save is TRUE
  if (save) {
    # Define the paths for the output and 'Tables' directory
    if (!dir.exists(path.out)) {
      dir.create(path.out, recursive = TRUE)
      message(paste("Output directory created:", path.out))
    }
    # Write the output to a CSV file in the 'Tables' directory
    output_file <- file.path(path.out, paste0(sensor, ".csv"))  # Ensure the correct file path
    write.table(srf, file = output_file, sep = ',', row.names = FALSE)

    # Confirm that the file has been written
    message(paste("Output written to:", output_file))
  }

  # Reshape the data to long format for easier plotting
  srf_long <- reshape2::melt(srf, id.vars = "wavelength",
                             variable.name = "spectra_band",
                             value.name = "Response")

  # Remove "X" prefix from the spectra_band column
  #srf_long$spectra_band <- gsub("^X", "", srf_long$spectra_band)

  # Conditional plotting based on the get.plot argument
  if (get.plot) {
    plot_response <- ggplot(srf_long, aes(x = wavelength, y = Response, color = spectra_band)) +
      geom_line(linewidth = 1) +
      scale_colour_viridis_d(option = "plasma") +
      labs(x = "Wavelength (nm)", y = "Spectral Response",
           title = paste(sensor, "Spectral Response Curves")) +
      theme_bw() +
      theme(legend.position = "none",  # This removes the legend
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9)
      )
    print(plot_response)
  }

  return(as.matrix(srf))
}
