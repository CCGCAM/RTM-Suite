#' Get Convolved Spectra
#'
#' This function performs spectral convolution using specified sensor spectral response functions (SRFs).
#' It can interpolate the input spectra to match the SRF resolution and plot the convolved and observed spectra.
#'
#' @param rfl A matrix or data frame of reflectance spectra, with columns representing wavelengths.
#' @param sensor A character string specifying the sensor to use for convolution. Options are "sentinel2a", "sentinel2b", or "prisma".
#' If NULL, "sentinel2a" is used by default.
#' @param plot.spectra A logical value indicating whether to plot the convolved and observed spectra. Default is TRUE.
#' @return A data frame of convolved spectra.
#' @export
#'
#' @examples
#' \dontrun{
#' #spectral response function downloaded from: https://oceancolor.gsfc.nasa.gov/resources/docs/rsr_tables/
#'
#' # 1. Generate a Lookup Table (LUT) for PROSAIL inputs
#'
#' inputs.prosail <- ToolsRTM::inputsPROSAIL
#' # Create a LUT with 1 entry using the defined PROSAIL inputs
#' df.LUT <- as.data.frame(ToolsRTM::getLUT(inputs = inputs.prosail, nLUT = 1, setseed = 1234))
#'
#' # 2. Simulate spectral reflectance using PROSAIL and apply BRDF effect
#'
#' # Simulate reflectance values using the PROSAIL model for the first LUT entry
#' reflectance_values <- ToolsRTM::foursail(inputLUT = df.LUT[1,], rsoil = rsoil, LeafModel = 'PROSPECT-PRO')
#' rdot <- reflectance_values[[1]]  # Direct reflectance
#' rsot <- reflectance_values[[2]]  # Diffuse reflectance
#'
#' # Compute the Bidirectional Reflectance Factor (BRF) using BRDF effects
#' reflectance_values <- ToolsRTM::Compute_BRF(rdot = rdot, rsot = rsot,
#'                                             tts = df.LUT[1, 'tts'],
#'                                             data.light = ToolsRTM::dataSpec_PDB)
#' # Convert the reflectance values to a matrix for further processing
#' sim.matrix <- as.matrix(t(reflectance_values))
#' colnames(sim.matrix)=paste0("X",wavelength)
#' # 3. Convolve the simulated reflectance with a sensor-specific SRF and plot the result
#'
#' # Here, we use the Sentinel-2A sensor for spectral resampling and plotting
#' convolved_spectra <- get.spectra.convolved(rfl = sim.matrix, sensor = "Sentinel2a", plot.spectra = TRUE)
#' }

#rfl=sim.canopy; sensor="prisma"
get.spectra.convolved <- function(rfl, sensor, plot.spectra=T){

  #Convert rfl to a matrix if it's not already one
  if (!is.matrix(rfl)) {
    rfl <- as.matrix(rfl)
  }


  # Check if rfl is valid
  if (is.null(rfl) || nrow(rfl) == 0) {
    stop("Input 'rfl' must be a non-empty matrix.")
  }

  # Initialize spectral response function (SRF) based on the sensor
  sensor_list <- c("Sentinel2a", "Sentinel2b", "PRISMA")
  if (missing(sensor) || is.null(sensor) || !sensor %in% sensor_list) {
    sensor <- "Sentinel2a"
  }

  if (sensor == "Sentinel2a" || sensor == "Sentinel2b") {
    srf <- ToolsRTM::srf.sentinel2a
    print(paste("Spectral resampling function to", toupper(sensor), "is being processed ..."))
    colnames(srf) <- gsub("S2A_SR_AV_", "", colnames(srf))
    colnames(srf)[colnames(srf) == "SR_WL"] <- "wavelength"
    bands <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B9", "B10", "B11", "B12")
  } else if (sensor == "PRISMA") {
    srf <- ToolsRTM::srf.prisma
    print('Spectral resampling function to PRISMA is being processed ...')
    bands <- grep("^X[0-9]+\\.?[0-9]*$", colnames(srf), value = TRUE)
  }

  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # 2.   Convolve with SRF table ----
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # Convolve with SRF table
  # Normalize SRF bands
  for (b in bands) {
    srf[, b] <- srf[, b] / sum(srf[, b], na.rm = TRUE)
  }

  # Create a new dataframe to store convolution results
  conv <- data.frame(id = 1:nrow(rfl), stringsAsFactors = FALSE)


  # Initialize the progress bar
  pb <- txtProgressBar(min = 0, max = length(bands), style = 3)
  # Convolve process
  for (idx in seq_along(bands)) {
    b <- bands[idx]
    srf_b <- as.data.frame(t(subset(srf, select = c("wavelength", b))))
    colnames(srf_b) <- paste0("X", srf_b[1, ])
    srf_b <- srf_b[-c(1), , drop = FALSE]
    srf_b <- srf_b[, srf_b[1, ] != 0, drop = FALSE]
    srf_b <- srf_b[, intersect(colnames(srf_b), colnames(rfl)), drop = FALSE]

    # Handle cases with one or more rows in rfl
    i <- rfl[, intersect(colnames(srf_b), colnames(rfl)), drop = FALSE]

    # If 'i' has one row, ensure it matches the structure for sweep
    if (nrow(i) == 1) {
      i_ <- sweep(as.matrix(i), 2, as.numeric(srf_b), FUN = "*")
    } else {
      i_ <- sweep(i, 2, as.numeric(srf_b), FUN = "*")
    }

    i_ <- as.data.frame(rowSums(i_, na.rm = TRUE))
    colnames(i_) <- c(paste(b))
    conv <- cbind(conv, i_)

    # Update progress bar
    setTxtProgressBar(pb, idx)
  }

  # Close the progress bar
  close(pb)

  # Band-name -> wavelength mapping, needed whether or not plot.spectra is
  # requested (previously only computed inside the plot.spectra==TRUE branch,
  # which meant plot.spectra=FALSE - the normal batch-processing case -
  # always crashed with "object 'db.sensor' not found").
  if (sensor == "Sentinel2a") {
    db.sensor <- subset(ToolsRTM::sensor.characteristics, Sensor == 'Sentinel2a')
    nm <- cbind(bands, as.numeric(db.sensor$average))
  } else if (sensor == "Sentinel2b") {
    db.sensor <- subset(ToolsRTM::sensor.characteristics, Sensor == 'Sentinel2b')
    nm <- cbind(bands, as.numeric(db.sensor$average))
  } else if (sensor == "PRISMA") {
    nm <- cbind(bands, as.numeric(ToolsRTM::fwhm.prisma$wavelength))
    db.sensor <- list(average = nm[, 2])
  }
  colnames(nm) <- c("band", "nm")

  #  interpolated_df <- cbind(rfl[non_x_columns], interpolated_df)
  ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # 3.   Plot
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  if (plot.spectra==T ) {
    if(sensor=="Sentinel2a" || sensor=="Sentinel2b" || sensor=="PRISMA") {
      for (b in colnames(conv)) {
        if (b %in% nm[, "band"]) {
          colnames(conv)[colnames(conv) == b] <- as.character(nm[nm[, "band"] == b, "nm"])
        }
      }
    }

      #Prepared convolved and observed reflectance for plotting
      #Convolved
      conv.plot=as.data.frame(conv)
      conv.plot$id=NULL
      conv.plot<- tidyr::gather(conv.plot, key = "band", value = "reflectance")

      summary_conv <- conv.plot |>
        dplyr::group_by(band) |>
        dplyr::summarise(
          average = mean(reflectance),
          median = median(reflectance),
          percentile_10 = quantile(reflectance, 0.10),
          percentile_50 = quantile(reflectance, 0.50),
          percentile_90 = quantile(reflectance, 0.90),
        )
      summary_conv$band <- as.numeric(as.character(summary_conv$band))

      #Observed
      rfl=as.data.frame(rfl)
      rfl$id=NULL
      rfl<- tidyr::gather(rfl, key = "band", value = "reflectance")

      summary_rfl <- rfl |>
        dplyr::group_by(band) |>
        dplyr::summarise(
          average = mean(reflectance),
          median = median(reflectance),
          percentile_10 = quantile(reflectance, 0.10),
          percentile_50 = quantile(reflectance, 0.50),
          percentile_90 = quantile(reflectance, 0.90),
        )

      summary_rfl$band <- as.numeric(as.character(gsub("X", "", summary_rfl$band)))

      plot_convoluted<-ggplot(summary_conv, aes(x = band)) +

        #Result from convolve
        geom_line(aes(y = average), color = "black", linewidth = 0.8) +
        geom_line(aes(y = median),  linetype = "dashed", color = "black", linewidth = 0.8) +
        geom_ribbon(aes(ymin = percentile_10, ymax = percentile_90),fill = "black", alpha = 0.3) +
        geom_point(aes(y = average), color = "black", fill = "darkgray", size = 2) +  # Points for convoluted data

        # Observed data
        geom_line(data = summary_rfl, aes(x = band, y = average), colour = "blue", linewidth = 0.6) +
        geom_line(data = summary_rfl, aes(x = band, y = median),  linetype = "dashed", color = "blue", linewidth = 0.6) +
        geom_ribbon(data = summary_rfl,aes(x = band, ymin = percentile_10, ymax = percentile_90), linetype = "dashed",fill = "blue", alpha = 0.3) +
        #geom_point(data = summary_rfl, aes(x = band, y = average), fill = "darkblue", size = 1.6) +
        labs(x = "wavelength (nm)",y = "Reflectance") +
        theme_bw() +  theme(legend.position="none",
                            plot.title = element_text(hjust = 0.5, size=12,face="bold"),
                            panel.background = element_rect(fill="white"),
                            axis.title = element_text(face="bold", size=14),
                            axis.text.y=element_text(color = "black", hjust = 0.5, size=12,face="bold"),
                            axis.text.x=element_text(color = "black", hjust = 0.5, size=12,face="bold"),
                            legend.title=element_blank())
      print(plot_convoluted)

      return(conv)
  } else {
    
     colnames(conv) <- c('id',db.sensor$average)
    
    return(conv)
      }




  }




