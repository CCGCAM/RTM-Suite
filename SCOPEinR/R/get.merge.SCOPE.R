#' Get.merge.Output generate the LUT table adding the apparent reflectance, radiance or fluorescence emission
#'
#' @param paths.sims character. Path to the folder of a SCOPE simulation run (containing the \code{Parameters/} subfolder and the output CSVs to merge).
#' @param inputs A vector specifying the input parameters to include in the final lookup table. If missing or NULL all values will be returned.
#' @param options A character vector specifying which additional data tables to merge into the LUT.
#' @param spectra.data A string specifying the type of spectral data to merge ("rfl.app", "rfl", "Lo", "Lo.sif", "sif", "sif.leaves", "sif.shaded", "sif.sunlit", "sif.scattered", "sif.hemis").
#'
#' @return The merged lookup table.
#' @export
#'
#' @examples
#' \dontrun{
#' paths.sims <- "path/to/sims/directory/"
#' inputs <- c('Cab','Car','Anth','LMA','EWT','Vcmax25','CBC','Prot','Cs',
#'              'Cbrown','Cx','LIDFa','LIDFb','LAI','Rin','Rli')
#' merged_LUT <- get.merge.SCOPE(paths.sims = paths.sims, inputs = inputs,
#'                                options = c("vegetation", "fluxes", "fluorescence", "aPAR"),
#'                                spectra.data = "rfl.app")
#' }

get.merge.SCOPE <- function(paths.sims=NULL, inputs, options = c("vegetation", "fluxes", "fluorescence", "aPAR"), spectra.data = 'rfl') {

  if (missing(options) | is.null(options)) {

    LUT <- data.table::fread(paste(paths.sims, 'Parameters/', 'inputLUT.csv', sep = ''), header = TRUE, sep = ',')
  } else {

    # Load necessary tables based on options
    loaded_tables <- list()
    if ("veg" %in% options) {
      loaded_tables[["veg"]] <- data.table::fread(paste(paths.sims, 'vegatation.csv', sep = ''), header = TRUE, sep = ',')
    }
    if ("fluxes" %in% options) {
      loaded_tables[["fluxes"]] <- data.table::fread(paste(paths.sims, 'fluxes.csv', sep = ''), header = TRUE, sep = ',')
    }
    if ("fluorescence" %in% options) {
      loaded_tables[["fluorescence"]] <- data.table::fread(paste(paths.sims, 'fluorescence_scalar.csv', sep = ''), header = TRUE, sep = ',')
    }
    if ("aPAR" %in% options) {
      loaded_tables[["aPAR"]] <- data.table::fread(paste(paths.sims, 'aPAR.csv', sep = ''), header = TRUE, sep = ',')
    }

    dataset <- do.call(cbind,loaded_tables)

    LUT <-data.table::fread(paste(paths.sims, 'Parameters/','inputLUT.csv',sep=''),header=T,sep=',')

    if (missing(inputs) | is.null(inputs)){
      inputs <- names(LUT)
    }
    LUT <- cbind(LUT[, ..inputs, with = FALSE], dataset)

  }


  # Perform the specified merge operation
  if (spectra.data == "rfl.app") {
    data.spectral<-SCOPEinR::get.spectra.SCOPE(getSpectral = T)
    wave <-data.spectral$wlS[1:2001]
    data_ <-data.table::fread(paste(paths.sims,'reflapp.csv',sep=''),header=T,sep=',')
    colnames(data_) <- paste0("R.",wave,sep='')

    LUT <- cbind(LUT,data_)

  } else if (spectra.data == "rfl") {

    data.spectral<-SCOPEinR::get.spectra.SCOPE(getSpectral = T)
    wave <-data.spectral$wlS[1:2001]
    data_ <-data.table::fread(paste(paths.sims,'refl.csv',sep=''),header=T,sep=',')
    colnames(data_) <- paste0("R.",wave,sep='')

    LUT <- cbind(LUT,data_)

  } else if (spectra.data == "Lo") {

    data.spectral<-SCOPEinR::get.spectra.SCOPE(getSpectral = T)
    wave <-data.spectral$wlS[1:2001]
    data_ <-data.table::fread(paste(paths.sims,'Lo_spectrum.csv',sep=''),header=T,sep=',')
    colnames(data_) <- paste0("Lo.",wave,sep='')

    LUT <- cbind(LUT,data_)

  } else if (spectra.data == "Lo.sif") {

    data.spectral<-SCOPEinR::get.spectra.SCOPE(getSpectral = T)
    wave <-data.spectral$wlS[1:2001]
    data_ <-data.table::fread(paste(paths.sims,'Lo_spectrum_includingF.csv',sep=''),header=T,sep=',')
    colnames(data_) <- paste0("Lo.",wave,sep='')

    LUT <- cbind(LUT,data_)

  } else if (spectra.data == "sif") {

    data.spectral<-SCOPEinR::get.spectra.SCOPE(getSpectral = T)
    wave <-data.spectral$wlF
    data_ <-data.table::fread(paste(paths.sims,'fluorescence.csv',sep=''),header=T,sep=',')
    colnames(data_) <- paste0("F.",wave,sep='')

    LUT <- cbind(LUT,data_)

  } else if (spectra.data == "sif.leaves") {

    data.spectral<-SCOPEinR::get.spectra.SCOPE(getSpectral = T)
    wave <-data.spectral$wlF
    data_ <-data.table::fread(paste(paths.sims,'fluorescence_All_Leaves.csv',sep=''),header=T,sep=',')
    colnames(data_) <- paste0("F.",wave,sep='')

    LUT <- cbind(LUT,data_)

  } else if (spectra.data == "sif.shaded") {

    data.spectral<-SCOPEinR::get.spectra.SCOPE(getSpectral = T)
    wave <-data.spectral$wlF
    data_ <-data.table::fread(paste(paths.sims,'fluorescence_shaded.csv',sep=''),header=T,sep=',')
    colnames(data_) <- paste0("F.",wave,sep='')

    LUT <- cbind(LUT,data_)

  } else if (spectra.data == "sif.sunlit") {

    data.spectral<-SCOPEinR::get.spectra.SCOPE(getSpectral = T)
    wave <-data.spectral$wlF
    data_ <-data.table::fread(paste(paths.sims,'fluorescence_sunlit.csv',sep=''),header=T,sep=',')
    colnames(data_) <- paste0("F.",wave,sep='')

    LUT <- cbind(LUT,data_)

  } else if (spectra.data == "sif.scattered") {

    data.spectral<-SCOPEinR::get.spectra.SCOPE(getSpectral = T)
    wave <-data.spectral$wlF
    data_ <-data.table::fread(paste(paths.sims,'fluorescence_scattered.csv',sep=''),header=T,sep=',')
    colnames(data_) <- paste0("F.",wave,sep='')

    LUT <- cbind(LUT,data_)

  } else if (spectra.data == "sif.hemis") {

    data.spectral<-SCOPEinR::get.spectra.SCOPE(getSpectral = T)
    wave <-data.spectral$wlF
    data_ <-data.table::fread(paste(paths.sims,'fluorescence_hemis.csv',sep=''),header=T,sep=',')
    colnames(data_) <- paste0("F.",wave,sep='')

    LUT <- cbind(LUT,data_)

  } else {
    stop("Invalid merge type. Please use 'rfl.app', 'rfl', 'Lo', 'Lo.sif','sif', 'sif.leaves', 'sif.shaded', 'sif.sunlit','sif.scattered' or 'sif.hemis'.")
  }


  return(LUT)
}
