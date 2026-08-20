#' Get CSV Data from Multiple Folders
#'
#' This function reads CSV files from multiple folders and combines them into a list.
#'
#' @param path.out The path to the root folder containing the subfolders with CSV files.
#' @param n.folders The number of folders from the end of the directory structure to consider.
#' @param files.names The type of files to look for. Options: "Fluorescence", "Reflectance", or "Radiance".
#'
#' @return A list with \code{data} (list, one element per file pattern in \code{files.names}, each the row-bound data across the last \code{n.folders} folders), \code{names} (the file patterns with the \code{.csv} suffix stripped), and \code{lut} (the row-bound \code{Parameters/inputLUT.csv} across those same folders).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' path.out <- "outs"
#' n.folders <- 4
#' files.names <- "Fluorescence"
#' combined_data <- getCSV(path.out, n.folders, files.names)
#' }

getCSV <- function(path.out=NULL, n.folders=4, files.names='All') {


  ## find names ofr each type
  if (files.names == 'Fluorescence') {
    names_pattern <- c('fluorescence_All_Leaves.csv','fluorescence_scattered.csv','fluorescence_scattered.csv','fluorescence_shaded.csv',
                       'fluorescence_sunlit.csv','fluorescence_soil.csv',
                       'fluorescence_hemis.csv','fluorescence_ReabsCorr.csv')
  } else if (files.names == 'Reflectance'){
    names_pattern <-  c('refl.csv','reflapp.csv','rdd.csv','rdo.csv','rsd.csv','rso.csv')
  } else if (files.names == 'Radiance'){
    names_pattern <- c('Lo_spectrum_includingF.csv','Lo_spectrum.csv', 'Eout_spectrum.csv')
  } else if (files.names == 'Irradiance'){
    names_pattern <- c('Esky.csv','Esun.csv')
  } else if (files.names == 'Others'){
    names_pattern <- c('aPAR.csv','fluxes.csv','fluorescence_scalar.csv','radiation.csv','resistance.csv','vegatation.csv')
  } else if(files.names == 'All'){
    names_Lo_ <- c('Lo_spectrum_includingF.csv','Lo_spectrum.csv')
    names_rfl <- c('refl.csv','reflapp.csv','rdd.csv','rdo.csv','rsd.csv','rso.csv')
    names_irrad <- c('Esky.csv','Esun.csv')
    names_sif <-  c('fluorescence_All_Leaves.csv','fluorescence_scattered.csv','fluorescence_scattered.csv','fluorescence_shaded.csv',
                    'fluorescence_sunlit.csv','fluorescence_soil.csv',
                    'fluorescence_hemis.csv','fluorescence_ReabsCorr.csv')
    names_o <- c('aPAR.csv','fluxes.csv','fluorescence_scalar.csv','radiation.csv','resistance.csv','vegatation.csv')
    names_pattern <-c(names_rfl,names_Lo_,names_sif,names_irrad,names_o)



  }



  n.lut <- NULL

  # one accumulator per file pattern (not per folder): each element collects
  # that pattern's rows across all folders, so different file types (which
  # can share the same number of columns, e.g. same wavelength grid) never
  # get merged into one another
  data.all <- vector("list", length(names_pattern))

  folder_list <- list.dirs(path.out, full.names = TRUE, recursive = FALSE)
  folder_list <- tail(folder_list, n.folders)


  for (folder in folder_list) {

    for (j in 1:length(names_pattern)){
      names_pattern.j <-names_pattern[j]
      #print(names_pattern.j)

      # Read file from Parameters folder
      csv_file <- file.path(folder, names_pattern.j)


      if (file.exists(csv_file)) {
        data <- read.csv(csv_file)

        if (is.null(data.all[[j]])) {
          data.all[[j]] <- data
        } else if (ncol(data) == ncol(data.all[[j]])) {
          data.all[[j]] <- rbind(data.all[[j]], data)
        }

      }

    }

    # Read inputLUT from Parameters folder
    inputLUT_file <- file.path(folder, "Parameters", "inputLUT.csv")
    if (file.exists(inputLUT_file)) {
      lut <- read.csv(inputLUT_file)
      n.lut  <- rbind(n.lut, lut)
    }
  }

  # Get the current system time
  current_time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  # Create a folder with the current system time as part of the folder name
  folder_out <- paste(path.out,'/', current_time, '-merge/',sep = "")
  if (!dir.exists(folder_out)) {
    dir.create(folder_out)
  }

    # save using data.table
    # BUG (fixed): names_pattern[k] is already a full filename ending in
    # ".csv" (e.g. "refl.csv") -- appending another ".csv" here produced
    # doubled extensions like "refl.csv.csv". Beyond being cosmetically
    # wrong, this broke get.SCOPE.plots()'s exact-name setdiff() filter
    # (e.g. excluding "fluorescence_scalar.csv" no longer matched
    # "fluorescence_scalar.csv.csv", so that scalar/summary file leaked
    # into the spectral file list and crashed pivot_longer() downstream
    # with "cols must select at least one column").
    for (k in c(1:length(names_pattern))){
      if (!is.null(data.all[[k]])) {
        data.table::fwrite(data.all[[k]], paste(folder_out,'/',names_pattern[k],sep=''))
      }
    }

    folder_lut <- paste(folder_out,'Parameters/',sep = "")
    if (!dir.exists(folder_lut)) {
      dir.create(folder_lut)
    }

    if (!is.null(n.lut)) {
      data.table::fwrite(n.lut, paste(folder_lut,'inputLUT.csv',sep=''))
    }


  # Extracting the names from the patterns
  extracted_names <- gsub("\\.csv$", "", names_pattern)


  cat("\n")
  cat("List of files:\n")
  cat(names_pattern, sep=", ")
  cat("\n")
  return(list(data = data.all, names = extracted_names, lut = n.lut))
}






