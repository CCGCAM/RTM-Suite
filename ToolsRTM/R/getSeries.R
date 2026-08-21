#' get series in a dataframe from rasters
#'
#' @param pathRaster folder with the raster in tif format
#' @param shapefile shapefile file (.shp)
#' @param band_names names of the bands of the raster
#' @param factorR factor for the reflectance bands
#' @param get.indices A boolean is True, get spectral indices pre-define in getIndicesSE2, is not (FALSE) provide only bands
#'
#' @return a dataframe with values
#' @export
#'
#'
getSeries<-function(pathRaster=NULL, shapefile=NULL, band_names=NULL,factorR=NULL,get.indices = T){
  options(warn=-1) ###avoid warnings

  if (is.null(factorR)){
    message('please provide factor scale for reflectance bands')
    stop()
  } else{
    factor = factorR
  }

  files = list.files(pathRaster,pattern="*.tif$", full.names=TRUE)
  dates = list.files(pathRaster,pattern="*.tif$", full.names=F)
  #dates=substr(dates,6,15)
  dates <- extract_dates_from_tiff_files(dates)

  if  (length(files) == 0){
    message('please provide raster in TIFF format')
    stop('TIF format is needed')
  } else {
      message('processing rasters ...')

      ##################################################################################

      progress_bar = txtProgressBar(min=0, max=length(files), style = 3, char="=")
      list.indices<-list()
      for (k in c(1:length(files))){
        setTxtProgressBar(progress_bar, k)
        #rs <- stack(files)
        ## by order
        # NB (raster->terra migration): bare brick() previously only resolved
        # if the caller happened to have library(raster) attached (raster was
        # never namespace-qualified here) -- see the ToolsRTM/vignettes/*
        # entries on this same bug class. terra::extract(cells=) replaces
        # raster::extract(cellnumbers=); terra::extract() already returns a
        # data.frame by default (no df= argument) and does not itself drop NA
        # rows (na.rm= has no terra equivalent).
        rs <- terra::rast(files[k])
        names(rs) <- band_names

        if (class(shapefile)[1] == 'SpatialPointsDataFrame'){
          r.extract <- terra::extract(rs, shapefile, cells = FALSE)
          data.write<-data.frame(ID=r.extract[,'ID'],Date= dates[k],r.extract[,c(2:dim(r.extract)[2])]* factor)
          #sensor.bands<-names(data.write[,grep(colnames(data.write),pattern="B",fixed = TRUE)])
        } else {
          r.extract <- terra::extract(rs, shapefile, cells = TRUE)
          data.write<-data.frame(ID=r.extract[,'ID'],cell = r.extract[,'cell'],Date= dates[k],r.extract[,c(3:dim(r.extract)[2])]* factor)
          #sensor.bands<-names(data.write[,grep(colnames(data.write),pattern="B",fixed = TRUE)])
        }

        #r.extract[is.nan(r.extract)] <- NA

        if (is.null(get.indices) | get.indices == F){

          list.indices[[k]]<- data.write
        } else {
          message(cat(' '))
          message(paste('adding spectral indices for the SE-2 image: ',k,'/',length(files),sep = ''))
          list.indices[[k]]<-ToolsRTM::getIndicesSE2(df=data.write,sensor='Sentinel-2a', df.data=data.write)

        }

      }
      ##################################################################################
  }
close(progress_bar)
data.export<-do.call(rbind.data.frame, list.indices)
return(data.export)

}


#' Extract acquisition dates from a set of TIFF file names
#'
#' @param tiff_files character vector. File paths or names containing a date in either \code{YYYY-MM-DD} or \code{YYYYMMDD} format, which will be extracted via regular expression.
#'
#' @return dates
#' @export
#'
extract_dates_from_tiff_files <- function(tiff_files) {
  # create an empty vector to store the dates
  dates <- vector("list", length(tiff_files))

  # iterate over each file name and extract the date using regular expressions
  for (i in seq_along(tiff_files)) {
    # define a regular expression to match both date formats
    regex <- "\\d{4}-\\d{2}-\\d{2}|\\d{8}"

    # extract the date from the file name using the regular expression
    date_str <- stringr::str_extract(tiff_files[i], regex)

    # print the file name and the extracted date for debugging
    #cat(sprintf("File: %s, Date: %s\n", tiff_files[i], date_str))

    # convert the date to Date format if it was found
    if (!is.na(date_str)) {
      dates[[i]] <- date_str #as.Date(date_str, format = "%Y-%m-%d", tryFormats = c("%Y%m%d"))
      #print(dates[[i]])
    } else {
      dates[[i]] <- 'missing'
      stop('files should contained date in this format: %Y-%m-%d or YYYYMMDD')
    }
  }

  # remove any missing dates from the resulting vector and return it
  return(unlist(dates[!is.na(dates)]))
}

