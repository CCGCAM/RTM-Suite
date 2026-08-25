



#' Extract all bands from SE images (adapted to NETCDF dataset)
#'
#' @param netCDFs  path with the Netcdf
#' @param bands  a vector with the names of the inputs of the NetCDF
#' @param output  path of the outputs
#' @param pattern  File pattern to extract the netcdf (Area or studied region)
#'
#' @return tiffs
#' @export
#'
#' 
getTIFFs<-function(netCDFs=NULL, bands=NULL, output=NULL, pattern=NULL){

  if (!requireNamespace("ncdf4", quietly = TRUE)) {
    stop("getTIFFs() requires the 'ncdf4' package. Install it with install.packages(\"ncdf4\") and try again.", call. = FALSE)
  }

  is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
  
  if (is.null(pattern)){
    files_nc = list.files(netCDFs,pattern="*.nc$", full.names=TRUE)
  } else {
    files_nc = list.files(netCDFs,pattern=pattern, full.names=TRUE)
  }
  
  message(' TIFFs conversion is processing ...')

  progress_bar = txtProgressBar(min=0, max=length(files_nc), style = 3, char="=")
  
  for (nc_file in c(1:length(files_nc))){
    setTxtProgressBar(progress_bar, nc_file)
    ##############################################################################################################################
    #	1.  Get variables from NetCDF    -----    
    ##############################################################################################################################
    
    nc_data <- ncdf4::nc_open(paste(files_nc[nc_file],sep=''))
    # Save the print(nc) dump to a text file
    {
      sink(paste(files_nc[nc_file],'.txt',sep=''))
      print(nc_data)
      sink()
    }
    
    lon <- ncdf4::ncvar_get(nc_data, "x")
    lat <- ncdf4::ncvar_get(nc_data, "y", verbose = F)
    t <- ncdf4::ncvar_get(nc_data, "time")
    tunits<-ncdf4::ncatt_get(nc_data,"time",attname="units")
    tustr<-strsplit(tunits$value, " ")
    dates<-as.Date(t,origin=unlist(tustr)[3])

    ##############################################################################################################################
    #	2.  Save SE-2A Images by Day and create a stack     -----    
    ##############################################################################################################################
    
    
    for (i in c(bands)){
      B.array <- ncdf4::ncvar_get(nc_data, i) # store the data in a 3-dimensional array
      dim(B.array) 
      
      fillvalue <- ncdf4::ncatt_get(nc_data, i, "_FillValue")
      fillvalue
      
      B.array[B.array == fillvalue$value] <- NA
      for (k in c(1:dim(B.array)[3])){ 
        B.slice <- B.array[, , k] 
        CRS_ncdf<-'+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '
        # terra::rast()'s crs= takes a plain proj4/WKT string directly, no
        # separate CRS wrapper needed.
        r <- terra::rast(t(B.slice), extent = terra::ext(min(lon), max(lon), min(lat), max(lat)), crs = CRS_ncdf)
        path_daily=paste(output,'/daily','',sep='')
        ifelse(!dir.exists(path_daily), dir.create(path_daily), FALSE)
        path_out=paste(path_daily,'/day-',dates[k],'',sep='')
        ifelse(!dir.exists(path_out), dir.create(path_out), FALSE)
        # terra::writeRaster() infers format from the filename extension
        # rather than a separate format= argument -- filename needs an
        # explicit extension for that inference to work.
        filename<-paste(path_out,'/SE2A-',i,'-',dates[k],'.tif',sep="")
        terra::writeRaster(r, filename, filetype="GTiff", overwrite=TRUE)
        
      }
      
    }
    
    ncdf4::nc_close(nc_data) 
    
  }  ##### end loop for netcdf
  close(progress_bar)
  m_final<-message('TIFFs files were generated sucessfully')
  return(m_final)
}


