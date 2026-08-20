



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
    
    nc_data <- nc_open(paste(files_nc[nc_file],sep=''))
    # Save the print(nc) dump to a text file
    {
      sink(paste(files_nc[nc_file],'.txt',sep=''))
      print(nc_data)
      sink()
    }
    
    lon <- ncvar_get(nc_data, "x")
    lat <- ncvar_get(nc_data, "y", verbose = F)
    t <- ncvar_get(nc_data, "time")
    tunits<-ncatt_get(nc_data,"time",attname="units")
    tustr<-strsplit(tunits$value, " ")
    dates<-as.Date(t,origin=unlist(tustr)[3])

    ##############################################################################################################################
    #	2.  Save SE-2A Images by Day and create a stack     -----    
    ##############################################################################################################################
    
    
    for (i in c(bands)){
      B.array <- ncvar_get(nc_data, i) # store the data in a 3-dimensional array
      dim(B.array) 
      
      fillvalue <- ncatt_get(nc_data, i, "_FillValue")
      fillvalue
      
      B.array[B.array == fillvalue$value] <- NA
      for (k in c(1:dim(B.array)[3])){ 
        B.slice <- B.array[, , k] 
        CRS_ncdf<-'+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '
        r <- raster(t(B.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS(CRS_ncdf))
        path_daily=paste(output,'/daily','',sep='')
        ifelse(!dir.exists(path_daily), dir.create(path_daily), FALSE)
        path_out=paste(path_daily,'/day-',dates[k],'',sep='') 
        ifelse(!dir.exists(path_out), dir.create(path_out), FALSE)
        filename<-paste(path_out,'/SE2A-',i,'-',dates[k],sep="")
        writeRaster(r, filename, "GTiff", overwrite=TRUE)
        
      }
      
    }
    
    nc_close(nc_data) 
    
  }  ##### end loop for netcdf
  close(progress_bar)
  m_final<-message('TIFFs files were generated sucessfully')
  return(m_final)
}


