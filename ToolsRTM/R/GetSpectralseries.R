#' Extract all bands from SE images (adapted to NETCDF dataset)
#'
#' @param netCDFs  List of NetCDF
#' @param bands  a vector with the names of the inputs of the NetCDF
#' @param shapefile  shapefile (point)
#' @param factorSE  factor to apply in reflectance files (point)
#' @param Indices  if TRUE or null, the function estimate the Indices for each date.
#' @return spectral indices
#' @export
#'
#'

GetSpectralseries<-function(netCDFs=NULL, bands=NULL, shapefile=NULL, factorSE = 1/10000, Indices = T){


  if (is.null(bands)) {
    message('indicate the bands of the sensor')
    close()
  }
  #bands<-c('B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B11','B12', 'SCL')
  #shapefile=shape.sb
  #netCDFs =path_netCDFs

  is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))
  #files_nc = list.files(netCDFs,pattern="*.nc$", full.names=TRUE)
  files_nc = netCDFs
  message(' Getting time series  ...')

  progress_bar = txtProgressBar(min=0, max=length(files_nc), style = 3, char="=")
  table.to.save<-list()

  for (nc_file in c(1:length(files_nc))){
    setTxtProgressBar(progress_bar, nc_file)
    nc_data <- ncdf4::nc_open(paste(files_nc[nc_file],sep=''))
    names_var = names(nc_data$var)
    lon <- ncdf4::ncvar_get(nc_data, "x")
    lat <- ncdf4::ncvar_get(nc_data, "y", verbose = F)
    t <- ncdf4::ncvar_get(nc_data, "time")
    tunits<-ncdf4::ncatt_get(nc_data,"time",attname="units")
    tustr<-strsplit(tunits$value, " ")
    dates<-as.Date(t,origin=unlist(tustr)[3])

    table.info<-list()
      for (i in c(bands)){
        B.array <- raster::brick(paste(files_nc[nc_file],sep=''), varname=i)
        if (class(shapefile)[1] == 'SpatialPointsDataFrame'){
          table.info[[i]]<- t(raster::extract(B.array, shapefile, df = F, na.rm = T, cellnumbers = F) * factorSE)

        }  else if (class(shapefile)[1] == 'sf'){
          shape.to<- as_Spatial(shapefile)
          table.extract<- (raster::extract(B.array, shape.to, df = T, na.rm = T, cellnumbers = T))
          levels(table.extract)<-levels(as.factor(shape.to[[names(shape.to)[1]]]))

          #table.extract<-cbind(ID_shape=levels(as.factor(shapefile[[names(shapefile)[1]]])),ID=table.extract[,'ID'],cell = table.extract[,'cell'],stack(table.extract[,c(3:dim(table.extract)[2])]* factorSE))
          table.extract <- cbind(ID=table.extract[,'ID'],cell = table.extract[,'cell'],stack(table.extract[,c(3:dim(table.extract)[2])]* factorSE))
          table.extract$ind<-substr(table.extract$ind,2,12)
          names(table.extract)<-c('ID','cell',i,'Date')
          #levels(table.extract$ID) = levels(as.factor(shapefile$id))

          table.info[[i]] <-table.extract
        } else {

          table.extract<- (raster::extract(B.array, shapefile, df = T, na.rm = T, cellnumbers = T))
          levels(table.extract)<-levels(as.factor(shapefile[[names(shapefile)[1]]]))

          #table.extract<-cbind(ID_shape=levels(as.factor(shapefile[[names(shapefile)[1]]])),ID=table.extract[,'ID'],cell = table.extract[,'cell'],stack(table.extract[,c(3:dim(table.extract)[2])]* factorSE))
          table.extract <- cbind(ID=table.extract[,'ID'],cell = table.extract[,'cell'],stack(table.extract[,c(3:dim(table.extract)[2])]* factorSE))
          table.extract$ind<-substr(table.extract$ind,2,12)
          names(table.extract)<-c('ID','cell',i,'Date')
          #levels(table.extract$ID) = levels(as.factor(shapefile$id))

          table.info[[i]] <-table.extract
        }

      }

    #### for points

    if (class(shapefile)[1] == 'SpatialPointsDataFrame'){
      table.in<-data.frame(dates,do.call(cbind, table.info))
      colnames(table.in)<-c('Date',bands)
      table.to.save[[nc_file]]<-table.in
      #### for polygons
    } else {


      Ids_ <- lapply(table.info, function(x) x[1])
      table.ids<-data.frame(do.call(cbind, Ids_[[1]]))

      #Ids_i <- lapply(table.info, function(x) x[2])
      #table.ids_i<-data.frame(do.call(cbind, Ids_i[[1]]))

      cell_ <- lapply(table.info, function(x) x[2])
      table.cell<-data.frame(do.call(cbind, cell_[[1]]))

      Date_ <- lapply(table.info, function(x) x[4])
      table.date<-data.frame(do.call(cbind, Date_[[1]]))

      bands_cols <- lapply(table.info, function(x) x[3])
      table.bands<-data.frame(do.call(cbind, bands_cols))

      table.in <-cbind(table.ids,table.cell,table.date,table.bands)

      #table.in <-cbind(ID=table.extract$ID,cell=table.extract$cell,Date=table.extract$Date,table.bands)
      table.to.save[[nc_file]]<-table.in
    }




  }

  table.f<-data.frame(do.call(rbind, table.to.save))

  if (Indices == T | is.null(Indices)) {

    #wavelengths.sentinel<-c(442.7,492.4,559.8,664.6,704.1,740.5,782.8,832.8,864.7,945.1,1613.7,2202.4)
    if (class(shapefile)[1] == 'SpatialPointsDataFrame'){
      table.f<-ToolsRTM::getIndicesSE2(df=table.f[,n_bands.points],sensor = "Sentinel-2a",df.data=table.f)
    } else {
      table.f<-ToolsRTM::getIndicesSE2(df=table.f[,n_bands.points],sensor = "Sentinel-2a",df.data=table.f)
    }
  }

  close(progress_bar)
  return(table.f)
  }
