# NB: this function predates the raster->terra migration and has pre-existing
# issues unrelated to that migration, left untouched here since fixing them
# needs domain context (ICOS flux-tower / SCOPE anomaly-mapping workflow) this
# pass doesn't have -- flagged, not fixed: (1) a hardcoded debug override of
# `rasterFiles`/`indice` right below shadows whatever the caller actually
# passed in; (2) `names.files` is read (line originally ~40, before the
# migration) before it is ever assigned (only set two lines later, from the
# PREVIOUS loop iteration's value on iteration 1); (3) the function `return()`s
# `rast.2020`, a variable that is never defined anywhere in this function;
# (4) `readOGR()` depends on the `rgdal` package, archived from CRAN since
# 2023 and no longer installable via `install.packages()`. Only the the raster package's 
# dependency itself is migrated to terra:: below; the function remains
# non-callable as written.
getAnomalies<-function(rasterFiles=NULL,indice=NULL, Sensor='Sentinel2a',
                       factorR=NULL, freqAnomalies='Monthly',
                       date.min='2017-01-01', date.max='2019-12-31',
                       year.to.compare=2020, output=NULL ){

  ############################################################
  ############################################################
  #to rm
  rasterFiles='examples/SEdata/TIFFs/Area1/Indices/'
  indice='NDVI'
  ############################################################
  ############################################################

  if (is.null(factorR)){
    message('factor apply will be 1')
    factorR=1
  }
  
  if (Sensor == 'Sentinel2a'){
    acron='SE2A-'
  } else{
    message('please, this function is only for Sentinel2a data')
    stop()
  }
  
  ## Prepare the variables
  nfiles = list.files(rasterFiles,pattern=indice, full.names=T)
  
  ########## get all rasters in a list
  # initiate progress bar
  message('step 1: processing complete serie in a list.')
  bar.progress <- progress::progress_bar$new(format = "Processing [:bar] :percent in :elapsedfull, estimated time remaining :eta",
          total = length(nfiles), clear = F, width= 100)
  raslist<-list()
  dtimes<-c()
  for (i in c(1:length(nfiles))){
    bar.progress$tick()
    
    rast<- terra::rast(nfiles[i])
    names(rast)<-names.files
    raslist[[i]]<-rast
    Split <- strsplit(nfiles[i], paste(indice,'-',acron,sep=''))
    names.files = Split[[1]][length(Split[[1]])]
    date_subs<-(sub('\\.tif$', '', names.files))
    dtimes[i]<-date_subs
  }
  message('step 2: processing complete serie in a stack')
  dtimes<-as.Date(dtimes,"%Y-%m-%d")
  rastStack <- terra::rast(raslist)
  terra::time(rastStack) <- dtimes

  rast.anomalies <- rastStack[[which(terra::time(rastStack) >= date.min & (terra::time(rastStack) <= date.max))]]
  rast.to.compare <- rastStack[[which(terra::time(rastStack) >= as.Date(ISOdate(year.to.compare, 1, 1)))]]

  v.anomalies <- as.numeric(format(as.Date(terra::time(rastStack),format = "%Y-%m-%d"), format = "%m"))
  v.anomalies.to <- as.numeric(format(as.Date(terra::time(rastStack),format = "%Y-%m-%d"), format = "%m"))
  #sum layers by months
  r.monthly<- terra::tapp(rast.anomalies, index = v.anomalies, fun = mean)
  names(r.monthly) <- paste('Month-',substr(names(r.monthly),7,9),sep='')

  r.monthly.to.compare<- terra::tapp(rast.to.compare, index = v.anomalies.to, fun = mean)
  names(r.monthly.to.compare) <- paste('Month-',substr(names(r.monthly.to.compare),7,9),sep='')

  path_out<-paste(output,'/Anomalies/',sep="")
  ifelse(!dir.exists(path_out), dir.create(path_out), FALSE)


  terra::writeRaster(r.monthly, filename=paste0(path_out,indice,'-SE2A-monthly-',date.min,'-',date.max,'.tif'), overwrite=TRUE)
  terra::writeRaster(r.monthly.to.compare, filename=paste0(path_out,indice,'-SE2A-monthly-',year.to.compare,'.tif'), overwrite=TRUE)


  ###get temper
  shape.ICOS<-terra::vect(paste('Shapefiles/ICOS/Areas_affected_healthy.shp',sep=""))
  
  
  
  
  
  bar.progress <- progress::progress_bar$new(format = "Processing [:bar] :percent in :elapsedfull, estimated time remaining :eta",
                                             total = length(nfiles), clear = F, width= 100)
    



    
    rast.anomalies <- rastStack[[which(terra::time(rastStack) >= '2017-01-01' & (terra::time(rastStack) <= '2019-12-31'))]]
    rast.anomalies2 <- rastStack[[which(terra::time(rastStack) >= '2017-01-01' & (terra::time(rastStack) <= '2018-12-31'))]]
 
    raslist[[i]]<-stack.image
    
    return(rast.2020)
    
  }


