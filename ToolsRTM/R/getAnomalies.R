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
    
    rast<- raster::raster(nfiles[i]) 
    names(rast)<-names.files
    raslist[[i]]<-rast
    Split <- strsplit(nfiles[i], paste(indice,'-',acron,sep=''))
    names.files = Split[[1]][length(Split[[1]])]
    date_subs<-(sub('\\.tif$', '', names.files))
    dtimes[i]<-date_subs
  }
  message('step 2: processing complete serie in a stack')
  dtimes<-as.Date(dtimes,"%Y-%m-%d")
  rastStack <- raster::stack(raslist)
  rastStack <- raster::setZ(rastStack, dtimes)

  rast.anomalies <- subset(rastStack, which(getZ(rastStack) >= date.min & (getZ(rastStack) <= date.max)))
  rast.to.compare <- subset(rastStack, which(getZ(rastStack) >= as.Date(ISOdate(year.to.compare, 1, 1))))

  v.anomalies <- as.numeric(format(as.Date(rastStack@z$time,format = "%Y-%m-%d"), format = "%m"))
  v.anomalies.to <- as.numeric(format(as.Date(rastStack@z$time,format = "%Y-%m-%d"), format = "%m"))
  #sum layers by months
  r.monthly<- stackApply(rast.anomalies, v.anomalies, fun = mean)
  names(r.monthly) <- paste('Month-',substr(names(r.monthly),7,9),sep='')
  
  r.monthly.to.compare<- stackApply(rast.to.compare, v.anomalies.to, fun = mean)
  names(r.monthly.to.compare) <- paste('Month-',substr(names(r.monthly.to.compare),7,9),sep='')
  
  path_out<-paste(output,'/Anomalies/',sep="")
  ifelse(!dir.exists(path_out), dir.create(path_out), FALSE)
  
  
  writeRaster(r.monthly, filename=paste(path_out,indice,'-SE2A-monthly-',date.min,'-',date.max,sep=""), "GTiff", overwrite=TRUE,bylayer=T)
  writeRaster(r.monthly.to.compare, filename=paste(path_out,indice,'-SE2A-monthly-',year.to.compare,sep=""), "GTiff", overwrite=TRUE,bylayer=T)
  
  
  ###get temper
  shape.ICOS<-readOGR(paste('Shapefiles/ICOS/Areas_affected_healthy.shp',sep=""))
  
  
  
  
  
  bar.progress <- progress::progress_bar$new(format = "Processing [:bar] :percent in :elapsedfull, estimated time remaining :eta",
                                             total = length(nfiles), clear = F, width= 100)
    



    
    rast.anomalies <- subset(rastStack, which(getZ(rastStack) >= '2017-01-01' & (getZ(rastStack) <= '2019-12-31')))
    rast.anomalies2 <- subset(rastStack, which(getZ(rastStack) >= '2017-01-01' & (getZ(rastStack) <= '2018-12-31')))
 
    raslist[[i]]<-stack.image
    
    return(rast.2020)
    
  }


