#' function for Spatial mapping of the main traits
#'
#' @param rasterFiles path with the image (in RasterBrick/RasterStack formats) or a single stackRaster
#' @param ForestLayer  path with the images of the forest mask. This code is adapted to from
#' COPERNICUS Land Monitoring Service
#' @param Sensor by default is 'Sentinel2a', is not needed
#' @param saveFile path with the path to save the spatial trait mapping
#' @param proj projection of the files (shapefile and raster should be in same projection)
#' @param shapeLayer path with the shapefile with your study area. This shape is use for cropping images
#' @param model.ML  Machine learning model
#' @param trait trais to estimate (actually for Cab and LAI)
#' @param factorR a numeric scaling factor applied to the raster reflectance values before index computation (e.g. to convert to 0-1 reflectance)
#' @return a list with: i) a scatter-plot between trait and best indicator, and ii) spatial trait mappping mask with the forest map
#' @export
#'
#' 


getSpatialTrait<-function(rasterFiles=NULL,ForestLayer=NULL,Sensor=NULL,
                          saveFile=NULL,proj=3035,shapeLayer=NULL,
                          model.ML =NULL, trait='Cab',
                          factorR=NULL){
  options(warn=-1) ###avoid warnings
  #message('Processing ....')
  
  shape.area<-sf::read_sf(shapeLayer,quiet = TRUE)
  files.sensor = rasterFiles#list.files(rasterFiles,pattern="*.tif$", full.names=TRUE)
  Split <- strsplit(files.sensor, "/")
  names.files = Split[[1]][length(Split[[1]])]#list.files(rasterFiles,pattern="*.tif$", full.names=FALSE)
  ############################################################
  ########### is projec is null apply proj =3035
  ############################################################
  if (is.null(proj)) {
    message('missing projection')
    message('EPSG:3035 will be used')
    proj=3035
  } 
  #proj_db <- system.file("proj/proj.db", package = "sf")
  #if (proj_db == "") proj_db <- proj_db_path
  #crs_table <- sf::read_sf(proj_db, "crs_view") # extracts the "crs_view" table
  #subset(crs_table, grepl("Belg|Ostend", name) & auth_name == "EPSG")[2:5]
  
  crs_files <- sf::st_crs(proj)
  crs_to_raster<-crs_files$proj4string
  ############################################################
  ### check bands for Sensor (only Sentinel-2 is implemented)
  ############################################################
  if (is.null(Sensor) | Sensor == 'Sentinel2a'){
    sensor.bands <- data.frame('B01'=442.7,'B02'=492.4, 'B03'=559.8, 'B04'=664.6, 'B05'=704.1, 'B06'=740.5,
                        'B07' = 782.8, 'B08' = 832.8, 'B8A' = 864.7, 'B09' = 945.1, 'B11' = 1613.7, 'B12' = 2202.4)
    
    sensor.rfl<-paste('R.',c(442.7,492.4,559.8,664.6,704.1,740.5,782.8,832.8,864.7,945.1,1613.7,2202.4),sep='')
  } else{
    message('Please provide information for bands')
    stop(' needs a sensor information')
  }
  ############################################################
  ## get Mosaic with FOrest Maks (Forest type 2018 later from COPERNICUS)
  ############################################################
  files.forest = ForestLayer
  # NB (raster->terra migration): bare brick() previously only resolved if
  # the caller happened to have library(raster) attached.
  forest_m<-terra::rast(files.forest[1])

  ##### Generate a mask for forest map in fishnet area
  forest_m<-terra::crop(forest_m,shape.area)
  # NB: neither this clamp() call nor the two later in this function pass
  # lower=/upper= -- both the raster package's clamp() and terra::clamp() default those to
  # -Inf/Inf, so with useValues=FALSE/values=FALSE and no bounds this was
  # already a no-op in the original the raster package's  code (nothing is ever outside
  # [-Inf, Inf]); preserved as-is here rather than guessing intended bounds.
  forest_m<-terra::clamp(forest_m, values=FALSE)
  forest_m[forest_m < 1] = NA
  #the raster package's values(forest_m)[the raster package's values(forest_m) < 1] = NA

  ##get_clouds of points
  ncases= round((terra::ncell(forest_m)-terra::freq(forest_m, value = NA)) / 100,0)
  if (ncases < 1500){
    ncases=1500
  }
  # NB: the raster package's sampleRandom(..., sp=TRUE) returned an sp Spatial*DataFrame;
  # terra's equivalent is spatSample(..., as.points=TRUE), which returns a
  # SpatVector instead -- a genuine type change, not just a renamed function.
  # $ID<- assignment, sf::st_as_sf(), and merge() below all have SpatVector
  # methods in terra/sf so this should carry through unchanged, but this
  # specific path could not be exercised against real forest-mask data here.
  shp_points <- terra::spatSample(forest_m, size = ncases, method = "random",
                                  xy = TRUE, as.points = TRUE, na.rm = TRUE)
  shp_points$ID<-c(1:dim(shp_points)[1])
  shp_sf<-sf::st_as_sf(shp_points,crs = st_crs(3035),quiet = TRUE)
  
  path_cloudPoints<-paste(saveFile,'points/',sep='')
  ifelse(!dir.exists(path_cloudPoints), dir.create(path_cloudPoints), FALSE)
  sf::st_write(shp_sf, paste(path_cloudPoints,"Forest_points.shp",sep=''),delete_layer = TRUE,quiet = TRUE) # overwrites
  
  # initiate progress bar for the files
  #progress_bar = txtProgressBar(min=0, max=length(files.sensor), style = 3, char="=")
  list.trait<-list()
  list.plot<-list()
  #for (i in c(1:length(files.sensor))){
  raster.to<- terra::rast(files.sensor[1]) * factorR
  names(raster.to)<-names(sensor.bands)
  ##get Reflectance for predict specific trait
  df.e <- terra::extract(raster.to, shp_points, cells = TRUE)
  df.e <- na.omit(df.e)
  colnames(df.e)<-c('ID', 'cell',sensor.rfl)
  df.e$pred<-predict(object = model.ML,df.e[,sensor.rfl])

  
  ############################################################
  ### filter by trait some outliers (negative and out of range)
  if (trait == 'Cab') {
    df.e$pred = ifelse(df.e$pred < 0.5, NA, df.e$pred)
    df.e$pred = ifelse(df.e$pred > 90, NA, df.e$pred)
    df.e<-na.omit(df.e)
    shape_togetIndices<-merge(shp_points[,'ID'],df.e, by.x='ID' ,by.y='ID', all.x=F)
  } else if (trait == 'LAI') {
    df.e$pred = ifelse(df.e$pred < 0.1, NA, df.e$pred)
    df.e$pred = ifelse(df.e$pred > 7, NA, df.e$pred)
    df.e<-na.omit(df.e)
    shape_togetIndices<-merge(shp_points[,'ID'],df.e, by.x='ID' ,by.y='ID', all.x=F)
  }
  
  ############################################################ 
  ### remove outliers
  ############################################################
  filter.column<-list()
  names_j<-colnames(df.e)[c(3:dim(df.e)[2])]
  for (j in c(1:length(names_j))){
    filter.column[[j]]<-ToolsRTM::filter_outliers(df=df.e, input=names_j[j])
  }
  ############################################################
  df<-cbind(df.e[,1:2],do.call(cbind, filter.column))
  colnames(df)<-names(df.e)
  df<-na.omit(df)
  wave=as.numeric( sub("R.", "", sensor.rfl, fixed = TRUE))
  df.indices<-ToolsRTM::getIndicesSE2(df=df, df.data=NULL)
  names_indices<-names(df.indices)
  names_df<-c('ID','cell',names(sensor.bands),'pred')
  colnames(df)<-names_df
  df.rfl.indices<-cbind(df,df.indices)
  
  ### change according to trait
  if (trait == 'Cab'){
    df.rfl.indices$pred_g<- cut(df.rfl.indices[,'pred'],breaks = seq(0,100,0.2))
    df.mean = aggregate(df.rfl.indices[,c(names(sensor.bands), names_indices, 'pred')],
                        by = list(TraitClass=df.rfl.indices$pred_g),
                        FUN=quantile, probs  = 0.5,  na.rm = TRUE)
    is.na(df.mean) <- sapply(df.mean, is.infinite)
  } else if (trait == 'LAI'){
    df.rfl.indices$pred_g<- cut(df.rfl.indices[,'pred'],breaks = seq(0,7,0.1))
    df.mean = aggregate(df.rfl.indices[,c(names(sensor.bands), names_indices, 'pred')],
                        by = list(TraitClass=df.rfl.indices$pred_g),
                        FUN=quantile, probs  = 0.5,  na.rm = TRUE)
    is.na(df.mean) <- sapply(df.mean, is.infinite)
  }

  r.corr<-list()
  intercept<-list()
  slope<-list()
  p.value<-list()
  
  ############################################################
  ### correlation datast
  ############################################################

  names.mean<-colnames(df.mean)
  
  for (k in 1:length(names.mean)) {
    k_=names.mean[k]
    
    df.sb<-df.mean[c('pred',k_)]
    names(df.sb)<-c('pred','index')
   
    #df.sb<-na.omit(df.sb)
    
    if ((k_ == 'pred') | (k_ == 'TraitClass') | (nrow(df.sb) == 0) ){
      r.corr[[k]]<-NA
      p.value[[k]] <-NA
      intercept[[k]]<-NA
      slope[[k]]<-NA
    } else {
          fmla <- as.formula(paste('pred'," ~ ", paste('index', collapse= "")))
          fit<-lm(fmla, data=df.sb)
          r.corr[[k]]<-round(cor(df.sb[,'pred'], df.sb[,'index'],use='pairwise.complete.obs')^2,5)
          p.value[[k]]<-summary(fit)$coefficients[2,4] 
          intercept[[k]]<-summary(fit)$coefficients[1,1] 
          slope[[k]]<-summary(fit)$coefficients[2,1]
          
    }
  }
  ############################################################
  ## select best index or band
  df.cor<-data.frame(Indices=names(df.mean),R=unlist(r.corr), R2=c(unlist(r.corr))**2, p.value=unlist(p.value),intercept=unlist(intercept),slope=unlist(slope))
  
  #sort by best R2
  df.cor<-na.omit(df.cor)
  df.cor <- df.cor[order(-df.cor$R2),] 

  path_corr<-paste(saveFile,'correlations/',sep='')
  ifelse(!dir.exists(path_corr), dir.create(path_corr), FALSE)
  
  write.table(df.cor, file=paste(path_corr,trait,'_',names.files,'_corr.csv',sep=''),
              sep = ",", append = T)
  message(df.cor[1,'Indices'])
  ############################################################
  ###  get plot using the best correlation
  ############################################################
  axis_x<-bquote(bold(.(df.cor[1,'Indices'])[''])) # axis x
  axis_y<-bquote(bold(.(trait)['predicted']))# axis y
  index=df.cor[1,'Indices']
  data.plot<-df.mean[,c('pred', index)]
  colnames(data.plot)<-c('y','index')

  mylabel.r.test = bquote(bold(r)^2 == .(format(df.cor[1,'R2'], digits = 3)))
  statsLabel = paste0("r2 = ", round(df.cor[1,'R2'],2))
  
  list.plot<-ggplot(data.plot, aes(y=y, x=index)) +
    geom_point(alpha=0.6) + geom_smooth(method=lm, formula = 'y ~ x') + theme_bw()+
    geom_abline(intercept = 0, slope = 1,linetype="dashed", size=0.5,color='gray')+
    xlab(axis_x) + ylab(axis_y) + ggtitle(statsLabel) 

  path_plots<-paste(saveFile,'plots/',sep='')
  ifelse(!dir.exists(path_plots), dir.create(path_plots), FALSE)
  print(list.plot)
  ggsave(paste(path_plots,trait,'_',names.files,'_using_',index,'.png',sep=''),width=6.45, height=7.47)
  
  ### Apply a glm Model
  getR2<-df.cor[1,3]
  if (getR2 <= 0.45) {
    
  }
  
  
  ############################################################
  ## get spatial index form Sentinel image
  ############################################################
  r.index<-ToolsRTM::getSpatial_index(rasterFiles=files.sensor,Sensor='Sentinel2a',
                                      SpectraltoCompute= index,
                              factorR=factorR)
  r.index<-unlist(r.index)

  ############################################################
  ### filter spatially outliers from trait map
  ############################################################
 
  r.trait<-df.cor[1,'slope']* r.index[[index]] +df.cor[1,'intercept']
  
  r.trait<-terra::clamp(r.trait, values=FALSE)
  if (trait == 'Cab'){
    r.trait[r.trait <= 0] = NA
    r.trait[r.trait >= 90] = NA
    #the raster package's values(r.trait)[the raster package's values(r.trait) <= 0] = NA
    #the raster package's values(r.trait)[the raster package's values(r.trait) >= 90] = NA
  } else if (trait == 'LAI') {
    r.trait[r.trait <= 0] = NA
    r.trait[r.trait >= 7] = NA
  }
  
  ############################################################
#  r.trait<-terra::crop(r.trait,ext(forest_m))
  #r.trait<-terra::crop(r.trait,forest_m)
  forest.resample <- terra::resample(forest_m,
                                     r.trait,
                              method='bilinear')
  ### save in a list the raster
  list.trait<-terra::mask(r.trait,forest.resample)
  list.trait<-terra::clamp(list.trait, values=FALSE)

  ############################################################
  ### write raster only if there is a path
  ############################################################

  if (!is.null(saveFile)){
    terra::writeRaster(list.trait, filename=paste(saveFile,trait,'_',names.files,sep=''),
                       filetype="GTiff", overwrite=TRUE)

  }
    ############################################################
    ############################################################  
  #} end for
  #close(progress_bar)
  ############################################################
  ############################################################
  to_export = list('raster'=list.trait,'plot'=list.plot)
  return(to_export)
}


#' remove outliers
#'
#' @param df  a dataframe with the column to remove outliers
#' @param input  string, name of the variable
#'
#' @return a vector 
#' @export
#'
#' 
filter_outliers<-function(df=NULL, input=NULL){
  qnt <- quantile(df[,input], probs=c(.25, .60), na.rm = T)
  outs <- quantile(df[,input], probs=c(.05, .95), na.rm = T)
  H <- 1.5 * IQR(df[,input], na.rm = T)
  df[,input][ df[,input] < (qnt[1] - H)] <- NA #outs[1]
  df[,input][ df[,input] > (qnt[2] + H)] <- NA #outs[2]
  vector<-df[,input]
  return(vector)
}


#' Get a Mosaic from TIFFs Tiles from (Forest type 2018 later from COPERNICUS)
#'
#' @param ForestLayers  path with the tiles
#' @param shapeLayer  a shapefile format with the extension of the study region
#' @param output path to write the forest mask
#' @param proj projection, by default 3035
#'
#' @return a mosaic
#' @export
#'
#' 
GetMosaics<-function(ForestLayers=NULL, shapeLayer=NULL,
                     output=NULL, proj=3035){
  ############################################################
  if (is.null(proj)) {
    message('missing projection')
    message('EPSG:3035 will be used')
    proj=3035
  } 
  crs_files <- sf::st_crs(proj)
  crs_to_raster<-crs_files$proj4string
  
  files.forest = list.files(ForestLayers,pattern="*.tif$", full.names=TRUE)
  shape.area<-sf::read_sf(shapeLayer,quiet = TRUE)
  Split <- strsplit(shapeLayer, "/")
  names.files = Split[[1]][length(Split[[1]])]#
  # NB (raster->terra migration): extent()/raster()/projection()<-/writeRaster()
  # were all called bare here too (no the raster package's  prefix) -- same "only worked if
  # the caller happened to have library(raster) attached" issue as elsewhere
  # in this file.
  e <- terra::ext(shape.area)
  template <- terra::rast(e)
  terra::crs(template) <- crs_to_raster
  terra::writeRaster(template, filename=paste(output,'Forest_mask_',names.files,'.tif',sep=''), filetype="GTiff", overwrite=TRUE)
  gdal.file.list <- lapply(files.forest, terra::rast)
  gdal.mosaic <- do.call(terra::mosaic, gdal.file.list)
  terra::writeRaster(gdal.mosaic, filename=paste(output,'Forest_mask_',names.files,'.tif',sep=''), overwrite=TRUE)
  ## clean memory in raster folder
  file.remove(list.files(paste(tempdir(),'raster/',sep=''),full.names = T))

  forest_m<-terra::rast(paste(output,'Forest_mask_',names.files,'.tif',sep=''))
  terra::crs(forest_m) <- crs_to_raster
  ##### Generate a mask for forest map in fishnet area
  forest.sb<-terra::crop(forest_m,shape.area)
  terra::values(forest.sb)[terra::values(forest.sb) < 1] = NA
  terra::writeRaster(forest.sb, filename=paste(output,'Forest_mask_',names.files,'_crop.tif',sep=''), filetype="GTiff", overwrite=TRUE)
  
  
  return(message('mosaic processing is done'))
}
  

