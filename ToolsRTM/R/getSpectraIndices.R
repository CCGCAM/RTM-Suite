
#'  this function produces Spatial maps of the spectral indices
#'
#' @param rasterFiles path with the images (in RasterBrick/RasterStack formats)
#' @param Sensor character with the sensor 'Sentinel2a'
#' @param SpecEq Character. Refer to an expression corresponding to the spectral index to compute
#' @param factorR numeric. multiplying factor used to write reflectance in image (==10000 for S2)
#' @param SpectraltoCompute List with the spectral index to compute. by default 'All'
#' @param path.export path_to save Bands
#' @param single.bands extract singles_bands, by deault is False
#' @return spectral indices
#' @export
#'
#' 
getSpectraIndices<-function(rasterFiles=NULL,Sensor='Sentinel2a',SpecEq=NULL,SpectraltoCompute= 'All',
                            factorR=NULL, path.export=NULL, single.bands=T){
  options(warn=-1) ###avoid warnings
  
  if (is.null(Sensor) | Sensor == 'Sentinel2a'){
    Bands <- data.frame('B01'=442.7,'B02'=492.4, 'B03'=559.8, 'B04'=664.6, 'B05'=704.1, 'B06'=740.5,
                          'B07' = 782.8, 'B08' = 832.8, 'B8A' = 864.7, 'B09' = 945.1, 'B11' = 1613.7, 'B12' = 2202.4)
  } else{
    message('not implemented yet')
    stop()
  } 
  if (is.null(SpectraltoCompute)) {
    message('please specific one spectral indicator: see ?getSpectraIndices')
    message('or use All to estimate all indices')
    stop()
  }
  
  if (is.null(factorR)) {
    factorR=1/10000
    message('factor scale used: factorR=1/10000')
  }

  files = list.files(rasterFiles,pattern="*.tif$", full.names=TRUE)
  names_files = list.files(rasterFiles,pattern="*.tif$", full.names=F)
  indices <- list()
  # initiate progress bar
  bar.progress <- progress::progress_bar$new(format = "Processing [:bar] :percent in :elapsedfull, estimated time remaining :eta",
                                             total = length(files), clear = F, width= 100)
  for (i in c(1:length(files))){
   
    bar.progress$tick()
    # NB (raster->terra migration): bare brick() previously only resolved if
    # the caller happened to have library(raster) attached.
    r<- terra::rast(files[i]) * factorR
    names(r)<-names(Bands)
  
      if (is.null(single.bands) | single.bands == F){
        message('only spectral indices will be computed')
      } else {
        indices[['B01']] <- r[['B01']]
        indices[['B02']] <- r[['B02']]
        indices[['B03']] <- r[['B03']]
        indices[['B04']] <- r[['B04']]
        indices[['B05']] <- r[['B05']]
        indices[['B06']] <- r[['B06']]
        indices[['B07']] <- r[['B07']]
        indices[['B08']] <- r[['B08']]
        indices[['B8A']] <- r[['B8A']]
        indices[['B09']] <- r[['B09']]
        indices[['B11']] <- r[['B11']]
        indices[['B12']] <- r[['B12']]
      }
  
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'NDVI'){
        # Rouse et al. (1974)
        indices[['NDVI']] <- (r[['B08']] - r[['B04']]) / (r[['B08']] + r[['B04']])
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'RDVI'){
        # Rougean and Breon (1995)
        indices[['RDVI']] <- (r[['B08']] - r[['B04']]) / (r[['B08']] + r[['B04']]) ** 0.5
      }
    
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'SR'){
        indices[['SR']] <- r[['B08']] / r[['B04']]
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'MSR'){
        # Chen (1996)
        indices[['MSR']] <- (r[['B08']] / r[['B04']] - 1) / ((r[['B08']] / r[['B04']]) ** 0.5 + 1)
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'OSAVI'){
        # Rondeaux et al. (1996)
        indices[['OSAVI']] <- ((1 + 0.16) * (r[['B08']] - r[['B04']]) / (r[['B08']] + r[['B04']] + 0.16))
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'MSAVI'){
        # Qi et al. (1994)    
        indices[['MSAVI']] <- 1 / 2 * (2 * r[['B08']] + 1 - sqrt(((2 * r[['B08']] + 1) ^ 2) - 8 * (r[['B08']] - r[['B04']])))
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'MTVI1'){
        # Broge & Leblanc (2000); Haboudane et al. (2004)
        indices[['MTVI1']] <- 1.2 * (1.2 * (r[['B08']] - r[['B03']]) - 2.5 * (r[['B04']] - r[['B03']]))
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'MTVI2'){
        # Haboudane et al. (2004)    
        indices[['MTVI2']] <- (1.5 * (1.2 * (r[['B08']] - r[['B03']]) - 2.5 * (r[['B04']] - r[['B03']]))) / sqrt((2 * r[['B08']] + 1) ^ 2 - (6 * r[['B08']] - 5 * sqrt(r[['B04']])) - 0.5)
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'MCARI'){
        #  Hermann et al. (2010)
        indices[['MCARI']] <- ((r[['B05']] - r[['B04']]) - 0.2 * (r[['B05']] - r[['B03']])) * (r[['B05']] / r[['B04']])
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'MCARI1'){
        # Haboudane et al. (2004)
        indices[['MCARI1']] <- 1.5 * (2.5 * (r[['B08']] - r[['B04']]) - 1.3 * (r[['B08']] - r[['B03']]))
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'MCARI2'){
        # Haboudane et al. (2004)  
        indices[['MCARI2']] <- (1.5 * (2.5 * (r[['B08']] - r[['B04']]) - 1.3 * (r[['B08']] - r[['B03']]))) / sqrt((2 * r[['B08']] + 1) ^ 2 - (6 * r[['B08']] - 5 * sqrt(r[['B04']])) - 0.5)
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'EVI'){
        # Huete et al. (2002)    
        indices[['EVI']] <- 2.5 * (r[['B08']] - r[['B04']]) / (r[['B08']] + 6 * r[['B04']] - 7.5 * r[['B02']] + 1)
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'GM1'){
        # Gitelson and Merzlyak (1997)
        indices[['GM1']] <- r[['B06']] / r[['B03']]
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'GM2'){
        # Gitelson and Merzlyak (1997)
        indices[['GM2']] <- r[['B06']] / r[['B05']]
      }
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'TCARI'){
        # Haboudane et al. (2002)
        indices[['TCARI']] <- 3 * ((r[['B05']] - r[['B04']]) - 0.2 * (r[['B05']] - r[['B03']]) * (r[['B05']] / r[['B04']]))
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'TCARI_OSAVI'){
        # Haboudane et al. (2002)
        iTCARI <- 3 * ((r[['B05']] - r[['B04']]) - 0.2 * (r[['B05']] - r[['B03']]) * (r[['B05']] / r[['B04']]))
        iOSAVI <- ((1 + 0.16) * (r[['B08']] - r[['B04']]) / (r[['B08']] + r[['B04']] + 0.16))
        indices[['TCARI_OSAVI']] <- iTCARI/iOSAVI
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'TVI'){
        # Broge and Leblanc (2000)
        indices[['TVI']] <- 0.5 * (120 * (r[['B06']] - r[['B03']]) - 200 * (r[['B04']] - r[['B03']]))
      }
     
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'SIPI'){
        # Pe??uelas et al. (1995)
        indices[['SIPI']] <- (r[['B08']] - r[['B01']]) / (r[['B08']] + r[['B04']])
      }
       
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'ARI'){
        #Anthocyanin reflectance index
        indices[['ARI']] <- (1/r[['B03']])- (1/r[['B05']]) 
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'ARI2'){
        indices[['ARI2']] <- (r[['B08']]/r[['B02']])- (r[['B08']]/r[['B03']]) 
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'BAI'){
        indices[['BAI']]  <- (1/((0.1-r[['B04']])**2+(0.06-r[['B08']])**2))
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'BAIS2'){
        indices[['BAIS2']]  <-  (1-((r[['B06']]*r[['B07']]*r[['B8A']])/r[['B04']])**0.5)*((r[['B12']]-r[['B8A']])/((r[['B12']]+r[['B8A']])**0.5)+1)
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'GNDVI'){
        indices[['GNDVI']] <- (r[['B08']] - r[['B03']]) / (r[['B08']] + r[['B03']])
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'CIg'){ 
        indices[['CIg']] <- r[['B08']] / r[['B03']] -1 
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'ARVI'){
          y = 0.069;
        indices[['ARVI']] <- (r[['B8A']] - r[['B08']] - y * (r[['B04']] -r[['B02']]) ) / (r[['B8A']] + r[['B04']] - y * (r[['B04']] -r[['B02']]) )
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'AVI'){
        indices[['AVI']] <- 2.0 * r[['B09']] - r[['B04']]
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'ARV2'){
        #Atmospherically Resistant Vegetation Index 2  (abbrv. ARVI2)
        indices[['ARV2']] <- -0.18 + 1.17 *(r[['B08']] - r[['B04']]) / (r[['B08']] + r[['B04']])
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'NBR'){
        #Normalized Difference NIR/SWIR Normalized Burn Ratio (abbrv. NBR)
        indices[['NBR']] <- (r[['B08']] - r[['B12']]) / (r[['B08']] + r[['B12']])
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'NBR-2'){
       indices[['NBR-2']] <- (r[['B11']] - r[['B12']]) / (r[['B11']] + r[['B12']])
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'NDRE'){
        #Normalized Difference NIR/Rededge Normalized Difference Red-Edge (abbrv. NDRE)
        indices[['NDRE']] <- (r[['B08']] - r[['B05']]) / (r[['B08']] + r[['B05']])
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'MNDVI'){
        #Normalized Difference NIR/MIR Modified Normalized Difference Vegetation Index (abbrv. MNDVI)
        indices[['MNDVI']] <- (r[['B08']] - r[['B11']]) / (r[['B08']] + r[['B11']])
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'RedEg1'){
        #Red edge 1  (abbrv. Rededge1)
        indices[['RedEg1']] <- r[['B05']] / r[['B04']]
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'RedEg2'){
        indices[['RedEg2']] <- (r[['B05']] - r[['B04']]) / (r[['B05']] + r[['B04']])
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'WDRVI'){
        #Wide Dynamic Range Vegetation Index  (abbrv. WDRVI)
        indices[['WDRVI']] <- (0.1 * r[['B08']] - r[['B04']]) / (0.1 * r[['B08']] + r[['B04']])
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'NDWI'){
        #Normalized Difference Water Index
        indices[['NDWI']] <- (r[['B8A']] - r[['B11']]) / (r[['B8A']] + r[['B11']])
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'NDWI2'){
        indices[['NDWI2']] <- (r[['B8A']] - r[['B12']]) / (r[['B8A']] + r[['B12']])
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'CR_SWIR'){
        #CR_SWIR from J.B.Feret
        indices[['CR_SWIR']] <- r[['B11']]/(r[['B8A']]+(Bands[['B11']]-Bands[['B8A']])*(r[['B12']]-r[['B8A']])/(Bands[['B12']]-Bands[['B8A']]))
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'CIre'){
        #CIre
        indices[['CIre']] <- (r[['B07']] / r[['B05']])-1
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'CIgreen'){
        indices[['CIgreen']] <- (r[['B08']] / r[['B03']])-1
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'IRECI'){
        indices[['IRECI']] <- (r[['B07']] - r[['B04']]) / (r[['B05']] / r[['B06']])
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'S2REP'){
        indices[['S2REP']] <- 700 + 35*( ((r[['B07']] - r[['B04']]/2) -  r[['B05']])/ (r[['B06']] - r[['B05']]))
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'RVI'){
        indices[['RVI']] <-  (r[['B08']] / r[['B04']])
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'PVI'){
        #Perpendicular Vegetation Index 
        #Initialize parameters
        a = 0.149
        ar = 0.374
        b = 0.735
        indices[['PVI']] <-  (1.0 /sqrt(a** 2.0+ 1.0)) * (r[['B08']] - ar - b)
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'REIP1'){
        #Red-Edge Inflection Point 1  (abbrv. REIP1)
        indices[['REIP1']] <- 700 + 405 * ( ((r[['B04']] - r[['B07']]/2) -  r[['B05']])/ (r[['B06']] - r[['B05']]))
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'REIP2'){
        indices[['REIP2']] <- 700 + 405 * ( ((r[['B04']] - r[['B07']]/2) -  r[['B05']])/ (r[['B06']] - r[['B05']]))
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'Greeness'){
        indices[['Greeness']] <- r[['B03']] / r[['B04']]
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'Redness'){
        # Gitelson et al. (2000)
        indices[['Redness']] <- r[['B05']] / r[['B04']]
      }
      
      if (SpectraltoCompute == 'All' | SpectraltoCompute == 'PSSRa'){
        # Blackburn (1998)    
        indices[['PSSRa']] <- r[['B08']] / r[['B04']]
      }
        

      ifelse(!dir.exists(path.export), dir.create(path.export), FALSE)
      # NB (raster->terra migration): raster::writeRaster(..., bylayer=TRUE)
      # wrote one file per layer from a single call using a vector of
      # filenames -- terra::writeRaster() has no bylayer= equivalent, so each
      # layer is written individually here instead. Filenames also need an
      # explicit .tif extension now: terra infers format from it rather than
      # from a separate format= argument.
      indices_stack <- terra::rast(indices)
      out_filenames <- paste(path.export,'/',names(indices),'-',names_files[i],'.tif',sep='')
      for (lyr in seq_len(terra::nlyr(indices_stack))) {
        terra::writeRaster(indices_stack[[lyr]], filename = out_filenames[lyr], filetype='GTiff', overwrite=TRUE)
      }
      
  } # end for files
  

  m_final<-message('indices files were generated sucessfully')
  return('process done!')
  
}