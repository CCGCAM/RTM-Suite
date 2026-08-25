#'  this function produces Spatial maps of the spectral index
#'
#' @param rasterFiles path with the image (in RasterBrick/RasterStack formats)
#' @param Sensor character with the sensor 'Sentinel2a'
#' @param factorR numeric. multiplying factor used to write reflectance in image (==10000 for S2)
#' @param SpectraltoCompute List with the spectral index to compute. by default 'All'
#' 
#' @return a spatial map.
#' @export
#' @examples here adding examples ....
#' 
#' 
getSpatial_index<-function(rasterFiles=NULL,Sensor='Sentinel2a',SpectraltoCompute= NULL,
                            factorR=NULL){
  options(warn=-1) ###avoid warnings
  
  if (is.null(Sensor) | Sensor == 'Sentinel2a'){
    Bands <- data.frame('B01'=442.7,'B02'=492.4, 'B03'=559.8, 'B04'=664.6, 'B05'=704.1, 'B06'=740.5,
                          'B07' = 782.8, 'B08' = 832.8, 'B8A' = 864.7, 'B09' = 945.1, 'B11' = 1613.7, 'B12' = 2202.4)
  } else{
    message('not implemented yet')
    stop()
  } 
  if (is.null(SpectraltoCompute)) {
    message('please specific one spectral indicator: see ?getSpectraindex')
    message('or use All to estimate all index')
    stop()
  }
  
  if (is.null(factorR)) {
    factorR=1/10000
    message('factor scale used: factorR=1/10000')
  }

    r<- brick(rasterFiles) * factorR
    names(r)<-names(Bands)
    index<-list()
    
    if (SpectraltoCompute == 'B01'){
        index[['B01']] <- r[['B01']]
    }  
    if (SpectraltoCompute == 'B02'){
      index[['B02']] <- r[['B02']]
    }  
    if (SpectraltoCompute == 'B03'){
      index[['B03']] <- r[['B03']]
    }  
    
    if (SpectraltoCompute == 'B04'){
      index[['B04']] <- r[['B04']]
    }  
    if (SpectraltoCompute == 'B05'){
      index[['B05']] <- r[['B05']]
    }  
    if (SpectraltoCompute == 'B06'){
      index[['B06']] <- r[['B06']]
    }  
    if (SpectraltoCompute == 'B07'){
      index[['B07']] <- r[['B07']]
    }  
    if (SpectraltoCompute == 'B08'){
      index[['B08']] <- r[['B08']]
    }  
    if (SpectraltoCompute == 'B8A'){
      index[['B8A']] <- r[['B8A']]
    }  
    if (SpectraltoCompute == 'B09'){
      index[['B09']] <- r[['B09']]
    }  
    if (SpectraltoCompute == 'B11'){
      index[['B11']] <- r[['B11']]
    }  
    if (SpectraltoCompute == 'B12'){
      index[['B12']] <- r[['B12']]
    }  

    if (SpectraltoCompute == 'NDVI'){
      # Rouse et al. (1974)
      index[['NDVI']] <- (r[['B08']] - r[['B04']]) / (r[['B08']] + r[['B04']])
    }
    
    if (SpectraltoCompute == 'RDVI'){
      # Rougean and Breon (1995)
      index[['RDVI']] <- (r[['B08']] - r[['B04']]) / (r[['B08']] + r[['B04']]) ** 0.5
    }
  
    if (SpectraltoCompute == 'SR'){
      index[['SR']] <- r[['B08']] / r[['B04']]
    }
    
    if (SpectraltoCompute == 'MSR'){
      # Chen (1996)
      index[['MSR']] <- (r[['B08']] / r[['B04']] - 1) / ((r[['B08']] / r[['B04']]) ** 0.5 + 1)
    }
    
    if (SpectraltoCompute == 'OSAVI'){
      # Rondeaux et al. (1996)
      index[['OSAVI']] <- ((1 + 0.16) * (r[['B08']] - r[['B04']]) / (r[['B08']] + r[['B04']] + 0.16))
    }
    
    if (SpectraltoCompute == 'MSAVI'){
      # Qi et al. (1994)    
      index[['MSAVI']] <- 1 / 2 * (2 * r[['B08']] + 1 - sqrt(((2 * r[['B08']] + 1) ^ 2) - 8 * (r[['B08']] - r[['B04']])))
    }
    
    if (SpectraltoCompute == 'MTVI1'){
      # Broge & Leblanc (2000); Haboudane et al. (2004)
      index[['MTVI1']] <- 1.2 * (1.2 * (r[['B08']] - r[['B03']]) - 2.5 * (r[['B04']] - r[['B03']]))
    }
    
    if (SpectraltoCompute == 'MTVI2'){
      # Haboudane et al. (2004)    
      index[['MTVI2']] <- (1.5 * (1.2 * (r[['B08']] - r[['B03']]) - 2.5 * (r[['B04']] - r[['B03']]))) / sqrt((2 * r[['B08']] + 1) ^ 2 - (6 * r[['B08']] - 5 * sqrt(r[['B04']])) - 0.5)
    }
    
    if (SpectraltoCompute == 'MCARI'){
      #  Hermann et al. (2010)
      index[['MCARI']] <- ((r[['B05']] - r[['B04']]) - 0.2 * (r[['B05']] - r[['B03']])) * (r[['B05']] / r[['B04']])
    }
    
    if (SpectraltoCompute == 'MCARI1'){
      # Haboudane et al. (2004)
      index[['MCARI1']] <- 1.5 * (2.5 * (r[['B08']] - r[['B04']]) - 1.3 * (r[['B08']] - r[['B03']]))
    }
    
    if (SpectraltoCompute == 'MCARI2'){
      # Haboudane et al. (2004)  
      index[['MCARI2']] <- (1.5 * (2.5 * (r[['B08']] - r[['B04']]) - 1.3 * (r[['B08']] - r[['B03']]))) / sqrt((2 * r[['B08']] + 1) ^ 2 - (6 * r[['B08']] - 5 * sqrt(r[['B04']])) - 0.5)
    }
    
    if (SpectraltoCompute == 'EVI'){
      # Huete et al. (2002)    
      index[['EVI']] <- 2.5 * (r[['B08']] - r[['B04']]) / (r[['B08']] + 6 * r[['B04']] - 7.5 * r[['B02']] + 1)
    }
    
    if (SpectraltoCompute == 'GM1'){
      # Gitelson and Merzlyak (1997)
      index[['GM1']] <- r[['B06']] / r[['B03']]
    }
    
    if (SpectraltoCompute == 'GM2'){
      # Gitelson and Merzlyak (1997)
      index[['GM2']] <- r[['B06']] / r[['B05']]
    }
    if (SpectraltoCompute == 'TCARI'){
      # Haboudane et al. (2002)
      index[['TCARI']] <- 3 * ((r[['B05']] - r[['B04']]) - 0.2 * (r[['B05']] - r[['B03']]) * (r[['B05']] / r[['B04']]))
    }
    
    if (SpectraltoCompute == 'TCARI_OSAVI'){
      # Haboudane et al. (2002)
      iTCARI <- 3 * ((r[['B05']] - r[['B04']]) - 0.2 * (r[['B05']] - r[['B03']]) * (r[['B05']] / r[['B04']]))
      iOSAVI <- ((1 + 0.16) * (r[['B08']] - r[['B04']]) / (r[['B08']] + r[['B04']] + 0.16))
      index[['TCARI_OSAVI']] <- iTCARI/iOSAVI
    }
    
    if (SpectraltoCompute == 'TVI'){
      # Broge and Leblanc (2000)
      index[['TVI']] <- 0.5 * (120 * (r[['B06']] - r[['B03']]) - 200 * (r[['B04']] - r[['B03']]))
    }
   
    if (SpectraltoCompute == 'SIPI'){
      # Pe??uelas et al. (1995)
      index[['SIPI']] <- (r[['B08']] - r[['B01']]) / (r[['B08']] + r[['B04']])
    }
     
    if (SpectraltoCompute == 'ARI'){
      #Anthocyanin reflectance index
      index[['ARI']] <- (1/r[['B03']])- (1/r[['B05']]) 
    }
    
    if (SpectraltoCompute == 'ARI2'){
      index[['ARI2']] <- (r[['B08']]/r[['B02']])- (r[['B08']]/r[['B03']]) 
    }
    
    if (SpectraltoCompute == 'BAI'){
      index[['BAI']]  <- (1/((0.1-r[['B04']])**2+(0.06-r[['B08']])**2))
    }
    
    if (SpectraltoCompute == 'BAIS2'){
      index[['BAIS2']]  <-  (1-((r[['B06']]*r[['B07']]*r[['B8A']])/r[['B04']])**0.5)*((r[['B12']]-r[['B8A']])/((r[['B12']]+r[['B8A']])**0.5)+1)
    }
    
    if (SpectraltoCompute == 'GNDVI'){
      index[['GNDVI']] <- (r[['B08']] - r[['B03']]) / (r[['B08']] + r[['B03']])
    }
    
    if (SpectraltoCompute == 'CIg'){ 
      index[['CIg']] <- r[['B08']] / r[['B03']] -1 
    }
    
    if (SpectraltoCompute == 'ARVI'){
        y = 0.069;
      index[['ARVI']] <- (r[['B8A']] - r[['B08']] - y * (r[['B04']] -r[['B02']]) ) / (r[['B8A']] + r[['B04']] - y * (r[['B04']] -r[['B02']]) )
    }
    
    if (SpectraltoCompute == 'AVI'){
      index[['AVI']] <- 2.0 * r[['B09']] - r[['B04']]
    }
    
    if (SpectraltoCompute == 'ARV2'){
      #Atmospherically Resistant Vegetation Index 2  (abbrv. ARVI2)
      index[['ARV2']] <- -0.18 + 1.17 *(r[['B08']] - r[['B04']]) / (r[['B08']] + r[['B04']])
    }
    
    if (SpectraltoCompute == 'NBR'){
      #Normalized Difference NIR/SWIR Normalized Burn Ratio (abbrv. NBR)
      index[['NBR']] <- (r[['B08']] - r[['B12']]) / (r[['B08']] + r[['B12']])
    }
    
    if (SpectraltoCompute == 'NBR-2'){
     index[['NBR-2']] <- (r[['B11']] - r[['B12']]) / (r[['B11']] + r[['B12']])
    }
    
    if (SpectraltoCompute == 'NDRE'){
      #Normalized Difference NIR/Rededge Normalized Difference Red-Edge (abbrv. NDRE)
      index[['NDRE']] <- (r[['B08']] - r[['B05']]) / (r[['B08']] + r[['B05']])
    }
    
    if (SpectraltoCompute == 'MNDVI'){
      #Normalized Difference NIR/MIR Modified Normalized Difference Vegetation Index (abbrv. MNDVI)
      index[['MNDVI']] <- (r[['B08']] - r[['B11']]) / (r[['B08']] + r[['B11']])
    }
    
    if (SpectraltoCompute == 'RedEg1'){
      #Red edge 1  (abbrv. Rededge1)
      index[['RedEg1']] <- r[['B05']] / r[['B04']]
    }
    
    if (SpectraltoCompute == 'RedEg2'){
      index[['RedEg2']] <- (r[['B05']] - r[['B04']]) / (r[['B05']] + r[['B04']])
    }
    
    if (SpectraltoCompute == 'WDRVI'){
      #Wide Dynamic Range Vegetation Index  (abbrv. WDRVI)
      index[['WDRVI']] <- (0.1 * r[['B08']] - r[['B04']]) / (0.1 * r[['B08']] + r[['B04']])
    }
    
    if (SpectraltoCompute == 'NDWI'){
      #Normalized Difference Water Index
      index[['NDWI']] <- (r[['B8A']] - r[['B11']]) / (r[['B8A']] + r[['B11']])
    }
    
    if (SpectraltoCompute == 'NDWI2'){
      index[['NDWI2']] <- (r[['B8A']] - r[['B12']]) / (r[['B8A']] + r[['B12']])
    }
    
    if (SpectraltoCompute == 'CR_SWIR'){
      #CR_SWIR from J.B.Feret
      index[['CR_SWIR']] <- r[['B11']]/(r[['B8A']]+(Bands[['B11']]-Bands[['B8A']])*(r[['B12']]-r[['B8A']])/(Bands[['B12']]-Bands[['B8A']]))
    }
    
    if (SpectraltoCompute == 'CIre'){
      #CIre
      index[['CIre']] <- (r[['B07']] / r[['B05']])-1
    }
    
    if (SpectraltoCompute == 'CIgreen'){
      index[['CIgreen']] <- (r[['B08']] / r[['B03']])-1
    }
    
    if (SpectraltoCompute == 'IRECI'){
      index[['IRECI']] <- (r[['B07']] - r[['B04']]) / (r[['B05']] / r[['B06']])
    }
    
    if (SpectraltoCompute == 'S2REP'){
      index[['S2REP']] <- 700 + 35*( ((r[['B07']] - r[['B04']]/2) -  r[['B05']])/ (r[['B06']] - r[['B05']]))
    }
    
    if (SpectraltoCompute == 'RVI'){
      index[['RVI']] <-  (r[['B08']] / r[['B04']])
    }
    
    if (SpectraltoCompute == 'PVI'){
      #Perpendicular Vegetation Index 
      #Initialize parameters
      a = 0.149
      ar = 0.374
      b = 0.735
      index[['PVI']] <-  (1.0 /sqrt(a** 2.0+ 1.0)) * (r[['B08']] - ar - b)
    }
    
    if (SpectraltoCompute == 'REIP1'){
      #Red-Edge Inflection Point 1  (abbrv. REIP1)
      index[['REIP1']] <- 700 + 405 * ( ((r[['B04']] - r[['B07']]/2) -  r[['B05']])/ (r[['B06']] - r[['B05']]))
    }
    
    if (SpectraltoCompute == 'REIP2'){
      index[['REIP2']] <- 700 + 405 * ( ((r[['B04']] - r[['B07']]/2) -  r[['B05']])/ (r[['B06']] - r[['B05']]))
    }
    
    if (SpectraltoCompute == 'Greeness'){
      index[['Greeness']] <- r[['B03']] / r[['B04']]
    }
    
    if (SpectraltoCompute == 'Redness'){
      # Gitelson et al. (2000)
      index[['Redness']] <- r[['B05']] / r[['B04']]
    }
    
    if (SpectraltoCompute == 'PSSRa'){
      # Blackburn (1998)    
      index[['PSSRa']] <- r[['B08']] / r[['B04']]
    }
      
  
  return(index)
  
}
