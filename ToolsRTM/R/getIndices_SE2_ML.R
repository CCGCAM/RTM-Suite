
#' extract spectral indices at Sentinel-2 resolution a sort number of indices
#'
#' @param df a dataframe with reflectance where each rows correspond with an spectrum
#' @param sensor Sensor options: 'Sentinel-2a', or 'Sentinel-2b'
#' @param df.data  dataset with IDs that corresponde with each spectrum, is null is also enable
#' @param fast.process  when the bands are ordered for SE2, please use fast.process = T, otherwise use False or nothing
#'
#'
#' @return a dataframe with indices and your dataset
#' @export
#'
#'
#'

getIndicesSE2.ML <- function(df,sensor='Sentinel-2a', df.data=NULL, fast.process=NULL) {

  if (class(data)[1] == "matrix"){
    df<-as.data.frame(df)
  }

  #wave.avalaible = sort(as.numeric(gsub(".*?([0-9]+).*", "\\1", wavelengths)))

  S2.provided.bands <- sort(names(df))
  S2.sort<-c('B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B11','B12')
  S2.provided.bands<-S2.provided.bands[order(match(S2.provided.bands,S2.sort))]
  #print(S2.provided.bands)
  ##original conf form J.B Feret for CR.SWIR
  bandset.SE2a <- c('B1'=442.7,'B2'=492.7, 'B3'=559.8, 'B4'=664.6, 'B5'=704.1, 'B6'=740.5,
                 'B7' = 782.8, 'B8' =832.8,'B8A' = 864.7,
                 'B9'= 945.1,'B10'=1373.5,'B11' = 1613.7 ,'B12' = 2202.4)

  bandset.SE2b <- c('B1'=442.2,'B2'=492.3, 'B3'=558.9, 'B4'=664.9, 'B5'=703.8, 'B6'=739.1,
                 'B7' =779.7, 'B8' =  832.9, 'B8A' = 864.0,
                 'B9'= 943.2,'B10'=1376.9,'B11'= 1610.4,'B12' =2185.7 )


  if (is.null(sensor) == TRUE){
    bandset.SE<-  bandset.SE2a
  } else if (sensor == 'Sentinel-2a') {
    bandset.SE<-  bandset.SE2a
  } else if (sensor == 'Sentinel-2b'){
    bandset.SE<-  bandset.SE2b
  }


  indices.list = list()
  # create progress bar
  total=dim(df)[1]

 barProgress <- txtProgressBar(min = 1, max = total, style = 3)

for (i in c(1:dim(df)[1])){

  values<-df[i,]
  r <- as.numeric(values[1,])
  names(r)<-S2.provided.bands
  list.r<-list(bandset.SE,r)

  if (fast.process == F | missing(fast.process) | is.null(fast.process)){
    # take our list and rbind it into a data.frame, filling in missing values with NA
    list.r<-plyr::ldply(list.r , rbind)
    r <- list.r[2,]
  } else {
    r <- r
  }


  indices = c()

  #if(header) indices['Structural'] <- 'Structural'

  # NDVI (R800-R670)/(R800+R670)
  # Rouse et al. (1974)
  indices['NDVI'] <- (r['B8'] - r['B4']) / (r['B8'] + r['B4'])
  # OSAVI [(1+0.16)*(R800-R670)/(R800+R670+0.16)]
  # Rondeaux et al. (1996)
  indices['OSAVI'] <- ((1 + 0.16) * (r['B8'] - r['B4']) / (r['B8'] + r['B4'] + 0.16))

  # MSAVI 1/2*[(2*R800+1-?((?(2*R800+1)?^2)-8*(R800-R670))]
  # Qi et al. (1994)
  indices['MSAVI'] <- 1 / 2 * (2 * r['B8'] + 1 - sqrt(((2 * r['B8'] + 1) ^ 2) - 8 * (r['B8'] - r['B4'])))

  # MCARI ((RB5-R670) - 0.2*(RB5-R550))*(RB5/R670)
  # Hermann et al. (2010)
  indices['MCARI'] <- ((r['B5'] - r['B4']) - 0.2 * (r['B5'] - r['B3'])) * (r['B5'] / r['B4'])

  # EVI 2.5*(R800-R670)/(R800+6*R670-7.5*R400+1)
  # Huete et al. (2002)
  indices['EVI'] <- 2.5 * (r['B8'] - r['B4']) / (r['B8'] + 6 * r['B4'] - 7.5 * r['B2'] + 1)

  # TCARI 3*[(RB5-R670)-0.2*(RB5-R550)*(RB5/R670)]
  # Haboudane et al. (2002)
  indices['TCARI'] <- 3 * ((r['B5'] - r['B4']) - 0.2 * (r['B5'] - r['B3']) * (r['B5'] / r['B4']))
  a.factor = 0.496

  # TCARI/OSAVI TCARI/OSAVI
  # Haboudane et al. (2002)
  indices['TCARI_OSAVI'] <- as.numeric(indices['TCARI']) / as.numeric(indices['OSAVI'])
  #Anthocyanin reflectance index

  indices['Datt1'] <- (r['B8'] - r['B5']) / (r['B8'] + r['B4'])

  # GM1 R750/R550
  # Gitelson and Merzlyak (1997)
  indices['GM1'] <- r['B6'] / r['B3']

  #Normalized Difference NIR/MIR Modified Normalized Difference Vegetation Index (abbrv. MNDVI)
  indices['MNDVI'] <- (r['B8'] - r['B11']) / (r['B8'] + r['B11'])

  #Normalized Difference NIR/Rededge Normalized Difference Red-Edge (abbrv. NDRE)
  indices['NDRE'] <- (r['B8'] - r['B5']) / (r['B8'] + r['B5'])

  #Red edge 1  (abbrv. Rededge1)
  indices['RedEg1'] <- r['B5'] / r['B4']

  #Wide Dynamic Range Vegetation Index  (abbrv. WDRVI)
  indices['WDRVI'] <- (0.1 * r['B8'] - r['B4']) / (0.1 * r['B8'] + r['B4'])
  #Normalized Difference Water Index
  indices['NDWI'] <- (r['B8A'] - r['B11']) / (r['B8A'] + r['B11'])
  indices['NDWI2'] <- (r['B8A'] - r['B12']) / (r['B8A'] + r['B12'])

  #CR_SWIR from J.B.Feret
  indices['CR.SWIR'] <- r['B11']/(r['B8']+(bandset.SE['B11']-bandset.SE['B8'])*(r['B12']-r['B8'])/(bandset.SE['B12']-bandset.SE['B8']))
  # PSSRa R800/R680
  # Blackburn (1998)
  indices['PSSRa'] <- r['B8'] / r['B4']

  #CIre
  indices['CIre'] <- (r['B7'] / r['B5'])-1
  indices['CIgreen'] <- (r['B8'] / r['B3'])-1
  indices['IRECI'] <- (r['B7'] - r['B4']) / (r['B5'] / r['B6'])

  #### Using Band 8
  indices['CR.red.nir.1'] <- r['B6']  / (r['B8'] + (bandset.SE['B6'] - bandset.SE['B8']) * (r['B7'] - r['B8']) / (bandset.SE['B7'] - bandset.SE['B8']))
  indices['CR.red.nir.6'] <- r['B7']  / (r['B8'] + (bandset.SE['B7'] - bandset.SE['B8']) * (r['B5'] - r['B8']) / (bandset.SE['B5'] - bandset.SE['B8']))


  indices['SBI'] <- 0.3037 * r['B2'] + 0.2793 * r['B3'] + 0.4743 * r['B4'] + 0.5585 * r['B8'] + 0.5082 * r['B11'] + 0.1863 * r['B2']
  indices['GVI'] <- -0.2848 * r['B2'] - 0.2435 * r['B3']  - 0.5436 * r['B4'] + 0.7243 * r['B8'] + 0.0840 * r['B11'] - 0.1800 * r['B12']
  indices['WET'] <- 0.1509 * r['B2'] + 0.1973 * r['B3'] + 0.3279 * r['B4'] + 0.3406 * r['B8A'] - 0.7112 * r['B11'] - 0.4572 * r['B12']

  indices['BF.Anth'] <- (r['B3'] - r['B2']) / (r['B3'] + r['B2'])
  indices['CR.red.nir']  <- r['B7']  / (r['B8'] + ( (r['B6'] - r['B8']) / (740.5- 833) ) * (782.8 - 832.8))
  indices['CR.Brown']  <- r['B3']  / (r['B7'] + ( (r['B2'] - r['B7']) / (492.4- 782.8) ) * (559.8 - 782.8))
  # Calculate NDVIv
  indices['NDVIv'] <- ((r['B8'] - r['B4']) / (r['B8'] + r['B4'])) * r['B8']
  sigma <- 0.15
  knr <- exp(-(r['B8'] - r['B4'])^2/(2*sigma^2))
  indices['kNDVI'] <- (1-knr) / (1+knr)

  indices.list[[i]] = indices
  setTxtProgressBar(barProgress, i)

}

  df.indices <- data.frame(matrix(unlist(indices.list), nrow=length(indices.list), byrow=T))
  colnames(df.indices)<-names(indices)
  if (is.null(df.data)){
    ## remove indices wih no data
    df.indices[sapply(df.indices, is.infinite)] <- NA
    df.indices = df.indices[, colSums(is.na(df.indices)) != nrow(df.indices)]
    #df.indices <- df.indices[!is.infinite(colSums(df.indices)),]


    df.indices_<-df.indices


  } else{
    ## remove indices wih no data
    df.indices[sapply(df.indices, is.infinite)] <- NA
    df.indices = df.indices[, colSums(is.na(df.indices)) != nrow(df.indices)]
    #df.indices <- df.indices[!is.infinite(colSums(df.indices)),]
    df.indices_<-cbind(df.data,df.indices)

  }

  close(barProgress)
  return(df.indices_)
}
