
#' extract spectral indices at Sentinel-2 resolution
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

getIndicesSE2 <- function(df,sensor='Sentinel-2a', df.data=NULL, fast.process=NULL) {

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

  # RDVI (R800-R670)/(R800+R670)^0.5
  # Rougean and Breon (1995)
  indices['RDVI'] <- (r['B8'] - r['B4']) / (r['B8'] + r['B4']) ** 0.5

  # SR R800/R670
  # Jordan (1969)
  indices['SR'] <- r['B8'] / r['B4']

  # MSR (R800/R670-1)/((R800/R670)^0.5+1)
  # Chen (1996)
  indices['MSR'] <- (r['B8'] / r['B4'] - 1) / ((r['B8'] / r['B4']) ** 0.5 + 1)

  # OSAVI [(1+0.16)*(R800-R670)/(R800+R670+0.16)]
  # Rondeaux et al. (1996)
  indices['OSAVI'] <- ((1 + 0.16) * (r['B8'] - r['B4']) / (r['B8'] + r['B4'] + 0.16))

  # MSAVI 1/2*[(2*R800+1-?((?(2*R800+1)?^2)-8*(R800-R670))]
  # Qi et al. (1994)
  indices['MSAVI'] <- 1 / 2 * (2 * r['B8'] + 1 - sqrt(((2 * r['B8'] + 1) ^ 2) - 8 * (r['B8'] - r['B4'])))

  # MTVI1 1.2*[1.2*(R800-R550)-2.5*(R670-R550)]
  # Broge & Leblanc (2000); Haboudane et al. (2004)
  indices['MTVI1'] <- 1.2 * (1.2 * (r['B8'] - r['B3']) - 2.5 * (r['B4'] - r['B3']))

  # MTVI2 (1.5*[1.2*(R800-R550)-2.5*(R670-R550)])/SQR((2*R800+1)^2-(6*R800-5*SQR(R670))-0.5)
  # Haboudane et al. (2004)
  indices['MTVI2'] <- (1.5 * (1.2 * (r['B8'] - r['B3']) - 2.5 * (r['B4'] - r['B3']))) / sqrt((2 * r['B8'] + 1) ^ 2 - (6 * r['B8'] - 5 * sqrt(r['B4'])) - 0.5)

  # MCARI ((RB5-R670) - 0.2*(RB5-R550))*(RB5/R670)
  # Hermann et al. (2010)
  indices['MCARI'] <- ((r['B5'] - r['B4']) - 0.2 * (r['B5'] - r['B3'])) * (r['B5'] / r['B4'])

  # MCARI1 1.2* [2.5* (R800-  R670)-  1.3* (R800 -R550) ]
  # Haboudane et al. (2004)
  indices['MCARI1'] <- 1.5 * (2.5 * (r['B8'] - r['B4']) - 1.3 * (r['B8'] - r['B3']))

  # MCARI2 (1.5*[2.5*(R800-R670)-1.3*(R800-R550) ])/SQRT((2*R800+1)^2-(6*R800-5*SQRT(R670))-0.5)
  # Haboudane et al. (2004)
  indices['MCARI2'] <- (1.5 * (2.5 * (r['B8'] - r['B4']) - 1.3 * (r['B8'] - r['B3']))) / sqrt((2 * r['B8'] + 1) ^ 2 - (6 * r['B8'] - 5 * sqrt(r['B4'])) - 0.5)

  # EVI 2.5*(R800-R670)/(R800+6*R670-7.5*R400+1)
  # Huete et al. (2002)
  indices['EVI'] <- 2.5 * (r['B8'] - r['B4']) / (r['B8'] + 6 * r['B4'] - 7.5 * r['B2'] + 1)

  ## Pigmentos
  #if(header) indices['Pigments'] <- 'Pigments'


  # GM1 R750/R550
  # Gitelson and Merzlyak (1997)
  indices['GM1'] <- r['B6'] / r['B3']

  # GM2 R750/RB5
  # Gitelson and Merzlyak (1997)
  indices['GM2'] <- r['B6'] / r['B5']

  # TCARI 3*[(RB5-R670)-0.2*(RB5-R550)*(RB5/R670)]
  # Haboudane et al. (2002)
  indices['TCARI'] <- 3 * ((r['B5'] - r['B4']) - 0.2 * (r['B5'] - r['B3']) * (r['B5'] / r['B4']))
  a.factor = 0.496
  indices['CARI2']  <-  abs((r['B5'] - r['B3']) / 150.0 * r['B4'] + r['B4'] + r['B3'] - (a.factor * r['B3'])) / ((a.factor ** 2) + 1.0) ** 0.5 * (r['B5'] / r['B4'])

  # TCARI/OSAVI TCARI/OSAVI
  # Haboudane et al. (2002)
  indices['TCARI_OSAVI'] <- as.numeric(indices['TCARI']) / as.numeric(indices['OSAVI'])

  # TVI 0.5*[120*(R750-R550)-200*(R670-R550) ]
  # Broge and Leblanc (2000)
  indices['TVI'] <- 0.5 * (120 * (r['B6'] - r['B3']) - 200 * (r['B4'] - r['B3']))
  # SRPI R430/R680

  #if(any(bandset.SE <= 445)){

  # SIPI (R800-R445)/(R800+R680)
  # Pe??uelas et al. (1995)
  indices['SIPI'] <- (r['B8'] - r['B1']) / (r['B8'] + r['B4'])


  ## SEntinel 2a
  #if(header) indices['SE2a'] <- 'SE2a Indices'

  #Anthocyanin reflectance index

  indices['Datt1'] <- (r['B8'] - r['B5']) / (r['B8'] + r['B4'])
  indices['Datt4'] <- (r['B4']) / (r['B3'] * r['B8'])
  indices['Datt6'] <- (r['B8A']) / (r['B3'] * r['B5'])
  # Canopy Chlorophyll Content Index  (abbrv. CCCI)
  indices['CCCI'] <- ((r['B8'] - r['B5']) / (r['B8'] + r['B5'])) / ((r['B8'] - r['B4'])/ (r['B8'] + r['B4']))

  indices['ARI'] <- (1/r['B3'])- (1/r['B5'])
  indices['GNDVI'] <- (r['B8'] - r['B3']) / (r['B8'] + r['B3'])
  indices['CIg'] <- r['B8'] / r['B3'] -1

  y = 0.069;
  indices['ARVI'] <- (r['B8A'] - r['B8'] - y * (r['B4'] -r['B2']) ) / (r['B8A'] + r['B4'] - y * (r['B4'] -r['B2']) )

  #if (any(bandset.SE == 945) | any(bandset.SE == 943)) {
  indices['AVI'] <- 2.0 * r['B9'] - r['B4']
  #}
  #Atmospherically Resistant Vegetation Index 2  (abbrv. ARVI2)
  indices['ARV2'] <- -0.18 + 1.17 *(r['B8'] - r['B4']) / (r['B8'] + r['B4'])
  #Normalized Difference NIR/SWIR Normalized Burn Ratio (abbrv. NBR)
  indices['NBR'] <- (r['B8'] - r['B12']) / (r['B8'] + r['B12'])
  indices['NBR.2'] <- (r['B11'] - r['B12']) / (r['B11'] + r['B12'])

  #Normalized Difference NIR/Rededge Normalized Difference Red-Edge (abbrv. NDRE)
  indices['NDRE'] <- (r['B8'] - r['B5']) / (r['B8'] + r['B5'])
  #Normalized Difference NIR/MIR Modified Normalized Difference Vegetation Index (abbrv. MNDVI)
  indices['MNDVI'] <- (r['B8'] - r['B11']) / (r['B8'] + r['B11'])

  #Red edge 1  (abbrv. Rededge1)
  indices['RedEg1'] <- r['B5'] / r['B4']
  indices['RedEg2'] <- (r['B5'] - r['B4']) / (r['B5'] + r['B4'])
  #Wide Dynamic Range Vegetation Index  (abbrv. WDRVI)
  indices['WDRVI'] <- (0.1 * r['B8'] - r['B4']) / (0.1 * r['B8'] + r['B4'])
  #Normalized Difference Water Index
  indices['NDWI'] <- (r['B8A'] - r['B11']) / (r['B8A'] + r['B11'])
  indices['NDWI2'] <- (r['B8A'] - r['B12']) / (r['B8A'] + r['B12'])

  #Leaf Water Content Index  (abbrv. LWCI)
  MIDIR = 0.101
  indices['LWCI'] <- log(1.0 -( r['B8'] - MIDIR)) / -log(1.0 * (r['B8'] - MIDIR))
  #CR_SWIR from J.B.Feret
  indices['CR.SWIR'] <- r['B11']/(r['B8A']+(bandset.SE['B11']-bandset.SE['B8A'])*(r['B12']-r['B8A'])/(bandset.SE['B12']-bandset.SE['B8A']))
  indices['CR.SWIR.2'] <- r['B11']/(r['B8']+((r['B12']-r['B8'])/(2185.7-832.8))*(1610.4-832.8))

  #CIre
  indices['CIre'] <- (r['B7'] / r['B5'])-1
  indices['CIrededge'] <- (r['B8'] / r['B5'])-1
  indices['CIgreen'] <- (r['B8'] / r['B3'])-1
  indices['Chlred.edge'] <- (r['B7'] / r['B5']) ** (-1)
  indices['CVI'] <- (r['B8'] * r['B4']) / r['B4']**2
  indices['IRECI'] <- (r['B7'] - r['B4']) / (r['B5'] / r['B6'])
  indices['REP'] <- 700 + 40* (((r['B4'] + r['B7'])/2) -  r['B5'])/ (r['B6'] - r['B5'])
  indices['RVI'] <-  (r['B8'] / r['B4'])
  #Perpendicular Vegetation Index
  #Initialize parameters
  a = 0.149
  ar = 0.374
  b = 0.735
  indices['PVI'] <-  (1.0 /sqrt(a** 2.0+ 1.0)) * (r['B8'] - ar - b)
  #Red-Edge Inflection Point 1  (abbrv. REIP1)
  indices['REIP1'] <- 700.0 + 40.0  * (( ((r['B4'] + r['B7'])/2) -  r['B5'])/ (r['B6'] - r['B5']))
  indices['REIP2'] <- 702.0 + 40.0  * (( ((r['B4'] + r['B7'])/2) -  r['B5'])/ (r['B6'] - r['B5']))
  ## BGR
  #if(header) indices['BGR'] <- 'BGR'

  # G R550/R670
  # -
  indices['Greeness'] <- r['B3'] / r['B4']

  # R RB5/R670
  # Gitelson et al. (2000)
  indices['Redness'] <- r['B5'] / r['B4']

  # RARS R746/R513

  ## NIR-VIS
  #if(header) indices['NIR-VIS'] <- 'NIR-VIS'

  # PSSRa R800/R680
  # Blackburn (1998)
  indices['PSSRa'] <- r['B8'] / r['B4']

  #if(header) indices['Red-edge'] <- 'Red-edge'

  if (is.na(r['B8A']) != TRUE){
    #### Using Band 8A
    indices['CR.red.nir.1']  <- r['B6']  / (r['B8A'] + ( (r['B7'] - r['B8A']) / (bandset.SE['B6'] - bandset.SE['B8A']) ) * (bandset.SE['B7'] - bandset.SE['B8A']))
    indices['CR.red.nir.2']  <- r['B5']  / (r['B7'] + ( (r['B6'] - r['B7']) / (bandset.SE['B6'] - bandset.SE['B7']) ) * (bandset.SE['B6'] - bandset.SE['B7']))
    indices['CR.red.nir.3']  <- r['B5']  / (r['B8A'] + ( (r['B6'] - r['B8A']) / (bandset.SE['B5'] - bandset.SE['B8A']) ) * (bandset.SE['B6'] - bandset.SE['B8A']))
    indices['CR.red.nir.4']  <- r['B4']  / (r['B8A'] + ( (r['B5'] - r['B8A']) / (bandset.SE['B4'] - bandset.SE['B8A']) ) * (bandset.SE['B5'] - bandset.SE['B8A']))
    indices['CR.red.nir.5']  <- r['B6']  / (r['B8A'] + ( (r['B5'] - r['B8A']) / (bandset.SE['B6'] - bandset.SE['B8A']) ) * (bandset.SE['B5'] - bandset.SE['B8A']))
    indices['CR.red.nir.6']  <- r['B7']  / (r['B8A'] + ( (r['B5'] - r['B8A']) / (bandset.SE['B5'] - bandset.SE['B8A']) ) * (bandset.SE['B7'] - bandset.SE['B8A']))


  } else{
    #### Using Band 8
    indices['CR.red.nir.1']  <- r['B6']  / (r['B8'] + ( (r['B7'] - r['B8']) / (bandset.SE['B6'] - bandset.SE['B8']) ) * (bandset.SE['B7'] - bandset.SE['B8']))
    indices['CR.red.nir.2']  <- r['B5']  / (r['B7'] + ( (r['B6'] - r['B7']) / (bandset.SE['B6'] - bandset.SE['B7']) ) * (bandset.SE['B6'] - bandset.SE['B7']))
    indices['CR.red.nir.3']  <- r['B5']  / (r['B8'] + ( (r['B6'] - r['B8']) / (bandset.SE['B5'] - bandset.SE['B8']) ) * (bandset.SE['B6'] - bandset.SE['B8']))
    indices['CR.red.nir.4']  <- r['B4']  / (r['B8'] + ( (r['B5'] - r['B8']) / (bandset.SE['B4'] - bandset.SE['B8']) ) * (bandset.SE['B5'] - bandset.SE['B8']))
    indices['CR.red.nir.5']  <- r['B6']  / (r['B8'] + ( (r['B5'] - r['B8']) / (bandset.SE['B6'] - bandset.SE['B8']) ) * (bandset.SE['B5'] - bandset.SE['B8']))
    indices['CR.red.nir.6']  <- r['B7']  / (r['B8'] + ( (r['B5'] - r['B8']) / (bandset.SE['B5'] - bandset.SE['B8']) ) * (bandset.SE['B7'] - bandset.SE['B8']))

  }


  indices['SBI'] <- 0.3037 * r['B2'] + 0.2793 * r['B3'] + 0.4743 * r['B4'] + 0.5585 * r['B8'] + 0.5082 * r['B11'] + 0.1863 * r['B2']
  indices['GVI'] <- -0.2848 * r['B2'] - 0.2435 * r['B3']  - 0.5436 * r['B4'] + 0.7243 * r['B8'] + 0.0840 * r['B11'] - 0.1800 * r['B12']
  indices['WET'] <- 0.1509 * r['B2'] + 0.1973 * r['B3'] + 0.3279 * r['B4'] + 0.3406 * r['B8A'] - 0.7112 * r['B11'] - 0.4572 * r['B12']

  indices['BF.Anth'] <- (r['B3'] - r['B2']) / (r['B3'] + r['B2'])
  indices['CR.red.nir']  <- r['B7']  / (r['B8'] + ( (r['B6'] - r['B8']) / (740.5- 833) ) * (782.8 - 832.8))
  

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
