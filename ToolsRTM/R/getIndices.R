#' GetIndices: A Function to Estimate Common Spectral Indices
#'
#' This function calculates common spectral vegetation indices from reflectance data,
#' using user-specified spectral domains and a scaling factor if needed.
#'
#' @param data A dataframe or matrix containing reflectance values in named columns.
#' @param pattern.rfl A character string indicating the prefix pattern for reflectance column names.
#'   Default is `'R.'`, which assumes column names like `'R.400'`, `'R.670'`, etc.
#' @param factor A numeric scaling factor applied to reflectance values (e.g., 1/10000 for scaling
#'   reflectance from digital numbers). Default is `NULL`, meaning no scaling is applied (factor = 1).
#' @param spectral.domain A character string specifying the spectral domain for index calculation.
#'   Available options are:
#'   - `"VNIR"` (default): Includes wavelengths between 400-850 nm.
#'   - `"SWIR"`: Includes wavelengths between 800-1750 nm.
#'   - `"VNIR-SWIR"`: Includes wavelengths between 400-1750 nm.
#'   If `NULL`, it defaults to `"VNIR"`.
#'
#' @return A dataframe containing the original data along with computed spectral indices.
#' @export
#'
#' @examples
#' # Example 1: Compute spectral indices using default VNIR domain
#' data_example <- data.frame(R.400 = runif(10, 0, 1), R.670 = runif(10, 0, 1), R.800 = runif(10, 0, 1))
#' result <- getIndices(data_example)
#'
#' # Example 2: Compute indices using SWIR domain with a scaling factor of 1/10000
#' data_example <- data.frame(R.850 = runif(10, 0, 1), R.1200 = runif(10, 0, 1), R.1600 = runif(10, 0, 1))
#' result <- getIndices(data_example, spectral.domain = "SWIR", factor = 1/10000)
#'
#' # Example 3: Compute indices for VNIR-SWIR domain with custom reflectance column pattern
#' data_example <- data.frame(Reflectance_450 = runif(10, 0, 1), Reflectance_900 = runif(10, 0, 1))
#' result <- getIndices(data_example, pattern.rfl = "Reflectance_", spectral.domain = "VNIR-SWIR")
#'


getIndices <- function(data, pattern.rfl='R.', factor = NULL, spectral.domain=NULL) {

  # If factor is NULL, set it to 1 (no scaling)
  if (is.null(factor)) {
    factor <- 1
  }

  # Ensure input data is a dataframe or matrix
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Error: 'data' must be a dataframe or matrix")
  }

  # Convert matrix to dataframe if needed
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }


  # Validate spectral domain
  valid_domains <- c("VNIR", "SWIR", "VNIR-SWIR")

  if (is.null(spectral.domain)) {
    spectral.domain <- "VNIR"  # Default to VNIR if NULL
  }

  if (!spectral.domain %in% valid_domains) {
    stop("Error: Please use a correct spectral domain. Options are 'VNIR', 'VNIR-SWIR', and 'SWIR'.")
  }

  print(paste("Estimating indices using", spectral.domain, "domain..."))


  # Extract reflectance column names matching the pattern
  reflectance_cols <- grep(paste0("^", pattern.rfl, "[0-9]+(\\.[0-9]+)?$"), names(data), value = TRUE)


  # Extract numeric wavelengths
  wavelengths <- as.numeric(gsub(paste0(pattern.rfl, "(\\d+\\.?\\d*)"), "\\1", reflectance_cols))

  # Define wavelength ranges for spectral domain selection
  if (spectral.domain == "VNIR") {
    selected_cols <- reflectance_cols[wavelengths >= 400 & wavelengths <= 850]
  } else if (spectral.domain == "SWIR") {
    selected_cols <- reflectance_cols[wavelengths >= 800 & wavelengths <= 2550]
  } else if (spectral.domain == "VNIR-SWIR") {
    selected_cols <- reflectance_cols[wavelengths >= 400 & wavelengths <= 2550]
  }

  if (length(selected_cols) == 0) {
    stop("Error: No reflectance columns fall within the specified spectral domain.")
  }

  # Define min and max range (rounded)
  min_ <- round(min(wavelengths)) + 10  # Ensure integer rounding works
  max_ <- round(max(wavelengths)) + 10


  range2interpo <-  c(min_:max_)

  indices.list = list()

  # create progress bar
  total=dim(data)[1]
  barProgress <- txtProgressBar(min = 0, max = total, style = 3)

  for (i in c(1:dim(data)[1])){
    # Must match `wavelengths` (derived from reflectance_cols, not the
    # domain-filtered selected_cols -- selected_cols is only used above to
    # validate the spectral domain actually has columns in range).
    values <- as.numeric(data[i, reflectance_cols]) * factor
    r = signal::interp1(wavelengths, values, range2interpo, extrap = T)
    names(r) <- range2interpo

    indices = c()

    if(spectral.domain %in% 'VNIR'){

      ##########################################
      ## Structurual Indices
      ##########################################

      # NDVI (R800-R670)/(R800+R670)
      # Rouse et al. (1974)
      indices['NDVI'] <- (r['800'] - r['670']) / (r['800'] + r['670'])

      # RDVI (R800-R670)/(R800+R670)^0.5
      # Rougean and Breon (1995)
      indices['RDVI'] <- (r['800'] - r['670']) / (r['800'] + r['670']) ** 0.5

      # SR R800/R670
      # Jordan (1969)
      indices['SR'] <- r['800'] / r['670']

      # MSR (R800/R670-1)/((R800/R670)^0.5+1)
      # Chen (1996)
      indices['MSR'] <- (r['800'] / r['670'] - 1) / ((r['800'] / r['670']) ** 0.5 + 1)

      # OSAVI [(1+0.16)*(R800-R670)/(R800+R670+0.16)]
      # Rondeaux et al. (1996)
      indices['OSAVI'] <- ((1 + 0.16) * (r['800'] - r['670']) / (r['800'] + r['670'] + 0.16))

      # MSAVI 1/2*[(2*R800+1-?((?(2*R800+1)?^2)-8*(R800-R670))]
      # Qi et al. (1994)
      indices['MSAVI'] <- 1 / 2 * (2 * r['800'] + 1 - sqrt(((2 * r['800'] + 1) ^ 2) - 8 * (r['800'] - r['670'])))

      # MTVI1 1.2*[1.2*(R800-R550)-2.5*(R670-R550)]
      # Broge & Leblanc (2000); Haboudane et al. (2004)
      indices['MTVI1'] <- 1.2 * (1.2 * (r['800'] - r['550']) - 2.5 * (r['670'] - r['550']))

      # MTVI2 (1.5*[1.2*(R800-R550)-2.5*(R670-R550)])/SQR((2*R800+1)^2-(6*R800-5*SQR(R670))-0.5)
      # Haboudane et al. (2004)
      indices['MTVI2'] <- (1.5 * (1.2 * (r['800'] - r['550']) - 2.5 * (r['670'] - r['550']))) / sqrt((2 * r['800'] + 1) ^ 2 - (6 * r['800'] - 5 * sqrt(r['670'])) - 0.5)

      # MCARI ((R700-R670) - 0.2*(R700-R550))*(R700/R670)
      # Hermann et al. (2010)
      indices['MCARI'] <- ((r['700'] - r['670']) - 0.2 * (r['700'] - r['550'])) * (r['700'] / r['670'])

      # MCARI1 1.2* [2.5* (R800-  R670)-  1.3* (R800 -R550) ]
      # Haboudane et al. (2004)
      indices['MCARI1'] <- 1.5 * (2.5 * (r['800'] - r['670']) - 1.3 * (r['800'] - r['550']))

      # MCARI2 (1.5*[2.5*(R800-R670)-1.3*(R800-R550) ])/SQRT((2*R800+1)^2-(6*R800-5*SQRT(R670))-0.5)
      # Haboudane et al. (2004)
      indices['MCARI2'] <- (1.5 * (2.5 * (r['800'] - r['670']) - 1.3 * (r['800'] - r['550']))) / sqrt((2 * r['800'] + 1) ^ 2 - (6 * r['800'] - 5 * sqrt(r['670'])) - 0.5)

      # EVI 2.5*(R800-R670)/(R800+6*R670-7.5*R400+1)
      # Huete et al. (2002)
      indices['EVI'] <- 2.5 * (r['800'] - r['670']) / (r['800'] + 6 * r['670'] - 7.5 * r['400'] + 1)

      # LIC1 (R800-R680)/(R800+R680)
      # Lichtenthaler et al. (1996)
      indices['LIC1'] <- (r['800'] - r['680']) / (r['800'] + r['680'])

      ##########################################
      ## Pigments Indices
      ##########################################

      # VOG R740/R720
      # Vogelmann et al. (1993)
      indices['VOG'] <- r['740'] / r['720']

      # VOG2 (R734-R747)/(R715+R726)
      # Vogelmann et al. (1993); Zarco-Tejada et al. (1999)
      indices['VOG2'] <- (r['734'] - r['747']) / (r['715'] + r['726'])

      # VOG3 (R734-R747)/(R715+R720)
      # Vogelmann et al. (1993); Zarco-Tejada et al. (1999)
      indices['VOG3'] <- (r['734'] - r['747']) / (r['715'] + r['720'])

      # GM1 R750/R550
      # Gitelson and Merzlyak (1997)
      indices['GM1'] <- r['750'] / r['550']

      # GM2 R750/R700
      # Gitelson and Merzlyak (1997)
      indices['GM2'] <- r['750'] / r['700']

      # TCARI 3*[(R700-R670)-0.2*(R700-R550)*(R700/R670)]
      # Haboudane et al. (2002)
      indices['TCARI'] <- 3 * ((r['700'] - r['670']) - 0.2 * (r['700'] - r['550']) * (r['700'] / r['670']))

      # TCARI/OSAVI TCARI/OSAVI
      # Haboudane et al. (2002)
      indices['T.O'] <- as.numeric(indices['TCARI']) / as.numeric(indices['OSAVI'])

      # CI R750/R710
      # Zarco-Tejada et al. (2001)
      indices['CI'] <- r['750'] / r['710']

      # TVI 0.5*[120*(R750-R550)-200*(R670-R550) ]
      # Broge and Leblanc (2000)
      indices['TVI'] <- 0.5 * (120 * (r['750'] - r['550']) - 200 * (r['670'] - r['550']))

      # SRPI R430/R680
      # Pe??uelas et al. (1995)
      indices['SRPI'] <- r['430'] / r['680']

      # NPQI (R415-R435)/(R415+R435)
      # Barnes (1992)
      indices['NPQI'] <- (r['415'] - r['435']) / (r['415'] + r['435'])

      # NPCI (R680-R430)/(R680+R430)
      # Pe??uelas et al. (1994)
      indices['NPCI'] <- (r['680'] - r['430']) / (r['680'] + r['430'])

      # CTR1 R695/R420
      # Carter (1994); Carter et al. (1996)
      indices['CTR1'] <- r['695'] / r['420']

      # CAR R515/R570
      # Hernandez-Clemente et al. (2012)
      indices['CAR'] <- r['515'] / r['570']

      # Datt-CabCx+c R672/((R550*(3*R708)))
      # Datt (1998)
      indices['DCabxc'] <- r['672'] / ((r['550'] * (3 * r['708'])))

      # DattNIRCabCx+c R860/((R550*R708))
      # Datt (1998)
      indices['DNCabxc'] <- r['860'] / ((r['550'] * r['708']))

      # SIPI (R800-R445)/(R800+R680)
      # Pe??uelas et al. (1995)
      indices['SIPI'] <- (r['800'] - r['445']) / (r['800'] + r['680'])

      # CRI550 (1/R510)-(1/R550)
      # Gitelson et al. (2003, 2006)
      indices['CRI550'] <-(1 / r['510']) - (1 / r['550'])

      # CRI700 (1/R510)-(1/R700)
      # Gitelson et al. (2003, 2006)
      indices['CRI700'] <- (1 / r['510']) - (1 / r['700'])

      # CRI550 (515 instead of 510) (1/R510)-(1/R550)
      # Gitelson et al. (2003, 2006)
      indices['CRI550m'] <- (1 / r['515']) - (1 / r['550'])

      # CRI700 (515 instead of 510) (1/R510)-(1/R700)
      # Gitelson et al. (2003, 2006)
      indices['CRI700m'] <- (1 / r['515']) - (1 / r['700'])

      # RNIR*CRI550 (1/R510)-(1/R550)*R770
      # Gitelson et al. (2003, 2006)
      indices['RCRI550'] <- (1 / r['510']) - (1 / r['550']) * r['770']

      # RNIR*CRI700 (1/R510)-(1/R700)*R770
      # Gitelson et al. (2003, 2006)
      indices['RCRI700'] <- (1 / r['510']) - (1 / r['700']) * r['770']

      # PSRI (R680-R500)/R750
      # Merzlyak et al. (1996)
      indices['PSRI'] <- (r['680'] - r['500']) / r['750']

      # LIC3 R440/R740
      # Lichtenhaler et al. (1996)
      indices['LIC3'] <- r['440'] / r['740']

      ## Red-edge Indices
      ##original conf form J.B Feret for CR_SWIR
      waves_ <- c('w02'=497, 'w03'=560.0, 'w04'=665, 'w05'=704, 'w06'=740,
                  'w07' = 782, 'w08' = 835, 'w8A' = 865, 'w11' = 1614, 'w12' = 2202,'w.800'=800,'w.762'=762)
      #
      indices['CIre'] <- (r['782'] / r['705'])-1
      indices['CIrededge'] <- (r['800'] / r['705'])-1
      indices['CIgreen'] <- (r['800'] / r['B3'])-1
      indices['Chlred.edge'] <- (r['780'] / r['705']) ** (-1)
      indices['CVI'] <- (r['800'] * r['665']) / r['665']**2
      indices['IRECI'] <- (r['780'] - r['665']) / (r['705'] / r['740'])
      indices['REP'] <- 700 + 40* (((r['665'] + r['780'])/2) -  r['705'])/ (r['740'] - r['705'])
      indices['RVI'] <-  (r['800'] / r['665'])

      indices['RedEg1'] <- r['705'] / r['665']
      indices['RedEg2'] <- (r['705'] - r['665']) / (r['705'] + r['665'])


      # ##########################################
      ## PRIs Indices
      ##########################################

      # PRI (R570-R530)/(R570+R530)
      # Gamon et al. (1992)
      indices['PRI'] <- (r['570'] - r['530']) / (r['570'] + r['530'])

      # PRI515 (R515-R530)/(R515+R530)
      # Hern??ndez-Clemente et al. (2011)
      indices['PRI515'] <- (r['515'] - r['530']) / (r['515'] + r['530'])

      # PRIM1 (R512-R531)/(R512+R531)
      # Gamon et al. (1993)
      indices['PRIM1'] <- (r['512'] - r['531']) / (r['512'] + r['531'])

      # PRIM2 ((R600-R531))/((R600+R531))
      # Gamon et al. (1993)
      indices['PRIM2'] <- (r['600'] - r['531']) / (r['600'] + r['531'])

      # PRIM3 (R670-R531)/(R670+R531)
      # Gamon et al. (1993)
      indices['PRIM3'] <- (r['670'] - r['531']) / (r['670'] + r['531'])

      # PRIM4 (R570-R531-R670)/(R571+ R531+ R670)
      # Gamon et al. (1993)
      indices['PRIM4'] <- (r['570'] - r['531'] - r['670']) / (r['570'] + r['531'] + r['670'])

      # PRIn PRI/(RDVI*R700/R670)
      # Zarco-Tejada et al. (2013)
      indices['PRIn'] <- as.numeric(indices['PRI']) / (as.numeric(indices['RDVI']) * r['700'] / r['670'])

      # PRI*CI ((R570-R530)/(R570+R530))*((R760/R700)-1)
      # Garrity et al. (2011)
      indices['PRI_CI'] <- as.numeric(indices['PRI']) * ((r['760'] / r['700']) - 1)

      ##########################################
      ## BGR Indices
      ##########################################

      # B R450/R490
      # -
      indices['B'] <- r['450'] / r['490']

      # G R550/R670
      # -
      indices['G'] <- r['550'] / r['670']

      # R R700/R670
      # Gitelson et al. (2000)
      indices['R'] <- r['700'] / r['670']

      # BGI1 R400/R550
      # Zarco-Tejada et al. (2005; 2012)
      indices['BGI1'] <- r['400'] / r['550']

      # BGI2 R450/R550
      # Zarco-Tejada et al. (2005; 2012)
      indices['BGI2'] <- r['450'] / r['550']

      # BF1: R400/R410
      # BF2: R400/R420
      # BF3: R400/R430
      # BF4: R400/R440
      # BF5: R400/R450
      indices['BF1'] <- r['400'] / r['410']
      indices['BF2'] <- r['400'] / r['420']
      indices['BF3'] <- r['400'] / r['430']
      indices['BF4'] <- r['400'] / r['440']
      indices['BF5'] <- r['400'] / r['450']

      # BRI1 R400/R690
      # Zarco-Tejada et al. (2012)
      indices['BRI1'] <- r['400'] / r['690']

      # BRI2 R450/R690
      # Zarco-Tejada et al. (2012)
      indices['BRI2'] <- r['450'] / r['690']

      # RGI R690/R550
      # -
      indices['RGI'] <- r['690'] / r['550']

      # RARS R746/R513
      # Chappelle et al. (1992)
      indices['RARS'] <- r['746'] / r['513']

      # LIC2 R440/R690
      # Lichtenthaler et al. (1996)
      indices['LIC2'] <-  r['440'] / r['690']

      # HI (R534-R698)/(R534+R698)-1/2*R704
      # Mahlein et al. (2013)
      indices['HI'] <- (r['534'] - r['698']) / (r['534'] + r['698']) - 1 / 2 * r['704']

      # CUR (R675*R690)/(R683)^2
      # Zarco-Tejada et al. (2000)
      indices['CUR'] <- (r['675'] * r['690']) / (r['683']) ** 2

      ##########################################
      ## NIR-VIS Indices
      ##########################################

      # PSSRa R800/R680
      # Blackburn (1998)
      indices['PSSRa'] <- r['800'] / r['680']

      # PSSRb R800/R635
      # Blackburn (1998)
      indices['PSSRb'] <- r['800'] / r['635']

      # PSSRc R800/R470
      # Blackburn (1998)
      indices['PSSRc'] <- r['800'] / r['470']

      # PSNDc (R800-R470)/(R800+R470)
      # Blackburn (1998)
      indices['PSNDc'] <- (r['800'] - r['470']) / (r['800'] + r['470'])



      ########################################################
      ## Adding spectral indices using  Red-edge channels
      ########################################################

      ##original conf form J.B Feret for CR_SWIR
      waves_ <- c('w02'=497, 'w03'=560.0, 'w04'=665, 'w05'=704, 'w06'=740,
                  'w07' = 782, 'w08' = 835, 'w8A' = 865, 'w11' = 1614, 'w12' = 2202,'w.800'=800,'w.762'=762)


      indices['CR.red.nir.1'] <- r['740']  / (r['800'] + (waves_['w06'] - waves_['w.800']) * (r['782'] - r['800']) / (waves_['w.762'] - waves_['w.800']))
      indices['CR.red.nir.2'] <- r['704']  / (r['782'] + (waves_['w05'] - waves_['w07']) * (r['740'] - r['782']) / (waves_['w06'] - waves_['w.762']))
      indices['CR.red.nir.3'] <- r['704']  / (r['800'] + (waves_['w05'] - waves_['w.800']) * (r['740'] - r['800']) / (waves_['w06'] - waves_['w.800']))
      indices['CR.red.nir.4'] <- r['665']  / (r['800'] + (waves_['w04'] - waves_['w.800']) * (r['704'] - r['800']) / (waves_['w06'] - waves_['w.800']))
      indices['CR.red.nir.5'] <- r['740']  / (r['800'] + (waves_['w06'] - waves_['w.800']) * (r['704'] - r['800']) / (waves_['w05'] - waves_['w.800']))
      indices['CR.red.nir'] <- r['782']  / (r['800'] + (waves_['w07'] - waves_['w.800']) * (r['704'] - r['800']) / (waves_['w05'] - waves_['w.800']))
      indices['CR.red.nir.7'] <- r['762']  / (r['800'] + (waves_['w.762'] - waves_['w.800']) * (r['704'] - r['800']) / (waves_['w05'] - waves_['w.800']))


      indices.list[[i]] = indices

    }

    if (spectral.domain %in% "SWIR" | spectral.domain %in% "VNIR-SWIR") {
      ##########################################
      ## SWIR Indices
      ##########################################

      # GnyLi ((R900_VNIR*R1050)-(R955*R1220))/((R900_VNIR*R1050)+(R955*R1220))
      # Gnyp et al. (2014)
      indices['GnyLi'] <- ((r['900'] * r['1050']) - (r['955'] * r['1220'])) / ((r['900'] * r['1050']) + (r['955'] * r['1220']))

      # GnyLi ((R850_VNIR*R1050)-(R955*R1220))/((R850_VNIR*R1050)+(R955*R1220))
      # Prev. mod. 850 nm instead of 900 nm
      indices['GnyLi.w850'] <- indices['GnyLi'] <- ((r['850'] * r['1050']) - (r['955'] * r['1220'])) / ((r['850'] * r['1050']) + (r['955'] * r['1220']))

      # GnyLi ((R850_VNIR*R1050)-(R955*R1220))/((R850_VNIR*R1050)+(R955*R1220))
      # Prev. mod. 950 nm instead of 900 nm
      indices['GnyLi.w950'] <- indices['GnyLi'] <- ((r['950'] * r['1050']) - (r['955'] * r['1220'])) / ((r['950'] * r['1050']) + (r['955'] * r['1220']))

      # CI1 ((R736-R735)/1)*(R990/R720)
      # Yansong et al. (2013)
      indices['CI1'] <- ((r['736'] - r['735']) / 1) * (r['990'] / r['720'])

      # CI2 (w/ VNIR) ((R736-R735)/1)*(R900/R720)
      # Yansong et al. (2013)
      indices['CI2'] <- ((r['736'] - r['735']) / 1) * (r['900'] / r['720'])

      # CI2 (w/ VNIR) ((R736-R735)/1)*(R900/R720)
      # Prev. mod. 850 nm instead of 900 nm
      #indices['CI2 (w/ VNIR) 850'] <- ((r['736'] - r['735']) / 1) * (r['850'] / r['720'])

      # CI2 (w/ VNIR) ((R736-R735)/1)*(R900/R720)
      # Prev. mod. 950 nm instead of 900 nm
      indices['CI2'] <- ((r['736'] - r['735']) / 1) * (r['950'] / r['720'])

      # MCARI_1510 ((R700-R1510) - 0.2*(R700-R550))*(R700/R1510)
      # Hermann et al. (2010)
      indices['MCARI.1510'] <- ((r['700'] - r['1510']) - 0.2 * (r['700'] - r['550'])) * (r['700'] / r['1510'])

      # TCARI_1510 3*((R700-R1510)-0.2*(R700-R550))*(R700/R1510)
      # Hermann et al. (2010)
      indices['TCARI.1510'] <- 3 * ((r['700'] - r['1510']) - 0.2 * (r['700'] - r['550']) * (r['700'] / r['1510']))

      # OSAVI_1510 (1+0.16)*(R800-R1510)/(R800+R1510+0.16)
      # Rondeaux et al. (1996)
      indices['OSAVI.1510'] <- ((1 + 0.16) * (r['800'] - r['1510']) / (r['800'] + r['1510'] + 0.16))

      # TCARI_OSAVI_1510 TCARI/OSAVI (1510 nm)
      # Hermann et al. (2010)
      indices['TCARI/OSAVI.1510'] <- as.numeric(indices['TCARI 1510']) / as.numeric(indices['OSAVI 1510'])

      # NRI_1510 (R1510-R660)/(R1510+R660)
      # Hermann et al. (2010)
      indices['NRI.1510'] <- (r['1510'] - r['660']) / (r['1510'] + r['660'])

      # RSI_990_720 R990/R720
      # Yao et al. (2010)
      indices['RSI.990.720'] <- r['990'] / r['720']

      # NRI_1770_693 (R1770-R693)/(R1770+R693)
      # Ferwerda et al. (2005)
      #indices['NRI 1770 693'] <- (r['1770'] - r['693']) / (r['1770'] + r['693'])

      # NDNI (log10(1/R1510)-log10(1/R1680))/(log10(1/R1510)+log10(1/R1680))
      # Serrano et al. (2002)
      indices['NDNI'] <- (log10(1 / r['1510']) - log10(1 / r['1680'])) / (log10(1 / r['1510']) + log10(1 / r['1680']))

      # S1080 (R1080-R660)/(R1080+R660)
      # Mahajan et al. (2014)
      indices['S1080'] <- (r['1080'] - r['660']) / (r['1080'] + r['660'])

      # S1260 (R1260-R660)/(R1260+R660)
      # Mahajan et al. (2014)
      indices['S1260'] <- (r['1260'] - r['660']) / (r['1260'] + r['660'])

      # N1645 (R1645-R1715)/(R1645+R1715)
      # Pimstein et al. (2011)
      indices['N1645'] <- (r['1645'] - r['1715']) / (r['1645'] + r['1715'])

      # N870 (R870-R1450)/(R870+R1450)
      # Pimstein et al. (2011)
      indices['N870'] <- (r['870'] - r['1450']) / (r['870'] + r['1450'])
      ## Camino et al. 2018
      indices['N850.1510'] <- (r['850'] - r['1510']) / (r['850'] + r['1510'])
      ## Camino et al. 2018
      indices['NN1510'] <- (r['1510']) / (r['850'])


      #CR_SWIR from J.B.Feret
      waves_swir <- c('w02'=497, 'w03'=560.0, 'w04'=665, 'w05'=704, 'w06'=740,
                      'w07' = 782, 'w08' = 835, 'w8A' = 865, 'w11' = 1614, 'w12' = 2202)
      indices['CR.SWIR'] <- r['1614']/(r['865']+(waves_swir['w11']-waves_swir['w8A'])*(r['2202']-r['865'])/(waves_swir['w12']-waves_swir['w8A']))



      indices.list[[i]] = indices

    }
    setTxtProgressBar(barProgress, i)

  }
  close(barProgress)
  # Combine the original dataset with calculated indices
  indices.df <- as.data.frame(do.call(rbind, indices.list))
  ## remove indices with no data (this used to compute df.indices then return
  ## the unfiltered indices.df instead -- fixed to actually apply the filter)
  indices.df <- indices.df[, colSums(is.na(indices.df)) != nrow(indices.df), drop = FALSE]
  ## merge back onto the original data, as documented in @return (this used
  ## to be commented out, so callers never actually got the original data --
  ## many downstream scripts work around this by cbind()-ing it back
  ## themselves and de-duplicating columns; that's still safe with this fix)
  df.indices <- cbind(data, indices.df)

  return(df.indices)
}
