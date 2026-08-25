#' Compute Top-of-Canopy and Top-of-Atmosphere Reflectance and Radiance with Simulations
#'
#' This function computes the top-of-canopy (TOC) reflectance, top-of-atmosphere (TOA)
#' reflectance, and TOA radiance using the Soil-Plant-Atmosphere Radiative Transfer (SPART)
#' model based on the provided input parameters and models. It also supports multiple
#' simulations to assess variability.
#'
#' @param inputLUT A lookup table (LUT) containing the necessary parameters for the radiative transfer model.
#' @param optipar Optional parameters for optical calculations (default is NULL).
#' @param CanopyModel A character string specifying the canopy model to use. Default is 'fourSAIL'.
#' @param LeafModel A character string specifying the leaf model to use. Default is 'PROSPECT-PRO'.
#' @param df.irradiance A data frame containing irradiance data for the calculations.
#' @param sensor.i A dataset with sensor information for which the calculations are performed. Options include "LANDSAT4.TM", "LANDSAT5.TM", 
#'               "LANDSAT5.TM, "LANDSAT7.ETM", "LANDSAT8.OLI","Sentinel2A.MSI", "Sentinel2B.MSI", "Sentinel3A.OLCI",
#'               "Sentinel3B.OLCI"and "TerraAqua.MODIS"
#' @param rsoil Optional numeric vector, length 2001 (400-2400nm, 1nm step),
#'   of soil reflectance to use instead of the built-in BSM soil model --
#'   e.g. from \code{\link{get.marmit.rsoil}}. See \code{\link{SPART}} for
#'   details; behaves identically here. Default \code{NULL} (use BSM).
#'
#' @return A list containing TOC reflectance, TOA reflectance, and TOA radiance for each simulation.
#' @export
#'
#' @examples
#' \dontrun{
#' # Requires a real LUT (see get.LUTfromRanges() / getLUT()) and real
#' # irradiance data — not run automatically as no example dataset ships
#' # with the package for this.
#' LUT <- get.LUTfromRanges(...)          # build/load a real LUT here
#' df.irradiance <- data.frame(...)        # real irradiance data here
#' results <- SPART.simN(inputLUT = LUT, optipar = NULL,
#'                        CanopyModel = 'fourSAIL', LeafModel = 'PROSPECT-PRO',
#'                        sensor.i = ToolsRTM::TerraAqua.MODIS, df.irradiance = df.irradiance,
#'                        N = 10)
#' }
SPART.simN <-  function(inputLUT, optipar=NULL, CanopyModel = 'fourSAIL',LeafModel='PROSPECT-PRO',
                   sensor.i = NULL , df.irradiance=NULL, rsoil = NULL){

  ## Soil-Plant-Atmosphere Radiative Transfer model for top-of-canopy and top-of-atmosphere reflectance
  # Developed by Peiqi Yang               (p.yang@utwente.nl)
  #              Christiaan van der Tol   (c.vandertol@utwente.nl)
  #              Wout Verhoef             (w.verhoef@utwente.nl)
  # ported to R C.Camino

  ##### Initial conditions
  # inputs = inputs.SPART
  # LUT<-as.data.frame(getLUT(inputs = inputs, nLUT=1, setseed = 1234))
  # inputLUT=LUT[1,]
  # optipar=optipar2021.Pro.CX;
  # CanopyModel = 'fourSAIL';LeafModel='PROSPECT-PRO';
  # sensor.i =LANDSAT7.ETM; get.plots=T;



  ##################################################################################
  ### 0.1 Check the model
  ##################################################################################

  if (missing(LeafModel)){
    LeafModel='PROSPECT-PRO'
  }

  if (missing(CanopyModel)){
    CanopyModel='fourSAIL'
  }

  if (missing(inputLUT)){
    inputLUT = inputs.SPART
    inputLUT<-as.data.frame(getLUT(inputs = inputLUT, nLUT=1, setseed = 1234))
  }
  ##################################################################################
  ### 0.2 Check the canopy model
  ##################################################################################

  if (missing(CanopyModel)){
    model = 'fourSAIL'

  } else {
    model = CanopyModel
  }

  ##################################################################################
  ### 0.0 Load optipar (optical leaf properties)
  ##################################################################################

  if (is.null(optipar)){
    # optipar2021.Pro.CX only exists in SCOPEinR, not ToolsRTM (see the same
    # fix already applied in SPART()) -- ToolsRTM::optipar is this package's
    # own equivalent dataset.
    optipar <- ToolsRTM::optipar
  }

  if (is.null(optipar$GSV) && all(c("GSV.1", "GSV.2", "GSV.3") %in% names(optipar))) {
    optipar$GSV <- cbind(optipar$GSV.1, optipar$GSV.2, optipar$GSV.3)
  }

  df.optipar <- data.frame(Kp = optipar$Kp,
                           Kcab = optipar$Kab,
                           Kca = optipar$Kca,
                           Kant = optipar$Kant,
                           Kcbc = optipar$Kcbc,
                           wave = optipar$wl)

  ref.optipar <- dataSpec_PRO
  df.ref <- data.frame(Kp = ref.optipar$spA_Protein,
                       Kcab = ref.optipar$spA_Cab,
                       Kca = ref.optipar$spA_Car,
                       Kant = ref.optipar$spA_Atn,
                       Kcbc = ref.optipar$spA_nonPro,
                       wave = ref.optipar$wave)

  ##################################################################################
  ### 0.1 Load spectral dataset (optical leaf properties)
  ##################################################################################

  data.spectral<-get.spectra.spart(getSpectral = T)

  IwlP <- data.spectral[['IwlP']]
  IwlT <- data.spectral[['IwlT']]
  nwlP <- data.spectral[['nwlP']]
  nwlT <- data.spectral[['nwlT']]

  ##################################################################################
  ### 1. Run the BSM model for soil reflectance
  ##################################################################################
  # Parameters for BSM model
  soilemp<-list()
  soilemp[['SMC']] = 25         # empirical parameter (fixed) for BSM (25, original value in SCOPE )
  soilemp[['film']] = 0.015    #empirical parameter (fixed) for BSM
  soilemp[['SMp']] = 15    #soil moisture volume percentage (5 - 55)

  # Was: data.soil <- getinputLUT(inputLUT = SCOPE.LUT.default[1,], dataset =
  # 'soil', ...), driven by data.opts/SCOPE.LUT.default -- but getinputLUT()
  # and both of those datasets only exist in SCOPEinR, not ToolsRTM (same
  # issue already fixed in SPART()). Build just what's actually needed.
  data.soil <- list(BSMBrightness = 0.5, BSMlat = 25, BSMlon = 45, rs_thermal = 0.06)

  wave.comp = c(data.spectral[['reg1']],data.spectral[['reg2']],data.spectral[['reg3']])

  if (!is.null(rsoil)) {
    if (length(rsoil) != length(data.spectral[['reg1']])) {
      stop("SPART.simN(): 'rsoil' must have length ", length(data.spectral[['reg1']]),
           " (", min(data.spectral[['reg1']]), "-", max(data.spectral[['reg1']]),
           "nm, 1nm step), got length ", length(rsoil), ".")
    }
    rfl.comp = c(rsoil, rep(data.soil[['rs_thermal']],length(data.spectral[['IwlT']])))
  } else {
    spec<-list()
    spec[['GSV']] <- optipar$GSV #optipar2020.prospectD.BSM2019$GSV
    spec[['Kw']] <-optipar$Kw # water absorption spectrum
    spec[['nw']] <- optipar$nw # water refraction index spectrum

    rsoil.BSM <- getBSM.toolsRTM(soilpar=data.soil,spec=spec,emp = soilemp);

    rfl.comp = c(c(rsoil.BSM), rep(data.soil[['rs_thermal']],length(data.spectral[['IwlT']])))
  }

  db.rsoil <- data.frame(wave=wave.comp, rfl.soil = rfl.comp)




  ##################################################################################
  ### 2. Run the model
  ##################################################################################


  model.sim<-foursail(inputLUT=inputLUT[1,],rsoil=db.rsoil[['rfl.soil']],
                                LeafModel = LeafModel, spectrum.all = F)

  rfl.canopy.brdf<-Compute_BRF(rdot=model.sim$rdot,rsot=model.sim$rsot,tts=inputLUT[1,'tts'],data.light=dataSpec_PDB,  short.waves=T)

  df.toc <- data.frame(wave= data.spectral[['wlO']], rfl.toc = rfl.canopy.brdf)

  ### get SMAC coeff. for specifi sensor
  sensors.properties = get.coef.SMAC(sensor = sensor.i)
  wlSensor =  sensors.properties[['wl.smac']]
  coefs.SMAC  = sensors.properties[['coefs.SMAC']]
  Sensor.name = sensors.properties[['Sensor.name']]
  center.wvl <- sensors.properties[['center.wvl']]
  wl.smac <- sensors.properties[['wl.smac']]

  df.sensor = data.frame(wlSensor= wlSensor,wl.smac=wl.smac  ) #center.wvl=center.wvl)
  ### Interpolation using sensor
  rv_so <- signal::interp1(data.spectral[['wlO']], model.sim[['rsot']] ,wlSensor,'spline',1E-4)
  rv_do  <- signal::interp1(data.spectral[['wlO']], model.sim[['rdot']] ,wlSensor,'spline',1E-4)
  rv_dd <- signal::interp1(data.spectral[['wlO']], model.sim[['rddt']] ,wlSensor,'spline',1E-4)
  rv_sd <- signal::interp1(data.spectral[['wlO']], model.sim[['rsdt']] ,wlSensor,'spline',1E-4)

  rfl.canopy.brdf.sensor <- signal::interp1(data.spectral[['wlO']], rfl.canopy.brdf ,wlSensor,'spline',1E-4)

   # Run the atmosphere radiative transfer model
  atmopt <- get.smac(inputLUT=inputLUT, sensor= sensor.i )

  # Extract results from the model output
  ta_ss <- atmopt[['Ta_ss']]
  ta_sd <- atmopt[['Ta_sd']]
  ta_oo <- atmopt[['Ta_oo']]
  ta_do <- atmopt[['Ta_do']]
  ra_dd <- atmopt[['Ra_dd']]
  ra_so <- atmopt[['Ra_so']]
  T_g   <- atmopt[['Tg']]

  # Upscale TOC to TOA
  R_TOC <- (ta_ss * rv_so + ta_sd * rv_do) / (ta_ss + ta_sd)

  # Calculate components for TOA reflectance
  rtoa1 <- (ta_sd * rv_do + ta_ss * rv_sd * ra_dd * rv_do) * ta_oo / (1 - rv_dd * ra_dd)
  rtoa2 <- (ta_ss * rv_sd + ta_sd * rv_dd) * ta_do / (1 - rv_dd * ra_dd)
  rtoa0 <- ra_so + ta_ss * rv_so * ta_oo

  # Calculate TOA reflectance
  R_TOA <- T_g * (rtoa0 + rtoa1 + rtoa2)



  # Default dataset when df.irradiance is NULL
  if (missing(df.irradiance) || is.null(df.irradiance)) {

    wl_Ea = Extraterrestrial_irradiance$wave # wl
    Ea0 = Extraterrestrial_irradiance$EIrrad # Ea0 = W M-2 nm-1
    df.irradiance <- data.frame(wave = wl_Ea, irrad = Ea0)
  }

  if (!is.null(df.irradiance) && is.data.frame(df.irradiance)) {

    df.irradiance <- data.frame(wave = df.irradiance[,1], irrad = df.irradiance[,2])
  } else {
   stop('please provide a dataframe with wave and Eo ...')
  }


  La_extra <-  get.spectral.convolution(df.irradiance = df.irradiance,sensor.i, get.plots=F)
  # Calculate TOA radiance
  L_TOA <- La_extra$Eo * R_TOA

  df.rad.toa<- data.frame(wave = La_extra$wave, rad.toa = L_TOA, rfl.toa = R_TOA)



  df.refl <- data.frame(wave = wlSensor,rfl.toa = R_TOA, rad.toa = L_TOA,rfl.toc = R_TOC,rfl.toc.brdf= rfl.canopy.brdf.sensor )

  ##::::::::::::::::::::::::::::::::::::::::::
  # Plot for Reflectance at TOA, TOC
  ##:::::::::::::::::::::::::::::::::::::::::



  df.all <- data.frame(wave = wlSensor, rad.toa = L_TOA, rfl.toa = R_TOA, rfl.toc = R_TOC, rfl.toc.BRDF= rfl.canopy.brdf.sensor)
  rownames(df.all) <-NULL
  outputs <- list(output = df.all, Irrad.toa = La_extra, rfl.toc.brdf = df.toc,
                  rf.soil = db.rsoil, sensor = df.sensor,
                  inputLUT=inputLUT)

  return(outputs)


}
