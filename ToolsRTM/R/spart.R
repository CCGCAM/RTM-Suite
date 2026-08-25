#' Compute Top-of-Canopy and Top-of-Atmosphere Reflectance and Radiance
#'
#' This function computes the top-of-canopy (TOC) reflectance, top-of-atmosphere (TOA)
#' reflectance, and TOA radiance using the Soil-Plant-Atmosphere Radiative Transfer (SPART)
#' model based on the provided input parameters and models.
#'
#' @param inputLUT A lookup table (LUT) containing the necessary parameters for the radiative transfer model.
#' @param optipar Optional parameters for optical calculations (default is NULL).
#' @param CanopyModel A character string specifying the canopy model to use. Default is 'fourSAIL'.
#' @param LeafModel A character string specifying the leaf model to use. Default is 'PROSPECT-PRO'.
#' @param df.irradiance A data frame containing irradiance data for the calculations.
#' @param sensor.i A dataset with sensor information for which the calculations are performed. Options include "LANDSAT4.TM", "LANDSAT5.TM",
#'               "LANDSAT5.TM, "LANDSAT7.ETM", "LANDSAT8.OLI","Sentinel2A.MSI", "Sentinel2B.MSI", "Sentinel3A.OLCI",
#'               "Sentinel3B.OLCI"and "TerraAqua.MODIS"
#' @param rsoil Optional numeric vector, length 2001 (400-2400nm, 1nm step --
#'   matching \code{get.spectra.spart(getSpectral = TRUE)$reg1}), of soil
#'   reflectance. When supplied, this is used directly as the soil spectrum
#'   instead of computing one from the built-in BSM (Brightness-Shape-Moisture)
#'   soil model -- e.g. to feed in a spectrum from \code{\link{get.marmit.rsoil}}
#'   (pass \code{wl.out = 400:2400} there to get a matching-length vector).
#'   The thermal region (2400nm+) is still padded the same way as the BSM
#'   path (a flat \code{0.06} reflectance). Default \code{NULL} (use BSM).
#'   When \code{rsoil} is \code{NULL}, the BSM model reads its 4 real soil
#'   parameters from \code{inputLUT} when present --
#'   \code{BSMBrightness} (overall soil brightness/albedo scaling),
#'   \code{BSMlat} (spectral shape "latitude", 20-40 deg), \code{BSMlon}
#'   (spectral shape "longitude", 45-65 deg), and \code{SMp} (soil moisture
#'   volume percentage, 5-55) -- falling back to \code{0.5}/\code{25}/\code{45}/
#'   \code{15} respectively for any of the four not present in \code{inputLUT}.
#'   Note \code{psoil}, if present in
#'   \code{inputLUT}, is \strong{not} one of these four and has no effect
#'   here -- it's a PROSAIL/4SAIL soil-brightness-mixing parameter (dry vs.
#'   wet reference spectrum interpolation), not a BSM parameter; BSM uses
#'   the 4 parameters above instead. \code{inputsSPART} ships all 4 as
#'   sampleable columns for exactly this reason.
#' @param get.plots A boolean indicating whether to generate plots of the results. Default is TRUE.
#'
#' @return A list containing TOC reflectance, TOA reflectance, and TOA radiance.
#' @export
#'
#' @examples
#' \dontrun{
#' # Requires a real LUT (see get.LUTfromRanges() / getLUT()) and real
#' # irradiance data — not run automatically as no example dataset ships
#' # with the package for this.
#' LUT <- get.LUTfromRanges(...)          # build/load a real LUT here
#' df.irradiance <- data.frame(...)        # real irradiance data here
#' results <- SPART(inputLUT = LUT, optipar = NULL,
#'                   CanopyModel = 'fourSAIL', LeafModel = 'PROSPECT-PRO',
#'                   sensor.i = ToolsRTM::TerraAqua.MODIS, df.irradiance = df.irradiance,
#'                   get.plots = TRUE)
#' }

SPART <-  function(inputLUT, optipar=NULL, CanopyModel = 'fourSAIL',LeafModel='PROSPECT-PRO',
                   sensor.i = NULL , df.irradiance=NULL, rsoil = NULL,
                   get.plots=T){

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
  ### 0.1 Check the plots (TRUE or FALSE)
  ##################################################################################

  # NB (fix): was `if (missing(get.plots)) get.plots <- FALSE else get.plots <- TRUE`
  # -- this ignored whatever value the caller actually passed and just checked
  # whether the argument was supplied at all: passing get.plots = FALSE
  # explicitly got silently overwritten to TRUE (any explicit value forced
  # TRUE), while omitting it entirely forced FALSE -- backwards from the
  # function's own declared default (get.plots = TRUE in the signature).
  # get.plots already holds the correct value via normal R argument binding
  # (the caller's value, or the signature default if omitted); no override
  # needed at all.

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
    inputLUT<-as.data.frame(getLUT(inputs = inputs, nLUT=1, setseed = 1234))
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
    # Was `if (missing(optipar)) optipar = optipar2021.Pro.CX`: (1) missing()
    # is FALSE when a caller explicitly passes optipar = NULL (this
    # function's own documented default, and what the @examples block does),
    # so the fallback never triggered for that -- the intended -- case; (2)
    # optipar2021.Pro.CX only exists in SCOPEinR, not ToolsRTM, so this failed
    # with "object not found" for anyone who hadn't also loaded SCOPEinR.
    # ToolsRTM::optipar is this package's own equivalent dataset.
    optipar <- ToolsRTM::optipar
  }

  # ToolsRTM::optipar stores the three GSV (Global Soil Vectors) components as
  # separate columns (GSV.1/.2/.3) rather than one GSV matrix column like
  # SCOPEinR's optipar variants -- reassemble it here so the BSM soil model
  # call below (which expects spec$GSV as an nwl x 3 matrix) works either way.
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
  if (get.plots ==  T){

    plot.optipar<- ggplot2::ggplot(data = df.optipar, aes(x = wave)) +
      labs(y = "specific absorption coefficient", x = "") + xlim(400, 2500) +
      geom_line(aes(y = Kp, color = "Proteins")) +
      geom_line(aes(y = Kcbc, color = "CBC")) +
      theme_bw() +
      guides(color = guide_legend(title = "Absorption coefficient"), linetype = guide_legend(title = "Absorption coefficient"), shape = guide_legend(title = "Absorption coefficient")) +

      theme(legend.position = "right")

    print(plot.optipar)

    plot.pigments<- ggplot2::ggplot(data = df.optipar, aes(x = wave)) +
      labs(y = "specific absorption coefficient", x = "") + xlim(400, 2500) +
      geom_line(aes(y = Kcab, color = "Cab")) +
      geom_line(aes(y = Kca, color = "Car")) +
      geom_line(aes(y = Kant, color = "Anth")) +
      theme_bw() +
      guides(color = guide_legend(title = "Absorption coefficient"), linetype = guide_legend(title = "Absorption coefficient"), shape = guide_legend(title = "Absorption coefficient")) +

      theme(legend.position = "right")

    print(plot.pigments)

  }

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
  # SMp is a real BSM soil parameter (soil moisture volume %, 5-55) and is
  # read from inputLUT$SMp when present -- SMC/film stay fixed
  # "empirical parameter (fixed)" constants, matching the original SCOPE
  # convention noted below (BSM's own calibration, not something a LUT
  # varies).
  soil_col <- function(name, default) {
    if (name %in% names(inputLUT) && !is.na(inputLUT[[name]][1])) inputLUT[[name]][1] else default
  }
  soilemp<-list()
  soilemp[['SMC']] = 25         # empirical parameter (fixed) for BSM (25, original value in SCOPE )
  soilemp[['film']] = 0.015    #empirical parameter (fixed) for BSM
  soilemp[['SMp']] = soil_col("SMp", 15)    #soil moisture volume percentage (5 - 55)

  # Was: data.soil <- getinputLUT(inputLUT = SCOPE.LUT.default[1,], dataset =
  # 'soil', calc.heat = ..., calc.rss_rbs = ...), driven by data.opts/
  # SCOPE.LUT.default -- but getinputLUT() and both of those datasets only
  # exist in SCOPEinR, not ToolsRTM (and ToolsRTM can't depend on SCOPEinR,
  # since SCOPEinR already depends on ToolsRTM -- that would be circular).
  # SPART only ever uses data.soil for two things below: the BSM soil model's
  # BSMBrightness/BSMlat/BSMlon, and rs_thermal for the thermal-region
  # padding -- the GAM/rss/rbs fields getinputLUT() also computes are for
  # SCOPE's thermal energy balance, not used anywhere in this (purely
  # optical) function. Build just what's actually needed, self-contained.
  data.soil <- list(BSMBrightness = soil_col("BSMBrightness", 0.5),
                     BSMlat = soil_col("BSMlat", 25),
                     BSMlon = soil_col("BSMlon", 45),
                     rs_thermal = 0.06)

  wave.comp = c(data.spectral[['reg1']],data.spectral[['reg2']],data.spectral[['reg3']])

  if (!is.null(rsoil)) {
    # Caller-supplied soil spectrum (e.g. from get.marmit.rsoil()) instead of
    # the built-in BSM model -- must cover the same optical region (reg1,
    # 400-2400nm, 1nm) that BSM would otherwise produce.
    if (length(rsoil) != length(data.spectral[['reg1']])) {
      stop("SPART(): 'rsoil' must have length ", length(data.spectral[['reg1']]),
           " (", min(data.spectral[['reg1']]), "-", max(data.spectral[['reg1']]),
           "nm, 1nm step), got length ", length(rsoil), ".")
    }
    rfl.comp = c(rsoil, rep(data.soil[['rs_thermal']],length(data.spectral[['IwlT']])))
  } else {
    spec<-list()
    spec[['GSV']] <- optipar$GSV #optipar2020.prospectD.BSM2019$GSV
    spec[['Kw']] <-optipar$Kw # water absorption spectrum
    ## this is for get soil with 2500nm
    #spec[['Kw']] <-ref.optipar$spA_Cw[1:10] # water absorption spectrumr
    ##optipar2020.prospectD.BSM2019$nw
    spec[['nw']] <- optipar$nw # water refraction index spectrum


    #print('getting  soil reflectance using a Brightness-Shape-Moisture soil model model (BSM) ...')

    rsoil.BSM <- getBSM.toolsRTM(soilpar=data.soil,spec=spec,emp = soilemp);

    rfl.comp = c(c(rsoil.BSM), rep(data.soil[['rs_thermal']],length(data.spectral[['IwlT']])))
  }

  db.rsoil <- data.frame(wave=wave.comp, rfl.soil = rfl.comp)


  ##::::::::::::::::::::::::::::::::::::::::::
  # Plot for Soil reflectance using BSM model
  ##:::::::::::::::::::::::::::::::::::::::::

  if (get.plots ==  T){

    plot.soil <- ggplot2::ggplot(data = db.rsoil, aes(x = wave, y = rfl.soil)) +
      labs(y= "soil reflectance", x = "")+ xlim(400,2400) +
      geom_line() + theme_bw()

    print(plot.soil)

  }

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

  ##::::::::::::::::::::::::::::::::::::::::::
  # Plot for Extraterrestrial irradiance
  ##:::::::::::::::::::::::::::::::::::::::::

  if (get.plots ==  T){

    plot.Irrad <- ggplot2::ggplot(data = df.irradiance, aes(x = wave, y = irrad)) +
      labs(y= " Extraterrestrial irradiance", x = "")+ xlim(400,2450) +
      geom_line() + theme_bw()

    print(plot.Irrad)

  }

  La_extra <-  get.spectral.convolution(df.irradiance = df.irradiance,sensor.i, get.plots=F)
  # Calculate TOA radiance
  L_TOA <- La_extra$Eo * R_TOA

  df.rad.toa<- data.frame(wave = La_extra$wave, rad.toa = L_TOA, rfl.toa = R_TOA)


  if (get.plots ==  T){

    plot.rad <- ggplot2::ggplot(data = df.rad.toa, aes(x = wave)) +
      labs(y= " TOA rfl-rad", x = "") +
      geom_line(aes(y = rad.toa, color = "Rad")) +
      geom_line(aes(y = rfl.toa, color = "RFL")) + theme_bw() +
      guides(color = guide_legend(title = "TOA"), linetype = guide_legend(title = "TOA"), shape = guide_legend(title = "TOA"))


    print(plot.rad)

  }

  df.refl <- data.frame(wave = wlSensor,rfl.toa = R_TOA, rad.toa = L_TOA,rfl.toc = R_TOC,rfl.toc.brdf= rfl.canopy.brdf.sensor )

  ##::::::::::::::::::::::::::::::::::::::::::
  # Plot for Reflectance at TOA, TOC
  ##:::::::::::::::::::::::::::::::::::::::::

  if (get.plots ==  T){

    plot.rfl.toa <- ggplot2::ggplot(data = df.refl, aes(x = wave)) +
      labs(y = "reflectance", x = "") +
      geom_point(aes(y = rfl.toa, color = "TOA rfl."), size = 1) +
      geom_line(aes(y = rfl.toa, color = "TOA rfl.")) +

      #geom_point(aes(y = rad.toa, color = "TOA (rad.)"), size = 1) +
      #geom_line(aes(y = rad.toa, color = "TOA (rad.)")) +

      geom_point(aes(y = rfl.toc, color = "TOC rfl. (SMAC)"), size = 1) +
      geom_line(aes(y = rfl.toc, color = "TOC rfl. (SMAC)")) +
      geom_point(aes(y = rfl.toc.brdf, color = "TOC rfl. (BRDF)"), size = 1) +
      geom_line(aes(y = rfl.toc.brdf, color = "TOC rfl. (BRDF)")) +
      theme_bw() +
      guides(color = guide_legend(title = "Reflectance:"), linetype = guide_legend(title = "Reflectance:"), shape = guide_legend(title = "Reflectance:")) +

      theme(legend.position = "right")

    print(plot.rfl.toa)

  }



  df.all <- data.frame(wave = wlSensor, rad.toa = L_TOA, rfl.toa = R_TOA, rfl.toc = R_TOC, rfl.toc.BRDF= rfl.canopy.brdf.sensor)
  rownames(df.all) <-NULL
  outputs <- list(output = df.all, Irrad.toa = La_extra, rfl.toc.brdf = df.toc,
                  rf.soil = db.rsoil, sensor = df.sensor,
                  inputLUT=inputLUT)

  return(outputs)


}
