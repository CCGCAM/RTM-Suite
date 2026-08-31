#' \code{get.SCOPE} get simulations based on SCOPE model
#'
#' @param LUT leaf, biochemistry, viewing angles,meteo, and  canopy properties needed for running SCOPE model.
#' @param n.LUT number of input used for simulating with SCOPE model
#' @param options.SCOPE simulation options controlling the SCOPE run, such as the atmospheric
#' correction and whether to correct the total irradiance if a specific value is provided instead of the usual Modtran output.
#' @param optipar leaf optical parameters used by the selected leaf model (e.g. Fluspect refractive index and specific absorption coefficients).
#' @param path.out folder for saving the SCOPE outputs.
#' @param canopy.model Not currently functional -- SCOPE always uses its own multi-layer canopy engine (\code{getRTMo}), which is not the same code as \code{ToolsRTM::foursail()}/\code{foursail2()}/\code{inform()} and isn't swappable without reimplementing those inside SCOPE's multi-layer/gap-probability framework. Passing anything other than \code{'fourSAIL'} triggers a warning. Kept as an argument for forward compatibility.
#' @param leaf.model Not currently functional -- SCOPE always uses its own multi-layer Fluspect-Cx leaf optics (\code{get.fluspect_mSCOPE}), required by the fluorescence pipeline (\code{RTMf}). Passing anything other than \code{'fluspect-CX'} triggers a warning. Kept as an argument for forward compatibility.
#' @param get.outputs  if get.outputs = 'ALL' all variables were retrieved; if get.outputs = 'Main' the main variables were retrieved.
#' By default SCOPE uses get.outputs = 'ALL', for processiong huge LUT is recommended get.outputs = 'Main'.
#' @param get.plots  is true plot the intermediate plots
#' @param rsoil Optional externally-supplied soil reflectance spectrum, length matching
#' \code{data.spectral$reg1} (400-2400nm, 1nm steps, 2001 values) -- e.g. from
#' \code{ToolsRTM::get.marmit.rsoil()} driven by \code{STEMMUS-R::get.stemmus()}'s simulated
#' soil moisture. When supplied, bypasses the file-lookup/BSM soil branches entirely (whichever
#' \code{options.SCOPE} would otherwise select) and uses this spectrum directly. Defaults to
#' \code{NULL}, preserving the original file/BSM-resolved behaviour.
#' @description
#' Main module for calling the SCOPE model
#'
#' @references
#' Yang, P., E. Prikaziuk, W. Verhoef, and C. van der Tol. 2020.
#' “SCOPE 2.0: A Model to Simulate Vegetated Land Surface Fluxes and Satellite Signals.”
#' Geoscientific Model Development Discussions 2020: 1–26. https://doi.org/10.5194/gmd-2020-251.
#'
#' Van der Tol, C., W. Verhoef, J Timmermans, A Verhoef, and Z Su. 2009.
#' “An Integrated Model of Soil-Canopy Spectral Radiances, Photosynthesis, Fluorescence, Temperature and Energy Balance.”
#' Biogeosciences 6 (12): 3109–29. https://doi.org/10.5194/bg-6-3109-2009.
#'
#' @author 	Wout Verhoef, Prikaziuk, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#' @export
#'
#' @examples
#' \dontrun{
#' out <- get.SCOPE(LUT, n.LUT = 10, options.SCOPE, optipar, path.out,
#'                   leaf.model = "fluspect-CX", canopy.model = "fourSAIL",
#'                   get.outputs = "Main", get.plots = FALSE)
#' }

get.SCOPE<-function(LUT, n.LUT=10,options.SCOPE,
                    optipar,path.out,
                    leaf.model='fluspect-CX',canopy.model='fourSAIL',
                    get.outputs,
                    get.plots = T, rsoil = NULL)  {

  # Alias immediately: the function body below reuses the name `rsoil` internally
  # (a local data.frame from the file-based soil branch) which would otherwise
  # shadow this parameter before the override check runs.
  rsoil.external <- rsoil

  message('SCOPE 2.1. version will be executed .............')

  ##################################################################################
  ### 0.0.1 leaf.model / canopy.model are NOT wired into the simulation (yet)
  ##################################################################################
  # SCOPE always runs its own multi-layer leaf optics (get.fluspect_mSCOPE, needed
  # for the Mb/Mf fluorescence matrices RTMf depends on) and its own multi-layer
  # canopy engine (getRTMo) -- neither is the same code as ToolsRTM::foursail()/
  # foursail2()/inform(), which have a different interface and don't produce the
  # vertical-profile/gap-probability structures the rest of the SCOPE pipeline
  # (ebal, RTMf, RTMz, RTMt.sb) requires. So these two parameters are currently
  # accepted but have no effect -- warn loudly instead of silently ignoring a
  # value the caller explicitly asked for, so nobody assumes a different leaf/
  # canopy model actually ran.
  if (!missing(canopy.model) && !identical(canopy.model, 'fourSAIL')) {
    warning("get.SCOPE(): canopy.model = '", canopy.model, "' was requested, but SCOPE's ",
            "canopy engine (getRTMo) is not swappable yet -- the simulation always uses ",
            "SCOPE's own multi-layer RTMo, not ToolsRTM::", canopy.model, "(). This argument ",
            "is currently a no-op; see get.SCOPE.R for details.", call. = FALSE)
  }
  if (!missing(leaf.model) && !(leaf.model %in% c('fluspect-CX', 'fluspect-Cx'))) {
    warning("get.SCOPE(): leaf.model = '", leaf.model, "' was requested, but SCOPE always ",
            "runs its own multi-layer Fluspect-Cx leaf optics (get.fluspect_mSCOPE), which ",
            "the fluorescence pipeline (RTMf) depends on -- other leaf models are not ",
            "swappable yet. This argument is currently a no-op; see get.SCOPE.R for details.",
            call. = FALSE)
  }

  ##################################################################################
  ### 0.0 Load optipar (optical leaf properties)
  ##################################################################################

  if (missing(optipar)){
    optipar =  SCOPEinR::optipar2017.ProspectD
  } else{

    optipar =  optipar

  }

  ##################################################################################
  ### 0.2 Check for plot or not
  ##################################################################################

  if (missing(get.plots)){
    get.plots = FALSE
  }

  ##################################################################################
  ### 0.3 Check the canopy model
  ##################################################################################


  if (missing(canopy.model)){
    model = 'fourSAIL'

  } else {
    model = canopy.model
  }

  ##################################################################################
  ### 0.4 Check the Leaf model
  ##################################################################################

  if (missing(leaf.model)){
    model = 'fluspect-CX'

  } else {
    model = leaf.model
  }


  ##################################################################################
  ### 1. Load main variables
  ##################################################################################


  # Get the path to the 'input' folder in your package's installation directory
  #path_input <- 'input/'
  path_input <- paste(system.file("input", package = "SCOPEinR"),'/',sep='')

  if (missing(path.out)) {
    path_outs <- 'outs/'
    message('folder outputs will be used for the SCOPE outputs')

    if (!dir.exists(path_outs)) {
      dir.create(path_outs)
    }
  } else {
    path_outs <- path.out
    if (!dir.exists(path_outs)) {
      dir.create(path_outs)
    }
  }


  constants <- SCOPEinR::constants

  data.spectral<-SCOPEinR::get.spectra.SCOPE(getSpectral = T)


  ##################################################################################
  ### 2. Get options for SCOPE model
  ##################################################################################

  if (missing(options.SCOPE)){
    stop('please use and options file for running the model ')
  }
  #Get options for running SCOPE mode
  data.opts =  options.SCOPE
  options.lite                <- data.opts[1,]   # Value 1 indicates the SCOPE will use the lite SCOPEversion
  options.calc_fluor          <-  data.opts[2,]  # Value 1 calculate chlorophyll fluorescence in observation direction
  options.calc_planck         <-  data.opts[3,]  # Value 1 calculate spectrum of thermal radiation
  options.calc_xanthophyllabs <-  data.opts[4,]  # Value 1 include simulation of reflectance dependence on de-epoxydation state
  options.soilspectrum        <-  data.opts[5,]  # Value0: use soil reflectance from file; 1: calculate soil reflectance with BSM
  options.Fluorescence_model  <-  data.opts[6,]  # Value 0: empirical, with sustained NPQ (fit to Flexas' data); 1: empirical, with sigmoid for Kn; 2: Magnani 2012 model
  options.apply_T_corr        <-  data.opts[7,]  # Value 1 indicates that SCOPE will correct Vcmax and rate constants for temperature in biochemical.m
  options.verify              <-  data.opts[8,] # Value 1 indicates that SCOPE will check with field data
  options.mSCOPE              <-  data.opts[9,] # Value 1 indicates that SCOPE will use the mSCOPE considering vertical variations in the vegetation canopies
  #options.simulation
  # 2: Lookup-Table (specify the values to be included) not implemeted
  # 3: Lookup-Table with random input (specify the ranges of values) not implemented
  options.simulation          <-  data.opts[10,] # Value 0 the SCOPE will execute by individual runs based on LUT table and Value 1 for time series (uses text files with meteo input as time series)
  options.calc_directional   <-  data.opts[11,]  # Value 1 calculate full BRDF for many angles
  options.calc_vert_profiles <-  data.opts[12,]  # Value 1  indicates that SCOPE will estimate vertical profiles
  options.soil_heat_method     <-  data.opts[13,]    # Value 0 will estimate the GAM parameters with Soil_Inertia0(lambdas); Value 1 will estimate the GAM with Soil_Inertia1(SMC); Value 2 will estimate  G by 0.35*Rn (where always in no TS)
  options.calc_rss_rbs          <-  data.opts[14,] # Value 0 is fixed and Value 1 is calculated

  calc.heat = options.soil_heat_method$Value
  calc.rss_rbs = options.calc_rss_rbs$Value

  options.MoninObukhov         <-  data.opts[15,]  # Value 1 indicates that SCOPE wil use the MoninObukhov correction

  ## News
  options.LIDF              <-  data.opts[16,]     # Value 1 SCOPE wil estimate LIDF from angles file and Value 0 will use LIDFa and LIDFb from LUT table for estimating LIDF
  options.irradiance          <-  data.opts[17,]   # Value 1 SCOPE wil estimate LIDF from angles file and Value 0 will use LIDFa and LIDFb from LUT table for estimating LIDF
  #options.calc_zo_d     <-  options[18,3]         #



  ##################################################################################
  ### 2.3 Get brdf angles
  ##################################################################################

  #  print("Step 1: Initializing variables")


  type_angles <-c('brdf.angles','angles.file')
  type_angles <- type_angles[2]
  if (options.calc_directional$Value == 1){


    if (type_angles == 'angles.file'){


      # Pattern to match the file names
      pattern <- "\\.dat$"
      # Get all files in the directory matching the pattern
      matching_files <- list.files(paste(path_input, "leafangles",sep=''), pattern = pattern)

      # Sort the files and select the first one
      name.file <- sort(matching_files)[1]

      brdf.file <- read.table(file.path(path_input, "leafangles",name.file), sep=',',skip = 3, header = F)  # read file
      angles_ <- c( 10,20,30,40,50,60,70,80,82,84,86,88,90)
      brdf.file <-cbind(angles_,brdf.file)

      data.directional <- list()
      data.directional[['tto']] <-brdf.file[,1]
      data.directional[['ps']] <-brdf.file[,2]
      data.directional[['noa']] <-length(data.directional[['tto']])

    } else if  (type_angles == 'brdf.angles') {

      brdf.file <- SCOPEinR::brdf_angles2

      data.directional <- list()
      data.directional[['tto']] <-brdf.file[,1]
      data.directional[['ps']] <-brdf.file[,1]
      data.directional[['noa']] <-length(data.directional[['tto']])

    } else {
      data.directional <- list()
    }

  }


  ##################################################################################
  ### 2.4 Get Irradiance
  ##################################################################################

  # print("Step 2: Getting Irradiance file ...")

  type_irradiance <-c('measurements WithE','MODTRAN atmospheric.file')

  if (options.irradiance[,3] >= 1) {


    if (options.irradiance[,3] == 1 ){

      type.irradiance = 'measurements WithE'
      # Pattern to match the file names
      pattern <- "\\csv$"
      # Get all files in the directory matching the pattern
      matching_files <- list.files(paste(path_input, "radiationdata",sep=''), pattern = pattern)
      name.file <- sort(matching_files)[1]
      #name.file='example_rad_file.csv'
      irrad_file <- read.table(file.path(path_input, "radiationdata",name.file), sep=',',skip = 0, header = F)  # read file
      atmo = data.frame(wave = data.spectral[['wlIrrad']],Esky_ = irrad_file[,1], Esun_ = irrad_file[,2])

      atmo.sb<-atmo[1:2000,]


      if (get.plots ==  T){

        ggplot(data = atmo.sb, aes(x = wave, y = Esky_)) +
          labs(y= "Irradiance", x = "")+
          geom_line(color='black') + theme_bw() +
          geom_line(aes(x = wave, y = Esun_),color='navyblue')
      }

    } else if (options.irradiance[,3] == 2 ) {
      # read the table
      type.irradiance <- 'MODTRAN atmospheric.file'
      # Pattern to match the file names
      pattern <- "\\.atm$"
      # Get all files in the directory matching the pattern
      matching_files <- list.files(paste(path_input, "radiationdata",sep=''), pattern = pattern)

      # Sort the files and select the first one
      name.file <- sort(matching_files)[1]
      atmo <- read.table(file.path(path_input, "radiationdata",name.file), header = F, skip = 2)


      # set column names
      colnames(atmo) <- c("WN", "WL", paste0("T", 1:18))
      names_ <-c("WN", "WL",'toasun','rso','rdd','tss','tsd','too','tdo','tsstoo','tsdtoo','tsstdo',
                 'tsdtdo','tssrdd','toordd','tssrddtoo','Lat','Lab','Labtoo','Labtdo')
      head(atmo)
      atmo.sb<-atmo[1:2000,]


      if (get.plots ==  T){

        ggplot(data = atmo.sb, aes(x = WN, y = T1)) +
          labs(y= "MODTRAN atmosferic profiles", x = "")+
          geom_line(color='black') + theme_bw() +
          geom_line(aes(x = WN, y = T15),color='navyblue')

      }

    }

  } else {
    Esky_ <-SCOPEinR::Esky_[[1]]
    Esun_ <-SCOPEinR::Esun_[[1]]
    atmo = data.frame(wave = data.spectral[['wlIrrad']],Esky_ = Esky_, Esun_ = Esun_)
    atmo.sb<-atmo[1:2000,]

    if (get.plots ==  T){

      ggplot(data = atmo.sb, aes(x = wave)) +
        labs(y= "Irradiance", x = "") +
        geom_line(aes(y = Esun_, color = "Direct solar"), linewidth = 0.5) +
        geom_line(aes(y = Esky_, color = "Diffuse sky"), linewidth = 0.5) +
        scale_color_manual(name = "Irradiance type",
                           values = c("Direct solar" = "black", "Diffuse sky" = "navyblue")) +
        theme_bw() + xlim(400,3000) +
        theme(legend.position = "top") +
        guides(color = guide_legend(title = NULL))
    }

  }


  ##################################################################################
  ### 2.5. Get Soil reflectance
  ##################################################################################

  # print("Step 3: Getting soil spectrum ...")

  if (options.soilspectrum[,3] == 0) {  # check if options.soilspectrum option is true


    ## adding a txt with rfl or selected spectra of soil from SCOPEinR (same data)
    # Pattern to match the file names
    pattern <- "\\.txt$"
    # Get all files in the directory matching the pattern
    matching_files <- list.files(paste(path_input, "soil_spectra",sep=''), pattern = pattern)

    # Sort the files and select the first one
    name.file <- sort(matching_files)[1]

    soil_file <- read.table(file.path(path_input, "soil_spectra",name.file), skip = 0, header = F)  # read file
    colnames(soil_file) <- c('wave','rfl.soil1','rfl.soil2','rfl.soil3')

    rfl.soil.formSCOPE <-LUT[1,'spectrum'] +1
    rsoil <- soil_file[,c(1,rfl.soil.formSCOPE)]
    colnames(rsoil) <- c('wave','rfl.soil')

    wave.comp = c(rsoil$wave,data.spectral[['reg2']],data.spectral[['reg3']])
    rfl.comp = c(rsoil$rfl.soil, rep(LUT[['rs_thermal']][1],length(data.spectral[['IwlT']])))


    if (get.plots ==  T){

      db.rsoil <- data.frame(wave=wave.comp, rfl.soil = rfl.comp)

      plot.soil<- ggplot(data = db.rsoil, aes(x = wave, y = rfl.soil)) +
        labs(y= "reflectance", x = "")+ xlim(400,2499) +
        geom_line() + theme_bw()
      print(plot.soil)
    }

  } else {

    #get soil spectra from BSM model
    # Parameters for BSM model
    soilemp<-list()
    soilemp[['SMC']] = 25         # empirical parameter (fixed) for BSM (25, original value in SCOPE )
    soilemp[['film']] = 0.015    #empirical parameter (fixed) for BMS
    soilemp[['SMp']] = 15    #soil moisture volume percentage (5 - 55)



    data.soil<-SCOPEinR::getinputLUT(inputLUT=LUT[1,], dataset='soil',
                                     calc.heat = calc.heat,
                                     calc.rss_rbs =  calc.rss_rbs)
    spec<-list()
    spec[['GSV']] <- optipar$GSV #optipar2020.prospectD.BSM2019$GSV

    ##optipar2020.prospectD.BSM2019$Kw
    spec[['Kw']] <-optipar$Kw # water absorption spectrum
    ##optipar2020.prospectD.BSM2019$nw
    spec[['nw']] <- optipar$nw # water refraction index spectrum

    #print('getting  soil reflectance using a Brightness-Shape-Moisture soil model model (BSM) ...')

    rsoil.BSM <- SCOPEinR::getBSM(soilpar=data.soil,spec=spec,emp = soilemp);

    wave.comp = c(data.spectral[['reg1']],data.spectral[['reg2']],data.spectral[['reg3']])
    rfl.comp = c(c(rsoil.BSM), rep(data.soil[['rs_thermal']],length(data.spectral[['IwlT']])))

    if (get.plots ==  T){


      db.rsoil <- data.frame(wave=wave.comp, rfl.soil = rfl.comp)

      plot.soil <- ggplot(data = db.rsoil, aes(x = wave, y = rfl.soil)) +
        labs(y= "reflectance", x = "")+ xlim(400,2499) +
        geom_line() + theme_bw()

      print(plot.soil)

    }

  }


  ##################################################################################
  ### 3. Run the model
  ##################################################################################

  if (missing(n.LUT)){
    nSamples = dim(LUT)[1]

  } else {
    nSamples = n.LUT
  }

  # print("Step 4: Running Simulations...")


  # Initialize progress bar
  progress_bar <- progress::progress_bar$new( format = "[:bar] :percent ETA: :eta", total = nSamples)


  db.sim <- list()
  for (sim.i in  c(1:nSamples)){

    # Update progress bar

    progress_bar$tick()

    # Get main inputs
    data.leafbio<-SCOPEinR::getinputLUT(inputLUT=LUT[sim.i,], dataset='leafbio')
    data.canopy<-SCOPEinR::getinputLUT(inputLUT=LUT[sim.i,], dataset='canopy')
    data.meteo<-SCOPEinR::getinputLUT(inputLUT=LUT[sim.i,], dataset='meteo')

    data.soil<-SCOPEinR::getinputLUT(inputLUT=LUT[sim.i,], dataset='soil',
                                     calc.heat = calc.heat,
                                     calc.rss_rbs =  calc.rss_rbs)

    data.timeseries<-SCOPEinR::getinputLUT(inputLUT=LUT[sim.i,], dataset='timeseries') ## this is xyt in SCOPE
    data.angles<-SCOPEinR::getinputLUT(inputLUT=LUT[sim.i,], dataset='angles')
    data.mly<-SCOPEinR::getinputLUT(inputLUT=LUT[sim.i,], dataset='mly')


    ##################################################################################
    ### 3.1 Checking LIDF angles
    ##################################################################################


    if (options.LIDF[,3] == 1) {  # check if LIDF_file option is true
      pattern <- "\\.dat$"
      # Get all files in the directory matching the pattern
      matching_files <- list.files(paste(path_input, "leafangles",sep=''), pattern = pattern)
      name.file <- sort(matching_files)[1]

      #name.file='soybean_wildtype_leafangles.dat'
      LIDF_file <- scan(file.path(path_input, "leafangles",name.file), what = numeric(), skip = 3, quiet = TRUE)  # read file
      data.canopy[['lidf']]  <- LIDF_file
    } else {
      canopy_angles <- leafangles(data.canopy[['LIDFa']],data.canopy[['LIDFb']])  # % This is 'ladgen' in the original SAIL model,
      data.canopy[['lidf']] <-canopy_angles$lidf
    }
    ##################################################################################
    ### 3.2. Define canopy structure and other 'fixed' parameters
    ##################################################################################

    data.canopy[['nlincl']] = 13
    data.canopy[['nlazi']] = 36
    data.canopy[['litab']] = c(seq(5,75,10),seq(81,89,2))   # a column, never change the angles unless 'ladgen' is also adapted
    data.canopy[['lazitab']]  = seq(5,355,10)            #a row


    ##################################################################################
    ### 3.3. Get zo and d from
    ##################################################################################

    ## Get zo and d from  zo_and_d when the zo and d is not in the input Table
    # zo and d will estimate them from LAI, CR, CD1, Psicor, and CSSOIL
    if (options.simulation[1,'Value'] == 0) {

    #if (!("d" %in% colnames(LUT)) & !("zo" %in% colnames(LUT))) {
      options.soil_heat_method     <-  data.opts[13,]    # Value 0 will estimate the GAM parameters with Soil_Inertia0(lambdas); Value 1 will estimate the GAM with Soil_Inertia1(SMC); Value 2 will estimate  G by 0.35*Rn (where always in no TS)
      options.calc_rss_rbs          <-  data.opts[14,] # Value 0 is fixed and Value 1 is calculated


      parameters.zo_d  <- get.zo_and_d(inputLUT=LUT[sim.i,],constants,
                                       calc.heat=calc.heat,calc.rss_rbs=calc.rss_rbs)


      data.canopy[['zo']]  <-    parameters.zo_d[['zom']]
      data.canopy[['d']]  <-    parameters.zo_d[['d']]
    }

    ##################################################################################
    ### 3.4. soil heat
    ##################################################################################

    if (options.simulation[1,'Value'] == 1) {

      if (options.soil_heat_method[1,'Value'] < 2) {

        if ((is.null(data.meteo[['Ta']]) || data.meteo[['Ta']] == -273)){

          data.meteo[['Ta']] = 20
        }

        data.soil[['Tsold']] = matrix(rep(data.meteo$Ta,12),nrow = 12,ncol = 2)

      }
    }
    ##################################################################################
    ### 3.5. temperature sensitivity of photosynthesis parameters
    ##################################################################################

    ## temperature sensitivity of photosynthesis parameters
    data.leafbio[['TDP']] = define_temp_response_biochem(getTDP = T) # temperature response C3 and C4 according to CLM4 model

    ##################################################################################
    ### 3.6.  Set the canopy nlayers and canopy xl
    ##################################################################################


    # Set canopy.nlayers and canopy.xl
    data.canopy[['nlayers']] <- ceiling(10 * data.canopy[['LAI']]) + ifelse((data.meteo[['Rin']]  < 200) & options.MoninObukhov$Value, 60, 0)
    data.canopy[['nlayers']] <- max(2,data.canopy[['nlayers']])
    nl <- data.canopy[['nlayers']]
    x <- seq(-1/nl, -1, length.out = nl) # a column vector
    data.canopy[['xl']] <- c(0, x) # add top level
    # canopy.xl[-nl] <- canopy.xl[-nl] + canopy.xl[-nl] - 1 / (2 * nl) # middle of the thin layer


    ##################################################################################
    ### 3.7.  Get models for temporal model (not implemented )
    ##################################################################################

    if (options.simulation$Value == 1) {
      ## this is the mly with time series (not implemented
      ## check  load_mSCOPE_ts.M; load_timeseries.M and SCOPE.m

      mly_ts <- list()
      data.mly = mly_ts ## this is the mly with time series (not implemented)
    }

    if (options.mSCOPE$Value == 1){
      data.mly[['nly']] <- 3
      data.mly[['pLAI']] <- c(0.5,	1,	1.5)
      data.mly[['pCab']] <- c(40,	20,	80)
      data.mly[['pCar']] <- c(20,	20,	20)
      data.mly[['pLMA']] <- c(0.012, 0.012,	0.012)
      data.mly[['pEWT']] <- c(0.009,	0.009,	0.009)
      data.mly[['pCs']] <-  c(0,	0,	0)
      data.mly[['pN']] <- rep(1.4,3)
      data.mly[['totLAI']] <- 3	#sum(  mly[['pLAI ']] )

    }


    ##################################################################################
    ### 3.8.  Get parameters for Lite configurations
    ##################################################################################

    if (options.lite$Value == 0) {
      ## For aggregating layers
      integrate.layer ='angles_and_layers'
    } else {

      integrate.layer ='layers'
    }


    ##############################################################
    ## 4. Get Simulations with RTMs
    ##############################################################

    ##############################################################
    # 4.1.1 leaf radiative transfer model FLUSPECT
    ##############################################################

    data.leafbio[['emis']] <- 1- data.leafbio[['rho_thermal']] - data.leafbio[['tau_thermal']]
    data.leafbio[['Cx']]   <- 0

    #leafopt <-getFluspect.Cx(inputsLeaf=leafbio,inputsOptipar=ToolsRTM::optipar, version = 'SCOPE')
    #print('getting  get.fluspect_mSCOPE ...')

    if (get.plots == T) {

      data.leafopt = get.fluspect_mSCOPE(mly=data.mly,spectral=data.spectral,leafbio=data.leafbio, soil=data.soil,
                                         optipar=optipar, nl,step=5,get.plots=T)
    } else {

      data.leafopt = get.fluspect_mSCOPE(mly=data.mly,spectral=data.spectral,leafbio=data.leafbio, soil=data.soil,
                                         optipar=optipar, nl,step=5,get.plots=F)
    }


    ######################################################################################
    # 4.1.2 leaf radiative transfer model FLUSPECT with options.calc_xanthophyllabs
    ######################################################################################


    if (options.calc_xanthophyllabs$Value == 1) {
      data.leafbio[['Cx']] <- 1
      #print('getting  get.fluspect_mSCOPE with Xantophyll ...')
      if (get.plots == T) {
        data.leafoptZ = get.fluspect_mSCOPE(data.mly,data.spectral,data.leafbio,soil=data.soil,optipar=optipar, nl, step=5,get.plots=T)
      } else{
        data.leafoptZ = get.fluspect_mSCOPE(data.mly,data.spectral,data.leafbio,soil=data.soil,optipar=optipar, nl, step=5,get.plots=F)

      }
      ## to save the values from 2163 = 0.01
      data.leafopt$reflZ <- data.leafoptZ$refl
      data.leafopt$tranZ <- data.leafoptZ$tran


      if (get.plots == T) {

        wave.fluspect = c(data.spectral[['reg1']],data.spectral[['reg2']],data.spectral[['reg3']])
        rfl.fuspectZ = c(data.leafopt$reflZ[1,], rep(data.leafbio[['rho_thermal']],length(data.spectral[['IwlT']])))
        df.rfl <- data.frame(wave.fluspect=wave.fluspect, rfl.fuspectZ = rfl.fuspectZ)

        p.z1 <- ggplot(data = df.rfl, aes(x = wave.fluspect, y = rfl.fuspectZ)) +
          labs(y= " leaf reflectance Z1 (fluspect-Cx)", x = "")+ xlim(400,2499) +
          geom_line() + theme_bw()

        print(p.z1)



      }

    } # endoptions.calc_xanthophyllabs

    ######################################################################################
    # 4.2. soil reflectance using a Brightness-Shape-Moisture soil model  (BSM)
    ######################################################################################


    if (options.soilspectrum[,3] == 0) {  # check if options.soilspectrum option is true


      data.soil[['rfl.soil.complete']] <- rfl.comp
      data.soil[['rfl.soil']] <-  (rsoil$rfl.soil)


      if (get.plots == T) {

        db.rsoil <- data.frame(wave=wave.comp, rfl.soil = rfl.comp)

        plot.soil <-ggplot(data = db.rsoil, aes(x = wave, y = rfl.soil)) +
          labs(y= "reflectance", x = "")+ xlim(400,2499) + ylim(0.05,0.45) +
          geom_line() + theme_bw()

        print(plot.soil)

      }

    } else {

      #get soil spectra from BSM model
      # Parameters for BSM model
      soilemp<-list()
      soilemp[['SMC']] = 25         # empirical parameter (fixed) for BSM (25, original value in SCOPE )
      soilemp[['film']] = 0.015    #empirical parameter (fixed) for BMS
      soilemp[['SMp']] = 15    #soil moisture volume percentage (5 - 55)

      data.soil<-SCOPEinR::getinputLUT(inputLUT=LUT[sim.i,], dataset='soil',
                                       calc.heat = calc.heat,
                                       calc.rss_rbs =  calc.rss_rbs)

      spec<-list()

      spec[['GSV']] <- optipar$GSV #optipar2020.prospectD.BSM2019$GSV

      ##optipar2020.prospectD.BSM2019$Kw
      spec[['Kw']] <- optipar$Kw # water absorption spectrum
      ##optipar2020.prospectD.BSM2019$nw
      spec[['nw']] <-  optipar$nw # water refraction index spectrum

      #print('getting  soil reflectance using a Brightness-Shape-Moisture soil model model (BSM) ...')

      rsoil.BSM <- getBSM(soilpar=data.soil,spec=spec,emp = soilemp)


      wave.comp = c(data.spectral[['reg1']],data.spectral[['reg2']],data.spectral[['reg3']])
      rfl.comp = c(c(rsoil.BSM), rep(data.leafbio[['rho_thermal']],length(data.spectral[['IwlT']])))


      data.soil[['rfl.soil.complete']] <- rfl.comp
      data.soil[['rfl.soil']] <- c(rsoil.BSM)

      if (get.plots == T) {

        db.rsoil <- data.frame(wave=wave.comp, rfl.soil = rfl.complete,   emisi.soil = 1-rfl.comp)


        p.soil <-ggplot(data = db.rsoil, aes(x = wave, y = rfl.soil)) +
          labs(y= "soil reflectance / soil emissivity", x = "")+ xlim(400,2499) +
          geom_line(color='navyblue')+ theme_bw() +
          geom_line(aes(x = wave, y = emisi.soil), color='forestgreen')

        print(p.soil)
      }

    }

    ######################################################################################
    # 4.2b. optional externally-supplied soil reflectance (e.g. STEMMUS+MARMIT), bypassing
    # the file-lookup/BSM branches above entirely when supplied
    ######################################################################################

    if (!is.null(rsoil.external)) {
      n_reg1 <- length(data.spectral[['reg1']])
      if (length(rsoil.external) != n_reg1) {
        stop("get.SCOPE(): rsoil must have length ", n_reg1,
             " (data.spectral$reg1, 400-2400nm at 1nm steps), got length ", length(rsoil.external), ".")
      }
      data.soil[['rfl.soil']] <- rsoil.external
    }



    #####################################################################################################
    # 4.3.  the four stream canopy radiative transfer model for incident radiation model  (fourSAIL)
    ################################################################################################

    ## here adding thermal profile to spectra between 400:2500

    # Calling the function with warnings suppressed
    original_warn <- options(warn = -1)  # Set warn option to -1 to suppress warnings for plots
    # Restore it however this function exits (including on error), so the
    # warning suppression doesn't leak into the rest of the caller's session.
    on.exit(options(original_warn), add = TRUE)

    #print('getting fourSAIL model (RTMo) ...')
    if (get.plots == T) {

      outs.RTMo<-getRTMo(data.spectral,atmo,data.soil,data.leafopt,data.canopy,data.leafbio,
                         data.angles,data.meteo,data.opts=data.opts,get.plots=T)
    } else {

      outs.RTMo<-getRTMo(data.spectral,atmo,data.soil,data.leafopt,data.canopy,data.leafbio,
                         data.angles,data.meteo,data.opts=data.opts,get.plots=F)
    }


    ######################################################################################
    # 4.4.  the balance model (ebal)
    ######################################################################################

    data.gap <-outs.RTMo[['data.gap']]
    data.rad <-outs.RTMo[['data.rad']]
    data.profiles <-outs.RTMo[['data.profiles']]
    data.canopy <-outs.RTMo[['data.canopy']]

    data.timeseries<-SCOPEinR::getinputLUT(inputLUT=LUT[sim.i,], dataset='timeseries')

    #print('getting balance model (ebal) ...')

    if (get.plots == T) {

      ## Here is not working the get.plot because is so slow, get.plots is force to FALSE
      outs.eba <- get.ebal(data.rad, data.gap, data.meteo, data.soil, data.canopy, data.leafbio,data.leafopt,
                           data.spectral, data.opts=data.opts,
                           integrate.layer=integrate.layer, k.maxit=100, get.plots=F)
    } else {

      outs.eba <- get.ebal(data.rad, data.gap, data.meteo, data.soil, data.canopy, data.leafbio,data.leafopt,
                           data.spectral, data.opts=data.opts,
                           integrate.layer=integrate.layer, k.maxit=100, get.plots=F)
    }


    iter <- outs.eba[['iter']]
    data.rad <- outs.eba[['data.rad']]
    data.thermal <- outs.eba[['data.thermal']]
    data.soil <- outs.eba[['data.soil']]
    data.bcu <- outs.eba[['data.bcu']]
    data.bch <- outs.eba[['data.bch']]
    data.fluxes <- outs.eba[['data.fluxes']]
    resist_out <- outs.eba[['resist_out']]
    data.meteo <- outs.eba[['data.meteo']]

    ######################################################################################
    # 4.5.  the fluorescence radiative transfer model (RMTf)
    ######################################################################################


    if (options.calc_fluor$Value == 1){
      #print('getting fluorescence emission (RTMf) ...')

      if (get.plots == T) {

        data.rad <- get.RTMf(data.spectral,data.rad,data.soil,data.leafopt,
                             data.canopy,data.gap,data.angles,
                             #relative fluorescence emission efficiency for sunlit leaves
                             data.etau=data.bcu[['eta']],
                             #relative fluorescence emission efficiency for shaded leaves
                             data.etah =data.bch[['eta']],get.plots=T)
      } else {

        data.rad <- get.RTMf(data.spectral,data.rad,data.soil,data.leafopt,
                             data.canopy,data.gap,data.angles,
                             #relative fluorescence emission efficiency for sunlit leaves
                             data.etau=data.bcu[['eta']],
                             #relative fluorescence emission efficiency for shaded leaves
                             data.etah =data.bch[['eta']],get.plots=F)
      }

    }


    ########################################################################################################
    # 4.6.  the radiative transfer model for the Violaxanthin into Zeaxanthin effect in leaves (RTMz)
    ########################################################################################################


    if (options.calc_xanthophyllabs$Value == 1){

      #print('running the RTMz model ...')

      if (get.plots == T) {

        data.rad <- get.RTMz(data.spectral,data.rad,data.soil,data.leafopt,
                             data.canopy,data.gap,data.angles,
                             data.Knu = data.bcu[['Kn']],data.Knh = data.bch[['Kn']],
                             get.plots=T)

      } else {


        data.rad <- get.RTMz(data.spectral,data.rad,data.soil,data.leafopt,
                             data.canopy,data.gap,data.angles,
                             data.Knu = data.bcu[['Kn']],data.Knh = data.bch[['Kn']],
                             get.plots=F)
      }

    }

    #print('running the RTMt.sb model ...')

    if (get.plots == T) {

      outs.RMTt.sb <- get.RTMt.sb( data.rad,data.soil,data.leafbio,data.canopy,data.leafopt,data.gap,
                                   Tcu = data.thermal[['Tcu']],  Tch = data.thermal[['Tch']],
                                   Tsu = data.thermal[['Tsu']],  Tsh = data.thermal[['Tsh']],
                                   obsdir=1,data.spectral,
                                   data.opts=data.opts,get.plots=T)

    } else {

      outs.RMTt.sb <- get.RTMt.sb( data.rad,data.soil,data.leafbio,data.canopy,data.leafopt,data.gap,
                                   Tcu = data.thermal[['Tcu']],  Tch = data.thermal[['Tch']],
                                   Tsu = data.thermal[['Tsu']],  Tsh = data.thermal[['Tsh']],
                                   obsdir=1,data.spectral,
                                   data.opts=data.opts,get.plots=F)
    }

    ######################################################################################
    # 4.7.  the radiative transfer model for PRI effects (RTMt_planck)
    ######################################################################################



    if (options.calc_planck$Value == 1){


      #print('getting PRI effects (RTMt.planck) ...')

      if (get.plots == T) {

        data.rad <- get.RTMt.planck(data.spectral,data.rad,data.soil,
                                    data.leafopt,
                                    data.canopy,data.gap,Tcu=data.thermal[['Tcu']],Tch=data.thermal[['Tch']],
                                    Tsu=data.thermal[['Tsu']],Tsh=data.thermal[['Tsh']],
                                    get.plots=T)
      } else {

        data.rad <- get.RTMt.planck(data.spectral,data.rad,data.soil,
                                    data.leafbio, data.leafopt,
                                    data.canopy,data.gap,Tcu=data.thermal[['Tcu']],Tch=data.thermal[['Tch']],
                                    Tsu=data.thermal[['Tsu']],Tsh=data.thermal[['Tsh']],
                                    get.plots=F)
      }

    }


    ######################################################################################
    # 5. Computation of data products
    #
    #    - aPAR, LST, NPQ, ETR, photosynthesis, SIF-reabsorption correction
    #    - aPAR [umol m-2 s-1, total canopy and total chlorphyll]th radiative transfer model for PRI effects (RTMt_planck)
    #
    ######################################################################################

    Ps = data.gap[['Ps']][1:nl]
    Ph = (1-Ps)


    data.canopy[['LAIsunlit']] <-  data.canopy[['LAI']] * mean(Ps, na.rm = T)
    data.canopy[['LAIshaded']] <-  data.canopy[['LAI']] - data.canopy[['LAIsunlit']]

    integrate.layer = 'layers'

    ######################################

    # net PAR Cab sunlit leaves (photons)
    data.canopy[['Pnsun_Cab']] <- data.canopy[['LAI']] * meanleaf.v2(data.canopy, data.rad[['Pnu_Cab']], canopy.choice=integrate.layer, Ps)

    # net PAR Cab shaded leaves (photons)
    data.canopy[['Pnsha_Cab']] <- data.canopy[['LAI']] * meanleaf.v2(data.canopy, data.rad[['Pnh_Cab']], canopy.choice='layers', Ph)

    # net PAR Cab leaves (photons)
    data.canopy[['Pntot_Cab']] <- data.canopy[['Pnsun_Cab']] + data.canopy[['Pnsha_Cab']]

    ######################################

    # net PAR Cab sunlit leaves (photons)
    data.canopy[['Pnsun_Car']] <- data.canopy[['LAI']] * meanleaf.v2(data.canopy, data.rad[['Pnu_Car']], canopy.choice=integrate.layer, Ps)

    # net PAR Cab shaded leaves (photons)
    data.canopy[['Pnsha_Car']] <- data.canopy[['LAI']] * meanleaf.v2(data.canopy, data.rad[['Pnh_Car']], canopy.choice='layers', Ph)

    # net PAR Cab leaves (photons)
    data.canopy[['Pntot_Car']] <-  data.canopy[['Pnsun_Car']] + data.canopy[['Pnsha_Car']]

    ######################################

    # net PAR sunlit leaves (photons)
    data.canopy[['Pnsun']] <- data.canopy[['LAI']] * meanleaf.v2(data.canopy, data.rad[['Pnu']], canopy.choice=integrate.layer, Ps)

    # net PAR shaded leaves (photons)
    data.canopy[['Pnsha']] <- data.canopy[['LAI']] * meanleaf.v2(data.canopy, data.rad[['Pnh']], canopy.choice='layers', Ph)

    # net PAR leaves (photons)
    data.canopy[['Pntot']] <-  data.canopy[['Pnsun']] + data.canopy[['Pnsha']]

    ######################################

    # net PAR Cab sunlit leaves  (radiance)
    data.canopy[['Rnsun_Cab']] <- data.canopy[['LAI']] * meanleaf.v2(data.canopy, data.rad[['Rnu_Cab']], canopy.choice=integrate.layer, Ps)

    # net PAR Cab shaded leaves  (radiance)
    data.canopy[['Rnsha_Cab']] <- data.canopy[['LAI']] * meanleaf.v2(data.canopy, data.rad[['Rnh_Cab']], canopy.choice='layers', Ph)

    # net PAR Cab leaves  (radiance)
    data.canopy[['Rntot_Cab']] <-  data.canopy[['Rnsun_Cab']] + data.canopy[['Rnsha_Cab']]

    ######################################

    # net PAR Car sunlit leaves  (radiance)
    data.canopy[['Rnsun_Car']] <- data.canopy[['LAI']] * meanleaf.v2(data.canopy, data.rad[['Rnu_Car']], canopy.choice=integrate.layer, Ps)

    # net PAR Cab shaded leaves  (radiance)
    data.canopy[['Rnsha_Car']] <- data.canopy[['LAI']] * meanleaf.v2(data.canopy, data.rad[['Rnh_Car']], canopy.choice='layers', Ph)

    # net PAR Car leaves  (radiance)
    data.canopy[['Rntot_Car']] <-  data.canopy[['Rnsun_Car']] + data.canopy[['Rnsha_Car']]



    # net PAR  sunlit leaves  (radiance)
    data.canopy[['Rnsun_PAR']] <- data.canopy[['LAI']] * meanleaf.v2(data.canopy, data.rad[['Rnu_PAR']], canopy.choice=integrate.layer, Ps)

    # net PAR  shaded leaves  (radiance)
    data.canopy[['Rnsha_PAR']] <- data.canopy[['LAI']] * meanleaf.v2(data.canopy, data.rad[['Rnh_PAR']], canopy.choice='layers', Ph)

    # net PAR  leaves  (radiance)
    data.canopy[['Rntot_PAR']] <-  data.canopy[['Rnsun_PAR']] + data.canopy[['Rnsha_PAR']]


    # LST [K] (directional, but assuming black-body surface!)
    sigmaSB <- subset(SCOPEinR::constants,constant == 'sigmaSB')[[2]] #
    data.canopy[['LST']] <- (pi * (data.rad[['Lot']] + data.rad[['Lote']]) / (sigmaSB * data.rad[['canopyemis']]))^0.25
    data.canopy[['emis']]  <- data.rad[['canopyemis']]

    # photosynthesis [mumol m-2 s-1]
    A_mean_bch.leaf <- meanleaf.v2(data.canopy, data.bch[['A']], canopy.choice='layers', Ph)
    A_mean_bcu.leaf <- meanleaf.v2(data.canopy, data.bcu[['A']], canopy.choice=integrate.layer, Ps)
    data.canopy[['A']] <- data.canopy[['LAI']] * (A_mean_bch.leaf + A_mean_bcu.leaf)

    # electron transport rate [mumol m-2 s-1]
    J_mean_bch.leaf <- meanleaf.v2(data.canopy, data.bch[['Ja']], canopy.choice='layers', Ph)
    J_mean_bcu.leaf <- meanleaf.v2(data.canopy, data.bcu[['Ja']], canopy.choice=integrate.layer, Ps)
    data.canopy[['Ja']] <- data.canopy[['LAI']] * (J_mean_bch.leaf + J_mean_bcu.leaf)


    # non-photochemical quenching (energy) [W m-2]

    #  ENPQ energy
    enpq_bch<- data.rad[['Rnh_Cab']] * data.bch[['Phi_N']]
    enpq_bcu<- data.rad[['Rnu_Cab']] * data.bcu[['Phi_N']]

    ENPQ_mean_bch.leaf <- meanleaf.v2(data.canopy, enpq_bch, canopy.choice='layers', Ph)
    ENPQ_mean_bcu.leaf <- meanleaf.v2(data.canopy, enpq_bcu, canopy.choice=integrate.layer, Ps)

    data.canopy[['ENPQ']] <- data.canopy[['LAI']] * (ENPQ_mean_bch.leaf + ENPQ_mean_bcu.leaf)

    # NPQ energy

    pnpq_bch<- data.rad[['Pnh_Cab']] * data.bch[['Phi_N']]
    pnpq_bcu<- data.rad[['Pnu_Cab']] * data.bcu[['Phi_N']]

    PNPQ_mean_bch.leaf <- meanleaf.v2(data.canopy, pnpq_bch, canopy.choice='layers', Ph)
    PNPQ_mean_bcu.leaf <- meanleaf.v2(data.canopy, pnpq_bcu, canopy.choice=integrate.layer, Ps)

    data.canopy[['PNPQ']] <- data.canopy[['LAI']] * (PNPQ_mean_bch.leaf + PNPQ_mean_bcu.leaf)


    # computation of re-absorption corrected fluorescence
    # Yang and Van der Tol (2019); Van der Tol et al. (2019)
    #
    # aPAR_Cab_eta <- data.canopy[['LAI']] *(meanleaf.v2(data.canopy,data.bch$eta * data.rad$Rnh_Cab,'layers',Ph) + meanleaf.v2(data.canopy, data.bcu$eta * data.rad$Rnu_Cab,integrate.layer,Ps))

    aPAR_Cab_eta <- data.canopy[['LAI']] *(meanleaf.v2(data.canopy,data.bch$eta * data.rad$Pnh_Cab,'layers',Ph) + meanleaf.v2(data.canopy,data.bcu$eta * data.rad$Pnu_Cab,integrate.layer,Ps))


    if (options.calc_fluor$Value == 1){
      A_constant <- subset(SCOPEinR::constants,constant == 'A')[[2]] #
      ep <- A_constant * get.ephoton(data.spectral[['wlF']] *1E-9,constants)

      data.rad[['PoutFrc']] <- data.leafbio[['fqe']] * aPAR_Cab_eta

      # 1E-6: umol2mol, 1E3: nm-1 to um-1

      #rad.PoutFrc     = leafbio.fqe*aPAR_Cab_eta;
      #rad.EoutFrc_    = 1E-3*ep.*(rad.PoutFrc*optipar.phi(spectral.IwlF)); %1E-6: umol2mol, 1E3: nm-1 to um-1
      #rad.EoutFrc     = 1E-3*Sint(rad.EoutFrc_,spectral.wlF);


      data.rad[['EoutFrc_']] <- 1E-3 * ep * optipar$phi[data.spectral[['IwlF']]]
      data.rad[['EoutFrc']] <- 1E-3 * Sint(data.rad[['EoutFrc_']],data.spectral[['wlF']])

      sigmaF <- pi * data.rad[['LoF_']] / data.rad[['EoutFrc_']]
      data.rad[['sigmaF']] <-  signal::interp1(data.spectral[['wlF']][1:4:length(data.spectral[['wlF']])],sigmaF[1:4:length(sigmaF)],xi=data.spectral[['wlF']])

      data.canopy[['fqe']] <-  data.rad[['EoutFrc']] / data.canopy[['Pntot_Cab']]

    } else {

      data.canopy[['fqe']] <- NA
    }


    data.rad[['Lotot_']] <- data.rad[['Lo_']]  +  data.rad[['Lot_']]
    data.rad[['Eout_']] <- data.rad[['Eout_']]  +  data.rad[['Eoutte_']]


    if (options.calc_fluor$Value == 1){

      data.rad[['Lototf_']] <- data.rad[['Lotot_']]
      data.rad[['Lototf_']][data.spectral[['IwlF']]] <-   data.rad[['Lototf_']][data.spectral[['IwlF']]] +  data.rad[['LoF_']]
      ## reflectance apparent
      data.rad[['reflapp']] <- data.rad[['refl']]
      data.rad[['reflapp']][data.spectral[['IwlF']]] <- pi * data.rad[['Lototf_']][data.spectral[['IwlF']]] / ( data.rad[['Esun_']][data.spectral[['IwlF']]] + data.rad[['Esky_']][data.spectral[['IwlF']]] )

    }


    if (options.calc_directional$Value == 1) {

      #print('getting brdf angles and effects ...')
      directional <- get.brdf(data.spectral=data.spectral,data.angles=data.angles,data.rad=data.rad,
                              data.directional=data.directional,
                              atmo=atmo,
                              data.soil=data.soil,data.leafopt=data.leafopt,
                              data.leafbio=data.leafbio,data.canopy=data.canopy,data.gap=data.gap,
                              data.meteo=data.meteo,data.thermal=data.thermal,
                              data.bcu=data.bcu,data.bch=data.bch,data.opts=data.opts,
                              get.plots=F)

    } else {

      directional <-NA
    }

    #### get rad Lo

    data.rad[['Lo']] <- 0.001 * Sint(data.rad[['Lo_']][data.spectral[['IwlP']]],data.spectral[['wlP']])




    if( missing(get.outputs) | get.outputs == 'ALL'){

      db.sim[[sim.i]] <- list( data.spectral =data.spectral,
                      data.angles=data.angles,
                      data.rad=data.rad,
                      atmo=atmo,
                      resist_out = resist_out,
                      data.fluxes = data.fluxes,
                      data.soil=data.soil,
                      data.leafopt=data.leafopt,
                      data.leafbio=data.leafbio,
                      data.canopy=data.canopy,
                      data.profiles = data.profiles,
                      data.gap=data.gap,
                      data.meteo=data.meteo,
                      data.thermal=data.thermal,
                      data.bcu=data.bcu,
                      data.bch=data.bch,
                      data.directional=directional,
                      iter.ebal =iter,
                      data.opts=data.opts)


    } else if ( get.outputs == 'Main') {

      db.sim[[sim.i]] <- list( data.spectral =data.spectral,
                      data.rad=data.rad,
                      resist_out = resist_out,
                      data.fluxes = data.fluxes,
                      data.canopy=data.canopy,
                      data.meteo=data.meteo,
                      data.opts=data.opts)




    }
  #progress_bar$tick()

  } #end loop for simulations
  # close bar
  #close(progress_bar)

  return(db.sim)

}






