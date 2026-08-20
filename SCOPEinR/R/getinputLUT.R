#' Get main input for SCOPE
#'
#' \code{getinputLUT} extracts one group of SCOPE model input parameters
#' (a single row, or a set of rows, of a LUT/input table) into the named
#' list structure expected by the rest of SCOPEinR (leaf biology, meteo,
#' canopy structure, geometry, soil, coordinates, time series, or multilayer
#' leaf biochemistry profile). For \code{dataset == 'soil'}, it also derives
#' the soil thermal inertia (\code{GAM}, via \code{\link{Soil_Inertia0}} or
#' \code{\link{Soil_Inertia1}}) and, optionally, the soil/canopy resistances
#' \code{rss}/\code{rbs} (via \code{\link{calc_rssrbs}}).
#'
#' @param inputLUT data.frame or matrix. LUT/input table with one column per named parameter (e.g. \code{N}, \code{Cab}, \code{LAI}, \code{tts}, \code{Ta}, ...) and one row per simulation.
#' @param dataset character. Which group of parameters to extract; one of \code{'leafbio'}, \code{'meteo'}, \code{'canopy'}, \code{'angles'}, \code{'soil'}, \code{'Coordinates'}, \code{'mly'} or \code{'timeseries'}.
#' @param calc.heat numeric. Soil heat flux method used only when \code{dataset == 'soil'}: \code{0} estimates \code{GAM} with \code{\link{Soil_Inertia0}} (from \code{cs}, \code{rhos}, \code{lambdas}), \code{1} estimates \code{GAM} with \code{\link{Soil_Inertia1}} (from \code{SMC}), any other value (default \code{2}) leaves \code{GAM} as \code{NA} (soil heat flux computed elsewhere as \code{0.35*Rn}).
#' @param calc.rss_rbs numeric. Flag used only when \code{dataset == 'soil'}: \code{0} (default) uses \code{rss}/\code{rbs} directly from \code{inputLUT}; any other value recomputes them from \code{SMC}/\code{LAI}/\code{rbs} via \code{\link{calc_rssrbs}}.
#'
#' @return A named list with the extracted parameters for the requested \code{dataset} (e.g. \code{leafbio}, \code{meteo}, \code{canopy}, \code{angles}, \code{soil}, \code{coord}, \code{timeseries}, or \code{mly}).
#' @export
#'
#'
#' @author Carlos Camino
#'
#' @examples
#' \dontrun{
#' leafbio <- getinputLUT(inputLUT, dataset = 'leafbio')
#' soil <- getinputLUT(inputLUT, dataset = 'soil', calc.heat = 2, calc.rss_rbs = 0)
#' }
getinputLUT<-function(inputLUT, dataset,
                      calc.heat, calc.rss_rbs){



  if (dataset == 'leafbio'){
    leafbio<-list()

    ## 1.1. PROSPECT parameters  ----
    leafbio$N=inputLUT[,'N']; leafbio$Cab=inputLUT[,'Cab'];

    if (is.null(inputLUT[,'Car'])){
      leafbio$Car = 0.25 * leafbio$Cab
    } else {
      leafbio$Car=inputLUT[,'Car'];
    }

    leafbio$Anth=inputLUT[,'Anth']; leafbio$Cbrown=inputLUT[,'Cbrown'];
    leafbio$Cs=inputLUT[,'Cs']
    leafbio$EWT=inputLUT[,'EWT']; leafbio$LMA=inputLUT[,'LMA']; leafbio$alpha=inputLUT[,'alpha']
    leafbio$Prot=inputLUT[,'Prot'];leafbio$CBC=inputLUT[,'CBC'];
    leafbio$Cx=inputLUT[,'Cx']; #
    leafbio$rho_thermal=inputLUT[,'rho_thermal']
    leafbio$tau_thermal=inputLUT[,'tau_thermal']
    ## 1.3. Leaf_Biochemical parameters ----
    leafbio$Vcmax25=inputLUT[,'Vcmax25']
    leafbio$BallBerrySlope=inputLUT[,'BallBerrySlope']
    leafbio$BallBerry0=inputLUT[,'BallBerry0']
    leafbio$Type=inputLUT[,'Type']
    leafbio$kV=inputLUT[,'kV']

    leafbio$Rdparam=inputLUT[,'Rdparam']
    leafbio$Tparam <- c('0.2','0.3','281','308','328')



    leafbio$Kn0 <- inputLUT[,'Kn0']
    leafbio$Knalpha <- inputLUT[,'Knalpha']
    leafbio$Knbeta <- inputLUT[,'Knbeta']

    ## 1.4.Fluorescence ----
    leafbio$fqe<- inputLUT[,'fqe']
    ## 1.3. change Leaf_Biochemical (mag0000i model) ----
    leafbio$Tyear <- inputLUT[,'Tyear']
    leafbio$beta <- inputLUT[,'beta']
    leafbio$kNPQs <- inputLUT[,'kNPQs']
    leafbio$qLs <- inputLUT[,'qLs']
    leafbio$stressfactor <- inputLUT[,'stressfactor']

    dataset<-leafbio

  } else if (dataset == 'meteo') {
    ## 1.7. Meteo (values in data files, in the time series option, can overrule these values) ----
    meteo <- list()

    meteo$z<- inputLUT[,'z']
    meteo$Rin<- inputLUT[,'Rin']
    meteo$Ta<- inputLUT[,'Ta']
    meteo$Rli<- inputLUT[,'Rli']
    meteo$p<- inputLUT[,'p']
    meteo$ea<- inputLUT[,'ea']
    meteo$u<- inputLUT[,'u']
    meteo$Ca<- inputLUT[,'Ca']
    meteo$Oa<- inputLUT[,'Oa']

    dataset<-meteo

  } else if (dataset == 'canopy') {

    ## parameters for fourSAIL (canopy)
    canopy<-list()

    canopy$LIDFa=inputLUT[,'LIDFa']; canopy$LIDFb=inputLUT[,'LIDFb'];
    canopy$TypeLidf=inputLUT[,'TypeLidf']; canopy$LAI=inputLUT[,'LAI']
    canopy$hotspot=inputLUT[,'hspot'];
    canopy$hc=inputLUT[,'hc']; canopy$leafwidth = inputLUT[,'leafwidth']
    canopy[['hot']]  = canopy$leafwidth / canopy$hc
    ## 1.8. Aerodynamic (values in data files, in the time series option, can overrule these values) ----
    canopy$zo<- inputLUT[,'zo']
    canopy$d<- inputLUT[,'d']
    canopy$Cd<- inputLUT[,'Cd']
    canopy$rb<- inputLUT[,'rb']
    canopy$CR<- inputLUT[,'CR']
    canopy$CD1<- inputLUT[,'CD1']
    canopy$Psicor <- inputLUT[,'Psicor']
    canopy$rbs<- inputLUT[,'rbs']
    canopy$rwc<- inputLUT[,'rwc']

    dataset<-canopy

  } else if (dataset == 'angles') {
    ## 1.10. Angles	 ----
    angles <- list()

    angles$tts=inputLUT[,'tts']; angles$tto=inputLUT[,'tto']; angles$psi=inputLUT[,'psi']

    dataset<-angles

  } else if (dataset == 'soil') {

    ## 1.5. Soil ----
    soil<-list()

    soil$spectrum<- inputLUT[,'spectrum'] ##class of soil in the file
    soil$rss<- inputLUT[,'rss']
    soil$rs_thermal<- inputLUT[,'rs_thermal']
    soil$tau_thermal<- inputLUT[,'tau_thermal']

    soil$cs <- inputLUT[,'cs']
    soil$rhos <- inputLUT[,'rhos']
    soil$lambdas <- inputLUT[,'lambdas']
    soil$SMC <- inputLUT[,'SMC']

    if ( (mean(soil$SMC, na.rm = T) > 1) == T){
      soil$SMC <-  soil$SMC / 100  # SMC from [0 100] to [0 1]
    }


    ## derived input
    # Value 0 will estimate the GAM parameters with Soil_Inertia0(lambdas);
    # Value 1 will estimate the GAM with Soil_Inertia1(SMC);
    # Value 2 will estimate  G by 0.35*Rn (where always in no TS)

    if (missing(calc.heat)){

      calc.heat = 2
    }

    if (is.null(calc.heat)){

      calc.heat = 2
    }


    if (calc.heat == 1) {

      soil$GAM <-  Soil_Inertia1(soil$SMC)

    } else if (calc.heat == 0) {

      soil$GAM  <- Soil_Inertia0(soil[['cs']],soil[['rhos']],soil[['lambdas']])

    } else {

      soil$GAM <- NA
    }

    soil$BSMBrightness<- inputLUT[,'BSMBrightness']

    soil$CSSOIL<- inputLUT[,'CSSOIL']
    soil$rbs  <- inputLUT[,'rbs']
    soil$SMC  <- inputLUT[,'SMC']
    soil$BSMlat <- inputLUT[,'BSMlat']
    soil$BSMlon <- inputLUT[,'BSMlon']

    if (missing(calc.rss_rbs)){
      calc.rss_rbs = 0
    }

    if (is.null(calc.rss_rbs)){

      calc.rss_rbs = 0
    }


    if (calc.rss_rbs == 0) {

      soil$rss  <-  inputLUT[,'rss']
      soil$rbs  <-  inputLUT[,'rbs']

    } else {

      #run this this is for options.calc_rss_rbs ==1
      canopy.LAI <- inputLUT[,'LAI']

      outputs_ <-calc_rssrbs(SMC=soil[['SMC']], LAI= canopy.LAI, rbs=soil[['rbs']])
      soil$rss  <- outputs_['rss']
      soil$rbs  <- outputs_['rbs']
    }

    dataset<-soil



  } else if (dataset == 'Coordinates') {

    ##Coordinates in Decimal degrees
    coord<-list()

    coord$BSMlat<- inputLUT[,'BSMlat']
    coord$BSMlon<- inputLUT[,'BSMlon']

    dataset<-coord

  } else if (dataset == 'timeseries') {

    ## 1.9. timeseries (this option is only for time series)	 ----
    timeseries<-list()

    timeseries$startDate<- inputLUT[,'startDate']
    timeseries$endDate<- inputLUT[,'endDate']
    timeseries$LAT<- inputLUT[,'LAT']
    timeseries$LON<- inputLUT[,'LON']
    timeseries$timezn<- inputLUT[,'timezn']

    dataset <-timeseries



  } else if (dataset == 'mly') {

    mly <- list()
    mly[['nly']] <- 1
    mly[['pLAI']]<- inputLUT[,'LAI']
    mly[['pCab']] <- inputLUT[,'Cab']
    mly[['pCar']] <- inputLUT[,'Car']
    mly[['pLMA']] <- inputLUT[,'LMA']
    mly[['pEWT']] <- inputLUT[,'EWT']
    mly[['pCs']] <- inputLUT[,'Cs']
    mly[['pN']]  <- inputLUT[,'N']
    totLAI <- sum(mly[['pLAI']])

    mly[['totLAI']] <- totLAI

    dataset <- mly
    return(dataset)

  }


}
