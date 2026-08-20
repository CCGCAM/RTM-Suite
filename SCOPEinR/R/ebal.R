#' @title get.ebal
#' \code{get.ebal} Calculates the energy balance of a vegetated surface
#'
#' @param data.rad incident radiation
#' @param data.gap  probabilities of direct light penetration and viewing
#' @param data.meteo meteo characterisitics
#' @param data.soil soil properties
#' @param data.canopy canopy properties
#' @param data.leafbio   leaf biochemical parameters
#' @param data.leafopt   leaf biochemical parameters
#' @param data.spectral spectral information for the model
#' @param k.maxit maximum number for iteration on the balance model
#' @param integrate.layer how estimate the sum of layer (angles, layers or angles_layers)
#' @param data.opts options for running the model
#' @param get.plots  is true plot the intermediate plots
#'
#' @description
#'
#' date   26 Nov 2007 (CvdT)
#' updates:
#' - 29 Jan 2008 (JT & CvdT)     converted into a function
#'  - 11 Feb 2008 (JT & CvdT)     improved soil heat flux and temperature calculation
#'  - 14 Feb 2008 (JT)            changed h in to hc (as h=Avogadro`s constant)
#'  - 31 Jul 2008 (CvdT)          Included Pntot in output
#'  - 19 Sep 2008 (CvdT)          Converted F0 and F1 from units per aPAR into units per iPAR
#'  - 07 Nov 2008 (CvdT)          Changed layout
#'  - 18 Sep 2012 (CvdT)          Changed Oc, Cc, ec
#'  - Feb 2012 (WV)            introduced structures for variables
#'  - Sep 2013 (JV, CvT)       introduced additional biochemical model
#'  - 10 Dec 2019 (CvdT)          made a light version (layer averaged fluxes)
#'
#' @return  a list with:
#'  - iter        numerical parameters used in the iteration for energy balance closure.
#'  - fluxes      energy balance, turbulent, and CO2 fluxes.
#'  - rad         radiation spectra.
#'  - thermal     temperatures, aerodynamic resistances and friction velocity.
#'  - bcu, bch    leaf biochemical outputs for sunlit and shaded leaves,
#               respectively.
#' @export

#' @author 	Christiaan van der Tol,  Joris Timmermans (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' out <- get.ebal(data.rad, data.gap, data.meteo, data.soil, data.canopy,
#'                  data.leafbio, data.leafopt, data.spectral, data.opts,
#'                  integrate.layer = "angles_layers", k.maxit = 100,
#'                  get.plots = FALSE)
#' }

get.ebal <- function(data.rad, data.gap, data.meteo, data.soil, data.canopy, data.leafbio,data.leafopt,
                     data.spectral, data.opts,integrate.layer, k.maxit,get.plots) { #k, xyt,

  # This function has the following dependencies:

  #  get.RTMt.sb.R, get.RTMt.planck.R (optional), get.RTMf.R (optional)
  #  get.resistances.R
  #  get.heatfluxes.R
  #  get.biochemical.R
  #  get.soil.respiration.R


  # Table of contents of the function
  #
  #   1. Initialization for the iteration loop
  #      1.1. initial values are attributed to variables
  #   2. Energy balance iteration loop
  #     2-1. iteration between thermal RTM and surface fluxes
  #   3. Write warnings whenever the energy balance did not close
  #   4. Calculate vertical profiles (optional)
  #   5. Calculate spectral integrated energy, water and CO2 fluxes



  # The energy balance iteration loop works as follows:
  #
  # RTMo: More or less the classic SAIL model for Radiative Transfer of sun and sky light (no emission by the vegetation)
  #    While continue	Here an iteration loop starts to close the energy balance, i.e. to match the micro-meteorological model and the radiative transfer model
  # 	    - get.RTMt.sb: A numerical Radiative Transfer Model for thermal radiation emitted by the vegetation
  # 	    - get.resistances: Calculates aerodynamic and boundary layer resistances of vegetation and soil (the micro-meteorological model)
  # 	    - get.biochemical: Calculates photosynthesis, fluorescence and stomatal resistance of leaves (or biochemical_MD12: alternative)
  # 	     - get.heatfluxes: Calculates sensible and latent heat flux of soil and vegetation Next soil heat flux is calculated, the energy balance is evaluated, and soil and leaf temperatures adjusted to force energy balance closure
  #    end {while continue}
  #
  #   - meanleaf: Integrates the fluxes over all leaf inclinations azimuth angles and layers, integrates over the spectrum
  #


  if (missing(data.opts)){
    stop('please use options for Lite, Calc_vert_profiles ')

  } else{

    options.Fluorescence_model  = data.opts[6,]   #0: empirical, with sustained NPQ (fit to Flexas' data); 1: empirical, with sigmoid for Kn; 2: Magnani 2012 model
    options.simulation          = data.opts[10,]
    options.MoninObukhov        = data.opts[15,]
    options.soil_heat_method     = data.opts[13,]    # 0 - GAM=Soil_Inertia0(lambdas), 1 - GAM=Soil_Inertia1(SMC), 2 - G=0.35*Rn (always in no TS)

  }


  if (missing(get.plots)){
    get.plots = FALSE
  }

  # 1. initialisations and other preparations for the iteration loop
  # parameters for the closure loop
  counter <- 0                  # iteration counter of ebal
  maxit <- k.maxit               # maximum number of iterations

  #################################################################################
  maxEBer <- 1                  # maximum energy balance error (any leaf) [Wm-2]
  #################################################################################
  Wc <- 1                       # update step (1 is nominal, [0,1] possible)
  CONT <- TRUE                  # boolean indicating whether iteration continues


  # Assign values to constants
  MH2O =  subset(SCOPEinR::constants,constant == 'MH2O')[[2]]
  Mair =  subset(SCOPEinR::constants,constant == 'Mair')[[2]]
  rhoa = subset(SCOPEinR::constants,constant == 'rhoa')[[2]]
  cp = subset(SCOPEinR::constants,constant == 'cp')[[2]]
  sigmaSB = subset(SCOPEinR::constants,constant == 'sigmaSB')[[2]]

  # input preparation
  nl <- data.canopy[['nlayers']]
  GAM <- data.soil[['GAM']]
  Ps <- data.gap[['Ps']]

  kV <- data.canopy[['kV']]
  xl <- data.canopy[['xl']]
  LAI <- data.canopy[['LAI']]
  rss <- data.soil[['rss']]


  # functions for saturated vapour pressure
  es_fun <- function(Temp) 6.107 * 10^(7.5 * Temp / (237.3 + Temp))
  s_fun <- function(es, Temp) es * 2.3026 * 7.5 * 237.3 / (237.3 + Temp)^2


  ## Available options for SoilHeat Method:
  # - 0: GAM=Soil_Inertia0(lambdas),
  # - 1: GAM=Soil_Inertia1(SMC),
  # - 2: G=0.35*Rn (always in no TS)

  SoilHeatMethod <- options.soil_heat_method$Value


  if (!(options.simulation$Value == 1)) {
    SoilHeatMethod <- 2
  }

  if (SoilHeatMethod < 2) {
    if (k > 1) {
      Deltat <- (as.numeric(data.timeseries$t[k]) - as.numeric(data.timeseries$t[k - 1])) * 86400     # Duration of the time interval (s)
    } else {
      Deltat <- 1 / 48 * 86400
    }
    x <- matrix(1:24, ncol = 2) * Deltat
    Tsold <- data.soil$Tsold
  }

  # meteo
  Ta <- data.meteo[['Ta']]
  ea <- data.meteo[['ea']]
  Ca <- data.meteo[['Ca']]
  p <- data.meteo[['p']]
  Rnuc <- data.rad[['Rnuc']]
  meteo_u.Q <- data.rad[['Pnu_Cab']];

  ech <- rep(ea, nl)          # Leaf boundary vapour pressure (shaded/sunlit leaves)
  Cch <- rep(Ca, nl)
  ecu <- ea + 0 * Rnuc
  Ccu <- Ca + 0 * Rnuc          # Leaf boundary CO2 (shaded/sunlit leaves)

  # other preparations
  e_to_q <- MH2O / Mair / p             # Conversion of vapour pressure [Pa] to absolute humidity [kg kg-1]
  Fc <- Ps[1:(nl)]
  Fs <- c(1 - Ps[nl], Ps[nl])      # Matrix containing values for 1-Ps and Ps of soil
  #Fc <- (1 - Ps[1:(nl-1)])/nl      # Matrix containing values for Ps of canopy
  fV <- exp(data.leafbio[['kV']] * xl[1:(nl)])      # Vertical profile of Vcmax (nl=59)

  # initial values for the loop
  Ts <- c(Ta+3, Ta+3)         # soil temperature (+3 for a head start of the iteration)
  Tch <- rep(Ta+0.1, nl)      # leaf temperature (shaded leaves)
  Tcu <- rep(Ta+0.3, length(Rnuc)) # leaf temperature (sunlit leaves)
  data.meteo[['L']] <- -1E6             # Monin-Obukhov length

  meteo_h <- data.meteo

  meteo_u <- data.meteo

  # this for is the exponential decline of Vcmax25. If 'lite' the dimensions
  # are [nl], otherwise [13,36,nl]:

  if (is.vector(Rnuc) == FALSE) {

    fVu <- array(1, dim = c(13, 36, nl))
    for (i in 1:nl) {
      fVu[,,i] <- fV[i]
      print('check this code part ...')
    }
  } else {
    fVu <- fV
  }

  ## 2.1 Energy balance iteration loop
  # Energy balance loop (Energy balance and radiative transfer)
  while (CONT) {
    # 2.1. Net radiation of the components
    # Thermal radiative transfer model for vegetation emission (with Stefan-Boltzman's equation)

    ## Here is not working the get.plot because is so slow, get.plots is force to FALSE

    outputs_RTMt.sb <- get.RTMt.sb( data.rad,data.soil,data.leafbio,data.canopy,data.leafopt,data.gap,
                                      Tcu, Tch, Tsu = Ts[2], Tsh=Ts[1], obsdir=0,data.spectral,
                                      data.opts=data.opts,get.plots=F)


    data.rad <- outputs_RTMt.sb[['data.rad']]

    ### data.rad[['Rnhc']] is so high
    Rnhc <- data.rad[['Rnhc']] + data.rad[['Rnhct']]     # Canopy (shaded) net radiation
    Rnuc <- data.rad[['Rnuc']] + data.rad[['Rnuct']]     # Canopy (sunlit) net radiation
    Rnhs <- data.rad[['Rnhs']] + data.rad[['Rnhst']]     # Soil (sun+sh) net radiation
    Rnus <- data.rad[['Rnus']] + data.rad[['Rnust']]
    Rns <- c(Rnhs, Rnus)

    # 2.3. Biochemical processes

    #shaded leaves
    meteo_h[['Temp']] <- Tch
    meteo_h[['eb']] <- ech
    meteo_h[['Cs']] <- Cch
    meteo_h[['Q']] <- data.rad[['Pnh_Cab']] #[60x1 double] net PAR absorbed by Cab (moles m-2 s-1) of shaded leaves

    ## sunlit leaves
    meteo_u[['Temp']]<- Tcu
    meteo_u[['eb']] <- ecu
    meteo_u[['Cs']] <- Ccu
    meteo_u[['Q']] <- data.rad[['Pnu_Cab']]  # [13x36x60 double] net PAR absorbed by Cab (moles m-2 s-1) of sunlit leaves

    ###################
    ##############

    if (options.Fluorescence_model$Value == 1) {

        ## Here remove the options for plotting to slow
        data.bch <- get.biochemical.MD12(data.leafbio=data.leafbio, data.meteo=meteo_h, fV,get.plots = F)
        data.bcu <- get.biochemical.MD12(data.leafbio=data.leafbio, data.meteo=meteo_u, fVu,get.plots =F)


    } else {
        ## Here remove the options for plotting to slow
        data.bch <- get.biochemical(data.leafbio=data.leafbio, data.meteo=meteo_h,data.opts= data.opts, fV,get.plots = F)
        data.bcu <- get.biochemical(data.leafbio=data.leafbio, data.meteo=meteo_u, data.opts= data.opts, fVu,get.plots = F)

    }

    # Aerodynamic roughness
    # calculate friction velocity [m s-1] and aerodynamic resistances [s m-1]
    resist_out <- get.resistances(data.soil, data.canopy, data.meteo)

    data.meteo[['ustar']] <- resist_out$ustar
    raa <- resist_out[['raa']]
    rawc <- resist_out[['rawc']]
    raws <- resist_out[['raws']]
    rac <- (LAI + 1) * (raa + rawc)
    ras <- (LAI + 1) * (raa + raws)

    # Fluxes (latent heat flux (lE), sensible heat flux (H) and soil heat flux G
    # in analogy to Ohm's law, for canopy (c) and soil (s). All in units of [W m-2]
    output_ch <- get.heatfluxes(ra=rac, rs=data.bch[['rcw']], Tc=Tch, ea=ea, Ta=Ta, e_to_q, Ca, Ci=data.bch[['Ci']])

    lEch <- output_ch[['lE']]
    Hch <- output_ch[['H']]
    ech <- output_ch[['ec']]
    Cch <- output_ch[['Cc']]
    lambdah <- output_ch[['lambda']]
    sh <- output_ch[['s']]

    output_cu <- get.heatfluxes(ra=rac, rs=data.bcu[['rcw']], Tc=Tcu, ea=ea, Ta=Ta, e_to_q, Ca, Ci=data.bcu[['Ci']])

    lEcu <- output_cu[['lE']]
    Hcu <- output_cu[['H']]
    ecu <- output_cu[['ec']]
    Ccu <- output_cu[['Cc']]
    lambdau <- output_cu[['lambda']]
    su <- output_cu[['s']]

    output_s <- get.heatfluxes(ra=ras, rs=rss, Tc=Ts, ea=ea, Ta=Ta, e_to_q, Ca, Ci=Ca)

    lEs <- output_s[['lE']]
    Hs <- output_s[['H']]
    lambdas <- output_s[['lambda']]
    ss <- output_s[['s']]



    # integration over the layers and sunlit and shaded fractions
    Hstot <- sum(Fs * Hs)

    ### Hctot is a vector of length  = nl when integrate.layer == 'angles'
    #integrate.layer = 'angles' //     integrate.layer = 'layers'
    Hctot <-get.aggregator.ebal(LAI, sunlit_flux=Hcu, shaded_flux=Hch, Fs=Ps[1:(length(Ps) - 1)],
                                                 data.canopy, canopy.choice=integrate.layer)

    #print(Hctot)
    Htot <- Hstot + Hctot

    if (options.MoninObukhov$Value == 1) {
      data.meteo[['L']] <- get.Monin.Obukhov(data.meteo, Htot)
    }

    # ground heat flux
    if (SoilHeatMethod == 2) {
      G <- 0.35 * Rns
      dG <- 4 * (1 - data.soil$rs_thermal) * sigmaSB * (Ts + 273.15)^3 * 0.35
    } else {
      G <- GAM / sqrt(pi) * 2 * sum((c(Ts, Tsold[1:(length(Tsold) - 1),]) - Tsold) / Deltat * (sqrt(x) - sqrt(x - Deltat)))
      G <- G[1]
      dG <- GAM / sqrt(pi) * 2 * ((sqrt(x[1]) - sqrt(x[1] - Deltat))) / Deltat * matrix(rep(1, times = 2))
    }


    # energy balance errors, continue criterion and iteration counter
    EBerch  <- Rnhc - lEch - Hch
    EBercu  <- Rnuc - lEcu - Hcu
    EBers   <- Rns - lEs - Hs - G

    counter     <- counter + 1  # Number of iterations
    ## Maximum energy balance error sunlit vegetation
    maxEBercu   <- max(abs(EBercu))
    ## Maximum energy balance error shaded vegetation
    maxEBerch   <- max(abs(EBerch))
    ## Energy balance error soil
    maxEBers    <- max(abs(EBers))

    CONT  <- (maxEBercu > maxEBer | maxEBerch > maxEBer | maxEBers > maxEBer) & (counter < maxit+1) # Continue iteration?

    #if ((counter < maxit + 1) & (maxEBercu > maxEBer | maxEBerch > maxEBer | maxEBers > maxEBer)) {
     # CONT <- FALSE
    #}



    if (!CONT) {
      if (any(is.nan(c(maxEBercu, maxEBerch, maxEBers)))) {
        cat(sprintf('WARNING: NaN in fluxes, counter = %i\n', counter))
      }
      break
    }

    if (counter == 10) Wc <- 0.8
    if (counter == 20) Wc <- 0.6

    # if (counter > 99) { plot(EBercu), hold on } # uncomment if needed

    # New estimates of soil (s) and leaf (c) temperatures, shaded (h) and sunlit (1)
    Tch <- Tch + Wc * EBerch / ((rhoa * cp) / rac +
                                  rhoa * lambdah * e_to_q * sh / (rac + data.bch[['rcw']]) +
                                  4 * data.leafbio[['emis']] * sigmaSB * (Tch + 273.15)^3)

    Tcu <- Tcu + Wc * EBercu / ((rhoa * cp) / rac +
                                  rhoa * lambdau * e_to_q * su / (rac + data.bcu[['rcw']]) +
                                  4 * data.leafbio[['emis']] * sigmaSB * (Tcu + 273.15)^3)

    # Was `Ts[1] + ...`: only the shaded-soil temperature ever fed the update
    # for BOTH components (shaded and sunlit), while the denominator already
    # correctly used the full 2-element Ts vector -- inconsistent with how
    # Tch/Tcu are updated just above (full vector on both sides). That mismatch
    # meant Ts[2] (sunlit soil) never actually updated from its own previous
    # value, which is a plausible root cause of the soil energy balance never
    # converging while the (correctly-vectorized) canopy terms did.
    Ts  <- Ts + Wc * EBers / (rhoa * cp / ras +
                                rhoa * lambdas * e_to_q * ss / (ras + rss) +
                                4 * (1 - data.soil[['rs_thermal']]) * sigmaSB * (Ts + 273.15)^3 + dG)

    Tch[abs(Tch) > 100] <- Ta
    Tcu[abs(Tcu) > 100] <- Ta

  }

  ## 2.2 emmissivity calculation
  ## Here remove the options for plotting to slow
  rad_ <- get.RTMt.sb( data.rad,data.soil,data.leafbio,data.canopy,data.leafopt,data.gap,
                         Tcu, Tch, Tsu = Ts[2], Tsh=Ts[1], obsdir=0,data.spectral,
                         data.opts=data.opts,get.plots=F)




  blackleaf<-list()
  blackleaf[['tau_thermal']] <- 0
  blackleaf[['rho_thermal']] <- 0
  blacksoil<-list()
  blacksoil[['rs_thermal']] <- 0

  ## Here remove the options for plotting to slow
  rad0_ <-  get.RTMt.sb( data.rad,data.soil=blacksoil,data.leafbio=blackleaf,data.canopy,data.leafopt,data.gap,
                           Tcu, Tch, Tsu = Ts[2], Tsh=Ts[1], obsdir=0,data.spectral,
                           data.opts=data.opts,get.plots=F)


  data.rad[['canopyemis']] <- rad_$data.rad[['Eoutte']] / rad0_$data.rad[['Eoutte']]



  ## 3. Print warnings whenever the energy balance could not be solved

  ### here not used when use the get.SCOPE function
  if (counter >= maxit) {
    cat("WARNING: maximum number of iterations exceeded\n")
    cat(sprintf("Maximum energy balance error sunlit vegetation = %4.2f W m-2\n", maxEBercu))
    cat(sprintf("Maximum energy balance error shaded vegetation = %4.2f W m-2\n", maxEBerch))
    cat(sprintf("Energy balance error soil = %4.2f W m-2\n", maxEBers))
    cat(sprintf("Mean error sunlit vegetation = %4.2f W m-2\n", mean(as.vector(EBercu))))
  }

  # 4. some more outputs

  iter <- list()
  iter[['counter']] <- counter
  iter[['maxit']] <- maxit
  iter[['maxEBercu']] <- maxEBercu
  iter[['maxEBerch']] <- maxEBerch
  iter[['maxEBers']] <- maxEBers


  data.thermal<-list()

  data.thermal[['Tcu']] <- Tcu
  data.thermal[['Tch']] <- Tch
  data.thermal[['Tsu']] <- Ts[2]
  data.thermal[['Tsh']] <- Ts[1]

  data.fluxes <- list()



  # net radiation leaves
  data.fluxes[['Rnctot']]<- sum(get.aggregator.ebal(LAI, sunlit_flux=Rnuc, shaded_flux=Rnhc, Fs=Fc, data.canopy, canopy.choice=integrate.layer))
  # latent heat leaves
  data.fluxes[['lEctot']] <- sum(get.aggregator.ebal(LAI, sunlit_flux=lEcu, shaded_flux=lEch, Fs=Fc, data.canopy, canopy.choice=integrate.layer))
  # sensible heat leaves
  data.fluxes[['Hctot']] <- sum(get.aggregator.ebal(LAI, sunlit_flux=Hcu, shaded_flux=Hch, Fs=Fc, data.canopy, canopy.choice=integrate.layer))
  # photosynthesis leaves
  data.fluxes[['Actot']] <- sum(get.aggregator.ebal(LAI, sunlit_flux=data.bcu[['A']], shaded_flux=data.bch[['A']], Fs=Fc, data.canopy, canopy.choice=integrate.layer))
  # mean leaf temperature
  data.fluxes[['Tcave']] <- sum(get.aggregator.ebal(1, sunlit_flux=Tcu, shaded_flux=Tch, Fs=Fc, data.canopy, canopy.choice=integrate.layer))
  # Net radiation soil
  data.fluxes[['Rnstot']] <- sum(Fs * Rns)
  # Latent heat soil
  data.fluxes[['lEstot']] <- sum(Fs * lEs)
  # Sensible heat soil
  data.fluxes[['Hstot']] <- sum(Fs * Hs)
  # Soil heat flux
  data.fluxes[['Gtot']] <- sum(Fs * G)
  # Soil temperature
  data.fluxes[['Tsave']] <- sum(Fs * Ts)

  # fluxes[['Resp']] <- Fs * equations.soil_respiration(Ts) #  Soil respiration = 0
  data.fluxes[['Rntot']] <- data.fluxes[['Rnctot']] + data.fluxes[['Rnstot']]
  data.fluxes[['lEtot']] <- data.fluxes[['lEctot']] + data.fluxes[['lEstot']]
  data.fluxes[['Htot']] <- data.fluxes[['Hctot']] + data.fluxes[['Hstot']]

  resist_out[['rss']] <- rss  # this is simply a copy of the input rss

  # update soil temperatures history
  if (SoilHeatMethod < 2) {
    Tsold[-1,] <- soil$Tsold[-nrow(soil$Tsold),]
    Tsold[1,] <- Ts
    if (is.na(Ts)) {
      Tsold[1,] <- Tsold[2,]
    }
    soil[['Tsold']] <- Tsold
  }

  return(list(iter=iter,
              data.rad = data.rad,

              data.thermal=data.thermal,data.soil=data.soil,

               data.bcu=data.bcu,
               data.bch=data.bch,

               data.fluxes=data.fluxes,

               resist_out=resist_out,

              data.meteo=data.meteo))


} # end function
