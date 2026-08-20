#' calculates the energy balance of a vegetated surface
#'
#' @param constants
#' @param data.options
#' @param rad
#' @param gap
#' @param meteo
#' @param soil
#' @param canopy
#' @param leafbio
#'
#' @return
#' @export
#'
#' date: 26 Nov 2007 (CvdT)
#' updates:
#'  - 29 Jan 2008 (JT & CvdT)     converted into a function
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
#' @author 	Christiaan van der Tol, Joris Timmermans  (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
ebal_sunshade <- function(constants, data.options, rad, gap, meteo, soil, canopy, leafbio) {

  # parent: SCOPE.R (script)
  # dependencies:
  # RTMt_sb.R, RTMt_planck.R (optional), RTMf.R (optional)
  # resistances.R
  # heatfluxes.R
  # biochemical.R
  # soil_respiration.R
  #
  # Table of contents of the function
  #
  #   1. Initialisations for the iteration loop
  #           intial values are attributed to variables
  #   2. Energy balance iteration loop
  #           iteration between thermal RTM and surface fluxes
  #   3. Write warnings whenever the energy balance did not close
  #   4. Calculate vertical profiles (optional)
  #   5. Calculate spectrally integrated energy, water and CO2 fluxes
  #
  # The energy balance iteration loop works as follows:
  #
  # RTMo              More or less the classic SAIL model for Radiative
  #                   Transfer of sun and sky light (no emission by the vegetation)
  # While continue	Here an iteration loop starts to close the energy
  #                   balance, i.e. to match the micro-meteorological model
  #                   and the radiative transfer model
  # 	RTMt_sb         A numerical Radiative Transfer Model for thermal
  #                   radiation emitted by the vegetation
  # 	resistances     Calculates aerodynamic and boundary layer resistances
  #                   of vegetation and soil (the micro-meteorological model)
  # 	biochemical     Calculates photosynthesis, fluorescence and stomatal
  #                   resistance of leaves (or biochemical_MD12: alternative)
  # 	heatfluxes      Calculates sensible and latent heat flux of soil and
  #                   vegetation
  #                   Next soil heat flux is calculated, the energy balance
  #                   is evaluated, and soil and leaf temperatures adjusted
  #                   to force energy balance closure
  # end {while continue}
  #
  # meanleaf          Integrates the fluxes over all leaf inclinations
  #                   azimuth angles and layers, integrates over the spectrum

  # The input and output are structures. These structures are further
  # specified in a readme file.
  #
  # Input:
  #
  #   iter        numerical parameters used in the iteration for energy balance closure
  #   options     calculation options
  #   spectral    spectral resolutions and wavelengths
  #   rad         incident radiation
  #   gap         probabilities of direct light penetration and viewing
  #   leafopt     leaf optical properties
  #   angles      viewing and observation angles
  #   soil        soil properties
  #   canopy      canopy properties
  #   leafbio     leaf biochemical parameters
  #
  # Output:
  #
  #   iter        numerical parameters used in the iteration for energy balance closure
  #   fluxes      energy balance, turbulent, and CO2 fluxes
  #   rad         radiation spectra
  #   thermal     temperatures, aerodynamic resistances and friction velocity
  #   bcu, bch    leaf biochemical outputs for sunlit and shaded leaves,
  #               respectively

  constants <- SCOPEinR::constants
  ##load input tables
  leafbio<-getinputLUT(inputLUT=inputLUT, dataset='leafbio')
  canopy<-getinputLUT(inputLUT=inputLUT, dataset='canopy')
  meteo<-getinputLUT(inputLUT=inputLUT, dataset='meteo')
  soil<-getinputLUT(inputLUT=inputLUT, dataset='soil',
                    calc.heat = calc.heat,
                    calc.rss_rbs =  calc.rss_rbs)
  timeseries<-getinputLUT(inputLUT=inputLUT, dataset='timeseries')

  # 1. initialisations and other preparations for the iteration loop
  # parameters for the closure loop
  counter <- 0                  # iteration counter of ebal
  maxit <- 100                  # maximum number of iterations
  maxEBer <- 1                  # maximum energy balance error (any leaf) [Wm-2]
  Wc <- 1                       # update step (1 is nominal, [0,1] possible)
  CONT <- TRUE                  # boolean indicating whether iteration continues

  # functions for saturated vapour pressure
  es_fun <- function(Temp) 6.107 * 10^(7.5 * Temp / (237.3 + Temp))
  s_fun <- function(es, Temp) es * 2.3026 * 7.5 * 237.3 / (237.3 + Temp)^2

  Ta <- meteo[['Ta']]
  ea <- meteo[['ea']]
  Ca <- meteo[['Ca']]
  Ts <- matrix(rep(meteo$Ta, 2), ncol = 1)
  p <- meteo[['p']]


  Rnuc <- rad[['Rnuc']]
  Tch <- (Ta + 0.1) * rep(1, nl)                 # Leaf temperature (shaded leaves)
  Tcu <- (Ta + 0.3) * rep(1, length(Rnuc))       # Leaf temperature (sunlit leaves)
  ech <- matrix(rep(ea, nl), ncol = 1)            # Leaf boundary vapour pressure (shaded/sunlit leaves)
  Cch <- matrix(rep(Ca, nl), ncol = 1)
  ecu <- ea + 0 * Rnuc
  Ccu <- Ca + 0 * Rnuc          # Leaf boundary CO2 (shaded/sunlit leaves)
  meteo[['L']] <- -1                           # Monin-Obukhov length

  # Assign values to constants
  MH2O =  subset(data.constant,constant == 'MH2O')[[2]]
  Mair =  subset(data.constant,constant == 'Mair')[[2]]
  rhoa = subset(data.constant,constant == 'rhoa')[[2]]
  cp = subset(data.constant,constant == 'cp')[[2]]
  sigmaSB = subset(data.constant,constant == 'sigmaSB')[[2]]

  # input preparation
  nl <- canopy[['nlayers']]
  GAM <- soil[['GAM']]
  Ps <- gap[['Ps']]

  kV <- canopy[['kV']]
  xl <- canopy[['xl']]
  LAI <- canopy[['LAI']]
  rss <- soil[['rss']]

  # other preparations
  e_to_q <- MH2O / Mair / p            # Conversion of vapour pressure [Pa] to absolute humidity [kg kg-1]
  Fs <- c(1 - Ps[length(Ps)], Ps[length(Ps)])      # Matrix containing values for 1-Ps and Ps of soil
  Fc <- colSums(1 - Ps[1:(length(Ps)-1)]) / nl      # Matrix containing values for Ps of canopy
  fV <- exp(kV * xl[1:(length(xl)-1)])       # Vertical profile of Vcmax

  if (dim(Rnuc)[2] > 1) {
    fVu <- array(1, dim = c(13, 36, nl))
    for (i in 1:nl) {
      fVu[, , i] <- fV[i]
    }
  } else {
    fVu <- fV
  }

  ## 2. Energy balance iteration loop

  # Energy balance loop (Energy balance and radiative transfer)
  while (CONT) {                          # while energy balance does not close
    # 2.1. Net radiation of the components
    # Thermal radiative transfer model for vegetation emission (with Stefan-Boltzman's equation)
    rad  <- RTMt_sb(constants, rad, soil, leafbio, canopy, gap, Tcu, Tch, Ts[2], Ts[1], 0)

    Rnhc <- rad[['Rnhc']] + rad[['Rnhct']]             # Canopy (shaded) net radiation
    Rnuc <- rad[['Rnuc']] + rad[['Rnuct']]             # Canopy (sunlit) net radiation
    Rnhs <- rad[['Rnhs']] + rad[['Rnhst']]             # Soil   (sun+sh) net radiation
    Rnus <- rad[['Rnus']] + rad[['Rnust']]
    Rns <- c(Rnhs, Rnus)

    # 2.2. Aerodynamic roughness
    # calculate friction velocity [m s-1] and aerodynamic resistances [s m-1]
    resist_out  <- resistances(constants, soil, canopy, meteo)

    meteo[['ustar']] <- resist_out[['ustar']]
    raa     <- resist_out[['raa']]
    rawc    <- resist_out[['rawc']]
    raws    <- resist_out[['raws']]

    # 2.3. Biochemical processes
    meteo_h <- meteo_u <- meteo
    meteo_h[['Temp']] <- mean(Tch)
    meteo_h[['eb']] <- mean(ech)
    meteo_h[['Cs']] <- mean(Cch)
    meteo_h[['Q']]  <- mean(rad[['Pnh_Cab']])
    meteo_u[['Temp']] <- mean(Tcu)
    meteo_u[['eb']] <- mean(ecu)
    meteo_u[['Cs']] <- mean(Ccu)
    meteo_u[['Q']] <- mean(rad[['Pnu_Cab']])
    if (options.Fluorescence_model$Value == 2) {
      b <- biochemical_MD12
    } else {
      b <- biochemical
    }

    #bch = b(leafbio, meteo_h, options, constants, fV);
    #bcu = b(leafbio, meteo_u, options, constants, fVu);

    bch <- b(leafbio, meteo_h, options, constants, 1)
    bcu <- b(leafbio, meteo_u, options, constants, 1)

    # 2.4. Fluxes (latent heat flux (lE), sensible heat flux (H) and soil heat flux G
    # in analogy to Ohm's law, for canopy (c) and soil (s). All in units of [W m-2]

    rss <- soil[['rss']]
    rac <- (LAI+1) * (raa+rawc)
    ras <- (LAI+1) * (raa+raws)

    lEch_Hch_ech_Cch_lambdah_sh <- heatfluxes(rac, bch$rcw, mean(Tch), ea, Ta, e_to_q, Ca, bch$Ci, constants, es_fun, s_fun)


    output_ch <- heatfluxes(rac, bch[['rcw']], Tch, ea, Ta, e_to_q, Ca, bch[['Ci']], constants, es_fun, s_fun)
    lEch <- output_ch[[1]]
    Hch <- output_ch[[2]]
    ech <- output_ch[[3]]
    Cch <- output_ch[[4]]
    lambdah <- output_ch[[5]]
    sh <- output_ch[[6]]


    output_cu <- heatfluxes(rac, bcu[['rcw']], Tcu, ea, Ta, e_to_q, Ca, bcu[['Ci']], constants, es_fun, s_fun)
    lEcu <- output_cu[[1]]
    Hcu <- output_cu[[2]]
    ecu <- output_cu[[3]]
    Ccu <- output_cu[[4]]
    lambdau <- output_cu[[5]]
    su <- output_cu[[6]]

    output_s <- heatfluxes(ras, rss, Ts, ea, Ta, e_to_q, Ca, Ca, constants, es_fun, s_fun)
    lEs <- output_s[[1]]
    Hs <- output_s[[2]]
    lambdas <- output_s[[5]]
    ss <- output_s[[6]]


    # integration over the layers and sunlit and shaded fractions
    Hstot <- Fs * Hs

    if (ncol(Hcu) > 1) {

      Hctot <- LAI * (Fc * Hch + meanleaf(canopy, Hcu, "angles_and_layers", Ps))
    } else {

      Hctot <- LAI * (Fc * Hch + (1 - Fc) * Hcu)
    }

    Htot <- Hstot + Hctot

    # ground heat flux
    G <- 0.35 * Rns

    # 2.5. Monin-Obukhov length L
    meteo[['L']] <- Monin_Obukhov(constants, meteo, Htot)  # [1]

    # 2.6. energy balance errors, continue criterion and iteration counter

    EBerch <- mean(Rnhc) - lEch - Hch
    EBercu <- mean(Rnuc) - lEcu - Hcu
    EBers <- Rns - lEs - Hs - G

    counter <- counter + 1 # Number of iterations
    maxEBercu <- max(abs(EBercu))
    maxEBerch <- max(abs(EBerch))
    maxEBers <- max(abs(EBers))

    CONT <- (maxEBercu > maxEBer | maxEBerch > maxEBer | maxEBers > maxEBer) & counter < maxit + 1 # Continue iteration?
    if (!CONT) {
      if (any(is.nan(c(maxEBercu, maxEBerch, maxEBers)))){
        cat(sprintf("WARNING: NaN in fluxes, counter = %i\n", counter))
      }
      break
    }
    if (counter == 10) { Wc <- 0.8 }
    if (counter == 60) { Wc <- 0.6 }
    if (counter == 100) { Wc <- 0.2 }

    Tch <- Tch + Wc * EBerch / ((rhoa * cp) / rac + rhoa * lambdah * e_to_q * sh / (rac + bch$rcw) + 4 * leafbio[['emis']] * sigmaSB * (meteo_h[['Temp']] + 273.15) ^ 3)
    Tcu <- Tcu + Wc * EBercu / ((rhoa * cp) / rac + rhoa * lambdau * e_to_q * su / (rac + bcu$rcw) + 4 * leafbio[['emis']] * sigmaSB * (meteo_u[['Temp']] + 273.15) ^ 3)
    Ts <- Ts + Wc * EBers / (rhoa * cp / ras + rhoa * lambdas * e_to_q * ss / (ras + rss) + 4 * (1 - soil[['rs_thermal']]) * sigmaSB * (Ts + 273.15) ^ 3)
    Tch[abs(Tch) > 100] <- Ta
    Tcu[abs(Tcu) > 100] <- Ta

  }

  ## 3. Print warnings whenever the energy balance could not be solved

  if(counter >= maxit){

    cat("WARNING: maximum number of iterations exceeded\n")
    cat(sprintf("Energy balance error sunlit vegetation = %4.2f W m-2\n", maxEBercu))
    cat(sprintf("Energy balance error shaded vegetation = %4.2f W m-2\n", maxEBerch))
    cat(sprintf("Energy balance error soil = %4.2f W m-2\n", maxEBers))

  }

  # 4. some more outputs
  iter[['counter']]    <- counter
  # this is not yet written to output but making it avg here breaks RTMt_sb
  # thermal.Tcu     = mean(Tcu(:));
  # thermal.Tch     = mean(Tch);
  thermal[['Tcu']] <- Tcu
  thermal[['Tch']] <- Tch
  thermal[['Tsu']] <- Ts[2]
  thermal[['Tsh']] <- Ts[1]

  fluxes[['Rnctot']] <- LAI * aggregator(mean(Rnhc), mean(Rnuc), Fc, Ps, canopy)  # net radiation leaves
  fluxes[['lEctot']] <- LAI * aggregator(lEch, lEcu, Fc, Ps, canopy)              # latent heat leaves
  fluxes[['Hctot']] <- LAI * aggregator(Hch, Hcu, Fc, Ps, canopy)                # sensible heat leaves
  fluxes[['Actot']] <- LAI * aggregator(bch$A, bcu$A, Fc, Ps, canopy)            # photosynthesis leaves
  fluxes[['Tcave']] <- aggregator(mean(Tch), mean(Tcu), Fc, Ps, canopy)          # mean leaf temperature

  fluxes[['Rnstot']] <- Fs * Rns           # Net radiation soil
  fluxes[['lEstot']] <- Fs * lEs           # Latent heat soil
  fluxes[['Hstot']]  <- Fs * Hs            # Sensible heat soil
  fluxes[['Gtot']]   <- Fs * G             # Soil heat flux
  fluxes[['Tsave']]  <- Fs * Ts            # Soil temperature
  # fluxes[['Resp']]   <- Fs * equations.soil_respiration(Ts)  # Soil respiration = 0

  fluxes[['Rntot']] <- fluxes[['Rnctot']] + fluxes[['Rnstot']]
  fluxes[['lEtot']] <- fluxes[['lEctot']] + fluxes[['lEstot']]
  fluxes[['Htot']] <- fluxes[['Hctot']] + fluxes[['Hstot']]
  fluxes[['rss']] <- rss

  ## faking back layer structure of SCOPE

  bch     <- b(leafbio, meteo_h, options, constants, fV)
  bcu     <- b(leafbio, meteo_u, options, constants, fVu)

  return(list(iter=iter,rad=rad,
              thermal=thermal,soil=soil,
              bcu=bcu,bch=bch,
              fluxes=fluxes,
              resist_out=resist_out,meteo=meteo))

}


#' #' Aggregator function  for fluxes based on LAI parameters
#'
#' @param shaded_flux
#' @param sunlit_flux
#' @param Fc
#' @param Ps
#' @param canopy
#'
#' @return
#' @export
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#'
aggregator <- function(shaded_flux, sunlit_flux, Fc, Ps, canopy) {

  if (ncol(sunlit_flux) > 1) {
    flux_tot <- Fc * shaded_flux + meanleaf(canopy, sunlit_flux, type = "angles_and_layers", Ps = Ps)
  } else {
    flux_tot <- Fc * shaded_flux + (1 - Fc) * sunlit_flux
  }

  return(flux_tot)
}
