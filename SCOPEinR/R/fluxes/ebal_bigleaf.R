#' calculates the energy balance of a vegetated surface
#'
#' @param constants
#' @param options
#' @param rad
#' @param gap
#'
#' @return
#' @export
#'
#' @author 	 Christiaan van der Tol,   Joris Timmermans (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' last updates
#' date          26 Nov 2007 (CvdT)
#' updates       29 Jan 2008 (JT & CvdT)     converted into a function
#'               11 Feb 2008 (JT & CvdT)     improved soil heat flux and temperature calculation
#'               14 Feb 2008 (JT)            changed h in to hc (as h=Avogadro`s constant)
#'               31 Jul 2008 (CvdT)          Included Pntot in output
#'               19 Sep 2008 (CvdT)          Converted F0 and F1 from units per aPAR into units per iPAR
#'               07 Nov 2008 (CvdT)          Changed layout
#'               18 Sep 2012 (CvdT)          Changed Oc, Cc, ec
#'                  Feb 2012 (WV)            introduced structures for variables
#'                  Sep 2013 (JV, CvT)       introduced additional biochemical model
#'              10 Dec 2019 (CvdT)          made a light version (layer averaged fluxes)
# parent: SCOPE.m (script)
#'
#' @references
#' @examples
#'
ebal_bigleaf<-function(constants,options,rad,gap){

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

  # dependencies functions: ebal.R
  # uses in these functions:
  #  RTMt_sb.R, RTMt_planck.R (optional), RTMf.R (optional)
  #  resistances.R
  #  heatfluxes.T
  #  biochemical.R
  #  soil_respiration.R

  # Output:
  #
  #   iter        numerical parameters used in the iteration for energy balance closure
  #   fluxes      energy balance, turbulent, and CO2 fluxes
  #   rad         radiation spectra
  #   thermal     temperatures, aerodynamic resistances and friction velocity
  #   bcu, bch    leaf biochemical outputs for sunlit and shaded leaves,
  #               respectively

  ## 1. initialisations and other preparations for the iteration loop




  # Set initial condition
  CONT <- TRUE  # energy balance closure
  counter <- 0  # iteration counter
  maxit = 400;
  maxEBer = 1;
  Wc = 1;

  leafbio<-getinputLUT(inputLUT=inputLUT, dataset='leafbio')
  canopy<-getinputLUT(inputLUT=inputLUT, dataset='canopy')
  meteo<-getinputLUT(inputLUT=inputLUT, dataset='meteo')
  soil<-getinputLUT(inputLUT=inputLUT, dataset='soil',
                    calc.heat = calc.heat,
                    calc.rss_rbs =  calc.rss_rbs)

  # functions for saturated vapour pressure

  es_fun <- function(Temp) (6.107 * 10 ^(7.5 * Temp /(237.3 + Temp)))
  s_fun <- function(ei,Temp) (ei * 2.3026 * 7.5 * 237.3 / (237.3 + Temp)^2)

  # Define functions for saturated vapour pressure
  es_fun <- function(Temp) {6.107 * 10^(7.5 * Temp / (237.3 + Temp))}
  s_fun <- function(es, Temp) {es * 2.3026 * 7.5 * 237.3 / (237.3 + Temp)^2}

  # Retrieve variables from the "meteo" and "canopy" objects
  Ta = meteo[['Ta']] # Air temperature
  ea = meteo[['ea']] # Water vapour pressure
  Ca = meteo[['Ca']] # CO2 concentration
  Ts <- matrix(meteo[['Ta']], nrow = 2, ncol = 1, byrow = TRUE)  ; # Initialize a matrix with two rows, one for each canopy layer, with the same values as Ta
  p = meteo[['p']]     # Air pressure

  nl = canopy[['nlayers']] # Number of canopy layers

  # Assign values to variables
  Rnuc <- rad$Rnuc

  Tch <- (meteo[['Ta']] + 0.1) * matrix(1, nl, 1)  # Leaf temperature (shaded leaves). Creates a nl x 1 matrix with all elements equal to 1, and multiplies it by (Ta + 0.1)
  Tcu <- (meteo[['Ta']] + 0.3) * matrix(1, nrow = nrow(Rnuc), ncol = ncol(Rnuc), byrow = TRUE)  # Leaf temperature (sunlit leaves). Creates a matrix with the same dimensions as Rnuc, filled with (Ta + 0.3)
  ech <- meteo[['ea']] * matrix(1, nl, 1)  # Leaf boundary vapour pressure (shaded/sunlit leaves). Creates a nl x 1 matrix with all elements equal to ea
  Cch <- meteo[['Ca']] * matrix(1, nl, 1)  # Leaf boundary CO2 (shaded leaves). Creates a nl x 1 matrix with all elements equal to Ca
  ecu <- meteo[['ea']] + 0 * Rnuc  # Leaf boundary vapour pressure (sunlit leaves). Creates a vector of length equal to the length of Rnuc, filled with 0s and adds it to ea
  Ccu <- meteo[['Ca']] + 0 * Rnuc  # Leaf boundary CO2 (sunlit leaves). Creates a vector of length equal to the length of Rnuc, filled with 0s and adds it to Ca
  meteo[['L']] <- -1  # Set Monin-Obukhov length to -1

  # Assign values to constants
  MH2O =  subset(data.constant,constant == 'MH2O')[[2]]
  Mair =  subset(data.constant,constant == 'Mair')[[2]]
  rhoa = subset(data.constant,constant == 'rhoa')[[2]]
  cp = subset(data.constant,constant == 'cp')[[2]]
  sigmaSB = subset(data.constant,constant == 'sigmaSB')[[2]]
  Ps = gap[['Ps']];

  # Assign values to canopy parameters
  nl = canopy[['nl']];
  kV = canopy[['kV']];
  xl = canopy[['xl']];
  LAI = canopy[['LAI']];


  # Assign values to variables and constants
  e_to_q <- MH2O/ Mair /  meteo[['p']]  # Conversion of vapour pressure [Pa] to absolute humidity [kg kg-1]
  Fs <- c(1 - gap[['Ps']][nrow(gap[['Ps']])], gap[['Ps']][nrow(gap[['Ps']])])  # Matrix containing values for 1-Ps and Ps of soil
  Fc <- sum(1 - gap[['Ps']][1:(nrow(gap[['Ps']]) - 1)]) / nl  # Matrix containing values for Ps of canopy
  fV <- exp(canopy[['kV']] * canopy[['xl']][1:(nl - 1)])  # Vertical profile of Vcmax

  # Create fVu matrix
  if (ncol(rad$Rnuc) > 1) {
    fVu <- array(1, dim = c(13, 36, nl))  # Create a 13 x 36 x nl array filled with 1s
    for (i in 1:nl) {
      fVu[, , i] <- fV[i]  # Replace values in each layer of the array with corresponding value of fV
    }
  } else {
    fVu <- fV  # If there's only one column in Rnuc, set fVu equal to fV
  }


  ## 2. Energy balance iteration loop

  # Energy balance loop (Energy balance and radiative transfer)
  while(CONT) {                         # while energy balance does not close
    # 2.1. Net radiation of the components
    # Thermal radiative transfer model for vegetation emission (with Stefan-Boltzman's equation)
    rad <- RTMt_sb(constants, rad, soil, leafbio, canopy, gap, Tcu, Tch, Ts[2], Ts[1], 0)
    Rnhc <- rad$Rnhc + rad$Rnhct             # Canopy (shaded) net radiation
    Rnuc <- rad$Rnuc + rad$Rnuct             # Canopy (sunlit) net radiation
    Rnhs <- rad$Rnhs + rad$Rnhst             # Soil (sun+sh) net radiation
    Rnus <- rad$Rnus + rad$Rnust
    Rns <- c(Rnhs, Rnus)
    #    rad2  = RTMt_sb(constants,rad,soil,leafbio,canopy,gap,Tcu,Tch,Ts(2),Ts(1),0);

    # 2.2. Aerodynamic roughness
    # calculate friction velocity [m s-1] and aerodynamic resistances [s m-1]
    resist_out  = resistances(constants,soil,canopy,meteo);

    meteo.ustar = resist_out[['ustar']];
    raa  = resist_out[['raa']];
    rawc = resist_out[['rawc']];
    raws    = resist_out[['raws']];


    # 2.3. Biochemical processes
    meteo_h <- meteo
    meteo_u <- meteo
    meteo_h[['T_']] <- mean(Tch) * (1 - Fc) + mean(Tcu) * Fc
    meteo_h[['eb']] <- mean(ech) * (1 - Fc) + mean(ecu) * Fc
    meteo_h[['Cs']] <- mean(Cch) * (1 - Fc) + mean(Cch) * Fc
    meteo_h[['Q']] <- mean(rad[['Pnh_Cab']]) * (1 - Fc) + mean(rad[['Pnu_Cab']]) * Fc
    #   meteo_u$T _<- mean(Tcu)
    #   meteo_u$eb <- mean(ecu)
    #   meteo_u$Cs <- mean(Ccu)
    #   meteo_u$Q <- mean(rad$Pnu_Cab)
    if (options$Fluorescence_model == 2) {
      b <- biochemical_MD12
    } else {
      b <- biochemical
    }
    #   bch <- b(leafbio, meteo_h, options, constants, fV)
    #   bcu <- b(leafbio, meteo_u, options, constants, fVu)
    bch <- b(leafbio, meteo_h, options, constants, 1)
    #   bcu <- b(leafbio, meteo_u, options, constants, 1)
    # 2.4. Fluxes (latent heat flux (lE), sensible heat flux (H) and soil heat flux G
    # in analogy to Ohm's law, for canopy (c) and soil (s). All in units of [W m-2]
    rss <- soil$rss
    rac <- (LAI + 1) * (raa + rawc)
    ras <- (LAI + 1) * (raa + raws)
    heatfluxes_output <- heatfluxes(rac, bch$rcw, meteo_h$T, ea, Ta, e_to_q, Ca, bch$Ci, constants, es_fun, s_fun)
    lEch <- heatfluxes_output[['lEch']]
    Hch <- heatfluxes_output[['Hch']]
    ech <- heatfluxes_output[['ech']]
    Cch <- heatfluxes_output[['Cch']]
    lambdah <- heatfluxes_output[['lambdah']]
    sh <- heatfluxes_output[['sh']]

    # [lEcu,Hcu,ecu,Ccu,lambdau,su]     = heatfluxes(rac,bcu.rcw,meteo_u.T,ea,Ta,e_to_q,Ca,bcu.Ci,constants, es_fun, s_fun);
    # Hctot <- LAI * (Fc * Hch + meanleaf(canopy, Hcu, 'angles_and_layers', Ps));
    Hctot <- LAI * (Fc * Hch)
    Hstot <- Fs * heatfluxes(ras, rss, Ts, ea, Ta, e_to_q, Ca, Ca, constants, es_fun, s_fun)$Hs
    Htot <- Hstot + Hctot

    # ground heat flux
    G <- 0.35 * Rns

    # 2.5. Monin-Obukhov length L
    meteo[['L']] <- Monin_Obukhov(constants, meteo, Htot)  # [1]

    # 2.6. energy balance errors, continue criterion and iteration counter
    Rnhc <- Rns * (1 - Fc) - G
    EBerch <- mean(Rnhc) - lEch - Hch
    EBers <- Rns - lEs - Hs - G
    counter <- counter + 1  # Number of iterations
    maxEBerch <- max(abs(EBerch))
    maxEBers <- max(abs(EBers))

    CONT <- (
      maxEBerch > maxEBer |
        maxEBers > maxEBer
    ) & counter < maxit + 1

    if (!CONT) {
      if (any(is.nan(c(maxEBerch, maxEBers)))){
        cat(sprintf('WARNING: NaN in fluxes, counter = %i\n', counter))
      }
      break
    }

    if (counter == 10) Wc <- 0.8
    if (counter == 60) Wc <- 0.6
    if (counter == 100) Wc <- 0.2

    # 2.7. New estimates of soil (s) and leaf (c) temperatures, shaded (h) and sunlit (1)
    Tch <- Tch + Wc * EBerch / ((rhoa * cp) / rac + rhoa * lambdah * e_to_q * sh / (rac + bch[['rcw']]) + 4 * leafbio[['emis']] * sigmaSB * (meteo_h[['T_']] + 273.15)^3)
    Tcu <- mean(Tch) * rep(1, length(Rnuc)) # + Wc * EBercu / ((rhoa * cp) / rac + rhoa * lambdau * e_to_q * su / (rac + bcu$rcw) + 4 * leafbio$emis * sigmaSB * (meteo_u$T + 273.15)^3)
    Ts <- Ts + Wc * EBers / (rhoa * cp / ras + rhoa * lambdas * e_to_q * ss / (ras + rss) + 4 * (1 - soil[['rs_thermal']]) * sigmaSB * (Ts + 273.15)^3)
    Tch[abs(Tch) > 100] <- Ta
    # Tcu[abs(Tcu) > 100] <- Ta

  }




## 3. Print warnings whenever the energy balance could not be solved
  if (counter >= maxit) {
    cat("WARNING: maximum number of iteratations exceeded\n")
    cat(sprintf("Energy balance error sunlit vegetation = %.2f W m-2\n", maxEBercu))
    cat(sprintf("Energy balance error shaded vegetation = %.2f W m-2\n", maxEBerch))
    cat(sprintf("Energy balance error soil              = %.2f W m-2\n", maxEBers))
  }


## 4. some more outputs


  iter.counter <- counter
  # this is not yet written to output but making it avg here breaks RTMt_sb
  # thermal.Tcu     = mean(Tcu(:));
  # thermal.Tch     = mean(Tch);
  #set thermal.Tcu to the mean of Tch, multiplied by a matrix of ones with the same size as Rnuc
  thermal.Tcu     <- mean(Tch) * array(1, dim = dim(Rnuc))
  thermal.Tch     <- Tch
  thermal.Tsu     <- Ts[2]
  thermal.Tsh     <- Ts[1]

  fluxes$Rnctot <- mean(Rnhc) * (1 - Fc) + mean(Rnuc) * Fc #LAI* aggregator(mean(Rnhc), mean(Rnuc(:)), Fc, Ps, canopy);     # net radiation leaves
  fluxes$lEctot <- LAI * (Fc * lEch) #aggregator(lEch, lEcu, Fc, Ps, canopy);     # latent heat leaves
  fluxes$Hctot  <- LAI * (Fc * Hch) #aggregator(Hch, Hcu, Fc, Ps, canopy);       # sensible heat leaves
  fluxes$Actot  <- LAI * (Fc * bch$A) #aggregator(bch.A, bcu.A, Fc, Ps, canopy);   # photosynthesis leaves
  fluxes$Tcave  <- aggregator(mean(Tch), mean(Tcu), Fc, Ps, canopy)            # mean leaf temperature

  fluxes$Rnstot <- Fs * Rns           # Net radiation soil
  fluxes$lEstot <- Fs * lEs           # Latent heat soil
  fluxes$Hstot  <- Fs * Hs            # Sensible heat soil
  fluxes$Gtot   <- Fs * G             # Soil heat flux
  fluxes$Tsave  <- Fs * Ts            # Soil temperature
  # fluxes$Resp   <- Fs * equations.soil_respiration(Ts) #  Soil respiration = 0

  # calculate the total sensible heat as the sum of the sensible heat leaves and the sensible heat soil
  fluxes$Rntot <- fluxes$Rnctot + fluxes$Rnstot
  # calculate the total latent heat as the sum of the latent heat leaves and the latent heat soil
  fluxes$lEtot <- fluxes$lEctot + fluxes$lEstot
  # calculate the total sensible heat as the sum of the sensible heat leaves and the sensible heat soil
  fluxes$Htot <- fluxes$Hctot + fluxes$Hstot
  fluxes$rss <- rss

  ## faking back layer structure of SCOPE

  bch     <- b(leafbio, meteo_h, options, constants, fV)
  bcu     <- b(leafbio, meteo_h, options, constants, fVu)  # notice meteo_h, not _u

}


#'  calculate the total flux (e.g., total radiation, latent heat, sensible heat, photosynthesis) in a vegetation canopy, considering both shaded and sunlit components. The variable Fc represents the fractional area of the canopy that is shaded, and Ps represents the angles and layers of the canopy.
#'
#' @param shaded_flux  vector of shaded fluxes
#' @param sunlit_flux  vector of  sunlit fluxes
#' @param Fc  is the fraction of canopy that is sunlit
#' @param Ps  is used in meanleaf() which calculates the average flux.
#' @param canopy  is used in meanleaf() which calculates the average flux.
#'
#' @return
#' which calculates the average flux. If the number of columns in sunlit_flux is greater than 1,
#' the mean leaf flux is calculated using meanleaf(),
#' otherwise, the flux is calculated using the weighted average of shaded_flux and sunlit_flux.
#' The result is returned as flux_tot.
#' @export
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
aggregator <- function(shaded_flux, sunlit_flux, Fc, Ps, canopy) {
  if (ncol(sunlit_flux) > 1) {
    flux_tot <- Fc*shaded_flux + meanleaf(canopy, sunlit_flux, 'angles_and_layers', Ps)
  } else {
    flux_tot <- Fc*shaded_flux + (1-Fc)*sunlit_flux
  }
  return(flux_tot)
}


