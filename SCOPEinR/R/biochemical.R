#' @title get.biochemical
#'
#' @description
#' \code{get.biochemical} Calculates net assimilation rate A, fluorescence F using biochemical model
#'
#'  This function calculates:
#'    - stomatal resistance of a leaf or needle (s m-1)
#'    - photosynthesis of a leaf or needle (umol m-2 s-1)
#'    - fluorescence of a leaf or needle (fraction of fluor. in the dark)
#'
#' @param data.leafbio list. Leaf biochemical parameters: \code{Type} ('C3' or 'C4'), \code{stressfactor}, \code{Vcmax25}, \code{BallBerry0}, \code{BallBerrySlope}, \code{Rdparam}, \code{Kn0}/\code{Knalpha}/\code{Knbeta} (NPQ model parameters), optional \code{g_m} (mesophyll conductance, mol m^-2 s^-1 bar^-1), and (when temperature correction is enabled) the \code{TDP} sub-list of temperature-dependence parameters.
#' @param data.meteo list. Meteorological/leaf micro-environment characteristics: \code{Q} (absorbed PAR, umol photons m^-2 s^-1), \code{Cs} (CO2 concentration in the leaf boundary layer, ppm), \code{Temp} (leaf temperature, deg C or K), \code{eb} (vapour pressure in the leaf boundary layer, hPa), \code{Oa} (O2 concentration, mmol/mol), \code{p} (air pressure, hPa).
#' @param data.opts data.frame. Simulation options; row 7 (\code{tempcor}) is a 0/1 flag for whether temperature correction of Vcmax and the rate constants is applied.
#' @param fV numeric. Scaling factor applied to \code{Vcmax25} (e.g. the canopy nitrogen/Vcmax profile factor estimated in the energy balance).
#' @param get.plots logical. If \code{TRUE}, returns diagnostic plots for stomatal conductance, assimilation, and Jmax/Vcmax rates.
#'
#' @return A list \code{biochem_out} with (among others): \code{A} (net assimilation rate, umol m^-2 s^-1), \code{Ci}/\code{Cc} (internal/chloroplast CO2 concentration, ppm), \code{rcw} (stomatal resistance, s m^-1), \code{gs} (stomatal conductance), \code{Vcmax}, \code{Rd}, \code{Ja} (actual electron transport rate), \code{ps}/\code{ps_rel} (photochemical yield and degree of light saturation), \code{eta} (fluorescence yield relative to dark-adapted), \code{qE}/\code{qQ} (non-photochemical/photochemical quenching), \code{fs}/\code{fo}/\code{fm} (steady-state/dark/light-saturated fluorescence yields), \code{Kn}/\code{NPQ}/\code{Kf}/\code{Kp}/\code{Kd} (rate constants), and \code{SIF} (solar-induced fluorescence, \code{fs * Q}).
#' @export
#'
#' @author 	Joe Berry and Christiaan van der Tol, Ari Kornfeld (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' last updates
#' Date: 	21 Sep 2012.
#'
#' Update:   20 Feb 2013.
#' Update:      Aug 2013: correction of L171: Ci <- Ci*1e6/ p * 1E3
#'
#' Update:   2016-10 - (JAK) major rewrite to accomodate an iterative solution to the Ball-Berry equation
#'                   - also allows for g_m to be specified for C3 plants, but only if Ci_input is provided.
#'
#' Update:   25 Feb 2021: Temperature reponse functions by Dutta et al. implemented
#'
#' @references
#'
#' Farquhar et al. 1980, Collatz et al (1991, 1992), and:
#
#' Dutta, D., Schimel, D. S., Sun, Y., Tol, C. V. D., & Frankenberg, C. (2019).
#' Optimal inverse estimation of ecosystem parameters from observations of carbon and energy fluxes.
#' Biogeosciences, 16(1), 77-103.
#'
#' Van der Tol, C., Berry, J. A., Campbell, P. K. E., & Rascher, U. (2014).
#' Models of fluorescence and photosynthesis for interpreting measurements
#' of solar induced chlorophyll fluorescence. Journal of Geophysical Research: Biogeosciences, 119(12), 2312-2327.
#'
#' Bonan, G. B., Lawrence, P. J., Oleson, K. W., Levis, S., Jung, M., Reichstein, M., ... & Swenson, S. C. (2011).
#' Improving canopy processes in the Community Land Model version 4 (CLM4)
#' using global flux fields empirically inferred from FLUXNET data.
#' Journal of Geophysical Research: Biogeosciences, 116(G2).
#'

#' @examples
#' \dontrun{
#' biochem_out <- get.biochemical(data.leafbio, data.meteo, data.opts, fV = 1, get.plots = FALSE)
#' }
get.biochemical<-function(data.leafbio,data.meteo,data.opts,fV, get.plots){

  # Calculates net assimilation rate A, fluorescence F using biochemical model
  # Input (units are important):

  #                                   example soil moisture, leaf age). Use 1 to "disable" (1 <- no stress)
  #  OPTIONAL
  # Kpep25 (kp)   # [umol/m2/s]         PEPcase activity at 25 deg C (defaults to Vcmax/56
  # atheta      # [0-1]                  smoothing parameter for transition between Vc and Ve (light- and carboxylation-limited photosynthesis)
  # useTLforC3  # boolean              whether to enable low-temperature attenuation of Vcmax in C3 plants (its always on for C4 plants)
  # po0         #  double            Kp,0 (Kp,max) <- Fv/Fm (for curve fitting)
  # g_m         # mol/m2/s/bar      Mesophyll conductance (default: Infinity, i.e. no effect of g_m)

  # Note: always use the prescribed units. Temperature can be either oC or K
  # Note: input can be single numbers, vectors, or n-dimensional matrices

  #fV estimate from from ebal
  #fV=1


  if (missing(data.opts)){
    stop('please use options for leaf temperature correction for Vcmax, apply_T_corr ')

  } else{
    # tempcor   # [0, 1]               boolean (0 or 1) whether or not temperature correction to Vcmax has to be applied.

    tempcor <- data.opts[7,]

  }


  ## input
  # atmospheric physics
  constants <- SCOPEinR::constants
  ## constants
  rhoa <-  subset(constants,constant == 'rhoa')[[2]]           # [kg m-3]       specific mass of air
  Mair <-  subset(constants,constant == 'Mair')[[2]]           # [g mol-1]      molecular mass of dry air
  R <- subset(constants,constant == 'R')[[2]]                 # [J mol-1K-1]   Molar gas constant

  #Q: absorbed PAR flux net radiation, PAR  in  [umol photons m-2 s-1]
  # fPAR  fraction of incident light that is absorbed by the leaf (default <- 1, for compatibility)    # [0-1]
  Q <- data.meteo[['Q']]

  #Cs: initial estimate of conc. of CO2 in the bounary layer of the leaf  in [ppmV or umol mol]
  Cs <- data.meteo[['Cs']]

  # T leaf temperature in [oC or K]
  # convert temperatures to K if not already
  if (data.meteo[['Temp']][1] < 200) {
    T_k <- data.meteo[['Temp']] + 273.15
  } else {
    T_k <- data.meteo[['Temp']]
  }


  # eb: intial estimate of the vapour pressure in leaf boundary layer  in [hPa <- mbar]
  eb <- data.meteo[['eb']]

  # O: concentration of O2 (in the boundary layer, but no problem to use ambient) in [mmol/mol]
  O <- data.meteo[['Oa']]

  # p:  air pressure in [hPa]
  p <- data.meteo[['p']]

  #biochemical parameters
  # Type      # ['C3', 'C4']          text parameter, either 'C3' for C3 or 'any'C4' for C4
  Type <- data.leafbio[['Type']]

  # stressfactor [0-1]               stress factor to reduce Vcmax (for
  stressfactor <- data.leafbio[['stressfactor']]


  # Vcmax25 (Vcmo): maximum carboxylation capacity @ 25 degC in [umol/m2/s]
  Vcmax25 <- fV * data.leafbio[['Vcmax25']]
  # BallBerry0 # []   (OPTIONAL) Ball-Berry intercept term 'b' (if present, an iterative solution is used)
  # setting this to zeo disables iteration. Default <- 0.01
  BallBerry0 <- data.leafbio[['BallBerry0']]

  # BallBerrySlope (m): Ball-Berry coefficient 'm' for stomatal regulation in []
  BallBerrySlope <- data.leafbio[['BallBerrySlope']] ###

  # RdPerVcmax25 (Rdparam); respiration as fraction of Vcmax25 in []
  RdPerVcmax25 <- data.leafbio[['Rdparam']] ### %   respiration as fraction of Vcmax25

  # effcon: number of CO2 per electrons - typically 1/5 for C3 and 1/6 for C4   [mol CO2/mol e-]
  if (Type == "C3") {
    ## for C3 plant
    effcon <- 1/5
  } else {
    ## for C4 plant
    effcon <- 1/6  # C4
  }


  Tref <- 25 + 273.15;        # [K]           absolute temperature at 25 oC

  #Jmax25      <- 1.97*Vcmax25;
  #TPU25           <- 0.06*Jmax25;              # Triose Phosphate utilization rate

  # All of the Following Values are Adopted from Bonan et. al 2011 paper and
  # pertaings to C3 photosythesis

  Kc25 <- 405                  # [umol mol-1]
  Ko25 <- 279                    # [mmol mol-1]
  spfy25 <- 2444                   # specificity (Computed from Bernacchhi et al 2001 paper)

  # convert all to bar: CO2 was supplied in ppm, O2 in permil, and pressure in mBar
  ppm2bar <-  1e-6 * (p *1E-3)

  # Cs: CO2 concentration in the boundary layer
  Cs <- Cs * ppm2bar
  O <- (O * 1e-3) * (p *1E-3) * ifelse(Type == 'C3', 1, 0)# force O to be zero for C4 vegetation (this is a trick to prevent oxygenase)
  Kc25 <- Kc25 * 1e-6
  Ko25 <- Ko25 * 1e-3
  Gamma_star25 <- 0.5  *O / spfy25      # [ppm] compensation point in absence of Rd
  Rd25 <- RdPerVcmax25 * Vcmax25

  atheta  <- 0.8

  # Mesophyll conductance: by default we ignore its effect
  #  so Cc <- Ci - A/gm <- Ci
  g_m <- Inf

  if ("g_m" %in% names(data.leafbio)) {
    g_m <-  data.leafbio[['g_m']] * 1e6  # convert from mol to umol
  }

  # fluorescence parameters

  # Knparams:  parameters for empirical Kn (NPQ) model: Kn <- Kno * (1+beta).*x.^alpha./(beta + x.^alpha); [Kno, Kn_alpha, Kn_beta]  # [], [], []
  #
  #  or, better, as individual fields:
  #   Kno                                     Kno - the maximum Kn value ("high light")
  #   Kn_alpha, Kn_beta                      alpha, beta: curvature parameters
  #
  #

  Knparams <- c(data.leafbio[['Kn0']], data.leafbio[['Knalpha']], data.leafbio[['Knbeta']])
  Kf <- 0.05             # []     rate constant for fluorescence
  Kd <- max(0.8738, 0.0301*(T_k - 273.15) + 0.0773)
  Kp <- 4.0              # []     rate constant for photochemistry

  ## temperature corrections
  fluorescence <- list('Vcmax' = 1, 'Rd' = 1, 'TPU' = 1, 'Kc' = 1, 'Ko' = 1, 'Gamma_star' = 1)


  ## physiological options

  if (tempcor[1,'Value'] == 1){
    #message("Correction of Vcmax and rate constants for temperature is processed")

    if (Type == "C4") {

      # RdPerVcmax25 <- 0.025;  # Rd25 for C4 is different than C3
      # Rd25 <- RdPerVcmax25 * Vcmax25;
      # Constant parameters for temperature correction of Vcmax
      Q10 <- data.leafbio[['TDP']]$Q10  # Unit is  []
      s1 <- data.leafbio[['TDP']]$s1   # Unit is [K]
      s2 <- data.leafbio[['TDP']]$s2   # Unit is [K^-1]
      s3 <- data.leafbio[['TDP']]$s3   # Unit is [K]
      s4 <- data.leafbio[['TDP']]$s4   # Unit is [K^-1]

      # Constant parameters for temperature correction of Rd
      s5 <- data.leafbio[['TDP']]$s5  # Unit is [K]
      s6 <- data.leafbio[['TDP']]$s6  # Unit is [K^-1]

      fHTv <- 1 + exp(s1 * (T_k - s2))
      fLTv <- 1 + exp(s3 * (s4 - T_k))
      Vcmax <- (Vcmax25 * Q10^(0.1 * (T_k - Tref))) / (fHTv * fLTv)  # Temp Corrected Vcmax

      # Temperature correction of Rd
      fHTv <- 1 + exp(s5 * (T_k - s6))
      Rd <- (Rd25 * Q10^(0.1 * (T_k - Tref))) / fHTv  # Temp Corrected Rd

      # Temperature correction of Ke
      Ke25 <- 20000 * Vcmax25  # Unit is  []
      Ke <- (Ke25 * Q10^(0.1 * (T_k - Tref)))  # Temp Corrected Ke
    } else if (Type == "C3") {

      # temperature correction of Vcmax
      deltaHa <- data.leafbio[['TDP']]$delHaV  # Unit is  [J K^-1]
      deltaS <- data.leafbio[['TDP']]$delSV  # unit is [J mol^-1 K^-1]
      deltaHd <- data.leafbio[['TDP']]$delHdV  # unit is [J mol^-1]
      fTv <- get.temperature.functionC3(Tref, R, T_k, deltaHa)
      fHTv <- get.high.temp.inhibtionC3(Tref, R, T_k, deltaS, deltaHd)
      fluorescence[['Vcmax']] <- fTv * fHTv

      ## temperature correction for TPU

      # deltaHa <- data.leafbio[['TDP]]$delHaP;                # Unit is  [J K^-1]
      # deltaS  <- data.leafbio[['TDP]]$delSP;                 # unit is [J mol^-1 K^-1]
      # deltaHd <- data.leafbio[['TDP]]$delHdP;                # unit is [J mol^-1]
      # fTv <- get.high.temp.inhibtionC3(Tref,R,T_k,deltaHa);
      # fHTv <- get.high.temp.inhibtionC3(Tref,R,T_k,deltaS,deltaHd);
      # f[['TPU']] <- fTv .* fHTv;


      # temperature correction for Rd
      deltaHa <- data.leafbio[['TDP']]$delHaR    # Unit is  [J K^-1]
      deltaS  <- data.leafbio[['TDP']]$delSR    # unit is [J mol^-1 K^-1]
      deltaHd <- data.leafbio[['TDP']]$delHdR   # unit is [J mol^-1]
      fTv <- get.temperature.functionC3(Tref,R,T_k,deltaHa);
      fHTv  <- get.high.temp.inhibtionC3(Tref,R,T_k,deltaS,deltaHd)
      fluorescence[['Rd']] <- fTv * fHTv;



      deltaHa <- data.leafbio[['TDP']]$delHaKc
      fTv <- get.temperature.functionC3(Tref, R, T_k, deltaHa)
      fluorescence[['Kc']] <- fTv

      deltaHa <- data.leafbio[['TDP']]$delHaKo
      fTv <- get.temperature.functionC3(Tref, R, T_k, deltaHa)
      fluorescence[['Ko']] <- fTv

      deltaHa <- data.leafbio[['TDP']]$delHaT
      fTv <- get.temperature.functionC3(Tref, R, T_k, deltaHa)
      fluorescence[['Gamma_star']] <- fTv

      Ke <- 1
    }
  } else {
   # message("Correction of Vcmax and rate constants for temperature is not processed")

    Ke <- 1
  }



  if (Type == "C3") {
    Vcmax <- Vcmax25 * fluorescence[['Vcmax']] * stressfactor
    Rd <- Rd25 * fluorescence[['Rd']] * stressfactor
    Kc <- Kc25 * fluorescence[['Kc']]
    Ko <- Ko25 * fluorescence[['Ko']]

  }
  Gamma_star <- Gamma_star25 * fluorescence[['Gamma_star']]




  # calculation of potential electron transport rate
  po0 <- Kp / (Kf + Kd + Kp)        # maximum dark photochemistry fraction, i.e. Kn <- 0 (Genty et al., 1989)
  Je <- 0.5 * po0 * Q            # potential electron transport rate (JAK: add fPAR);
  # Gamma_star  <- 0.5 * O / spfy; %[bar]       compensation point in absence of Rd (i.e. gamma*) [bar]

  if (Type == 'C3') {
    MM_consts <- (Kc * (1 + O / Ko))  # Michaelis-Menten constants
    Vs_C3 <- Vcmax/2
    # minimum Ci (as fraction of Cs) for BallBerry Ci. (If Ci_input is present we need this only as a placeholder for the function call)
    minCi <- 0.3
  } else {


    # C4
    MM_consts <- 0               # just for formality, so MM_consts is initialized
    Vs_C3 <- 0                   # the same
    minCi <- 0.1                 # C4
  }

  ## calculation of Ci (internal CO2 concentration)
  RH <- pmin(1, eb/satvap(T_k-273.15))

  ## get CO2 per electron and Assimilation

  # Definir la función computeA_fun
  computeA_fun <- function(x) {
    get.computeA(x, Type, g_m, Vs_C3, MM_consts, Rd, Vcmax, Gamma_star, Je, effcon, atheta, Ke)

    }
  # Definir la función for estimating Ci_out
  fun_Ci <- function(x) {
    Ci_values <-  get.Ci.next(x, Cs, RH, minCi, BallBerrySlope, BallBerry0, computeA_fun, ppm2bar)
    element <- Ci_values$err[i]  # Obtener el elemento específico correspondiente a 'i'
    return(element)
  }


  if (all(BallBerry0 == 0)) {
    # b <- 0: no need to iterate:
    Ci <- get.BallBerry(Cs, RH, NULL, BallBerrySlope, BallBerry0, minCi)[[2]]

    # A <-  computeA_fun(Ci)
  } else {
    # compute Ci using iteration (JAK)
    # it would be nice to use a built-in root-seeking function but fzero requires scalar inputs and outputs,
    # Here I use a fully vectorized method based on Brent's method (like fzero) with some optimizations.
    maxiter <- 1000
    tol <- 1e-7  # 0.1 ppm more-or-less


    Ci = c()
    for (i in 1:length(Cs)) {

      lower <- Cs[i] - 0.001
      upper <- Cs[i] + 0.001
      result <- pracma::brent(fun_Ci, a=lower, b=upper, maxiter = maxiter, tol = tol)
      #result <- pracma::brentDekker(fun_Ci, Cs[i]-10, b=Cs[i] + 100000, maxiter = maxiter, tol = tol)
      Ci[i] <- result$root

    }


    # NOTE: A is computed in Ci_next on the final returned Ci. fixedp_brent_ari() guarantees that it was done on the returned values.
    # A <-  computeA_fun(Ci)
  }

  ## Compute A, electron transport rate, and stomatal resistance
  params.biochem <- computeA_fun(Ci)

  # A: net assimilation rate of the leaves in (umol/m2/s)
  A <- params.biochem$A

  Ag <- params.biochem$Ag

  CO2_per_electron <- params.biochem$CO2_per_electron


  # stomatal conductance
  gs <- pmax(0, 1.6 * A * ppm2bar / (Cs - Ci))

  # actual electron transport rate
  Ja <- Ag / CO2_per_electron

  # stomatal resistance
  rcw <- (rhoa / (Mair * 1E-3)) / gs

  ps <- po0 * Ja / Je  # this is the photochemical yield

  nanPs <- is.na(ps)
  if (any(nanPs)) {
    if (length(po0) == 1) {
      ps[nanPs] <- po0
    } else {
      ps[nanPs] <- po0[nanPs]  # happens when Q <- 0, so ps <- po0 (other cases of NaN have been resolved)
    }
  }

  ps_rel <- pmax(0, 1 - ps / po0)  # degree of light saturation: 'x' (van der Tol e.a. 2014)


  db.Fluoresce  <- get.Fluorescence.model(ps, ps_rel, Kp,Kf,Kd,Knparams);
  # fluorescence as fraction of dark adapted (fs/fo)
  eta <- db.Fluoresce[['eta']]

  # qE: Non photochemical quenching
  qE <- db.Fluoresce[['qE']]

  # qQ photochemical quenching
  qQ <- db.Fluoresce[['qQ']]


  fs <- db.Fluoresce[['fs']]

  fo <- db.Fluoresce[['fo']]

  fm <- db.Fluoresce[['fm']]

  fm0 <- db.Fluoresce[['fm0']]

  fo0 <- db.Fluoresce[['fo0']]

  Kn <- db.Fluoresce[['Kn']]


  Kpa         <- ps /fs * Kf


  Cc <- NULL
  if (!is.null(g_m)) {
    Cc <- (Ci - A/g_m) / ppm2bar
  }
  Ci <- Ci / ppm2bar


  # Collect outputs

  # Output:

  # structure 'biochem_out' with the following elements:
  # A         # [umol/m2/s]           net assimilation rate of the leaves
  # Cs        # [umol/m3]             CO2 concentration in the boundary layer
  # eta0      # []                    fluorescence as fraction of dark
  #                                   ...adapted (fs/fo)
  # rcw       # [s m-1]               stomatal resistance
  # qE        # []                    non photochemical quenching
  # fs        # []                    fluorescence as fraction of PAR
  # Ci        # [umol/m3]             internal CO2 concentration
  # Kn        # []                    rate constant for excess heat
  # fo        # []                    dark adapted fluorescence (fraction of aPAR)
  # fm        # []                    light saturated fluorescence (fraction of aPAR)
  # qQ        # []                    photochemical quenching
  # Vcmax     # [umol/m2/s]           carboxylation capacity after
  #                                   ... temperature correction

  biochem_out <- list(A = A, Ci = Ci)

  if (get.plots ==  T){
    check.rate <- data.frame(profile=c(1:length(A)),A=A,Gs=gs,Vcmax=Vcmax,Jmax=Ja)


    p.rate <- ggplot(data = check.rate, aes(y = profile)) +
      labs(y= "nlayers", x = "")+
      geom_line(aes(x = Jmax, color = "Jmax"), linewidth = 0.5) +
      geom_line(aes(x = A, color = "A"), linewidth = 0.5) +
      geom_line(aes(x = Gs, color = "Gs"), linewidth = 0.5) +
      geom_line(aes(x = Vcmax, color = "Vcmax"), linewidth = 0.5) +
      scale_color_manual(name = "Biochemistry",
                         values = c("Jmax" = "forestgreen", 'A' = 'red','Gs' = 'grey',
                                    "Vcmax" = "navyblue")) +
      theme_bw()  +
      theme(legend.position = "top") +
      guides(color = guide_legend(title = NULL))
    print(p.rate)
  }




  if (!is.null(Cc)) {
    biochem_out[['Cc']] <- Cc
  }

  # output results
  biochem_out[['rcw']] <- rcw
  biochem_out[['gs']] <-  gs

  biochem_out[['RH']] <-  RH

  # Vcmax : carboxylation capacity after temperature correction
  biochem_out[['Vcmax']] <- Vcmax
  biochem_out[['Rd']] <- Rd
  biochem_out[['Ja']] <- Ja
  biochem_out[['ps']] <- ps# photochemical yield
  biochem_out[['ps_rel']] <- ps_rel  # degree of ETR saturation 'x' (van der Tol e.a. 2014)
  biochem_out[['Kd']]  <- Kd  # K_dark(T)
  biochem_out[['Kn']] <- Kn  # K_n(x);  x <- 1 - ps/p00 == 1 - Ja/Je
  biochem_out[['NPQ']] <- Kn / (Kf + Kd)#
  biochem_out[['Kf']] <- Kf  # Kf <- 0.05 (const)
  biochem_out[['Kp0']] <- Kp # Kp <- 4.0 (const): Kp, max
  biochem_out[['Kp']]  <- Kpa # Kp,actual


  biochem_out[['eta']] <- eta
  biochem_out[['qE']] <- qE
  biochem_out[['fs']] <- fs  # keep this for compatibility with SCOPE
  biochem_out[['ft']] <- fs  # keep this
  biochem_out[['SIF']] <- fs * Q
  biochem_out[['fo0']] <- fo0
  biochem_out[['fm0']] <- fm0
  biochem_out[['fo']] <- fo
  biochem_out[['fm']] <- fm
  biochem_out[['Fm_Fo']]  <- fm / fo  # parameters used for curve fitting
  biochem_out[['Ft_Fo']] <- fs / fo  # parameters used for curve fitting
  biochem_out[['qQ']] <- qQ
  biochem_out[['Phi_N']]  <- Kn/ (Kn + Kp + Kf + Kd)

  return(biochem_out)

}  # end of function biochemical




