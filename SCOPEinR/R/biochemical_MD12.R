#' @title get.biochemical.MD12
#' @description \code{get.biochemical.MD12} Calculates:
#' - CO2 concentration in intercellular spaces (umol/mol == ppmv)
#' - leaf net photosynthesis (umol/m2/s) of C3 or C4 species
#' - fluorescence yield of a leaf (fraction of reference fluorescence yield in dark-adapted and un-stressed leaf)
#'
#' Note: always use the prescribed units. Temperature can be either oC or K
#' Note: input can be single numbers, vectors, or n-dimensional matrices
#' Note: For consistency reasons, in C4 photosynthesis electron transport rates under CO2-limited conditions are computed by inverting the equation
#' applied for light-limited conditions(Ubierna et al 2013). A discontinuity would result when computing J from ATP requirements of Vp and Vco, as a
#' fixed electron transport partitioning is assumed for light-limited conditions
#'
#' @param data.leafbio  LUT table
#' @param data.meteo meteo characterisitics with L (Monin-Obukhov length), also carries the absorbed PAR (Q) and leaf temperature used internally.
#' @param fV fraction of Vcmax25 downregulated as function of cumulative absorbed PAR through the canopy (relative activity profile).
#' @param get.plots return plots for gs, assimilation and Jmax, Vcmax rate
#'
#' @return the following  parameters at eaf level:
#'  - A in umol/m2/s which is the net assimilation rate of the leaves
#'  - Ci in umol/mol which is the CO2 concentration in intercellular spaces (assumed to be the same as at carboxylation sites in C3 species)
#' - eta in (-) which is the amplification factor to be applied to PSII fluorescence yield spectrum
#'   relative to the dark-adapted, un-stressed yield calculated with either Fluspect or FluorMODleaf
#' @export
#' @description
#'
#' Date: 21 Sep 2012
#' Update:
#'  - 28 Jun 2013 Adaptation for use of Farquhar model of C3 photosynthesis (Farquhar et al 1980).
#'  - 18 Jul 2013 Inclusion of von Caemmerer model of C4 photosynthesis (von Caemmerer 2000, 2013).
#'  - 15 Aug 2013 Modified computation of CO2-limited electron transport in C4 species for consistency with light-limited value.
#'  - 22 Oct 2013 Included effect of qLs on Jmax and electron transport value of kNPQs re-scaled in input as NPQs.
#'  - 08 Jan 2019 (CvdT): minor modification to adjust to SCOPE_lite.

#
#' @author Federico Magnani, with contributions from Christiaan van der Tol  (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' out <- get.biochemical.MD12(data.leafbio, data.meteo, fV = 1, get.plots = FALSE)
#' }

get.biochemical.MD12<-function(data.leafbio,data.meteo,fV, get.plots){


  # Input (units are important when not otherwise specified, mol refers to mol C):
  constants <- SCOPEinR::constants
  # p  air pressure in [Pa]
  p = data.meteo[['p']] * 1e2

  # m: Ball-Berry coefficient 'm' for stomatal regulation in [mol/mol]

  BallBerrySlope = data.leafbio[['BallBerrySlope']]
  BallBerry0 = data.leafbio[['BallBerry0']]
  # O: ambient O2 concentration in [mmol/mol]
  O = data.meteo[['Oa']]

  # Cs: CO2 concentration at leaf surface in [umol/mol]
  Cs =data.meteo[['Cs']]

  # Type:   text parameter, either 'C3' for C3 or any other text for C4
  Type = data.leafbio[['Type']]

  # Tyear mean annual temperature in [°C]
  Tyear = data.leafbio[['Tyear']]

  # beta: fraction of photons partitioned to PSII (0.507 for C3, 0.4 for C4 Yin et al. 2006 Yin and Struik 2012) in []
  beta = data.leafbio[['beta']]

  # qLs:  fraction of functional reaction centres (Porcar-Castell 2011) in []
  qLs = data.leafbio[['qLs']]

  # NPQs: rate constant of sustained thermal dissipation, normalized to (kf+kD) (=kNPQs' Porcar-Castell 2011) in [s-1]
  NPQs = data.leafbio[['kNPQs']]

  # stress: optional input: stress factor to reduce Vcmax (for example soil moisture, leaf age). Default value = 1 (no stress). # []
  stress=data.leafbio[['stressfactor']]



  ##data.meteo = meteo_h

  # Q: photochemically active radiation absorbed by the leaf in [uE/m2/s]
  # Multiply the net PAR absorbed by Cab (in mol/m²/s) by the conversion factor of 2.3 μmol/J.
  # Photochemically active radiation absorbed by the leaf (in μmol/m²/s) = Net PAR absorbed by Cab (in mol/m²/s) * 2.3 μmol/J

  #Q = mean(data.meteo[['Q']]) *  2.3
  Q = data.meteo[['Q']]
  # Temp:  leaf temperature in [°C or K]
  Temp = data.meteo[['Temp']] #[1] ##nl

  # eb: vapour pressure in leaf boundary layer in [hPa]
  eb      = data.meteo[['eb']]

  # Vcmax25:  maximum carboxylation capacity at reference leaf temp in [umol/m2/s]
  Vcmax25    = fV * data.leafbio[['Vcmax25']]

  # Rdparam: respiration at reference temperature as fraction of Vcmax in[mol/mol]
  Rdparam =  data.leafbio[['Rdparam']] ##Check if its correct, original indicate RdPerVcmax25

  Q[Q == 0] <- 1e-9


  ## Global and site-specific constants
  R =  subset(SCOPEinR::constants,constant == 'R')[[2]] # [J/K/mol]     universal gas constant

  #---------------------------------------------------------------------------------------------------------
  ## Unit conversion and computation of environmental variables
  if (Temp[1] < 100){
    Temp = Temp + 273.15
  } else {
    Temp = data.meteo['Temp']
  }

  # [K]           convert temperatures to K if not already
  RH = pmin(1, eb / SCOPEinR::satvap(Temp-273.15))                       # []            relative humidity (decimal)
  Cs = Cs * p * 1.0E-11                            # [bar]         1E-6 to convert from ppm to fraction, 1E-5 to convert from Pa to bar
  O  = O * p *1.0e-08                            # [bar]         1E-3 to convert from mmol/mol to fraction, 1E-5 to convert from Pa to bar

  #---------------------------------------------------------------------------------------------------------
  ## Define photosynthetic parameters (at reference temperature)
  SCOOP     = 2862                                    # [mol/mol]     Relative Rubisco specificity for CO2 vs O2 at ref temp (Cousins et al. 2010)
  Rdopt     = Rdparam * Vcmax25                           # [umol/m2/s]   dark respiration at ref temperature from correlation with Vcmax25
  if (Type == 'C3'){
    Jmo =  Vcmax25 * 2.68                            # [umol/m2/s]   potential e-transport at ref temp from correlation with Vcmax25 (Leuning 1997)
  } else {    # C4 species
    Jmo =  Vcmax25 * 40 / 6                           # [umole-/m2/s] maximum electron transport rate (ratio as in von Caemmerer 2000)
    Vpmo  =  Vcmax25 * 2.33                             # [umol/m2/s]   maximum PEP carboxylase activity (Yin et al. 2011)
    Vpr =  80                                     # [umol/m2/s]   PEP regeneration rate, constant (von Caemmerer 2000)
    gbs =  (0.0207 * Vcmax25 + 0.4806) * 1000           # [umol/m2/s]   bundle sheath conductance to CO2 (Yin et al. 2011)
    x =  0.4                                     # []            partitioning of electron transport to mesophyll (von Caemmerer 2013)
    alpha =  0                                      # []            bundle sheath PSII activity (=0 in maize/sorghum >=0.5 in other cases von Caemmerer 2000)
  }                                           # C3 species

  #---------------------------------------------------------------------------------------------------------
  ## Parameters for temperature corrections
  TREF = 25 + 273.15                             # [K]            reference temperature for photosynthetic processes
  HARD = 46.39                                 # [kJ/mol]       activation energy of Rd
  CRD = 1000 * HARD / (R * TREF)                   # []             scaling factor in RD response to temperature

  HAGSTAR = 37.83                                 # [kJ/mol]       activation energy of Gamma_star
  CGSTAR = 1000 * HAGSTAR / (R * TREF)                # []             scaling factor in GSTAR response to temperature



  if (Type == 'C3'){                        # C3 species

    HAJ  = 49.88                                 # [kJ/mol]       activation energy of Jm (Kattge & Knorr 2007)
    HDJ  = 200                                   # [kJ/mol]       deactivation energy of Jm (Kattge & Knorr 2007)
    DELTASJ = (-0.75 * Tyear + 660) / 1000                # [kJ/mol/K]     entropy term for J  (Kattge and Knorr 2007)

    HAVCM = 71.51                                 # [kJ/mol]       activation energy of Vcm (Kattge and Knorr 2007)
    HDVC = 200                                   # [kJ/mol]       deactivation energy of Vcm (Kattge & Knorr 2007)
    DELTASVC= (-1.07 * Tyear + 668) / 1000                # [kJ/mol/K]     entropy term for Vcmax (Kattge and Knorr 2007)

    KCOP = 404.9                                 # [umol/mol]     Michaelis-Menten constant for CO2 at ref temp (Bernacchi et al 2001)
    HAKC = 79.43                                 # [kJ/mol]       activation energy of Kc (Bernacchi et al 2001)

    KOOP = 278.4                                 # [mmol/mol]     Michaelis-Menten constant for O2  at ref temp (Bernacchi et al 2001)
    HAKO = 36.38                                 # [kJ/mol]       activation energy of Ko (Bernacchi et al 2001)



  } else {
    # C4 species (values can be different as noted by von Caemmerer 2000)
    HAJ  = 77.9                                  # [kJ/mol]       activation energy of Jm  (Massad et al 2007)
    HDJ  = 191.9                              # [kJ/mol]       deactivation energy of Jm (Massad et al 2007)
    DELTASJ = 0.627                              # [kJ/mol/K]     entropy term for Jm (Massad et al 2007). No data available on acclimation to temperature.

    HAVCM = 67.29                           # [kJ/mol]       activation energy of Vcm (Massad et al 2007)
    HDVC = 144.57                                # [kJ/mol]       deactivation energy of Vcm (Massad et al 2007)
    DELTASVC = 0.472                                 # [kJ/mol/K]     entropy term for Vcm (Massad et al 2007). No data available on acclimation to temperature.

    HAVPM = 70.37                                 # [kJ/mol]       activation energy of Vpm  (Massad et al 2007)
    HDVP = 117.93                                # [kJ/mol]       deactivation energy of Vpm (Massad et al 2007)
    DELTASVP= 0.376                                 # [kJ/mol/K]     entropy term for Vpm (Massad et al 2007). No data available on acclimation to temperature.

    KCOP = 944                               # [umol/mol]     Michaelis-Menten constant for CO2 at ref temp (Chen et al 1994 Massad et al 2007)
    Q10KC = 2.1                                   # []             Q10 for temperature response of Kc (Chen et al 1994 Massad et al 2007)

    KOOP = 633                                  # [mmol/mol]     Michaelis-Menten constant for O2 at ref temp (Chen et al 1994 Massad et al 2007)
    Q10KO = 1.2                                   # []             Q10 for temperature response of Ko (Chen et al 1994 Massad et al 2007)

    KPOP = 82                                   # [umol/mol]     Michaelis-Menten constant of PEP carboxylase at ref temp (Chen et al 1994 Massad et al 2007)
    Q10KP = 2.1                                   # []             Q10 for temperature response of Kp (Chen et al 1994 Massad et al 2007)

  }


  #---------------------------------------------------------------------------------------------------------
  ## Corrections for effects of temperature and non-stomatal limitations
  dum1 = R / 1000 * Temp                                  # [kJ/mol]
  dum2 = R / 1000 * TREF                               # [kJ/mol]

  Rd = Rdopt * exp(CRD - HARD / dum1)                  # [umol/m2/s]    mitochondrial respiration rates adjusted for temperature (Bernacchi et al. 2001)
  SCO = SCOOP / exp(CGSTAR - HAGSTAR / dum1)            # []             Rubisco specificity for CO2 adjusted for temperature (Bernacchi et al. 2001)

  Jmax = Jmo * exp(HAJ *(Temp - TREF) /(TREF * dum1))
  Jmax = Jmax *(1 + exp((TREF * DELTASJ - HDJ) / dum2))
  Jmax = Jmax /(1 + exp((Temp * DELTASJ - HDJ) / dum1))     # [umol e-/m2/s] max electron transport rate at leaf temperature (Kattge and Knorr 2007 Massad et al. 2007)

  Vcmax = Vcmax25  * exp(HAVCM *(Temp - TREF) / (TREF * dum1))
  Vcmax = Vcmax * (1 + exp((TREF * DELTASVC - HDVC) / dum2))
  Vcmax = Vcmax / (1 + exp((Temp *DELTASVC - HDVC) / dum1 ))    # [umol/m2/s]    max carboxylation rate at leaf temperature (Kattge and Knorr 2007 Massad et al. 2007)

  if (Type == 'C3'){                                             # C3 species

    CKC = 1000 * HAKC / (R * TREF)                     # []             scaling factor in KC response to temperature
    Kc = KCOP * exp(CKC - HAKC / dum1) * 1e-11 * p     # [bar]          Michaelis constant of carboxylation adjusted for temperature (Bernacchi et al. 2001)

    CKO = 1000 * HAKO / (R * TREF)                     # []             scaling factor in KO response to temperature
    Ko = KOOP * exp(CKO - HAKO / dum1) * 1e-8 * p      # [bar]          Michaelis constant of oxygenation adjusted for temperature (Bernacchi et al. 2001)

  } else {

    # C4 species
    Vpmax  = Vpmo  * exp(HAVPM *(Temp - TREF) / (TREF * dum1))
    Vpmax  = Vpmax * (1 + exp((TREF * DELTASVP - HDVP) / dum2))
    Vpmax  = Vpmax / (1 + exp((Temp * DELTASVP - HDVP) / dum1))# [umol/m2/s]    max carboxylation rate at leaf temperature (Massad et al. 2007)

    Kc = KCOP * Q10KC ^ ((Temp - TREF)/10) * 1e-11 * p    # [bar]          Michaelis constant of carboxylation temperature corrected (Chen et al 1994 Massad et al 2007)

    Ko = KOOP * Q10KO  ^ ((Temp - TREF)/10) * 1e-8 * p     # [bar]          Michaelis constant of oxygenation  temperature corrected (Chen et al 1994 Massad et al 2007)

    Kp = KPOP * Q10KP ^ ((Temp - TREF ) /10) * 1e-11 * p    # [bar]          Michaelis constant of PEP carboxyl temperature corrected (Chen et al 1994 Massad et al 2007)

  }

  #---------------------------------------------------------------------------------------------------------
  ## Define electron transport and fluorescence parameters
  kf        = 3.e7                                    # [s-1]         rate constant for fluorescence
  kD        = 1.e8                                   # [s-1]         rate constant for thermal deactivation at Fm
  kd        = 1.95e8                                    # [s-1]         rate constant of energy dissipation in closed RCs (for theta=0.7 under un-stressed conditions)
  po0max    = 0.88                                     # [mol e-/E]    maximum PSII quantum yield, dark-acclimated in the absence of stress (Pfundel 1998)
  kPSII     = (kD + kf) * po0max/(1 - po0max)             # [s-1]         rate constant for photochemisty (Genty et al. 1989)
  fo0       = kf /(kf + kPSII + kD)                        # [E/E]         reference dark-adapted PSII fluorescence yield under un-stressed conditions

  kps       = kPSII * qLs                              # [s-1]         rate constant for photochemisty under stressed conditions (Porcar-Castell 2011)
  kNPQs     = NPQs * (kf + kD)                           # [s-1]         rate constant of sustained thermal dissipation (Porcar-Castell 2011)
  kds       = kd * qLs
  kDs       = kD + kNPQs
  Jms       = Jmax * qLs                               # [umol e-/m2/s] potential e-transport rate reduced for PSII photodamage
  po0       = kps / (kps + kf + kDs)                       # [mol e-/E]    maximum PSII quantum yield, dark-acclimated in the presence of stress
  THETA     = (kps - kds)/(kps + kf + kDs)                  # []            convexity factor in J response to PAR

  #---------------------------------------------------------------------------------------------------------
  ## Calculation of electron transport rate
  Q2     = beta * Q * po0
  J      = (Q2 + Jms - sqrt((Q2 + Jms) ^2 - 4 * THETA * Q2 *Jms)) / (2 * THETA) # [umol e-/m2/s]    electron transport rate under light-limiting conditions

  #---------------------------------------------------------------------------------------------------------
  ## Calculation of net photosynthesis
  if (Type == 'C3'){
    minCi = 0.3
  } else {
    minCi = 0.1
  }

  ###
  Ci = get.BallBerry(Cs, RH, A=NULL, BallBerrySlope=8, BallBerry0, minCi)

  if (Type == 'C3'){                                            # C3 species, based on Farquhar model (Farquhar et al. 1980)

    GSTAR = 0.5 * O / SCO                             # [bar]             CO2 compensation point in the absence of mitochondrial respiration
    Cc = Ci[[2]]                                        # [bar]             CO2 concentration at carboxylation sites (neglecting mesophyll resistance)

    Wc = Vcmax * Cc  / (Cc + Kc  * (1 + O / Ko))     # [umol/m2/s]       RuBP-limited carboxylation
    Wj = J  * Cc / (4.5 * Cc + 10.5 * GSTAR)            # [umol/m2/s]       electr transp-limited carboxyl

    W = pmin(Wc,Wj)                                # [umol/m2/s]       carboxylation rate
    Ag = (1 - GSTAR / Cc) * W                       # [umol/m2/s]       gross photosynthesis rate
    A = Ag - Rd                                   # [umol/m2/s]       net photosynthesis rate
    Ja = J * W / Wj                                 # [umole-/m2/s]     actual linear electron transport rate

  } else {  # C4 species, based on von Caemmerer model (von Caemmerer 2000)
    #Ci    =  max(9.9e-6*(p*1e-5),Cs.*(1-1.6./(m.*RH*stress)))
    # [bar]             intercellular CO2 concentration from Ball-Berry model (Ball et al. 1987)

    minCi = 0.1
    Ci = get.BallBerry(Cs, RH, A=NULL, BallBerrySlope, 0, minCi)
    Cm    =  Ci[[2]]                                     # [bar]             mesophyll CO2 concentration (neglecting mesophyll resistance)
    Rs    =  0.5 * Rd                               # [umol/m2/s]       bundle sheath mitochondrial respiration (von Caemmerer 2000)
    Rm    =  Rs                                     # [umol/m2/s]       mesophyll mitochondrial respiration
    gam   =  0.5 / SCO                               # []                half the reciprocal of Rubisco specificity for CO2

    Vpc   = Vpmax * Cm / (Cm + Kp)                     # [umol/m2/s]       PEP carboxylation rate under limiting CO2 (saturating PEP)
    Vp    = pmin(Vpc,Vpr)                            # [umol/m2/s]       PEP carboxylation rate

    # Complete model proposed by von Caemmerer (2000)
    dum1  =  alpha / 0.047                           # dummy variables, to reduce computation time
    dum2  =  Kc / Ko
    dum3  =  Vp - Rm + gbs * Cm
    dum4  =  Vcmax - Rd
    dum5  =  gbs * Kc *(1 + O / Ko)
    dum6  =  gam * Vcmax
    dum7  =  x * J / 2  - Rm + gbs * Cm
    dum8  =  (1 - x) * J /3
    dum9  =  dum8 - Rd
    dum10 =  dum8 + Rd * 7 / 3

    a =  1 - dum1 * dum2
    b  =  -(dum3 + dum4 + dum5 + dum1 * (dum6 + Rd * dum2))
    c  =  dum4 * dum3 - dum6 * gbs * O + Rd * dum5
    Ac =  (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)           # [umol/m2/s]       CO2-limited net photosynthesis

    a  =  1- 7 / 3 * gam * dum1
    b  =  -(dum7 + dum9 + gbs * gam * O *7 /3  + dum1 * gam  * dum10)
    c  =  dum7 * dum9 - gbs * gam * O * dum10
    Aj =  (-b - sqrt(b^2 - 4* a * c)) / (2 * a)           # [umol/m2/s]       light-limited net photosynthesis (assuming that an obligatory Q cycle operates)

    A     =  pmin(Ac,Aj)                             # [umol/m2/s]       net photosynthesis

    Ja    =  J                                      # [umole-/m2/s]     actual electron transport rate, CO2-limited

    if (any(A == Ac)) {

      ind <- which(A == Ac)

      x_ <- rep(x, length(A))

      a[ind] <- x_[ind] * (1 - x_[ind]) / 6 / A[ind]
      b[ind] <- (1 - x_[ind]) / 3 * (gbs[ind] / A[ind] * (Cm[ind] - Rm[ind] / gbs[ind] - gam[ind] * O) - 1 - alpha * gam[ind] / 0.047) - x_[ind] / 2 * (1 + Rd[ind] / A[ind])
      c[ind] <- (1 + Rd[ind] / A[ind]) * (Rm[ind] - gbs[ind] * Cm[ind] - 7 * gbs[ind] * gam[ind] * O / 3) + (Rd[ind] + A[ind]) * (1 - 7 * alpha * gam[ind] / 3 / 0.047)
      Ja[ind] <- (-b[ind] + sqrt(b[ind]^2 - 4 * a[ind] * c[ind])) / (2 * a[ind])
    }

  }


  #---------------------------------------------------------------------------------------------------------
  ## Calculation of PSII quantum yield and fluorescence
  ps     = Ja / (beta * Q)                            # [mol e-/E]    PSII photochemical quantum yield
  fs   =  get.MD12(ps,Ja,Jms,kps,kf,kds,kDs)            # [E/E]         PSII fluorescence yield
  eta    = fs / fo0                                   # []            scaled PSII fluorescence yield

  ## JP add
  rhoa =  subset(constants,constant == 'rhoa')[[2]]           # [kg m-3]       specific mass of air
  Mair =  subset(constants,constant == 'Mair')[[2]]          # [g mol-1]      molecular mass of dry air
  #rcw         = 0.625*(Cs-Ci)./A *rhoa/Mair*1E3    * 1e6 ./ p .* 1E5


  ppm2bar <- 1e-6 * (p * 1E-5)
  gs <- 1.6 * A * ppm2bar / (Cs - Ci[[2]])
  rcw <- (rhoa / (Mair * 1E-3)) / gs
  rcw[A <= 0 & rcw != 0] <- 0.625 * 1.0E6


  ## convert back to ppm
  Ci = Ci[[2]] * 1e6 / p * 1e5

  ##
  #
  # Output:
  # A         # [umol/m2/s]           net assimilation rate of the leaves
  # Ci        # [umol/mol]            CO2 concentration in intercellular spaces (assumed to be the same as at carboxylation sites in C3 species)
  # eta       # []                    amplification factor to be applied to PSII fluorescence yield spectrum
  #                                   relative to the dark-adapted, un-stressed yield calculated with either Fluspect or FluorMODleaf



  if (get.plots ==  T){
    check.rate <- data.frame(profile=c(1:length(A)),A=A,Gs=gs,Vcmax=Vcmax,Jmax=Jmax)


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


  biochem_out<-list()

  biochem_out$A = A
  biochem_out$Ci = Ci
  biochem_out$ps = ps
  biochem_out$eta = eta
  biochem_out$fs = fs
  biochem_out$rcw = rcw
  biochem_out$qE = rcw * NaN # dummy output, to be consistent with SCOPE
  biochem_out$Kn = NPQs + 0 * rcw #
  biochem_out$Phi_N  = kNPQs / (kNPQs + kD + kf + kps) + 0 * rcw
  biochem_out$Ja = Ja

  return(biochem_out)

}
