#' @title getRTMt.sb
#' \code{get.RTMt.sb} Calculates total outgoing radiation in hemispherical
#' direction and total absorbed radiation per leaf and soil component.
#' Radiation is integrated over the whole thermal spectrum with Stefan-Boltzman's equation. This function is a simplified version of
#' 'getRTMt.planck', and is less time consuming since it does not do the calculation for each wavelength separately.
#'
#' @param data.rad a large number of radiative fluxes: spectrally distributed and integrated, and canopy radiative transfer coefficients
#' @param data.soil soil properties
#' @param data.leafbio leaf  properties
#' @param data.canopy canopy properties (such as LAI and height)
#' @param data.leafopt leaf optical properties
#' @param data.gap probabilities of direct light penetration and viewing
#' @param Tcu Temperature of sunlit leaves    (oC), (13x36x60)
#' @param Tch Temperature of shaded leaves    (oC), (13x36x60)
#' @param Tsu  Temperature of sunlit soil      (oC), (1)
#' @param Tsh  Temperature of shaded soil      (oC), (1)
#' @param obsdir logical, TRUE if directional (observation direction) radiative fluxes should be calculated.
#' @param data.spectral information about wavelengths and resolutions
#' @param data.opts simulation options. Here, the options need for the RT model
#' @param get.plots is true plot the intermediate plots
#' @description
#' date: 5  Nov 2007
#' update:
#'  - 13 Nov 2007
#'  - 16 Nov 2007 CvdT    improved calculation of net radiation
#'  - 27 Mar 2008 JT      added directional calculation of radiation
#'  - 24 Apr 2008 JT      Introduced dx as thickness of layer (see parameters)
#'  - 31 Oct 2008 JT      introduced optional directional calculation
#'  - 31 Oct 2008 JT      changed initialisation of F1 and F2 -> zeros
#'  - 07 Nov 2008 CvdT    changed layout
#'  - 16 Mar 2009 CvdT    removed Tbright calculation
#'  - Feb 2013 WV      introduces structures for version 1.40
#'  - 04 Dec 2019 CvdT    adapted for SCOPE-lite
#'  - 17 Mar 2020 CvdT    mSCOPE representation
#'
#' @author 	 Wout Verhoef, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @return a large number of radiative fluxes: spectrally distributed and integrated, and canopy radiative transfer coefficients.  Here, thermal fluxes are added
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data.rad <- get.RTMt.sb(data.rad, data.soil, data.leafbio, data.canopy,
#'                          data.leafopt, data.gap, Tcu, Tch, Tsu, Tsh,
#'                          obsdir = TRUE, data.spectral, data.opts,
#'                          get.plots = FALSE)
#' }
get.RTMt.sb<-function(data.rad,data.soil,data.leafbio,data.canopy,data.leafopt,
                      data.gap,Tcu,Tch,Tsu,Tsh,obsdir,data.spectral,data.opts=data.opts,get.plots){

  ######################################################################
  # Table of contents of the function
  #   0       preparations
  #       0.0     globals
  #       0.1     initialisations
  #       0.2     parameters
  #       0.3     geometric factors of Observer
  #       0.4     geometric factors associated with extinction and scattering
  #       0.5     geometric factors to be used later with rho and tau
  #       0.6     fo for all leaf angle/azumith classes
  #   1       calculation of upward and downward fluxes
  #   2       total net fluxes
  ######################################################################

  constants <- SCOPEinR::constants
  
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

  # 0 Preparations

  nl <- data.canopy[['nlayers']]
  lidf<-data.canopy[['lidf']]
  Ps <- data.gap[['Ps']]
  #Leaf/needle reflection
  rho <-data.leafbio[['rho_thermal']]    # [1]
  #Leaf/needle transmission
  tau <-data.leafbio[['tau_thermal']]      # [1]
  #Soil reflectance
  rs <-data.soil[['rs_thermal']]        # [1]

  rho_spec <- data.leafopt[['refl']] #S pectral reflectance of the leaf, 400 to 2400 nm
  tau_spec <- data.leafopt[['tran']]  # Spectral transmittance of the leaf, 400 to 2400 nm
  rsoil_spec <- data.soil[['rfl.soil']]

  rsoil <- data.soil[['rfl.soil']] # [nwl,nsoils] soil reflectance spectra
  epsc <- 1-rho-tau;            # [nwl]               Emissivity vegetation
  epss <- 1-rs                  # [nwl]               Emissivity soil

  LAI <- data.canopy[['LAI']] # leaf area index
  lidf <- data.canopy[['lidf']] # leaf area index
  dx <- 1/nl
  iLAI <-LAI * dx

  ## get last column for Xdd,Xsd, R_dd, R_sd,rho_dd and tau_dd
  Xdd <- data.rad[['Xdd']]
  Xdd<-Xdd[,ncol(Xdd)]
  Xsd <- data.rad[['Xsd']]
  Xsd<-Xsd[,ncol(Xsd)]

  Xss <- rep(data.rad[['Xss']][1], each = data.canopy[['nlayers']])

  R_dd  <- data.rad[['R_dd']]
  R_dd<-R_dd[,ncol(R_dd)]

  R_sd <- data.rad[['R_sd']]
  R_sd<-R_sd[,ncol(R_sd)]


  rho_dd <-  data.rad[['rho_dd']]
  rho_dd <-  rho_dd[,ncol(rho_dd)]

  tau_dd <-  data.rad[['tau_dd']]
  tau_dd <-  tau_dd[,ncol(tau_dd)]

  ## 1. calculation of upward and downward fluxes pag 305

  # 1.1 radiance by components

  # get Radiance by sunlit leaves
  Hcsu3  <-  epsc * get.Stefan_Boltzmann(Tcu)

  # get Radiance by shaded leaves
  Tch_s <-get.Stefan_Boltzmann(Tch)
  Hcsh  <-  epsc * Tch_s#matrix(get.Stefan_Boltzmann(Tch),ncol=1)

  # get Radiance by sunlit soil
  Hssu  <-  epss * get.Stefan_Boltzmann(Tsu)

  # get Radiance by shaded soil
  Hssh  <-  epss * get.Stefan_Boltzmann(Tsh)


  if (is.matrix(Hcsu3)== T) {
    if (dim(Hcsu3)[2] > 1){
      v1 <- rep(1/ncol(Hcsu3), ncol(Hcsu3))   # vector for computing the mean
      Hcsu2 <- Hcsu3

      lidf_m <- t(matrix(rep(lidf,nl), ncol = 59, nrow=13))
      ## get dim (2001 x 13)
      Hcsu2_lidf_v1 <- v1 * t(Hcsu2) %*% lidf_m
      # compute column means for each level
      Hcsu <-colMeans(Hcsu2_lidf_v1)
      Hcsu <-matrix(rep(Hcsu), ncol = 1, nrow=13)
      #### check it
    }

  } else if (is.matrix(Hcsu3) == F ) {
    Hcsu <- Hcsu3
  }

  Hcsu_<-Hcsu * Ps[1:nl]
  Hcsh_<-Hcsh * (1 - Ps[1:nl])
  Hc <- matrix(Hcsu_ + Hcsh_,nrow=1)          # hemispherical emittance by leaf layer
  Hs <-Hssu * Ps[nl+1] + Hssh * (1-Ps[nl+1])  # hemispherical emittance by soil surface



  #1.3 Diffuse radiation

  # direct, up and down diff. rad.
  U <- rep(0,nl+1)
  U[nl+1] <-   Hs
  Es_ <- rep(0, nl+1)  # initialize Es_ as a vector of zeros
  Emin <- rep(0, nl+1) # initialize Emin as a vector of zeros
  Eplu <- rep(0, nl+1)   # initialize Eplu as a vector of zeros

  Y <-rep(0, nl) ####
  for (j in nl:1) {

    Y[j] <- (rho_dd[j]*U[j+1] + Hc[j]*iLAI) / (1 - rho_dd[j]*R_dd[j+1])
    U[j] <- tau_dd[j]*(R_dd[j+1]*Y[j] + U[j+1]) + Hc[j]*iLAI
  }

  for (j in 1:nl) {
    Es_[j+1] <- Xss[j] * Es_[j]
    Emin[j+1] <- Xsd[j] * Es_[j] + Xdd[j] * Emin[j] + Y[j]
    Eplu[j] <- R_sd[j] * Es_[j] + R_dd[j] * Emin[j] + U[j]
  }



  Eplu[nl+1] <- R_sd[nl] * Es_[nl] + R_dd[nl] * Emin[nl] + Hs
  Eoutte <- Eplu[1]




  if ( missing(obsdir) == FALSE) {
    K <- data.gap[['K']]


    vb <- (data.rad[['vb']])[1,nl]


    vf <-  (data.rad[['vf']])[1,nl]
    # piLov is a scalar value that represents the directional radiation emitted
    # by vegetation and scattered radiation by vegetation for diffuse incidence

    piLov_term1 <- iLAI * (K * c(Hcsh) * (data.gap[['Po']][1:nl] - data.gap[['Pso']][1:nl]))
    piLov_term2<- iLAI * (K * Hcsu * data.gap[['Pso']][1:nl])
    piLov_term3<-  iLAI *  ((vb * Emin[1:nl] + vf * Eplu[1:nl]) * data.gap[['Po']][1:nl])
    piLov <- sum(piLov_term1 + piLov_term2 + piLov_term3)

    #  piLos is a scalar that represent the directional emitted radiation by soil.
    piLos <- Hssh * (data.gap[['Po']][nl+1] - data.gap[['Pso']][nl+1]) + Hssu * data.gap[['Pso']][nl+1]

    #piLot is the total directional radiation emitted and scattered by the surface (vegetation + soil)
    piLot <- piLov + piLos


    sigmaSB <- subset(SCOPEinR::constants, constant == 'sigmaSB')[[2]]
    Tbr <- (piLot / sigmaSB) ^ 0.25
    ##S pectral directional emissivity of the surface (vegetation + soil)
    data.rad[['Lote']] <- piLot / pi

    #Lot_ is the spectral directional radiance of the surface (vegetation + soil) in each wavelength band
    data.rad[['Lot_']]  <- get.Planck(wl=data.spectral$wlS, Tb=Tbr)  # Note that this is the directional blackbody radiance!

    Tbr2 <- (Eoutte / sigmaSB) ^ 0.25
    #Eoutte_ is the upward radiation flux density at the top of the atmosphere and
    data.rad[['Eoutte_']] <- get.Planck(wl=data.spectral$wlS, Tb=Tbr2)
  }

  ## 2. total net fluxes
  #net radiation per component, in W m-2 (leaf or soil surface)



  if (is.matrix(Hcsu3)== T) {

    if (dim(Hcsu3)[2] > 1){
      Rnuc <- array(0, dim = c(nrow(Hcsu3), ncol(Hcsu3), nl))

      for (j in 1:nl) {
        Rnuc[, , j] <- (Emin[j] + Eplu[j+1] - 2 * Hcsu3[,,j])    # sunlit leaf
      }

    }
  } else {

    Rnuc <- (Emin[1:(nl)] + Eplu[2:(nl+1)] - 2 * Hcsu)
  }


  #Rnhc <- (Emin[1:(nl)] + Eplu[2:(nl+1)] - 2 * Hcsh)
  Rnhc <- (Emin[-length(Emin)] + Eplu[-1] - 2 * Hcsh)

  Rnus <- (Emin[(nl+1)] - Hssu) # sunlit soil
  Rnhs <- (Emin[(nl+1)] - Hssh) # shaded soil


  ## 3. Write the output to the rad structure

  data.rad[['Emint']] <- Emin
  data.rad[['Eplut']] <- Eplu
  data.rad[['Eoutte']] <- Eoutte
  data.rad[['Rnuct']] <- Rnuc
  data.rad[['Rnhct']] <- Rnhc
  data.rad[['Rnust']] <- Rnus
  data.rad[['Rnhst']] <- Rnhs


  #########################################################################
  ###### Get some plots
  #########################################################################

  if (get.plots ==  T){


    #########################################################################
    ###### Get the hemispherical emittance by leaf layers plots
    #########################################################################

    check.rad<- data.frame(wave=data.spectral[['wlS']],Eoutte_= data.rad[['Eoutte_']],Lot_=data.rad[['Lot_']])


    p.rad <- ggplot(data = check.rad , aes(x = wave)) +
      labs(y= "radiation flux ", x = "") +
      geom_line(aes(y = Lot_, color = "Lot_"), linewidth = 0.5) +
      geom_line(aes(y = Eoutte_, color = "Eoutte_"), linewidth = 0.5) +

      scale_color_manual(name = "r",
                         values = c("Lot_" = "orange", "Eoutte_" ='brown' )) +
      theme_bw() + #xlim(400,2500) +
      theme(legend.position = "right") +
      guides(color = guide_legend(title = NULL))

    print(p.rad)


  }


  LRT<- list(data.rad =data.rad)

  return(LRT)


}











