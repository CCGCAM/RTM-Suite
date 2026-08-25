#' get.RTMt.planck
#' \code{get.RTMt.planck}analogue to get.RTMt.sb, this function calculates total outgoing radiation in hemispherical
#' direction and total absorbed radiation per leaf and soil component.
#' Radiation is integrated over the whole thermal spectrum with Stefan-Boltzman's equation. This function is a simplified version of
#' 'get.RTMt.planck', and is less time consuming since it does not do the calculation for each wavelength separately.
#'
#' @param data.spectral information about wavelengths and resolutions
#' @param data.rad a large number of radiative fluxes: spectrally distributed and integrated, and canopy radiative transfer coefficients.
#' @param data.soil  soil properties
#' @param data.leafbio leaf properties (Cab....)
#' @param data.leafopt leaf optical properties
#' @param data.canopy canopy properties (such as LAI and height)
#' @param data.gap probabilities of direct light penetration and viewing
#' @param Tcu Temperature of sunlit leaves    (oC), (13x36x60)
#' @param Tch Temperature of shaded leaves    (oC), (13x36x60)
#' @param Tsu  Temperature of sunlit soil      (oC), (1)
#' @param Tsh  Temperature of shaded soil      (oC), (1)
#' @param get.plots is true plot the intermediate plots
#'
#' @return a large number of radiative fluxes: spectrally distributed and integrated, and canopy radiative transfer coefficients., Here
#' thermal fluxes are added.
#' @export
#'
#' @description
#'
#' Date:   05  Nov 2007
#' Update:
#'  - 13 Nov 2007
#'  - 16 Nov 2007 CvdT: improved calculation of net radiation.
#'  - 27 Mar 2008 JT: added directional calculation of radiation.
#'  - 24 Apr 2008 JT: introduced dx as thickness of layer (see parameters).
#'  - 31 Oct 2008 JT: introduced optional directional calculation.
#'  - 31 Oct 2008 JT: changed initialisation of F1 and F2 -> zeros.
#'  - 07 Nov 2008 CvdT: changed layout.
#'  - 16 Mar 2009 CvdT: removed Tbright calculation.
#'  - Feb 2013 WV: introduces structures for version 1.40.
#'  - 04 Dec 2019 CvdT: adapted for SCOPE-lite.
#'  - 17 Mar 2020 CvdT: mSCOPE representation.

#' @author Wout Verhoef and Christiaan van der Tol  (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#' @examples
#' \dontrun{
#' data.rad <- get.RTMt.planck(data.spectral, data.rad, data.soil, data.leafbio,
#'                              data.leafopt, data.canopy, data.gap,
#'                              Tcu, Tch, Tsu, Tsh, get.plots = FALSE)
#' }
get.RTMt.planck <- function(data.spectral, data.rad, data.soil, data.leafbio,data.leafopt, data.canopy, data.gap,
                            Tcu,Tch,Tsu,Tsh, get.plots =T){

  if (missing(get.plots)){
    get.plots = FALSE
  }

  ###############################################################################
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
  ###############################################################################

  ## 0.2     parameters

  IT  <- data.spectral[['IwlT']]
  wlt <- data.spectral[['wlT']]

  deg2rad <- subset(SCOPEinR::constants,constant == 'deg2rad')[[2]] #

  nl <- data.canopy[['nlayers']]

  lidf<-data.canopy[['lidf']]
  Ps <- data.gap[['Ps']]
  adding_rho_thermal <- matrix(data.leafbio[['rho_thermal']], nrow = nl, ncol = 161)
  # Leaf/needle reflection
  rho <-data.leafopt[['refl']]
  rho <- cbind(rho, adding_rho_thermal)[,IT]

  #Leaf/needle transmission
  adding_tau_thermal <- matrix(data.leafbio[['tau_thermal']], nrow = nl, ncol = 161)
  tau <-data.leafopt[['tran']]
  tau <- cbind(tau, adding_tau_thermal)[,IT]
  # soil reflectance spectra ([nwl,nsoils])
  rsoil <- data.soil[['rfl.soil']]
  rsoil <- c(rsoil, rep(data.leafbio[['rho_thermal']],161))
 # Emissivity vegetation
  epsc <-1-rho-tau
  #Emissivity soil
  epss <- 1-rsoil

  LAI <- data.canopy[['LAI']] # leaf area index
  lidf <- data.canopy[['lidf']] # leaf area index
  dx <- 1/nl
  iLAI <-LAI * dx

  ## get last column for Xdd,Xsd, R_dd, R_sd,rho_dd and tau_dd
  Xdd <- data.rad[['Xdd']][,IT]
  Xdd<-Xdd[,ncol(Xdd)]
  Xsd <- data.rad[['Xsd']][,IT]
  Xsd<-Xsd[,ncol(Xsd)]

  Xss <- rep(data.rad[['Xss']][1], each = data.canopy[['nlayers']])

  R_dd  <- data.rad[['R_dd']][,IT]
  R_dd<-R_dd[,ncol(R_dd)]

  R_sd <- data.rad[['R_sd']][,IT]
  R_sd<-R_sd[,ncol(R_sd)]


  rho_dd <-  data.rad[['rho_dd']][,IT]
  rho_dd <-  rho_dd[,ncol(rho_dd)]

  tau_dd <-  data.rad[['tau_dd']][,IT]
  tau_dd <-  tau_dd[,ncol(tau_dd)]

  ## 0.2  initialization of output variables

  piLot_ <- Eoutte_ <- numeric(length(IT))
  Emin_ <- Eplu_ <- matrix(0, nrow = nl + 1, ncol = length(IT))


  ## 1. calculation of upward and downward fluxes pag 305
  Hcsu3 <- Hcsh <- Hssu <- Hssh <-Hc <- Hs <-c()

  ### loop for estimating radiance at several direction

  for (i in c(1:length(IT))) {
    # 1.1 radiance by components
    # get Radiance by sunlit leaves
    Hcsu3 <- pi * get.Planck(wlt[i],Tcu+273.15,epsc[1,i])

    # get Radiance by shaded leaves
    Hcsh  <- pi * get.Planck(wlt[i],Tch+273.15,epsc[1,i])

    # get Radiance by sunlit soil

    Hssu  <- pi * get.Planck(wlt[i],Tsu+273.15,epss[i])

    # get Radiance by shaded soil

    Hssh  <- pi * get.Planck(wlt[i],Tsh+273.15,epss[i])

    # 1.2 radiance by leaf layers Hv and by soil Hs (modified by JAK 2015-01)

    if (is.vector(Hcsu3) == F ) {

      if (is.matrix(Hcsu3) == T){
        Hcsu3<-matrix(Hcsu3,ncol=1)
        v1 <- rep(1/ncol(Hcsu3), ncol(Hcsu3))   # vector for computing the mean
        Hcsu2 = matrix(Hcsu3,nrow=dim(Hcsu3)[2], ncol=dim(Hcsu3)[1])   # create a block matrix from the 3D array

        lidf_m <- matrix(lidf, ncol =13, nrow=1)

        Hcsu2_lidf_v1 <- v1 * t(Hcsu2) %*% lidf_m
        Hcsu = Hcsu2_lidf_v1
      } else {
        # Hcsu3 here is [leaf angle (13), azimuth (36), layer (nl)] --
        # weighted mean over leaf angle (weight = lidf, the leaf
        # inclination distribution) and uniform mean over azimuth,
        # returning one value per layer.
        nlazi_here <- dim(Hcsu3)[2]
        Hcsu <- sapply(seq_len(dim(Hcsu3)[3]), function(k) {
          sum(Hcsu3[, , k] * lidf) / nlazi_here
        })
      }

    } else {
      Hcsu <- Hcsu3
    }


    # hemispherical emittance by leaf layers
    Hc <- Hcsu * Ps[1:nl] + Hcsh * (1 - Ps[1:nl])
    # hemispherical emittance by soil surface

    Hs <- Hssu * Ps[nl+1] + Hssh * (1 - Ps[nl+1])

    # 1.3 Diffuse radiation

    # direct, up and down diff. rad.
    U <- rep(0,nl+1)
    U[nl+1] <-   Hs
    Es_ <- rep(0, nl+1)  # initialize Es_ as a vector of zeros
    Emin <- rep(0, nl+1) # initialize Emin as a vector of zeros
    Eplu <- rep(0, nl+1)   # initialize Eplu as a vector of zeros

    Y <-rep(0, nl)
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

    # downwelling diffuse radiance per layer
    Emin_[,i] <-Emin
    # upwelling   diffuse radiance
    Eplu_[,i] <-Eplu

    Eoutte_[i] <- Eplu[1]
    # 1.4 Directional radiation and brightness temperature


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
    piLot_[i] <- piLov + piLos


  } # end for

  # Lot_ is the spectral directional radiance of the surface
  Lot_ <- piLot_ / pi



  ## 3. Write the output to the rad structure

  data.rad[['Lot_']] <- rep(0, length(data.spectral$wlS))
  data.rad[['Eoutte_']]<- rep(0, length(data.spectral$wlS))


  data.rad[['Lot_']][IT] <- Lot_
  #emitted  diffuse radiance at top
  data.rad[['Eoutte_']][IT] <- Eoutte_
  data.rad[['Eplut_']] <- Eplu_
  data.rad[['Emint_']] <- Emin_

  #########################################################################
  ###### Get some plots
  #########################################################################


  if (get.plots ==  T){


    #########################################################################
    ###### Get the hemispherical emittance by leaf layers plots
    #########################################################################

    check.profiles <- data.frame(layer=c(1:nl),Hc=Hc)
 # hemispherical emittance by leaf layers
    p.hc <- ggplot(data = check.profiles, aes(x = layer)) +
      labs(x= "n layers", x = "") +
      geom_line(aes(y = Hc, color = "Hc"), linewidth = 0.5) +


      scale_color_manual(name = "r",
                         values = c("Hc" = "orange")) +
      theme_bw() + #xlim(400,2500) +
      theme(legend.position = "right") +
      guides(color = guide_legend(title = NULL))

    print(p.hc)

    #########################################################################
    ###### Get the hemispherical emittance by leaf layers plots
    #########################################################################


    check.s <- data.frame(wave=wlt,Lot_=Lot_,piLot_=piLot_,
                          Eoutte_=Eoutte_,Eplu_=Eplu_[1,],Emin_=Emin_[2,])

    p.lot <- ggplot(data = check.s, aes(x = wave)) +
      labs(x= "long wavelengths", y = "spectral directional radiance") +
      geom_line(aes(y = Lot_, color = "Lot_"), linewidth = 0.5) +
      geom_line(aes(y = piLot_, color = "piLot_"), linewidth = 0.5) +
      geom_line(aes(y = Eoutte_, color = "Eoutte_"), linewidth = 0.5) +
      geom_line(aes(y = Eplu_, color = "Eplu_"), linewidth = 0.5) +
      geom_line(aes(y = Emin_, color = "Emin_"), linewidth = 0.5) +
      scale_color_manual(name = "r",
                         values = c("piLot_" = "navyblue",
                                    'Eoutte_' = 'black',
                                   'Lot_' = 'red',
                                   'Eplu_' = 'orange',
                                   'Emin_'= 'blue')) +
      theme_bw() + #xlim(400,2500) +
      theme(legend.position = "right") +
      guides(color = guide_legend(title = NULL))

    print(p.lot)

  }





  ######################################################################
  ##### Comments: by CvdT, 11 December 2015.
  #####
  ##### We subtract Emin(1), because ALL incident (thermal) radiation from Modtran
  ##### has been taken care of in RTMo. Not ideal but otherwise radiation budget will not close!
  ######################################################################

  return(data.rad)

}

