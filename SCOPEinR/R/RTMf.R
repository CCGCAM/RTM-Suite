#' get.RTMf
#' \code{get.RTMf} Calculates the spectrum of fluorescent radiance in the observer's direction and also
#' the TOC spectral hemispherical upward Fs flux.
#'
#' @param data.spectral information about wavelengths and resolutions
#' @param data.rad a large number of radiative fluxes: spectrally distributed and integrated, and canopy radiative transfer coefficients.
#' @param data.soil  soil properties
#' @param data.leafopt leaf optical properties
#' @param data.canopy canopy properties (such as LAI and height)
#' @param data.gap probabilities of direct light penetration and viewing
#' @param data.angles viewing and observation angles
#' @param data.etau  relative fluorescence emission efficiency for sunlit leaves
#
#' @param data.etah relative fluorescence emission efficiency for shaded leaves
#' @param get.plots logical. If \code{TRUE}, intermediate diagnostic plots are produced. Default \code{TRUE}.
#'
#' @return a rad object with a large number of radiative fluxes: spectrally distributed
#' and integrated, and canopy radiative transfer coefficients.
#' @export
#' @details
#' Date:  12 Dec 2007.
#' Update history:
#'  - 26 Aug 2008 CvdT        small correction to matrices
#'  - 07 Nov 2008 CvdT        changed layout
#'  - 19 Mar 2009 CvdT        major corrections: lines 95-96,101-107, and 119-120.
#'  - 07 Apr 2009 WV & CvdT   major correction: lines 89-90, azimuth
#'  dependence was not there in previous verions (implicit assumption of  azimuth(solar-viewing) = 0). This has been corrected
#'
#'  - May-June 2012 WV & CvdT Add calculation of hemispherical Fs fluxes
#'  - Jan-Feb 2013 WV         Inputs and outputs via structures for SCOPE Version 1.40
#'  - Aug-Oct 2016 PY         Re-write the calculation of emitted SIF of each layer. It doesnt use loop at all. with the function bsxfun, the  calculation is much faster
#'
#'  - Oct 2017-Feb 2018 PY    Re-write the RTM of fluorescence
#'  - Jan 2020 CvdT           Modified to include 'lite' option, mSCOPE representation
#'  - 25 Jun 2020 PY          Po, Ps, Pso. fix the problem we have with the oblique angles above 80 degrees
#'
#' @author 	Wout Verhoef and Christiaan van der Tol  (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' get.RTMf(data.spectral, data.rad, data.soil, data.leafopt, data.canopy,
#'          data.gap, data.angles, data.etau, data.etah)
#' }
get.RTMf <- function(data.spectral, data.rad, data.soil, data.leafopt, data.canopy, data.gap,
                     data.angles, data.etau, data.etah,
                     get.plots =T){

  if (missing(get.plots)){
    get.plots = FALSE
  }

  # Table of contents of the function:
  #   0.preparations
  #   0.0     globals
  #   0.1     initialisations
  #   0.2     geometric quantities
  #   0.3     solar irradiance factor and ext. in obs dir for all leaf angle/azumith classes
  #   1.0     calculation of fluorescence flux in observation direction
  #
  ## 0.1 initialisations

  constants <- SCOPEinR::constants

  wlS <- data.spectral[['wlS']] # SCOPE wavelengths, make column vectors
  wlF <- seq(640, 850, by = 4) # Fluorescence wavelengths
  wlE <- seq(400, 750, by = 5) # Excitation wavelengths
  iwlfi <- match(wlE, wlS)
  iwlfo <- match(wlF, wlS)

  nf <- length(iwlfo)
  nl <- data.canopy[['nlayers']]
  LAI <- data.canopy[['LAI']]
  litab <- data.canopy[['litab']]
  lazitab <- data.canopy[['lazitab']]
  lidf <- data.canopy[['lidf']]

  nlazi <- length(lazitab) # azumith angle
  nlinc <- length(litab) # inclination
  nlori <- nlinc * nlazi # total number of leaf orientations
  layers <- 1:nl


  Ps <- data.gap[['Ps']]
  Po <- data.gap[['Po']]
  Pso <- data.gap[['Pso']]

  Qs <- Ps[1:(length(Ps)-1)]


  # speed-up the calculation by only using wavelength i and wavelength o part of the spectrum
  ## Esun_f : Direct solar in the fluorescence emission range
  Esunf_ <- data.rad[['Esun_']][iwlfi]

  ## Eminf_: downward diffuse radiation in the canopy in the fluorescence emission range
  Eminf_ <- t(data.rad[['Emin_']][,iwlfi])  # transpose into [nwlfo,nl] matrix

  ##Epluf_: upward diffuse radiation in the canopy (mW m-2 um-1) in the fluorescence emission range
  Epluf_ <- t(data.rad[['Eplu_']][,iwlfi])

  iLAI <- LAI/nl  # LAI of a layer [1]

  ### check it and all is correct

  Xdd <- data.rad[['Xdd']][,iwlfo]
  # rho_dd: diffuse-diffuse reflectance for the thin layers
  rho_dd <- data.rad[['rho_dd']][,iwlfo]
  R_dd <- data.rad[['R_dd']][,iwlfo]
  # tau_dd:  diffuse-diffuse transmittance for the thin layers
  tau_dd <- data.rad[['tau_dd']][,iwlfo]
  #vb: directional backscatter scattering coefficient for diffuse incidence
  vb <- data.rad[['vb']][,iwlfo]
  #vf: directional forward scattering coefficient for diffuse incidence
  vf <- data.rad[['vf']][,iwlfo]
  ### check it and all is correct
  #plot(Xdd[1,],tpye='l')

  # geometric quantities
  Mb <- data.leafopt[['Mb']]
  Mf <- data.leafopt[['Mf']]

  ### check it and all is errors in firest elemets
  #plot(Mf[5,,1],tpye='l')

  #plot(Mb[,2,1], type='l')
  # geometric factors
  deg2rad <- subset(SCOPEinR::constants,constant == 'deg2rad')[[2]]
  tto <- data.angles[['tto']]
  tts <- data.angles[['tts']]
  psi <- data.angles[['psi']]

  rs <- data.soil[['rfl.soil']][iwlfo]  # [nwlfo] soil reflectance
  cos_tto <- cos(tto*deg2rad)  # cos observation zenith angle
  sin_tto <- sin(tto*deg2rad)  # sin observation zenith angle
  cos_tts <- cos(tts*deg2rad)  # cos solar angle
  sin_tts <- sin(tts*deg2rad)  # sin solar angle
  cos_ttli <- cos(litab*deg2rad)  # cos leaf inclinaation angles
  sin_ttli <- sin(litab*deg2rad)  # sin leaf inclinaation angles
  cos_phils <- cos(lazitab*deg2rad)  # cos leaf azimuth angles rel. to sun azi
  cos_philo <- cos((lazitab-psi)*deg2rad)  # cos leaf azimuth angles rel. to viewing azi

  # geometric factors for all leaf angle/azumith classes

  cds <- matrix(0, nrow = nlinc, ncol = nlazi)  # Initialize cds matrix
  cdo <- matrix(0, nrow = nlinc, ncol = nlazi)  # Initialize cdo matrix

  for (i in 1:nlinc) {
    cds[i, ] <- cos_ttli[i] * rep(cos_tts,nlazi) + sin_ttli[i] * sin_tts * cos_phils  # [nlinc,nlazi]
    cdo[i, ] <- cos_ttli[i] * rep(cos_tto,nlazi) + sin_ttli[i] * sin_tto * cos_philo  # [nlinc,nlazi]
  }


  ### check it and all is correct
  #plot(cdo[,1],tpye='l')


  fs <- cds/cos_tts  # [nlinc,nlazi]
  absfs <- abs(fs)  # [nlinc,nlazi]
  fo <- cdo/cos_tto  # [nlinc,nlazi]
  absfo <- abs(fo)  # [nlinc,nlazi]
  fsfo <- fs*fo  # [nlinc,nlazi]
  absfsfo <- abs(fsfo)  # [nlinc,nlazi]
  foctl <- fo*(cos_ttli)# %*% matrix(1,1,36))  # [nlinc,nlazi]
  fsctl <- fs*(cos_ttli)# %*% matrix(1,1,36))  # [nlinc,nlazi]
  ctl2 <- matrix(cos_ttli^2,nrow=nlinc, ncol=nlazi)  # [nlinc,nlazi]



  #plot(ctl2[,1], type='l')
  # reshape all the variables with dimension of nlori

  # Reshape all the variables
  absfs <- matrix(absfs, nlori, 1, byrow=F)   # [nlori,1]
  absfo <- matrix(absfo, nlori, 1, byrow=F)  # [nlori,1]
  fsfo <- matrix(fsfo, nlori, 1, byrow=F)     # [nlori,1]
  absfsfo <- matrix(absfsfo, nlori, 1, byrow=F)  # [nlori,1]
  foctl <- matrix(foctl, nlori, 1)  # [nlori,1]
  fsctl <- matrix(fsctl, nlori, 1, byrow=F)   # [nlori,1]
  ctl2 <- matrix(ctl2, nlori, 1, byrow=F)    # [nlori,1]


  #1.0 calculation of fluorescence flux in observation direction

  #fluorescence matrices and efficiencies

  U <- matrix(0, nrow = nl+1, ncol = dim(data.leafopt[['Mb']])[1] )
  Fmin_ <- matrix(0, nrow = nl+1, ncol = dim(data.leafopt[['Mb']])[1] )
  Fplu_ <- matrix(0, nrow = nl+1, ncol = dim(data.leafopt[['Mb']])[1] )

  Mplu <- 0.5 * (data.leafopt[['Mb']] + data.leafopt[['Mf']]) # [nwlfo,nwlfi]
  Mmin <- 0.5 * (data.leafopt[['Mb']] - data.leafopt[['Mf']]) # [nwlfo,nwlfi]


  #plot(Mplu[10,,1],tpye='l')
  # in-products: we convert incoming radiation to a fluorescence spectrum using the matrices.


  MpluEmin <- matrix(0, nf, nl)
  MpluEplu <- matrix(0, nf, nl)
  MminEmin <- matrix(0, nf, nl)
  MminEplu <- matrix(0, nf, nl)
  MpluEsun <- matrix(0, nf, nl)
  MminEsun <- matrix(0, nf, nl)

  for (j in 1:nl) {

    ep <- subset(SCOPEinR::constants,constant == 'A')[[2]] * get.ephoton(wlF * 1E-9, SCOPEinR::constants)
    Emin_photo <-get.e2phot(wlE * 1E-9, Eminf_[, j], SCOPEinR::constants)
    #plot(Emin_photo,type = 'l',col='forestgreen')
    #aa <-ep * (Mplu[, , j] %*% matrix(Emin_photo, nrow = length(Emin_photo), ncol = 30, byrow = TRUE))

    MpluEmin[, j] <- ep * (Mplu[, , j] %*% get.e2phot(wlE * 1E-9, Eminf_[, j], SCOPEinR::constants))
    #plot(Mplu[,,j],type = 'l',col='forestgreen')
    MpluEplu[, j] <- ep * (Mplu[, , j] %*% get.e2phot(wlE * 1E-9, Epluf_[, j], SCOPEinR::constants))
    MminEmin[, j] <- ep * (Mmin[, , j] %*% get.e2phot(wlE * 1E-9, Eminf_[, j], SCOPEinR::constants))
    MminEplu[, j] <- ep * (Mmin[, , j] %*% get.e2phot(wlE * 1E-9, Epluf_[, j], SCOPEinR::constants))
    MpluEsun[, j] <- ep * (Mplu[, , j] %*% get.e2phot(wlE * 1E-9, Esunf_, SCOPEinR::constants))
    MminEsun[, j] <- ep * (Mmin[, , j] %*% get.e2phot(wlE * 1E-9, Esunf_, SCOPEinR::constants))
  }



  #plot(MpluEmin[2,],tpye='l')
  #plot(Mplu[,1,2],tpye='l')

  #df.plot <-MpluEsun

  #plot(NA,NA,type = 'l',col='forestgreen',xlab='', ylab = 'Mplu',ylim=c(0, max(df.plot)),xlim=c(640,850))
  #for (i in 1:nl){
  #  par(new=T)
  #  plot(wlF,df.plot[1:53,i],xlab='',ylab='Mplu',
  #       type = 'l',col='forestgreen',ylim=c(0, max(df.plot)),xlim=c(640,850))
  #}


  # Assuming your lidf, etau, nlori, nl, laz, and etah are defined
  laz <- 1/36

  etau <- matrix((data.etau), nrow=length(data.etau),ncol=T)
  # Repetir el vector 13 veces en filas y combinarlo en una matriz
  #etau_ <- do.call(rbind, replicate(13, data.bcu$eta, simplify = FALSE))

  if (dim(etau)[2] < 2) {
    etau <- array(etau, dim = c(nl, 13, 36))
    etau <- aperm(etau, perm = c(2, 3, 1))
  }

  etah_reshape <- matrix(data.etah, nrow = nl, ncol = nlori, byrow = F)
  etau_reshape <- matrix(etau, nrow = nlori, ncol = nl, byrow = F)
  lidf_laz <-lidf * laz
  lidf_laz <-rep(lidf_laz,36)
  #plot(lidf_laz,type='l',lwd=1)
  #plot(etau_reshape[1,],type='l',lwd=1)
  #plot(etah_reshape[,2],type='l',lwd=1)


  etau_lidf <- etau_reshape * lidf_laz
  etah_lidf <- t(etah_reshape) * lidf_laz

  #plot(etah_lidf[1,],type='l',lwd=1)
  #plot(etau_lidf[1,],type='l',lwd=1)

  absfsfo_nl<- c(absfsfo)
  #plot(absfsfo_nl, type='l')
  fsfo_nl<- c(fsfo)
  #plot(fsfo_nl, type='l')

  ## NB (fix): each colSums(...) term below is a length-nl per-layer scalar
  ## that must scale the corresponding COLUMN of an (nf x nl) matrix. Plain
  ## `v * M` in R recycles v over M's column-major-flattened elements, which
  ## only gives correct per-column scaling when nrow(M) == length(v) -- here
  ## nrow(M) == nf (53) while length(v) == nl, so for essentially every real
  ## canopy (nl != 53) the naive form silently produced scrambled values.
  ## Fixed to use sweep(M, 2, v, '*') for proper per-column scaling.
  wfEs_1 <- sweep(MpluEsun, 2, colSums(etau_lidf  * absfsfo_nl), '*') ## dim nf x nl
  wfEs_2 <- sweep(MminEsun, 2, colSums(etau_lidf * fsfo_nl), '*') ## dim nf x nl

  wfEs <- wfEs_1 + wfEs_2
  #plot(wfEs[,1], type='l')

  absfs_nl <- c(absfs)  ## NB (fix): was `c(absfsfo)`, an uncorrected copy-paste from the wfEs block above
  fsctl_nl <- c(fsctl)
  #plot(absfs_nl, type='l')

  ### same values
  sfEs <-  sweep(MpluEsun, 2, colSums(etau_lidf * absfs_nl), '*')  - sweep(MminEsun, 2, colSums(etau_lidf * fsctl_nl), '*')
  sbEs <-  sweep(MpluEsun, 2, colSums(etau_lidf * absfs_nl), '*')  + sweep(MminEsun, 2, colSums(etau_lidf * fsctl_nl), '*')

  absfo_nl <- c(absfo)
  foctl_nl <- c(foctl)

  vfEplu_h <-  sweep(MpluEplu, 2, colSums(etah_lidf * absfo_nl), '*')  - sweep(MminEplu, 2, colSums(etah_lidf * foctl_nl), '*')
  vfEplu_u <-  sweep(MpluEplu, 2, colSums(etau_lidf * absfo_nl), '*')  - sweep(MminEplu, 2, colSums(etau_lidf * foctl_nl), '*')

  vbEmin_h <-  sweep(MpluEmin, 2, colSums(etah_lidf * absfo_nl), '*')  + sweep(MminEmin, 2, colSums(etah_lidf * foctl_nl), '*')
  vbEmin_u <-  sweep(MpluEmin, 2, colSums(etau_lidf * absfo_nl), '*')  + sweep(MminEmin, 2, colSums(etau_lidf * foctl_nl), '*')

  ctl2_nl <- c(ctl2)

  sigfEmin_h <-  sweep(MpluEmin, 2, colSums(etah_lidf), '*')  - sweep(MminEmin, 2, colSums(etah_lidf * ctl2_nl), '*')
  sigfEmin_u <-  sweep(MpluEmin, 2, colSums(etau_lidf), '*')  - sweep(MminEmin, 2, colSums(etau_lidf * ctl2_nl), '*')


  sigbEmin_h <-  sweep(MpluEmin, 2, colSums(etah_lidf), '*')  + sweep(MminEmin, 2, colSums(etah_lidf * ctl2_nl), '*')
  sigbEmin_u <-  sweep(MpluEmin, 2, colSums(etau_lidf), '*')  + sweep(MminEmin, 2, colSums(etau_lidf * ctl2_nl), '*')

  #plot(colSums(etah_lidf),lwd=2,type='l')

  sigfEplu_h <-  sweep(MpluEplu, 2, colSums(etah_lidf), '*')  - sweep(MminEplu, 2, colSums(etah_lidf * ctl2_nl), '*')
  sigfEplu_u <-  sweep(MpluEplu, 2, colSums(etau_lidf), '*')  - sweep(MminEplu, 2, colSums(etau_lidf * ctl2_nl), '*')

  sigbEplu_h <-  sweep(MpluEplu, 2, colSums(etah_lidf), '*')  + sweep(MminEplu, 2, colSums(etah_lidf * ctl2_nl), '*')
  sigbEplu_u <-  sweep(MpluEplu, 2, colSums(etau_lidf), '*')  + sweep(MminEplu, 2, colSums(etau_lidf * ctl2_nl), '*')

  ################################################################
  ##   Emitted fluorescence
  ################################################################
  #  sunlit for each layer
  piLs <- (wfEs + vfEplu_u + vbEmin_u) # sunlit for each layer
  #plot(piLs[,1], type='l')

  # shade leaf for each layer
  piLd <- vbEmin_h + vfEplu_h
  #plot(piLd[,1], type='l')


  #Eq. 29a for sunlit leaf
  Fsmin <- sfEs + sigfEmin_u + sigbEplu_u
  # Eq. 29b for sunlit leaf
  Fsplu <- sbEs + sigbEmin_u + sigfEplu_u
  #Eq. 29a for shade leaf
  Fdmin <- sigfEmin_h + sigbEplu_h

  #Eq. 29a for shade leaf
  Fdplu <- sigbEmin_h + sigfEplu_h

  ## NB (fix): same column-recycling issue as the ~16 expressions above --
  ## `Qs`/`1-Qs` are length-nl per-layer scalars that must scale Fsmin/
  ## Fdmin's (nf x nl) columns; plain `Qs * Fsmin` only does that correctly
  ## when nrow(Fsmin) == length(Qs) (nf, not nl). Fixed with sweep().
  Femmin <- iLAI * sweep(Fsmin, 2, Qs, '*') + iLAI * sweep(Fdmin, 2, 1 - Qs, '*')
  Femplu <- iLAI * sweep(Fsplu, 2, Qs, '*') + iLAI * sweep(Fdplu, 2, 1 - Qs, '*')
  #plot(Femmin[,1], type='l')



  Y <- matrix(0, nrow = nl, ncol = dim(data.leafopt[['Mb']])[1] )
  # from bottom to top

  for (j in nl:1) {
    Y[j, ]  <- (rho_dd[j, ] * U[j+1, ] + Femmin[, j]) / (1 - rho_dd[j, ] * R_dd[j+1, ])
    U[j, ] <- tau_dd[j, ] * (R_dd[j+1, ] * Y[j, ] + U[j+1, ]) + Femplu[, j]
  }

  #plot(U[1,], type='l')

  # from top to bottom

  for (j in 1:nl) {
    Fmin_[j+1, ] <- Xdd[j, ] * Fmin_[j, ] + Y[j, ]
    Fplu_[j, ] <- R_dd[j, ] * Fmin_[j, ] + U[j, ]
  }


  piLo1 <- iLAI * piLs[,1:nl] %*% Pso[1:nl]
  #plot(piLo1, type='l')
  piLo2 <- iLAI * piLd  %*% (Po[1:nl] - Pso[1:nl])
  #plot(piLo2, type='l')
  piLo3 <- iLAI * t(vb * Fmin_[layers, ] + vf * Fplu_[layers, ]) %*% Po[1:nl]
  #plot(piLo3, type='l')

  piLo4 <- rs *  (Fmin_[nl+1, ]) *  (Po[nl+1])

  #plot(piLo4, type='l')

  piLtot <- piLo1 + piLo2 + piLo3 + piLo4
  #plot(piLtot, type='l')


  LoF_ <- piLtot / pi
  #plot(LoF_,type='l')

  Fhem_ <- Fplu_[1, ]
  #plot(Fhem_, type='l')

  ###############################################################################################
  # Output: a rad object with a large number of radiative fluxes: spectrally distributed
  #         and integrated, and canopy radiative transfer coefficients.
  #         Here, fluorescence fluxes are added
  ###############################################################################################

  data.rad[['LoF_']] <- signal::interp1(wlF, c(LoF_), xi = data.spectral[['wlF']], method = 'spline')
  data.rad[['EoutF_']] <- signal::interp1(wlF, c(Fhem_), xi = data.spectral[['wlF']], method = 'spline')

  data.rad[['LoF_sunlit']] <- signal::interp1(wlF, c(piLo1/pi), xi = data.spectral[['wlF']], method = 'spline')
  data.rad[['LoF_shaded']] <- signal::interp1(wlF, c(piLo2/pi), xi = data.spectral[['wlF']], method = 'spline')
  data.rad[['LoF_scattered']] <- signal::interp1(wlF, c(piLo3/pi), xi = data.spectral[['wlF']], method = 'spline')
  data.rad[['LoF_soil']] <- signal::interp1(wlF, c(piLo4/pi), xi = data.spectral[['wlF']], method = 'spline')

  data.rad[['EoutF']] <- 0.001 * Sint(Fhem_, wlF)
  data.rad[['LoutF']] <- 0.001 * Sint(LoF_, wlF)

  ### new
  data.rad[['Femliave_']] <- signal::interp1(wlF, c( rowSums(Femmin + Femplu)), xi = data.spectral[['wlF']], method = 'spline')



  data.rad[['F685']] <- max(data.rad[['LoF_']][1:55], na.rm=T) ### 55 is 69
  iwl685 <- which.max(data.rad[['LoF_']][1:55])
  data.rad[['wl685']] <- data.spectral[['wlF']][iwl685] ### 55 is 69

  if (iwl685 == 55) {
    data.rad[['F685']] <- NA
    data.rad[['wl685']] <- NA
  }


  data.rad[['F740']] <- max(data.rad[['LoF_']][70:length(data.rad[['LoF_']])],na.rm=T)
  data.rad[['wl740']] <- data.spectral[['wlF']][which.max(data.rad[['LoF_']][70:length(data.rad[['LoF_']])]) + 69]
  data.rad[['F684']] <- data.rad[['LoF_']][685 - data.spectral[['wlF']][1]]
  data.rad[['F761']] <- data.rad[['LoF_']][762 - data.spectral[['wlF']][1]]

  return(data.rad)


  #####################################################################################################
  ###### PLots
  ####################################################################################################

  if (get.plots ==  T){


    #####################################################################################################
    ###### Esun_f : Direct solar in the fluorescence emission range
    ####################################################################################################


    df.plot_ <- Esunf_
    # Create a data frame for the spectra
    df.plot_ <- data.frame(wavelength = wlE,
                           Esunf_= Esunf_)

    # Plot the spectra
    plot.Esunf <-ggplot(df.plot_, aes(x = wavelength, y = Esunf_)) +
      geom_line() +
      labs(x = "Wavelength", y = "Direct solar") +
      scale_color_manual(values = c("blue"),
                         labels = c("Esum")) +

      theme_bw()  +
      theme(legend.position = "top") +
      guides(color = guide_legend(title = NULL))

    print(plot.Esunf)
    #####################################################################################################
    ###### Upward diffuse radiation in the canopy (mW m-2 um-1) in the fluorescence emission range
    ####################################################################################################

    df.plot_ <- Epluf_

    sp.top <- df.plot_[,1]
    sp.bottom <- df.plot_[,nl]


    # Create a data frame for the spectra
    df <- data.frame(wavelength = wlE,
                     spectrum1 = sp.top,
                     spectrum2 = sp.bottom)
    # Reshape the data frame to long format
    df_long <- tidyr::gather(df, spectrum, value, -wavelength)

    # Plot the spectra
    plot.Epluf <-ggplot(df_long, aes(x = wavelength, y = value, color = spectrum)) +
      geom_line() +
      labs(x = "Wavelength", y = "upward diffuse radiation in the canopy") +
      scale_color_manual(values = c("blue", "red"),
                         labels = c("Spectrum top", "Spectrum bottom")) +

      theme_bw()  +
      theme(legend.position = "top") +
      guides(color = guide_legend(title = NULL))

    print(plot.Epluf)


    #####################################################################################################
    ######  Eminf_: downward diffuse radiation in the canopy in the fluorescence emission range
    ####################################################################################################
    ##

    df.plot_ <- Eminf_

    sp.top <- df.plot_[,1]
    sp.bottom <- df.plot_[,nl]


    # Create a data frame for the spectra
    df <- data.frame(wavelength = wlE,
                     spectrum1 = sp.top,
                     spectrum2 = sp.bottom)
    # Reshape the data frame to long format
    df_long <- tidyr::gather(df, spectrum, value, -wavelength)

    # Plot the spectra
    plot.Eminf <-ggplot(df_long, aes(x = wavelength, y = value, color = spectrum)) +
      geom_line() +
      labs(x = "Wavelength", y = "downward diffuse radiation in the canopy") +
      scale_color_manual(values = c("blue", "red"),
                         labels = c("top layer", "bottom layer")) +

      theme_bw()  +
      theme(legend.position = "top") +
      guides(color = guide_legend(title = NULL))

    print(plot.Eminf)




    #####################################################################################################
    ######  directional forward/backscatter scattering coefficients for diffuse incidence
    ####################################################################################################


    df.plot_ <- data.frame(wavelength = wlF, vf=vf[1,], vb=vb[1,])

    # Plot the spectra
    plot_ <-ggplot(df.plot_, aes(x = wavelength, y = vf)) +
      geom_line(aes(y = vf, color = "vf"), linewidth = 0.5) +
      geom_line(aes(y = vb, color = "vb"), linewidth = 0.5) +
      labs(x = "Wavelength", y = "directional forward/backscatter scattering coef.") +

      scale_color_manual(values = c("blue",'red'),
                         labels = c("vf",'vb')) +

      theme_bw()  +
      theme(legend.position = "top") +
      guides(color = guide_legend(title = NULL))

    print(plot_)


    #####################################################################################################
    ######  diffuse-diffuse parameters
    ####################################################################################################


    df.plot_ <- data.frame(wavelength = wlF, tau_dd=tau_dd[1,], Xdd=Xdd[1,],
                           rho_dd=rho_dd[1,],R_dd=R_dd[1,])

    # Plot the spectra
    plot_ <-ggplot(df.plot_, aes(x = wavelength, y = tau_dd)) +
      geom_line(aes(y = tau_dd, color = "tau_dd"), linewidth = 0.5) +
      geom_line(aes(y = Xdd, color = "Xdd"), linewidth = 0.5) +
      geom_line(aes(y = rho_dd, color = "rho_dd"), linewidth = 0.5) +
      geom_line(aes(y = R_dd, color = "R_dd"), linewidth = 0.5) +
      labs(x = "Wavelength", y = 'diffuse-diffuse parameters') +

      scale_color_manual(values = c("blue",'red','forestgreen','black'),
                         labels = c("tau_dd",'Xdd','rho_dd','R_dd')) +

      theme_bw()  +
      theme(legend.position = "top") +
      guides(color = guide_legend(title = NULL))

    print(plot_)



    #####################################################################################################
    ######  Mb: backward scattering fluorescence matrix, I for PSI and II for PSII
    #####   Mf: forward scattering fluorescence matrix,  I for PSI and II for PSII
    ####################################################################################################

    spectrum.Mf <- Mf[,1,1]
    spectrum.Mb <- Mb[,1,1]


    # Create a data frame for the spectra
    df_MB_Mf <- data.frame(wavelength = wlF,
                           spectrum.Mf = spectrum.Mf,
                           spectrum.Mb = spectrum.Mb)
    # Reshape the data frame to long format
    df_plot <- tidyr::gather(df_MB_Mf,  spectrum, value, -wavelength)


    # Plot the spectra
    plot_ <-ggplot(df_plot, aes(x = wavelength, y = value, color = spectrum)) +
      geom_line() +
      labs(x = "Wavelength", y = "fluorescence matrix") +
      scale_color_manual(values = c("blue",'red'),
                         labels = c("Mf",'Mb')) +

      theme_bw()  +
      theme(legend.position = "top") +
      guides(color = guide_legend(title = NULL))

    print(plot_)


    #####################################################################################################
    ######  PLots for emitted  fluorescence
    ####################################################################################################


    check.fluo <- data.frame(wave=wlF,piLs=piLs[,1],piLd=piLd[,1],
                             Fsmin=Fsmin[,1],
                             Fsplu=Fsplu[,1],
                             Fdmin=Fdmin[,1], Fdplu=Fdplu[,1],
                             Femmin = Femmin[,1], Femplu=Femplu[,1])


    p.fluo <- ggplot(data = check.fluo, aes(x = wave)) +
      labs(y= "emitted  fluorescence", x = "")+
      geom_line(aes(y = piLs, color = "piLds"), linewidth = 0.5) +
      geom_line(aes(y = piLd, color = "piLd"), linewidth = 0.5) +
      geom_line(aes(y = Fsmin, color = "Fsmin"), linewidth = 0.5) +
      geom_line(aes(y = Fsplu, color = "Fsplu"), linewidth = 0.5) +
      geom_line(aes(y = Fdmin, color = "Fdmin"), linewidth = 0.5) +
      geom_line(aes(y = Fdplu, color = "Fdplu"), linewidth = 0.5) +
      geom_line(aes(y = Femmin, color = "Femmin"), linewidth = 0.5) +
      geom_line(aes(y = Femplu, color = "Femplu"), linewidth = 0.5) +

      theme_bw()  +
      theme(legend.position = "top") +

      guides(color = guide_legend(title = NULL))
    print(p.fluo)

  }
  #####################################################################################################
  ###### end get.plots
  ####################################################################################################

}

