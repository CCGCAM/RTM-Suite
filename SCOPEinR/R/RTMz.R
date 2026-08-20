#' get.RTMz
#' \code{get.RTMz} Calculates the small modification of TOC outgoing radiance
#' due to the conversion of Violaxanthin into Zeaxanthin in leaves
#'
#' @param data.spectral information about wavelengths and resolutions
#' @param data.rad a large number of radiative fluxes: spectrally distributed and integrated, and canopy radiative transfer coefficients.
#' @param data.soil  soil properties
#' @param data.leafopt leaf optical properties
#' @param data.canopy canopy properties (such as LAI and height)
#' @param data.gap probabilities of direct light penetration and viewing
#' @param data.angles viewing and observation angles
#' @param data.Knu  ...kn for  sunlit leaves (data.bcu$Kn)
#' @param data.Knh .. kn for shaded leavess (data.bch$Kn)
#' @param get.plots is true plot the intermediate plots
#'
#' @return a rad object with a large number of radiative fluxes: spectrally distributed
#' and integrated, and canopy radiative transfer coefficients.
#' @export
#'
#' @description
#'
#' Date:   08 Dec 2016
#' Update:
#'  - 17 Mar 2020 CvdT: added cluming, mSCOPE representation
#'  - 25 Jun 2020 CvdT: Po, Ps, Pso. fix the problem we have with the oblique angles above 80 degrees

#' @author Christiaan van der Tol  (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#' @examples
#' \dontrun{
#' data.rad <- get.RTMz(data.spectral, data.rad, data.soil, data.leafopt,
#'                       data.canopy, data.gap, data.angles, data.Knu, data.Knh,
#'                       get.plots = FALSE)
#' }
get.RTMz <- function(data.spectral, data.rad, data.soil, data.leafopt, data.canopy, data.gap,
                     data.angles, data.Knu, data.Knh,
                     get.plots =T){

  if (missing(get.plots)){
    get.plots = FALSE
  }

  ## 0.1 initialisations

  wlS <- data.spectral[['wlS']] # SCOPE wavelengths, make column vectors
  wlZ <- data.spectral[['wlZ']]# Excitation wavelengths

  dummy <- c(wlZ[1]: wlZ[2])
  iwlfi <- match(dummy, wlS)


  nwlZ <- length(dummy)
  nl <- data.canopy[['nlayers']]
  LAI <- data.canopy[['LAI']]

  iLAI <- LAI/nl # LAI of a layer        [1]

  litab <- data.canopy[['litab']]
  lazitab <- data.canopy[['lazitab']]
  lidf <- data.canopy[['lidf']]

  nlazi <- length(lazitab) # azumith angle
  nlinc <- length(litab) # inclination
  nlori <- nlinc * nlazi # total number of leaf orientations
  layers <- 1:nl

  RZ <- t(data.leafopt[['reflZ']][,iwlfi] - data.leafopt[['refl']][,iwlfi])
  TZ <- t(data.leafopt[['tranZ']][,iwlfi] - data.leafopt[['tran']][,iwlfi])
  #plot(RZ[,1],type='l')


  Ps <- data.gap[['Ps']]
  Po <- data.gap[['Po']]
  Pso <- data.gap[['Pso']]

  Qs <- Ps[1:(length(Ps)-1)]


  # speed-up the calculation by only using wavelength i and wavelength o part of the spectrum
  ## Esun_f : Direct solar in the fluorescence emission range
  Esunf_ <- data.rad[['Esun_']][iwlfi]

  ## Eminf_: downward diffuse radiation in the canopy in the fluorescence emission range
  Eminf_ <- array(0, dim = c(nwlZ, nl+1, 2))
  Eminf_[,,1] <- t(data.rad[['Emins_']][,iwlfi])
  Eminf_[,,2] <- t(data.rad[['Emind_']][,iwlfi])

  ##Epluf_: upward diffuse radiation in the canopy (mW m-2 um-1) in the fluorescence emission range
  Epluf_ <- array(0, dim = c(nwlZ, nl+1, 2))
  Epluf_[,,1] <- t(data.rad[['Eplus_']][,iwlfi])
  Epluf_[,,2] <- t(data.rad[['Eplud_']][,iwlfi])


  Xdd <- data.rad[['Xdd']][,iwlfi]
  # rho_dd: diffuse-diffuse reflectance for the thin layers
  rho_dd <- data.rad[['rho_dd']][,iwlfi]
  R_dd <- data.rad[['R_dd']][,iwlfi]
  # tau_dd:  diffuse-diffuse transmittance for the thin layers
  tau_dd <- data.rad[['tau_dd']][,iwlfi]
  #vb: directional backscatter scattering coefficient for diffuse incidence
  vb <- data.rad[['vb']][,iwlfi]
  #vf: directional forward scattering coefficient for diffuse incidence
  vf <- data.rad[['vf']][,iwlfi]





  #plot(Mb[,2,1], type='l')
  # geometric factors
  deg2rad <- subset(SCOPEinR::constants,constant == 'deg2rad')[[2]]
  tto <- data.angles[['tto']]
  tts <- data.angles[['tts']]
  psi <- data.angles[['psi']]

  rs <- data.soil[['rfl.soil']][iwlfi]  # [nwlfo] soil reflectance
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
  absfs <- matrix(absfs, nlori, 1)  # [nlori,1]
  absfo <- matrix(absfo, nlori, 1)  # [nlori,1]
  fsfo <- matrix(fsfo, nlori, 1)    # [nlori,1]
  absfsfo <- matrix(absfsfo, nlori, 1)  # [nlori,1]
  foctl <- matrix(foctl, nlori, 1)  # [nlori,1]
  fsctl <- matrix(fsctl, nlori, 1)  # [nlori,1]
  ctl2 <- matrix(ctl2, nlori, 1)    # [nlori,1]

  #1.0 calculation of 'flux' in observation direction

  Fmin_ <- array(0, dim = c(nl+1, nwlZ, 2))
  Fplu_ <- array(0, dim = c(nl+1, nwlZ, 2))
  LoF_ <- matrix(0, nwlZ, 2)
  laz <- 1/36


  Kn2Cx <- function(Kn){
    #Cx = 0.70*Kn      # empirical fit by N Vilfan
    Cx = 0.3187 * Kn  # empirical fit by N Vilfan (Vilfan et al, 2018, 2019)
    return(Cx)
  }
  etah <- Kn2Cx(data.Knh)

  if (is.vector(data.Knu) == T) {
    ## NB (fix): was hardcoded `dim = c(30, 13, 36)` -- silently wrong for
    ## any canopy with nl != 30 (a leftover debug/example value). Use `nl`.
    etau <- array(Kn2Cx(data.Knu), dim = c(nl, 13, 36))
  } else {
    etau <- array((data.Knu), dim = c(length(data.Knu), 13, 36)) #make dimensions [nl,nlinc,nlazi]
    etau <- aperm(etau, c(3, 2, 1))  # Permute dimensions to match [nl, 13, 30]
    etau <- array(etau, dim = c(nl, 468))  # Reshape to [nl, nlori]
  }



  etau_lidf <- matrix(0, nrow = nlori, ncol = nl)  # Initialize etau_lidf matrix
  etah_lidf <- matrix(0, nrow = nlori, ncol = nl)  # Initialize etah_lidf matrix

  etau_reshape <- matrix(etau, nrow = nlori, ncol = nl, byrow = F)

  lidf_laz <- matrix(lidf * laz, nrow = nlori, ncol = 1)


  for (i in 1:nlori) {
    etau_lidf[i, ] <- etau_reshape[i,] * lidf_laz[i, ]
    etah_lidf[i, ] <- etah * lidf_laz[i, ]
  }

 for (k in 1:2){
   U <- Y <- matrix(0, nrow = nl + 1, ncol = nwlZ)

   if (k < 2) {
     MpluEsun <- RZ * Esunf_
     MminEsun <- TZ * Esunf_
   }

   MpluEmin <- RZ * Eminf_[, 1:nl, k]
   MpluEplu <- RZ * Epluf_[, 1:nl, k]
   MminEmin <- TZ * Eminf_[, 1:nl, k]
   MminEplu <- TZ * Epluf_[, 1:nl, k]

   absfsfo_nl<-matrix(c(absfsfo),nrow=dim(etau_lidf)[1],ncol=dim(etau_lidf)[2])
   fsfo_nl<-matrix(c(fsfo),nrow=dim(etau_lidf)[1],ncol=dim(etau_lidf)[2])

   ## NB (fix): same column-recycling issue as SCOPEinR::get.RTMf -- each
   ## colSums(...) term is a length-nl per-layer scalar that must scale the
   ## matching COLUMN of an (nwlZ x nl) matrix; plain `v * M` only does that
   ## correctly when nrow(M) == length(v) (here nrow(M) == nwlZ, not nl).
   ## Fixed to use sweep(M, 2, v, '*').
   wfEs_1 <- sweep(MpluEsun, 2, colSums(etau_lidf * absfsfo_nl), '*') ## dim nf x nl
   wfEs_2 <- sweep(MminEsun, 2, colSums(etau_lidf * fsfo_nl), '*') ## dim nf x nl
   wfEs <- wfEs_1 + wfEs_2
   #plot(wfEs[,1], type='l')


   #plot(etah_lidf[1,], type='l')

   #plot(wfEs[,1], type='l')

   absfo_nl<-matrix(c(absfo),nrow=dim(etau_lidf)[1],ncol=dim(etau_lidf)[2])
   absfs_nl<-matrix(c(absfs),nrow=dim(etau_lidf)[1],ncol=dim(etau_lidf)[2])
   fsctl_nl<-matrix(c(fsctl),nrow=dim(etau_lidf)[1],ncol=dim(etau_lidf)[2])

   ### same values
   sfEs <-  sweep(MpluEsun, 2, colSums(etau_lidf * absfs_nl), '*')  + sweep(MminEsun, 2, colSums(etau_lidf * fsctl_nl), '*')
   sbEs <-  sweep(MpluEsun, 2, colSums(etau_lidf * absfs_nl), '*')  + sweep(MminEsun, 2, colSums(etau_lidf * fsctl_nl), '*')

   foctl_nl<-matrix(c(foctl),nrow=dim(etau_lidf)[1],ncol=dim(etau_lidf)[2])

   vfEplu_h <-  sweep(MpluEplu, 2, colSums(etah_lidf * absfo_nl), '*')  + sweep(MminEplu, 2, colSums(etah_lidf * foctl_nl), '*')
   vfEplu_u <-  sweep(MpluEplu, 2, colSums(etau_lidf * absfo_nl), '*')  + sweep(MminEplu, 2, colSums(etau_lidf * fsctl_nl), '*')

   vbEmin_h <-  sweep(MpluEmin, 2, colSums(etah_lidf * absfo_nl), '*')  + sweep(MminEmin, 2, colSums(etah_lidf * foctl_nl), '*')
   vbEmin_u <-  sweep(MpluEmin, 2, colSums(etau_lidf * absfo_nl), '*')  + sweep(MminEmin, 2, colSums(etau_lidf * foctl_nl), '*')

   ctl2_nl<-matrix(c(ctl2),nrow=dim(etau_lidf)[1],ncol=dim(etau_lidf)[2])

   sigfEmin_h <-  sweep(MpluEmin, 2, colSums(etah_lidf), '*')  + sweep(MminEmin, 2, colSums(etah_lidf * ctl2_nl), '*')
   sigfEmin_u <-  sweep(MpluEmin, 2, colSums(etau_lidf), '*')  + sweep(MminEmin, 2, colSums(etau_lidf * ctl2_nl), '*')


   sigbEmin_h <-  sweep(MpluEmin, 2, colSums(etah_lidf), '*')  + sweep(MminEmin, 2, colSums(etah_lidf * ctl2_nl), '*')
   sigbEmin_u <-  sweep(MpluEmin, 2, colSums(etau_lidf), '*')  + sweep(MminEmin, 2, colSums(etau_lidf * ctl2_nl), '*')


   sigbEplu_h <-  sweep(MpluEplu, 2, colSums(etah_lidf), '*')  + sweep(MminEplu, 2, colSums(etah_lidf * ctl2_nl), '*')
   sigfEplu_u <-  sweep(MpluEplu, 2, colSums(etau_lidf), '*')  + sweep(MminEplu, 2, colSums(etau_lidf * ctl2_nl), '*')

   sigfEplu_h <-  sweep(MpluEplu, 2, colSums(etah_lidf), '*')  + sweep(MminEplu, 2, colSums(etah_lidf * ctl2_nl), '*')
   sigbEplu_u <-  sweep(MpluEplu, 2, colSums(etau_lidf), '*')  + sweep(MminEplu, 2, colSums(etau_lidf * ctl2_nl), '*')

   ################################################################
   ##   Emitted fluorescence
   ################################################################

   piLs <- (wfEs + vfEplu_u + vbEmin_u) # sunlit for each layer

   piLd <- vbEmin_h + vfEplu_h       # shade leaf for each layer
   Fsmin <- sfEs + sigfEmin_u + sigbEplu_u #Eq. 29a for sunlit leaf
   Fsplu <- sbEs + sigbEmin_u + sigfEplu_u #Eq. 29b for sunlit leaf
   Fdmin <- sigfEmin_h + sigbEplu_h  #Eq. 29a for shade leaf
   Fdplu <- sigbEmin_h + sigfEplu_h  #Eq. 29a for shade leaf

   Qs_byiLAI = iLAI * (Qs)
   Qs_byiLAI_dif = iLAI * (1-Qs)
   ## NB (fix): same column-recycling issue as the sweep()-fixed expressions
   ## above -- Qs_byiLAI/Qs_byiLAI_dif are length-nl per-layer scalars that
   ## must scale Fsmin/Fdmin's (nwlZ x nl) columns; plain `v * M` only does
   ## that correctly when nrow(M) == length(v) (here nwlZ, not nl).
   Femmin <- sweep(Fsmin, 2, Qs_byiLAI, '*') + sweep(Fdmin, 2, Qs_byiLAI_dif, '*')
   Femplu <- sweep(Fsplu, 2, Qs_byiLAI, '*') + sweep(Fdplu, 2, Qs_byiLAI_dif, '*')

   ################################################################
   ## NB (fix): the layer recursion (Y/U), the outgoing-radiance terms
   ## (piLo1..4) and LoF_[,k] used to live inside `if (get.plots==T)`,
   ## AFTER this k-loop had already ended -- since the loop's `k`/piLs/etc.
   ## variables are not k-indexed there, that block only ever ran (when it
   ## ran at all) using the LAST iteration's (k=2) values, and only with
   ## get.plots=TRUE; with the normal get.plots=FALSE call, LoF_ (and
   ## everything downstream of it) was never computed at all, so
   ## get.RTMz() silently returned its input unmodified. Moved inside the
   ## loop so it runs unconditionally, for both k=1 (sun) and k=2 (sky).
   ################################################################

   # from bottom to top
   for (j in nl:1) {
     Y[j, ]  <- (rho_dd[j, ] * U[j+1, ] + Femmin[, j]) / (1 - rho_dd[j, ] * R_dd[j+1, ])
     U[j, ] <- tau_dd[j, ] * (R_dd[j+1, ] * Y[j, ] + U[j+1, ]) + Femplu[, j]
   }

   # from top to bottom
   for (j in 1:nl) {
     Fmin_[j+1,,k ] <- Xdd[j, ] * Fmin_[j,,k ] + Y[j, ]
     Fplu_[j,,k ] <- R_dd[j, ] * Fmin_[j,,k ] + U[j, ]
   }

   piLo1 <- iLAI * (piLs[,1:nl]) %*% Pso[1:nl]
   piLo2 <- iLAI * piLd  %*% (Po[1:nl] - Pso[1:nl])
   piLo3 <- iLAI * t(vb * Fmin_[layers,,k] + vf * Fplu_[layers,,k ]) %*% Po[1:nl]
   ## NB (fix): was `Po[1:nl+1]` -- R parses this as `(1:nl)+1`, a length-nl
   ## vector, not the intended scalar soil-boundary gap probability
   ## `Po[nl+1]` (matches the equivalent, correctly-written line in
   ## SCOPEinR::get.RTMf).
   piLo4 <- rs * Fmin_[nl+1,,k ] * Po[nl+1]

   piLtot <- piLo1 + piLo2 + piLo3 + piLo4

   LoF_[,k] <- piLtot / pi

 }






  if (get.plots ==  T){
    check.fluo <- data.frame(wave=dummy,piLs=piLs[,1],piLd=piLd[,1],
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

  Fhem_ <- colSums(Fplu_[1, , , drop = FALSE])
  Fhem_ <- rowSums(Fhem_)
  #plot(Fhem_, type='l')

  ###############################################################################################
  # Output: a rad object with a large number of radiative fluxes: spectrally distributed
  #         and integrated, and canopy radiative transfer coefficients.
  #         Here, fluorescence fluxes are added
  ###############################################################################################

  ## NB (fix): was `sum(LoF_,2)` -- a MATLAB-to-R mistranslation. In MATLAB,
  ## `sum(A,2)` sums along the 2nd dimension (row sums); R's `sum()` has no
  ## such dimension argument -- it just treats `2` as another value to add,
  ## collapsing the whole (nwlZ,2) matrix to one scalar (total sum + 2) and
  ## recycling it across every wavelength in `iwlfi`, instead of adding the
  ## per-wavelength sun+sky total. Confirmed via direct repro
  ## (`sum(matrix(1:6,3,2),2)` == 23, not `rowSums(...)` == c(5,7,9)).
  data.rad[['Lo_']][iwlfi] <-   data.rad[['Lo_']][iwlfi] + rowSums(LoF_)
  data.rad[['rso']][iwlfi] <- data.rad[['rso']][iwlfi]+ LoF_[,1] / (data.rad[['Esun_']][iwlfi] )
  data.rad[['rdo']][iwlfi] <- data.rad[['rdo']][iwlfi]+ LoF_[,2] / (data.rad[['Esky_']][iwlfi] )
  data.rad[['rfl']][iwlfi] <- pi * data.rad[['Lo_']][iwlfi] / ( data.rad[['Esky_']][iwlfi] + data.rad[['Esun_']][iwlfi] )     #[nwl]
  data.rad[['Eout_']][iwlfi] <- data.rad[['Eout_']][iwlfi] + Fhem_

  return(data.rad)


}

