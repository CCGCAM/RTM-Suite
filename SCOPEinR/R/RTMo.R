

#' \code{getRTMo}
#' Calculates the spectra of hemisperical and directional observed visible
#' and thermal radiation (fluxes E and radiances L), as well as the single
#' and bi-directional gap probabilities
#'
#' @param data.spectral information about wavelengths and resolutions
#' @param atmo MODTRAN atmospheric parameters
#' @param data.soil soil properties
#' @param data.leafopt leaf optical properties
#' @param data.canopy canopy properties (such as LAI and height)
#' @param data.leafbio leaf biochemical parameters (Cab, Car...)
#' @param data.angles viewing and observation angles
#' @param data.meteo has the meteorological variables. Is only used to correct,
#' the total irradiance if a specific value is provided instead of the usual Modtran output.
#' @param data.opts simulation options. Here, the option
#' @param canopy.model Selection of canopy model. THe canopy models available are 'fourSAIL', 'fourSAIL' and 'INFORM'. By default fourSAIL model will be used.
#' @param get.plots  is true plot the intermediate plots
#' @description
#'
#' updates:
#'  - 10 Sep 2007 (CvdT)  - calculation of Rn
#'  - 5 Nov 2007          - included observation direction
#'  - 12 Nov 2007         - included abs. PAR spectrum output
#'                        - improved calculation efficiency
#'
#' - 13 Nov 2007          - written readme lines
#'
#' - 11 Feb 2008 (WV&JT)  - changed Volscat
#'                        - small change in calculation Po,Ps,Pso  (author:JT)
#'                        - introduced parameter 'lazitab'
#'                        - changed nomenclature
#'                        - Appendix IV: cosine rule
#'
#' - 04 Aug 2008 (JT)     - Corrections for Hotspot effect in the probabilities
#'
#' - 05 Nov 2008 (CvdT)   - Changed layout
#'
#' - 04 Jan 2011          - Included Pso function (Appendix IV) (JT&CvdT)
#'                        - removed the analytical function (for checking)(JT&CvdT)
#' - 02 Oct 2012 (CvdT)   - included incident PAR in output
#'
#' - Jan/Feb 2013 (WV)    - Major revision towards SCOPE version 1.40:
#'                        - Parameters passed using structures
#'                        - Improved interface with MODTRAN atmospheric data
#'                        - Now also calculates 4-stream
#'                        - reflectances rso, rdo, rsd and rdd analytically
#'
#' - Apri 2013 (CvT)         - improvements in variable names and descriptions
#'
#' - Dec 2019 CvdT        mSCOPE representation, lite option
#'
#' @references
#' Verhoef (1998), 'Theory of radiative transfer models applied in
#' optical remote sensing of vegetation canopies'. PhD Thesis Univ. Wageninegn.
#' Verhoef, W., Jia, L., Xiao, Q. and Su, Z. (2007) Unified optical -
#' thermal four - stream radiative transfer theory for homogeneous
#' vegetation canopies. IEEE Transactions on geoscience and remote sensing, 45,6.
#' Verhoef (1985), 'Earth Observation Modeling based on Layer Scattering Matrices',
#' Remote sensing of Environment, 17:167-175.
#'
#' @author 	 Wout Verhoef, Christiaan van der Tol, Joris Timmermans (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' rad <- getRTMo(data.spectral, atmo, data.soil, data.leafopt, data.canopy,
#'                 data.leafbio, data.angles, data.meteo, data.opts,
#'                 canopy.model = "fourSAIL", get.plots = FALSE)
#' }


getRTMo<-function(data.spectral,atmo,
                  data.soil,
                  data.leafopt,data.canopy,data.leafbio,
                  data.angles, data.meteo,
                  data.opts,canopy.model,get.plots = T)  {


  ###############################################
  ###############################################
  # Table of contents of the function
  #
  #   0.      Preparations
  #       0.1     parameters
  #       0.2     initialisations
  #   1.      Geometric quantities
  #       1.1     general geometric quantities
  #       1.2     geometric factors associated with extinction and scattering
  #       1.3     geometric factors to be used later with rho and tau
  #       1.4     solar irradiance factor for all leaf orientations
  #       1.5     probabilities Ps, Po, Pso
  #   2.      Calculation of upward and downward fluxes
  #   3.      Outgoing fluxes, hemispherical and in viewing direction, spectrum
  #   4.      Net fluxes, spectral and total, and incoming fluxes
  #   A1      functions J1 and J2 (introduced for stable solutions)
  #   A2      function volscat
  #   A3      function e2phot
  #   A4      function Pso
  ###############################################
  ###############################################


  if (missing(canopy.model)){
    model = 'fourSAIL'

  } else {
    model = canopy.model
  }

  if (missing(data.opts)){
   stop('please use options for Lite, Calc_vert_profiles ')
  } else{
    options.lite                = data.opts[1,]   # lite version
    options.simulation          = data.opts[10,]
    options.calc_vert_profiles   = data.opts[12,]

  }


  constants <- SCOPEinR::constants

  if (missing(get.plots)){
    get.plots = FALSE
  }
  ## 0. Preparations
  deg2rad <- subset(SCOPEinR::constants,constant == 'deg2rad')[[2]] # degree to rad

  #wlS is all regions reg1,re2,r3
  wl <- data.spectral[['wlS']] # SCOPE wavelengths as a column-vector
  nwl <- length(wl) #

  wlPAR <- data.spectral[['wlPAR']] # PAR wavelength range
  minPAR <- min(wlPAR) # min PAR
  maxPAR <- max(wlPAR) # max PAR
  Ipar <- which(wl>=minPAR & wl<=maxPAR) # Indices for PAR wavelenghts within wl
  tts <- data.angles[['tts']] # solar zenith angle
  tto <- data.angles[['tto']] # observer zenith angle
  psi <- data.angles[['psi']] # relative azimuth anglee

  nl <- data.canopy[['nlayers']] # number of canopy layers (nl)
  litab <- data.canopy[['litab']] # SAIL leaf inclibation angles # leaf inclination angles PY
  lazitab <- data.canopy[['lazitab']] # leaf azimuth angles relative to the sun
  nlazi <- data.canopy[['nlazi']] # number of azimuth angles (36)
  LAI <- data.canopy[['LAI']] # leaf area index
  lidf <- data.canopy[['lidf']] # leaf inclination distribution function
  xl <- data.canopy[['xl']] # all levels except for the top
  dx <- 1/nl

  ############################################################################
  #### Here you can use the Leaf model for estimating the canopy reflectance
  ############################################################################

  # Add 161 values of 0.01 to each row

  adding_rho_thermal <- matrix(data.leafbio[['rho_thermal']], nrow = nl, ncol = 161)
  adding_tau_thermal <- matrix(data.leafbio[['tau_thermal']], nrow = nl, ncol = 161)



  rho <- data.leafopt[['refl']] #S pectral reflectance of the leaf, 400 to 2400 nm
  rho <- cbind(rho, adding_rho_thermal)

  tau <- data.leafopt[['tran']]  # Spectral transmittance of the leaf, 400 to 2400 nm

  tau <- cbind(tau, adding_tau_thermal)

  kChlrel <- data.leafopt[['kChlrel']] #relative portion of chlorophyll contribution to reflectance
  kChlrel <- cbind(kChlrel, adding_rho_thermal)

  kCarrel <- data.leafopt[['kCarrel']] #relative portion of carotenoids  contribution to reflectance
  kCarrel <- cbind(kCarrel, adding_rho_thermal)

  rsoil <- data.soil[['rfl.soil']] # [nwl,nsoils] soil reflectance spectra
  rsoil <- c(rsoil, rep(data.leafbio[['rho_thermal']],161))

  if (get.plots ==  T){

    check.leaf<- data.frame(wave=data.spectral$wlS,rfl=rho[1,],trans=tau[1,], rsoil= rsoil,
                            kCarrel=kCarrel[1,],
                            kChlrel=kChlrel[1,])


    p.leaf <- ggplot(data = check.leaf, aes(x = wave)) +
      labs(y= "Reflectance", x = "") +
      geom_line(aes(y = rfl, color = "leaf reflectance"), linewidth = 0.5) +
      geom_line(aes(y = rsoil, color = "soil reflectance"), linewidth = 0.5) +
      scale_color_manual(name = "r",
                         values = c("leaf reflectance" = "navyblue", "soil reflectance" ='black')) +
      theme_bw() + xlim(400,2499) +
      theme(legend.position = "right") +
      guides(color = guide_legend(title = NULL))

    print(p.leaf)


    p.leaf.proportions <- ggplot(data = check.leaf, aes(x = wave)) +
      labs(y= "contribution of main pigments", x = "") +
      geom_line(aes(y = kChlrel, color = "chlorophylls"), linewidth = 0.5) +
      geom_line(aes(y = kCarrel, color = "carotenoids"), linewidth = 0.5) +
      scale_color_manual(name = "r",
                         values = c("chlorophylls" = "forestgreen","carotenoids" ='brown')) +
      theme_bw() + xlim(400,800) +
      theme(legend.position = "right") +
      guides(color = guide_legend(title = NULL))

    print(p.leaf.proportions)
  }


  ############################################################################
  ############################################################################

  epsc <- 1-rho-tau # [nl,nwl] emissivity of leaves

  epss <- 1-rsoil # [nwl] emissivity of soil

  iLAI <- LAI/nl # [1] LAI of elementary layer


  ###########################################################################
  # Initializations  variable for the sum
  ###########################################################################

  Rndif <- rep(0, nl)      # [nl]         abs. diffuse rad soil+veg


  Pdif <- rep(0, nl)       # [nl]           incident PAR flux of veg layers
  Pndif <- rep(0, nl)        # [nl]           net PAR veg
  Pndif_Cab <- rep(0, nl)    # [nl]           net PAR veg due to Cab
  Rndif_Cab <- rep(0, nl)    # [nl]           abs. diffuse rad soil+veg due to Cab

  Pndif_Car <- rep(0, nl)    # [nl]           net PAR veg due to Car
  Rndif_Car <- rep(0, nl)    # [nl]           abs. diffuse rad soil+veg due to Car

  Rndif_PAR <- rep(0, nl)      # [nl,nwl]       abs diff and PAR veg.


  Rndif_ <- matrix(0, nl, nwl)          # [nl,nwl]       abs diff and PAR veg.
  Pndif_ <- matrix(0, nl, ncol=length(wlPAR))      # [nl,nwlPAR]    net PAR veg

  Rndif_PAR_ <- matrix(0, nl, ncol=length(wlPAR))  # [nl,nwlPAR]    abs. diffuse rad soil+veg

  ##    IwlP is (400:2500)
  Pndif_Cab_ <- matrix(0, nl, length(data.spectral[['IwlP']]))        # [nl,length(IwlP)] net PAR veg due to Cab
  Rndif_Cab_ <- matrix(0, nl, length(data.spectral[['IwlP']]))        # [nl,length(IwlP)] abs. diffuse rad soil+veg due to Cab

  Pndif_Car_ <- matrix(0, nl, length(data.spectral[['IwlP']]))        # [nl,length(IwlP)] net PAR veg due to Car
  Rndif_Car_ <- matrix(0, nl, length(data.spectral[['IwlP']]))        # [nl,length(IwlP)] abs. diffuse rad soil+veg due to Car

  ###########################################################################


  # 1. Geometric quantities
  # 1.1 general geometric quantities. these variables are scalars
  #deg2rad <- pi/180   # degree to rad
  cos_tts <- cos(tts*deg2rad)  # cos solar angle
  tan_tto <- tan(tto*deg2rad)  # tan observation angle

  cos_tto <- cos(tto*deg2rad)  # cos observation angle
  sin_tts <- sin(tts*deg2rad)  # sin solar angle
  tan_tts <- tan(tts*deg2rad)  # tan observation angle

  psi <- abs(psi-360*round(psi/360))  # to ensure that volscatt is symmetric for psi=90 and psi=270
  dso <- sqrt(tan_tts^2 + tan_tto^2 - 2*tan_tts*tan_tto*cos(psi*deg2rad))

  # 1. Geometric quantities
  # 1.1 general geometric quantities. these variables are scalars
  cos_tts <- cos(tts * pi / 180)             # cos solar angle
  tan_tto <- tan(tto * pi / 180)             # tan observation angle
  cos_tto <- cos(tto * pi / 180)             # cos observation angle
  sin_tts <- sin(tts * pi / 180)             # sin solar angle
  tan_tts <- tan(tts * pi / 180)             # tan observation angle
  psi <- abs(psi - 360 * round(psi / 360))   # to ensure that volscatt is symmetric for psi=90 and psi=270
  dso <- sqrt(tan_tts^2 + tan_tto^2 - 2 * tan_tts * tan_tto * cos(psi * pi / 180))

  # 1.2 geometric factors associated with extinction and scattering
  volscatt_geo <-  get.volscatt.scope(tts,tto,psi,litab)

  chi_s <- volscatt_geo[[1]]
  chi_o<-volscatt_geo[[2]]
  frho<-volscatt_geo[[3]]
  ftau<-volscatt_geo[[4]]

  cos_ttlo <- cos(lazitab * pi / 180)         # cos leaf azimuth angles
  cos_ttli <- cos(litab * pi / 180)           # cos leaf angles
  sin_ttli <- sin(litab * pi / 180)           # sinus leaf angles

  ksli <- chi_s / cos_tts                     # extinction coefficient in direction of sun per leaf angle
  koli <- chi_o / cos_tto                     # extinction coefficient in direction of observer per leaf angle
  sobli <- frho * pi / (cos_tts * cos_tto)    # area scattering coefficient fractions
  sofli <- ftau * pi / (cos_tts * cos_tto)
  bfli <- cos_ttli^2

  # integration over angles (using a vector inproduct) -> scalars

  k <- sum(ksli * lidf) # extinction coefficient in direction of sun
  K <- sum(koli * lidf)  # extinction coefficient in direction of observer
  bf <- sum(bfli * lidf)
  sob <-sum(sobli * lidf)  # weight of specular2directional back scatter coefficient
  sof <- sum(sofli * lidf)  # weight of specular2directional forward scatter coefficient

  # 1.3 geometric factors to be used later with rho and tau, f1 f2 of pag 304:
  # these variables are scalars
  sdb <- 0.5 * (k + bf)                   # fs*f1
  sdf <- 0.5 * (k - bf)                   # fs*f2     weight of specular2diffuse     foward  scatter coefficient

  ddb <- 0.5 * (1 + bf)                   # f1^2+f2^2 weight of diffuse2diffuse      back    scatter coefficient
  ddf <- 0.5 * (1 - bf)                   # 2*f1*f2   weight of diffuse2diffuse      forward scatter coefficient

  dob <- 0.5 * (K + bf)                   # fo*f1     weight of diffuse2directional  back    scatter coefficient
  dof <- 0.5 * (K - bf)                   # fo*f2     weight of diffuse2directional  forward scatter coefficient

  # 1.4 solar irradiance factor for all leaf orientations
  Css <- cos_ttli * cos_tts
  Ss <- sin_ttli * sin_tts
  ###

  cos_deltas <- Css %*% matrix(1, nrow = 1, ncol = nlazi) + Ss %*% matrix(cos_ttlo, ncol = nlazi)
  fs <- abs(cos_deltas / cos_tts)

  # 2. Calculation of reflectance

  # 2.1 reflectance, transmittance factors in a thin layer
  # the following are vectors with length [nl,nwl]
  sigb <- ddb[1] * rho + ddf[1] * tau          # [nl,nwl] sigmab, p305{1} diffuse backscatter scattering coefficient for diffuse incidence
  sigf <- ddf[1] * rho + ddb[1] * tau          # [nl,nwl] sigmaf, p305{1} diffuse forward scattering coefficient for forward incidence
  sb <- sdb[1] * rho + sdf[1] * tau            # [nl,nwl] sb, p305{1} diffuse backscatter scattering coefficient for specular incidence
  sf <- sdf[1] *rho + sdb[1] * tau            # [nl,nwl] sf, p305{1} diffuse forward scattering coefficient for specular incidence
  vb <- dob[1] * rho + dof[1] * tau            # [nl,nwl] vb, p305{1} directional backscatter scattering coefficient for diffuse incidence
  vf <- dof[1] * rho + dob[1] * tau            # [nl,nwl] vf, p305{1} directional forward scattering coefficient for diffuse incidence
  w <- sob[1] * rho + sof[1] * tau             # [nl,nwl] w, p309{1} bidirectional scattering coefficient (directional-directional)
  a <- 1 - sigf                      # [nl,nwl] attenuation

  ## 3. Flux calculation
  # diffuse fluxes within the vegetation covered part

  #iLAI <- LAI/nl
  tau_ss <- matrix(1-k*iLAI, nrow=nl, ncol=nwl)  # REPLACE when LIDF profile ready.
  #tau_ss <- matrix(rep(1-k*LAI, nl), ncol=13, byrow=TRUE)   # direct-direct transmittance for the thin layers
  #tau_ss <- matrix(1 - k * iLAI, nrow = nl, ncol = length(k), byrow = TRUE)
  # Create tau_ss using matrix replication
  #tau_ss <- matrix(rep(1 - k * iLAI, nl), nrow = nl, ncol = length(k), byrow = TRUE)
  tau_dd <- (1-a * iLAI)   # diffuse-diffuse transmittance for the thin layers
  tau_sd <- sf*iLAI       # direct-diffuse transmittance for the thin layers
  rho_sd <- sb*iLAI       # direct-diffuse reflectance for the thin layers
  rho_dd <- sigb*iLAI     # diffuse-diffuse reflectance for the thin layers

  ### Get reflectance from n layers for the reg1 region (400:2500)
  refl <- get.reflectances(tau_ss, tau_sd, tau_dd, rho_dd, rho_sd, rsoil, nl, nwl=2162)

  ## adding the thermal reflectance for
  R_sd <- refl[['R_sd']]
  R_dd <- refl[['R_dd']]
  Xss <- refl[['Xss']][1]
  Xsd <- refl[['Xsd']] #
  Xdd <- refl[['Xdd']]
  ### Here the 4-stream reflectances rsd and rdd
  rdd <- R_dd[1,]        # TOC hemispherical-hemispherical reflectance
  rsd <- R_sd[1,]       # TOC directional-hemispherical reflectance


  ## give errors due dim of rsoil is up to 2400 nm
  #rfl.fourSA <- m4SAIL(inputLUT=inputLUT.nsamples[100,],rsoil=rsoil, LeafModel='PRO')
  #data.foursail_pro<-ToolsRTM::Compute_BRF(rdot=rfl.fourSA$rdot,rsot=rfl.fourSA$rsot,tts=inputLUT.nsamples[100,'tts'],data.light =ToolsRTM::dataSpec_PDB)

  #data.four = data.frame(wave= c(400:2500),rfl.4 = data.foursail_pro)

  if (get.plots ==  T){

    check.rfl <- data.frame(wave=data.spectral$wlS,rdd=rdd,rsd=rsd,Xsd=Xsd[1,])

    p.toc <- ggplot(data = check.rfl, aes(x = wave, y = rdd)) +
      labs(y= "TOC reflectance", x = "")+
      geom_line(color='black') + theme_bw() + xlim(400,2499) +
      geom_line(aes(x = wave, y = rsd),color='navyblue') +
      geom_line(aes(x = wave, y = Xsd),color='forestgreen') #+
    # geom_line(data=data.four,aes(x = wave, y = rfl.4),color='forestgreen')
    print(p.toc)
  }

  ## when model == 'fourSAIL' save some values
  #data.rad <- list() dont forget to remove this in section input at # 6.2. get rad outputs

  #data.rad[['rsd']] <- rsd# TOC directional-hemispherical reflectance
  #data.rad[['rdd']] <- rdd# TOC hemispherical-hemispherical reflectance


  ### using the complete spectra
  ## Get  direct solar radiance and diffuse sky radiance

  TOCirr <- get.calcTOCirr(atmo, data.meteo, rdd, rsd, wl, nwl = 2162)


  Esun_ <- TOCirr[['Esun_']] # the direct solar radiance at the top layer.
  Esky_ <- TOCirr[['Esky_']] #  the diffuse sky radiance at the top layer.

  if (get.plots ==  T){
    check.irrad <- data.frame(wave=data.spectral$wlS,Esun_=Esun_,Esky_=Esky_)

    p.ir <- ggplot(data = check.irrad, aes(x = wave)) +
      labs(y= "Irradiance", x = "") +
      geom_line(aes(y = Esun_, color = "Direct solar"), linewidth = 0.5) +
      geom_line(aes(y = Esky_, color = "Diffuse sky"), linewidth = 0.5) +
      scale_color_manual(name = "Irradiance type",
                         values = c("Direct solar" = "black", "Diffuse sky" = "navyblue")) +
      theme_bw() + xlim(400,3000) +
      theme(legend.position = "top") +
      guides(color = guide_legend(title = NULL))
    print(p.ir)
  }


  ### Estimate the radiation flux profiles at different stages of the calculation

  # Es_ represents the direct solar irradiance at each layer of the atmosphere.
  # Emin_ represents the downward diffuse irradiance at each layer of the atmosphere.
  # Eplu_ represents the upward diffuse irradiance at each layer of the atmosphere.

  Eflux1 <- get.fluxprofile(Esun_, 0*Esky_, rsoil=rsoil, Xss, Xsd, Xdd, R_sd, R_dd, nl, nwl,rs.thermal=0.06)
  #Emins_ and Emind_ are the downward diffuse irradiance flux profiles for the two sets of input parameters.

  Emins_ <- Eflux1[['Emin_']]  #Emins_ is the downward diffuse irradiance flux profiles for the sets of input parameters.

  Eplus_ <- Eflux1[['Eplu_']] #Eplus_ is the upward diffuse irradiance flux profiles for the sets of input parameters.


  Eflux2 <- get.fluxprofile(0*Esun_, Esky_, rsoil=rsoil, Xss, Xsd, Xdd, R_sd, R_dd, nl, nwl,rs.thermal=0.06)


  Emind_ <- Eflux2[['Emin_']] #Emind_ is the downward diffuse irradiance flux profiles for the sets of input parameters.

  Eplud_ <- Eflux2[['Eplu_']] #Eplud_ is the upward diffuse irradiance flux profiles for the sets of input parameters.


  #Emin_ and Eplu_ are the net downward and upward diffuse irradiance flux profiles, respectively.

  Emin_ <- Emins_ + Emind_
  Eplu_ <- Eplus_ + Eplud_



  if (get.plots ==  T){
    check.emi <- data.frame(wave=data.spectral$wlS,Emin_1=Emin_[1,],Emin_5=Emin_[5,],Eplu_1=Eplu_[1,],Eplu_5=Eplu_[5,])


    p.emi.1 <- ggplot(data = check.emi, aes(x = wave)) +
      labs(y= "downward diffuse irradiance flux ", x = "")+
      geom_line(aes(y = Emin_1, color = "canopy layer 1"), linewidth = 0.5) +
      geom_line(aes(y = Emin_5, color = "canopy layer 5"), linewidth = 0.5) +
      scale_color_manual(name = "Irradiance type",
                         values = c("canopy layer 1" = "forestgreen", "canopy layer 5" = "navyblue")) +
      theme_bw() + xlim(400,800) +
      theme(legend.position = "top") +
      guides(color = guide_legend(title = NULL))
    print(p.emi.1)

    p.emi.2 <- ggplot(data = check.emi, aes(x = wave)) +
      labs(y= "upward diffuse irradiance flux ", x = "")+
      geom_line(aes(y = Eplu_1, color = "canopy layer 1"), linewidth = 0.5) +
      geom_line(aes(y = Eplu_5, color = "canopy layer 5"), linewidth = 0.5) +
      scale_color_manual(name = "Irradiance type",
                         values = c("canopy layer 1" = "forestgreen", "canopy layer 5" = "navyblue")) +
      theme_bw() + xlim(400,800) +
      theme(legend.position = "top") +
      guides(color = guide_legend(title = NULL))
    print(p.emi.2)

  }





  ########################################
  # 1.5 probabilities Ps, Po, Pso
  ########################################

  Ps <- exp(k[1] * xl * LAI) #probability of viewing a leaf in solar dir
  Po <- exp(K[1] * xl * LAI)
  Ps[1:nl] <- Ps[1:nl] * (1 - exp(-k[1] * LAI * dx)) / (k[1] * LAI * dx) #Correct Ps/Po for finite dx
  Po[1:nl] <- Po[1:nl] * (1 - exp(-K[1] * LAI * dx)) / (K[1] * LAI * dx) #Correct Ps/Po for finite dx


  q <- data.canopy[['hot']]

  #Pso is a factor that represents the correlation between  PS and Po
  #where Ps is the probability of viewing a leaf in the solar directionand
  #Po is the probability of viewing a leaf in the opposite direction.

  Pso <- matrix(0, nrow = length(Po), ncol = 1)

  for (j in 1:length(xl)) {

    Pso[j, ] <- integrate(Vectorize(function(y) get.Pso(K, k, LAI, q, dso, y)),
                          lower = xl[j] - dx, upper = xl[j])$value / dx
  }

  Pso[Pso > Po] <- apply(cbind(Po[Pso > Po], Ps[Pso > Po]), 1, min)
  Pso[Pso > Ps] <- apply(cbind(Po[Pso > Ps], Ps[Pso > Ps]), 1, min)

  data.gap <- list()
  data.gap[['Pso']] <- Pso


  # 3.3 outgoing fluxes, hemispherical and in viewing direction, spectrum
  # in viewing direction, spectral due to diffuse light

  # vegetation contribution
  if (dim(vb)[2] == 2001){

    rs_thermal <- data.soil[['rs_thermal']]###  before I use this value 0.01
    rs_vec <- 2162 - 2001

    r_thermal <-matrix(rs_thermal, ncol = rs_vec, nrow= nrow(vb))
    vb_rs_thermal <- cbind(vb,r_thermal)
    vf_rs_thermal <- cbind(vf,r_thermal)
    w_rs_thermal <- cbind(w,r_thermal)

    epsc_rs_thermal <- cbind(epsc,r_thermal)
    epss_rs_thermal <- c(epss,r_thermal[1,])
    rsoil.thermal = c(rsoil, rep(rs_thermal,rs_vec))
  } else {

    vb_rs_thermal <- vb
    vf_rs_thermal <- vf
    w_rs_thermal <- w

    epsc_rs_thermal <- epsc
    epss_rs_thermal <- epss
    rsoil.thermal = rsoil
  }


  #####   ######   ######   ######   ######   ######   ######   ######
  #piLocd_ is the calculated flux density from the vegetation layer (in W/m^2)
  # Vegetation contribution

  vb_Po_Emin_ <- vb_rs_thermal * Po[1:nl] * Emind_[1:nl, ]


  vf_Po_Eplud_ <- vf_rs_thermal * Po[1:nl] * Eplud_[1:nl, ]
  #Fixed error here, iLAI multiply the sum final
  #piLocd_ <- colSums(vb_Po_Emin_) * iLAI  + colSums(vf_Po_Eplud_) * iLAI
  # give same
  piLocd_ <- (colSums(vb_Po_Emin_) + colSums(vf_Po_Eplud_)) * iLAI

  # soil contribution
  # piLosd_ is the downward flux density (W m^-2 nm^-1) received by the soil surface
  # due to diffuse radiation reflected downwards by the canopy

  # nl + 1 = 60 return zero, for this reason I choose nl=59

  Emin_Po <-Emind_[nl,] * Po[nl]

  piLosd_ <- rsoil.thermal * Emin_Po

  if (get.plots ==  T){
    check.flux.d <- data.frame(wave=data.spectral$wlS,piLosd_=piLosd_, piLocd_=piLocd_)


    p.flux.piLosd <- ggplot(data = check.flux.d, aes(x = wave)) +
      labs(y= "downward flux density (W m^-2 nm^-1) ", x = "")+
      geom_line(aes(y = piLosd_, ),color = "forestgreen", linewidth = 0.5)  +
      theme_bw() + xlim(400,2500)

    print(p.flux.piLosd)


    p.flux.piLocsd <- ggplot(data = check.flux.d, aes(x = wave)) +
      labs(y= "flux density from the vegetation layer (in W/m^2)", x = "")+
      geom_line(aes(y = piLocd_, ),color = "forestgreen", linewidth = 0.5)  +
      theme_bw() + xlim(400,2500)

    print(p.flux.piLocsd)

  }

  # in viewing direction, spectral due to direct solar light
  # vegetation contribution

  vb_Po_Emins_<- vb_rs_thermal * Po[1:nl] * Emins_[1:nl, ]
  vf_Po_Eplus_ <- vf_rs_thermal * Po[1:nl] * Eplus_[1:nl, ]
  w_Po_Esun  <- colSums(w_rs_thermal * Pso[1:nl]) * Esun_


  #Fixed error here, iLAI multiply the sum final
  #piLocu_ <-  colSums(vb_Po_Emins_) * iLAI +colSums(vf_Po_Eplus_) * iLAI
  piLocu_ <-  (colSums(vb_Po_Emins_) + colSums(vf_Po_Eplus_) + w_Po_Esun) * iLAI

  # soil contribution
  Emins_Po_Esun_ <-Emins_[nl,] * Po[nl] +  Esun_*Pso[nl]
  piLosu_ <- rsoil.thermal*(Emins_Po_Esun_)



  piLod_ <- piLocd_ + piLosd_ # [nwl] piRad in obsdir from Esky
  piLou_ <- piLocu_ + piLosu_ # [nwl] piRad in obsdir from Eskun
  piLoc_ <- piLocu_ + piLocd_ # [nwl] piRad in obsdir from vegetation
  piLos_ <- piLosu_ + piLosd_ # [nwl] piRad in obsdir from soil


  piLo_ <- piLoc_ + piLos_ # [nwl] piRad in obsdir
  Lo_ <- piLo_/pi  # [nwl] Rad in obsdir

  if (get.plots ==  T){

    check.PiLo_ <- data.frame(wave=data.spectral$wlS,piLod_=piLod_, piLou_=piLou_,
                              piLoc_=piLoc_,piLos_=piLos_, Lo_=Lo_)


    p.piLo_ <- ggplot(data = check.PiLo_, aes(x = wave)) +
      labs(y= " piRad in obsdir (Esky,Esksun,veg and soil)", x = "")+
      geom_line(aes(y = piLod_, ),color = "grey", linewidth = 0.5)  +
      geom_line(aes(y = piLou_, ),color = "navyblue", linewidth = 0.5)  +
      geom_line(aes(y = piLoc_, ),color = "forestgreen", linewidth = 0.5)  +
      geom_line(aes(y = piLos_, ),color = "brown", linewidth = 0.5)  +
      theme_bw() + xlim(400,2500)

    print(p.piLo_)

    p.Lo_ <- ggplot(data = check.PiLo_, aes(x = wave)) +
      labs(y= "Rad in obsdir", x = "")+
      geom_line(aes(y = Lo_, ),color = "navyblue", linewidth = 0.5)  +

      theme_bw() + xlim(400,2500)

    print(p.Lo_)

  }



  ##  4-stream reflectances rso, rdo
  #rso is TOC directional-directional reflectance
  rso <- piLou_/Esun_# [nwl] obsdir reflectance of solar beam
  #plot(rso)
  #rdo is TOC hemispherical-directional reflectance
  rdo <- piLod_/Esky_ # [nlw] obsir reflectance of sky irradiance
  #plot(rdo)
  #Refl is TOC reflectance
  Refl <- piLo_/(Esky_+Esun_) # [nwl]
  #plot(Refl)
  # prevents numerical instability in absorption windows
  Refl[Esky_ < 1E-4] <- rso[Esky_ < 1E-4]
  I <- which(Esky_ < 2E-4 * max(Esky_))

  Refl[I] <- rso[I] # prevents numerical instability in absorption windows
  #plot(Refl)

  if (get.plots ==  T){
    check.rfl_ <- data.frame(wave=data.spectral$wlS,rfl=Refl, rdo=rdo,
                             rso=rso)


    p.rfl_ <- ggplot(data = check.rfl_, aes(x = wave)) +
      labs(y= "reflectance", x = "")+
      geom_line(aes(y = rfl, ),color = "forestgreen", linewidth = 0.5)  +
      geom_line(aes(y = rdo, ),color = "navyblue", linewidth = 0.5)  +
      geom_line(aes(y = rso, ),color = "brown", linewidth = 0.5)  +
      theme_bw() + xlim(400,2500)

    print(p.rfl_)
  }



  # 4. net fluxes, spectral and total, and incoming fluxes
  # 4.1 incident PAR at the top of canopy, spectral and spectrally integrated
  P_ <- get.e2phot(lambda=wl[Ipar]*1E-9,E=(Esun_[Ipar]+Esky_[Ipar]),constants)
  P <- 0.001*Sint(P_,wlPAR) # mol m-2s-1


  EPAR_ <- Esun_[Ipar] + Esky_[Ipar]
  EPAR <- 0.001*Sint(EPAR_, wlPAR)

  if (get.plots ==  T){

    check.par <- data.frame(wave=wl[Ipar],P_=P_,
                            EPAR_=EPAR_)

    p.epar_ <- ggplot(data = check.par, aes(x = wave)) +
      labs(y= "incident PAR", x = "")+
      geom_line(aes(y = EPAR_, ),color = "forestgreen", linewidth = 0.5)  +
      theme_bw() + xlim(400,700)

    print(p.epar_)

    p.photons <- ggplot(data = check.par, aes(x = wave)) +
      labs(y= "photons energy", x = "")+
      geom_line(aes(y = P_, ),color = "orange", linewidth = 0.5)  +
      theme_bw() + xlim(400,700)

    print(p.photons)
  }

  #Psun <- 0.001*Sint(e2phot(wlPAR*1E-9,Esun_[Ipar],constants),wlPAR) # Incident solar PAR in PAR units
  # Incident and absorbed solar radiation


  ###############################################################
  # 4.2 Absorbed radiation
  #    absorbed radiation in Wm-2         (Asun)
  #    absorbed PAR in mol m-2s-1         (Pnsun)
  #    absorbed PAR in Wm-2               (Rnsun_PAR)
  #    absorbed PAR by Chl in mol m-2s-1  (Pnsun_Cab)
  ###############################################################


  # initialize variables
  Asun <-  rep(0, nl)
  Pnsun <- rep(0, nl)
  Rnsun_PAR <- rep(0, nl)
  Pnsun_Cab <- rep(0, nl)
  Rnsun_Cab <- rep(0, nl)
  Pnsun_Car <- rep(0, nl)
  Rnsun_Car <- rep(0, nl)


  #########



  for (j in 1:nl) {
    #as a list element
    total_ <-Sint(Esun_ * epsc_rs_thermal[j,], wl)
    Asun[j] <- total_ * 0.001  # Total absorbed solar radiation
    ## here change 0.001 by 1
    factor_ <- 0.001
    Pnsun[j] <- factor_ * Sint(get.e2phot(wlPAR * 1E-9, Esun_[Ipar] * epsc[j, Ipar], constants), wlPAR)  # Absorbed solar radiation in PAR range in moles m-2 s-1

    Rnsun_PAR[j] <- factor_ * Sint(Esun_[Ipar] * epsc[j, Ipar], wlPAR)

    Rnsun_Cab[j] <- factor_ * Sint(Esun_[data.spectral[['IwlP']]] * epsc[j, data.spectral[['IwlP']]] * kChlrel[j,data.spectral[['IwlP']]], data.spectral[['wlP']])
    Pnsun_Cab[j] <- factor_ * Sint(get.e2phot(data.spectral[['wlP']] * 1E-9, kChlrel[j,data.spectral[['IwlP']]] * Esun_[data.spectral[['IwlP']]] * epsc[j, data.spectral[['IwlP']]], constants), data.spectral[['wlP']])

    Rnsun_Car[j] <- factor_ * Sint(Esun_[data.spectral[['IwlP']]] * epsc[j, data.spectral[['IwlP']]] * kCarrel[j,data.spectral[['IwlP']]], data.spectral[['wlP']])
    Pnsun_Car[j] <- factor_ * Sint(get.e2phot(data.spectral[['wlP']] * 1E-9, kCarrel[j,data.spectral[['IwlP']]] * Esun_[data.spectral[['IwlP']]] * epsc[j, data.spectral[['IwlP']]], constants), data.spectral[['wlP']])

  }


  #4.3 total direct radiation (incident and net) per leaf area (W m-2 leaf)
  # total direct radiation (incident and net) per leaf area (W m-2 leaf)
  #Pdir = fs * Psun                        # [13 x 36]   incident


  if (options.lite$Value == 1) {

    # where j is nl
    #fs <- abs(cos_deltas / cos_tts)
    #fs <- matrix(lidf,ncol=1) %*% apply(fs, 2, mean)
    fs_ <- crossprod(lidf, rowMeans(fs))

    #Rndir <- Pndir <- Pndir_Cab <- Rndir_Cab <- Rndir_PAR <- Pndir_Car <- Rndir_Car <- array(0, dim = c(13, 36, nl))
    Rndir <- fs_ * Asun[j]
    Pndir <- fs_ * Pnsun[j]

    Pndir_Cab <- fs_ * Pnsun_Cab[j]
    Rndir_Cab <- fs_ * Rnsun_Cab[j]

    Pndir_Car <- fs_ * Pnsun_Car[j]
    Rndir_Car <- fs_ * Rnsun_Car[j]

    Rndir_PAR <- fs_ * Rnsun_PAR[j]

  } else {
    ## get an array with dims (13,36,nl)
    Rndir <- Pndir <- Pndir_Cab <- Rndir_Cab <- Rndir_PAR <- Pndir_Car <- Rndir_Car <- array(0, dim = c(13, 36, nl))
    for (j in 1:nl) {
      Rndir[, , j] <- fs * Asun[j]
      Pndir[, , j] <- fs * Pnsun[j]
      Pndir_Cab[, , j] <- fs * Pnsun_Cab[j]
      Rndir_Cab[, , j] <- fs * Rnsun_Cab[j]
      Pndir_Car[, , j] <- fs * Pnsun_Car[j]
      Rndir_Car[, , j] <- fs * Rnsun_Car[j]
      Rndir_PAR[, , j] <- fs * Rnsun_PAR[j]
    }
  }


  # 4.4 total diffuse radiation (net) per leaf area (W m-2 leaf)

  # Initializations  variable for each layer

  E_ <- list() # diffuse incident radiation for each layer 1:nl (59)



  #### 1 is top nl is bottom

  for (j in 1:nl) {
    # diffuse incident radiation for the present layer 'j' (mW m-2 um-1)
    E_[[j]] <- 0.5 * (Emin_[j,] + Emin_[j+1,] + Eplu_[j,] + Eplu_[j+1,])

    # incident PAR flux, integrated over all wavelengths (moles m-2 s-1)
    Pdif[j] <- 0.001 * Sint(get.e2phot(lambda=wlPAR * 1E-9, E_[[j]][Ipar], constants), wlPAR)  # [nl], including conversion mW >> W

    # net radiation (mW m-2 um-1) and net PAR (moles m-2 s-1 um-1), per wavelength
    Rndif_[j,] <- E_[[j]] * epsc_rs_thermal[j,]  # [nl,nwl] Net (absorbed) radiation by leaves

    #  Net (absorbed) as PAR photons
    Pndif_[j,] <- 0.001 * (get.e2phot(lambda=wlPAR * 1E-9, Rndif_[j, Ipar], constants))  # [nl,nwl] Net (absorbed) as PAR photons

    Rndif_Cab_[j,] <- (kChlrel[j,data.spectral[['IwlP']]] * Rndif_[j, data.spectral$IwlP])  # [nl,nwl] Net (absorbed) as PAR photons by Cab
    Pndif_Cab_[j,] <- 0.001 * (get.e2phot(data.spectral$wlP * 1E-9, (kChlrel[j,data.spectral[['IwlP']]] * Rndif_[j, data.spectral$IwlP]), constants))  # [nl,nwl] Net (absorbed) as PAR photons by Cab

    Rndif_Car_[j,] <- (kCarrel[j,data.spectral[['IwlP']]] * Rndif_[j, data.spectral$IwlP])  # [nl,nwl] Net (absorbed) as PAR photons by Car
    Pndif_Car_[j,] <- 0.001 * (get.e2phot(data.spectral$wlP * 1E-9, (kCarrel[j,data.spectral[['IwlP']]] * Rndif_[j, data.spectral$IwlP]), constants))  # [nl,nwl] Net (absorbed) as PAR photons by Car

    #Net (absorbed) as PAR energy
    Rndif_PAR_[j,] <- Rndif_[j, Ipar]  # [nl,nwlPAR] Net (absorbed) as PAR energy

    # net radiation (W m-2) and net PAR (moles m-2 s-1), integrated over all wavelengths
    ### Integral of values

    Pndif[j] <- Sint(Pndif_[j, Ipar], wlPAR)  # [nl] Absorbed PAR
    Rndif[j] <- 0.001 * Sint(Rndif_[j,], wl)  # [nl] Full spectrum net diffuse flux

    Pndif_Cab[j] <- Sint(Pndif_Cab_[j,], data.spectral$wlP)  # [nl] Absorbed PAR by Cab integrated
    Rndif_Cab[j] <- 0.001 * Sint(Rndif_Cab_[j,], data.spectral$wlP)  # [nl] Absorbed PAR by Cab integrated

    Pndif_Car[j] <- Sint(Pndif_Car_[j,], data.spectral$wlP) #Absorbed PAR by Car integrated
    Rndif_Car[j] <- 0.001* Sint(Rndif_Car_[j,], data.spectral$wlP) #Absorbed PAR by Car integrated

    Rndif_PAR[j] <- 0.001 * Sint(Rndif_PAR_[j,Ipar],wlPAR)    # [nl]  Absorbed PAR by Cab integrated


  }


  # soil layer, direct and diffuse radiation

  #Absorbed solar flux by the soil
  Rndirsoil <- 0.001 * Sint(Esun_* epss_rs_thermal, wl)
  # Absorbed diffuse downward flux by the soil (W m-2)
  Rndifsoil <- 0.001 * Sint(Emin_[nl+1,] * t(epss_rs_thermal),wl)

  # net (n) radiation R and net PAR P per component:
  # - sunlit (u)
  # - shaded (h)
  # - soil(s)
  # - canopy (c)
  # Units: W m-2 leaf or soil surface um-1

  Rnhc <- Rndif            # [nl] shaded leaves or needles
  Pnhc <- Pndif            # [nl] shaded leaves or needles
  Pnhc_Cab <- Pndif_Cab    # [nl] shaded leaves or needles
  Rnhc_Cab <- Rndif_Cab    # [nl] shaded leaves or needles
  Pnhc_Car <- Pndif_Car    # [nl] shaded leaves or needles
  Rnhc_Car <- Rndif_Car    # [nl] shaded leaves or needles
  Rnhc_PAR <- Rndif_PAR    # [nl] shaded leaves or needles

  if (options.lite$Value != 1) {
    #Normal SCOPE execution with [13 x 36 x nlayers] sunlit leaves
    Rnuc <- array(0, dim = c(13, 36, nl))
    Pnuc <- array(0, dim = c(13, 36, nl))
    Pnuc_Cab <- array(0, dim = c(13, 36, nl))
    Rnuc_PAR <- array(0, dim = c(13, 36, nl))
    Rnuc_Cab <- array(0, dim = c(13, 36, nl))
    Rnuc_Car <- array(0, dim = c(13, 36, nl))
    Pnuc_Car <- array(0, dim = c(13, 36, nl))


    for (j in 1:nl) {
      #Puc[, , j] <- Pdir[, , j] + Pdif[j] # [13,36,nl] Total fluxes on sunlit leaves or needles
      Rnuc[, , j] <- Rndir[, , j] + Rndif[j] # [13,36,nl] Total fluxes on sunlit leaves or needles
      Pnuc[, , j] <- Pndir[, , j] + Pndif[j] # [13,36,nl] Total fluxes on sunlit leaves or needles
      Pnuc_Cab[, , j] <- Pndir_Cab[, , j] + Pndif_Cab[j] # [13,36,nl] Total fluxes on sunlit leaves or needles
      Pnuc_Car[, , j] <- Pndir_Car[, , j] + Pndif_Car[j] # [13,36,nl] Total fluxes on sunlit leaves or needles
      Rnuc_PAR[, , j] <- Rndir_PAR[, , j] + Rndif_PAR[j] # [13,36,nl] Total fluxes on sunlit leaves or needles
      Rnuc_Cab[, , j] <- Rndir_Cab[, , j] + Rndif_Cab[j] # [13,36,nl] Total fluxes on sunlit leaves or needles
      Rnuc_Car[, , j] <- Rndir_Car[, , j] + Rndif_Car[j] # [13,36,nl] Total fluxes on sunlit leaves or needles
    }


  } else {

    #Lite SCOPE execution with [nlayers x 1] sunlit leaves, sunlit leaf inclinations are not accounted for

    Rnuc <- Rndir[1] + Rndif           # [nl] Total fluxes on sunlit leaves or needles
    Pnuc <- Pndir[1] + Pndif          # [nl] Total fluxes on sunlit leaves or needles
    Pnuc_Cab <- Pndir_Cab[1] + Pndif_Cab # [nl] Total fluxes on sunlit leaves or needles
    Rnuc_Cab <- Rndir_Cab[1] + Rndif_Cab # [nl] Total fluxes on sunlit leaves or needles
    Pnuc_Car <- Pndir_Car[1] + Pndif_Car # [nl] Total fluxes on sunlit leaves or needles
    Rnuc_Car <- Rndir_Car[1] + Rndif_Car # [nl] Total fluxes on sunlit leaves or needles
    Rnuc_PAR   <- Rndir_PAR[1] + Rndif_PAR # [nl] Total fluxes on sunlit leaves or needles


    if (get.plots ==  T){

      check.RN.vertical <- data.frame(pvertical=1:nl,Rnuc=Rnuc, Rnhc=Rnhc)


      p.vertical <-ggplot(data = check.RN.vertical, aes(y = pvertical)) +
        labs(y= "Net radiation ", x = "") +
        geom_line(aes(x = Rnuc, color = "Rnuc"), linewidth = 0.5)  +
        geom_line(aes(x = Rnhc, color = "Rnhc"), linewidth = 0.5)  +
        scale_color_manual(name = "Net radiation",
                           values = c("Rnuc" = "gold3",
                                      "Rnhc" = "black")) +
        theme_bw() +
        theme(legend.position = "top") +
        guides(color = guide_legend(title = NULL))

      print(p.vertical)
    }

  }

  # sunlit soil (1)
  Rnus <- Rndifsoil + Rndirsoil
  # shaded soil (1)
  Rnhs <- Rndifsoil


  if (options.calc_vert_profiles$Value == 1  ) {

    if (options.lite$Value != 1){

      Pnu1d <- meanleaf(data.canopy, F_= Pnuc, canopy.choice='angles',Ps=Ps)       # [nli, nlo, nl] mean net radiation sunlit leaves
      Pnu1d_Cab <- meanleaf(data.canopy, Pnuc_Cab, canopy.choice='angles',Ps=Ps)     # [nli, nlo, nl] mean net radiation sunlit leaves
      Pnu1d_Car <- meanleaf(data.canopy, Pnuc_Car, canopy.choice='angles',Ps=Ps)   # [nli, nlo, nl] mean net radiation sunlit leaves
      profiles <- list()
      profiles[['Pn1d']] <- ((1-Ps[1:nl]) * Pnhc + Ps[1:nl] * Pnu1d)     # [nl] mean photos leaves, per layer
      profiles[['Pn1d_Cab']] <- ((1-Ps[1:nl]) * Pnhc_Cab + Ps[1:nl] * Pnu1d_Cab) # [nl] mean photos leaves, per layer
      profiles[['Pn1d_Car']] <- ((1-Ps[1:nl]) * Pnhc_Car + Ps[1:nl] * Pnu1d_Car) # [nl] mean photos leaves, per layer

      if (get.plots ==  T){
        check.profile.vertical <- data.frame(pvertical=1:nl,
                                             Pn1d= profiles[['Pn1d']],
                                             Pn1d_CabP= profiles[['Pn1d_Cab']],
                                             Pn1d_Car= profiles[['Pn1d_Car']])


        p.vertical.2 <-ggplot(data = check.profile.vertical, aes(y = pvertical)) +
          labs(y= "net radiation sunlit leaves", x = "") +
          geom_line(aes(x = Pn1d, color = "Pn1d"), linewidth = 0.5)  +
          geom_line(aes(x = Pn1d_CabP, color = "Pn1d_CabP"), linewidth = 0.5)  +
          scale_color_manual(name = "Net radiation",
                             values = c("Pn1d" = "gold3",
                                        "Pn1d_CabP" = "black")) +
          theme_bw() +
          theme(legend.position = "top") +
          guides(color = guide_legend(title = NULL))

        print(p.vertical.2)
      }



    }
 } else {
    profiles <- list()
  }


  ####################################################################################
  # 5 Model output
  # up and down and hemispherical out, cumulative over wavelenght

  ####################################################################################

  Eout_ <- Eplu_[1,]
  # get hemispherical out, in optical range (W m-2) return a single value

  Eouto <- 0.001 * Sint(Eout_[data.spectral[['IwlP']]],data.spectral[['wlP']])
  # get hemispherical out, in thermal range (W m-2) return a single value
  Eoutt <- 0.001 * Sint(Eout_[data.spectral[['IwlT']]],data.spectral[['wlT']])
  # get hemispherical out, in thermal range (W m-2) return a single value

  Lot <- 0.001 * Sint(Lo_[data.spectral[['IwlT']]],data.spectral[['wlT']])

  ###################################################
  #### 6. get output list files
  ##################################################


  # 6.1. get gap outputs

  data.gap[['k']] <- k  # extinction cofficient in the solar direction
  data.gap[['K']] <- K  # extinction cofficient in the viewing direction
  data.gap[['Ps']] <- Ps # gap fraction in the solar direction
  data.gap[['Po']] <- Po # gap fraction in the viewing direction


  # 6.2. get rad outputs
  data.rad <- list()
  # 6.2.1. get TOC reflectance

  data.rad[['rsd']] <- rsd# TOC directional-hemispherical reflectance
  data.rad[['rdd']] <- rdd# TOC hemispherical-hemispherical reflectance
  data.rad[['rdo']] <- rdo# TOC hemispherical-directional reflectance
  data.rad[['rso']]  <- rso# TOC directional-directional reflectance
  data.rad[['refl']]  <- Refl# TOC reflectance

  # 6.2.2. get diffuse transmittance

  data.rad[['rho_dd']] <- rho_dd   # diffuse-diffuse reflectance for the thin layers
  data.rad[['tau_dd']] <- tau_dd   # diffuse-diffuse transmittance for the thin layers
  data.rad[['rho_sd']] <- rho_sd   # direct-diffuse reflectance for the thin layers
  data.rad[['tau_ss']] <- tau_ss   # direct-direct transmittance for the thin layers
  data.rad[['tau_sd']] <- tau_sd   # direct-diffuse transmittance for the thin layers

  data.rad[['R_sd']] <- R_sd ## adding and extra layer  (nl=60)
  data.rad[['R_dd']] <- R_dd ## adding and extra layer  (nl=60)

  data.rad[['Xdd']] <-   Xdd
  data.rad[['Xsd']] <- Xsd
  data.rad[['Xss']] <- Xss


  # 6.2.3. get coefficient for diffuse incidence

  data.rad[['vb']] <- vb# directional backscatter coefficient for diffuse incidence
  data.rad[['vf']] <- vf# directional forward scatter coefficient for diffuse incidence

  # 6.2.4. get coefficient for specular flux

  data.rad[['sigf']] <- sigf# forward scatter coefficient for specular flux
  data.rad[['sigb']] <- sigb# backscatter coefficient for specular flux

  # 6.2.5. get incident solar and sky spectrum

  data.rad[['Esun_']] <- Esun_    # [nwlx1 double]   incident solar spectrum (mW m-2 um-1)
  data.rad[['Esky_']] <- Esky_    # [nwlx1 double]   incident sky spectrum (mW m-2 um-1)

  # 6.2.6. get incident PAR

  data.rad[['PAR']] <- P * 1E6    # [1 double]       incident spectrally integrated PAR (micromoles m-2 s-1)
  data.rad[['EPAR']]    <- EPAR# [1 double]       incident PAR in energy units (W m-2)

  # 6.2.7. get diffuse radiation in the canopy to sevral directions

  data.rad[['Eplu_']] <- Eplu_    # [nlxnwl double]  upward diffuse radiation in the canopy (mW m-2 um-1)
  data.rad[['Emin_']] <- Emin_    # [nlxnwl double]  downward diffuse radiation in the canopy (mW m-2 um-1)
  data.rad[['Emins_']] <- Emins_   # [nlxnwl double]  downward diffuse radiation in the canopy due to direct solar rad (mW m-2 um-1)
  data.rad[['Emind_']] <- Emind_   # [nlxnwl double]  downward diffuse radiation in the canopy due to sky rad (mW m-2 um-1)
  data.rad[['Eplus_']] <- Eplus_   # [nlxnwl double]  upward diffuse radiation in the canopy due to direct solar rad (mW m-2 um-1)
  data.rad[['Eplud_']] <- Eplud_   # [nlxnwl double]  upward diffuse radiation in the canopy due to sky rad (mW m-2 um-1)

  # 6.2.7. get the TOC radiance in observation direction

  data.rad[['Lo_']] <- Lo_# [nwlx1 double]   TOC radiance in observation direction (mW m-2 um-1 sr-1)

  # 6.2.8. get the  TOC radiations

  data.rad[['Eout_']] <- Eout_    # [nwlx1 double]   TOC upward radiation (mW m-2 um-1)
  data.rad[['Eouto']] <- Eouto    # [1 double]        TOC spectrally integrated upward optical ratiation (W m-2)
  data.rad[['Eoutt']] <- Eoutt    # [1 double]        TOC spectrally integrated upward thermal ratiation (W m-2)
  data.rad[['Lot']] <- Lot

  # 6.2.9. getthe net radiation (shaded- sunlit leaves)

  data.rad[['Rnhs']] <- Rnhs# [1 double]        net radiation (W m-2) of shaded soil
  data.rad[['Rnus']] <- Rnus# [1 double]        net radiation (W m-2) of sunlit soil
  data.rad[['Rnhc']] <- Rnhc# [60x1 double]     net radiation (W m-2) of shaded leaves
  data.rad[['Rnuc']] <- Rnuc# [13x36x60 double] net radiation (W m-2) of sunlit leaves

  # 6.2.10. get net PAR  shaded and sunit leaves

  data.rad[['Pnh']] <- 1E6*Pnhc# [60x1 double]     net PAR (moles m-2 s-1) of shaded leaves
  data.rad[['Pnu']] <- 1E6*Pnuc# [13x36x60 double] net PAR (moles m-2 s-1) of sunlit leaves
  data.rad[['Pnh_Cab']] <- 1E6*Pnhc_Cab# [60x1 double]      net PAR absorbed by Cab (moles m-2 s-1) of shaded leaves
  data.rad[['Pnu_Cab']] <- 1E6*Pnuc_Cab # [13x36x60 double] net PAR absorbed by Cab (moles m-2 s-1) of sunlit leaves
  data.rad[['Rnh_Cab']] <- Rnhc_Cab # [60x1 double]    net PAR absorbed by Cab (W m-2) of shaded leaves
  data.rad[['Rnu_Cab']] <- Rnuc_Cab # [13x36x60 double] net PAR absorbed by Cab (W m-2) of sunlit leaves
  data.rad[['Pnh_Car']] <- 1E6*Pnhc_Car# [60x1 double]      net PAR absorbed by Cab (moles m-2 s-1) of shaded leaves
  data.rad[['Pnu_Car']] <- 1E6*Pnuc_Car # [13x36x60 double] net PAR absorbed by Cab (moles m-2 s-1) of sunlit leaves
  data.rad[['Rnh_Car']] <- Rnhc_Car # [60x1 double]    net PAR absorbed by Cab (W m-2) of shaded leaves
  data.rad[['Rnu_Car']] <- Rnuc_Car # [13x36x60 double] net PAR absorbed by Cab (W m-2) of sunlit leaves
  data.rad[['Rnh_PAR']] <- Rnhc_PAR # [60x1 double]     net PAR absorbed by Cab (W m-2) of shaded leaves
  data.rad[['Rnu_PAR']] <- Rnuc_PAR # [13x36x60 double] net PAR absorbed (W m-2) of sunlit

  LRT<- list(data.rad=data.rad, data.gap=data.gap,
             data.profiles=profiles, data.canopy=data.canopy)
  return(LRT)

  # 7. get some plots

  if (get.plots ==  T){

    check.Rndif.top <- data.frame(wave=wl,E_=E_[[1]],Rndif_=Rndif_[1,])

    p.Rndif_ <-ggplot(data = check.Rndif.top, aes(x = wave)) +
      labs(y= "Net radiation ", x = "") +
      geom_line(aes(y = E_, color = "diffuse incident radiation"), linewidth = 0.5)  +
      geom_line(aes(y = Rndif_, color = "radiation absorbed by leaves"), linewidth = 0.5)  +
      scale_color_manual(name = "Net radiation",
                         values = c("radiation absorbed by leaves" = "forestgreen",
                                    "diffuse incident radiation" = "navyblue")) +
      theme_bw() + xlim(400,2500) +
      theme(legend.position = "right") +
      guides(color = guide_legend(title = NULL))

    print(p.Rndif_)


    check.pigms.top <- data.frame(wave=data.spectral[['wlP']],Pndif_Cab_=Pndif_Cab_[1,],
                                  Rndif_Cab_=Rndif_Cab_[1,],
                                  Pndif_Car_=Pndif_Car_[1,],
                                  Rndif_Car_=Rndif_Car_[1,])


    p.RP.cab <-ggplot(data = check.pigms.top, aes(x = wave)) +
      labs(y= "Net PAR at top layer", x = "") +
      geom_line(aes(y = Pndif_Cab_, color = "Chorophylls"), linewidth = 0.5)  +
      geom_line(aes(y = Pndif_Car_, color = "Carotenoids"), linewidth = 0.5)  +
      scale_color_manual(name = "Net radiation",
                         values = c("Chorophylls" = "forestgreen",
                                    "Carotenoids" = "navyblue")) +
      theme_bw() + xlim(400,800) +
      theme(legend.position = "right") +
      guides(color = guide_legend(title = NULL))

    print(p.RP.cab)
  }


} ## end getRTMo

