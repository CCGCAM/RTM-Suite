#' @title get.fluspect_mSCOPE
#'
#' @description
#' \code{get.fluspect_mSCOPE} is an adaptation of the Fluspect leaf model for
#' the multi-layer (mSCOPE) canopy representation: leaf optical properties
#' (reflectance, transmittance, fluorescence excitation-emission matrices
#' \code{Mb}/\code{Mf}, and pigment contribution factors) are computed once
#' per distinct leaf biochemistry profile layer (\code{mly$nly} layers), and
#' then replicated across the \code{nl} canopy layers that fall within each
#' profile layer, weighted by \code{mly$pLAI}.
#'
#' @param mly list. Multi-layer leaf biochemistry profile, with elements \code{nly} (number of distinct biochemistry layers), \code{pLAI} (LAI fraction of each layer), and per-layer biochemistry vectors \code{pCab}, \code{pEWT}, \code{pCar}, \code{pLMA}, \code{pCs}, \code{pN} (one value per layer, same units as the corresponding \code{leafbio} entries).
#' @param spectral list. Spectral configuration, including wavelength vectors \code{wlP}, \code{wlE}, \code{wlF}, \code{wlS} and region boundaries \code{reg1}/\code{reg2}/\code{reg3}/\code{IwlT}.
#' @param leafbio list. Baseline leaf biochemical/structural properties; the fields \code{Cab}, \code{EWT}, \code{Car}, \code{LMA}, \code{Cs}, \code{N} are overwritten per layer from \code{mly}.
#' @param soil list. Soil properties; only \code{rs_thermal} (thermal-region soil reflectance) is used, and only for the optional diagnostic plots.
#' @param optipar list. Leaf-level optical parameters used by the underlying Fluspect model.
#' @param nl integer. Number of canopy layers over which the per-profile-layer leaf optics are replicated.
#' @param step numeric, optional. Wavelength step (nm) used to compute the \code{Mb}/\code{Mf} excitation-emission matrices; if missing, a default step of 5 nm is used (matrix of 53 x 71), a step of 1 nm gives a 211 x 351 matrix. Must be lower than 8 nm.
#' @param get.plots logical. If \code{TRUE}, produces diagnostic plots of leaf reflectance, pigment contribution factors, PSI/PSII quantum yield fractions and the \code{Mb}/\code{Mf} matrices. Default \code{TRUE}.
#'
#' @return A list \code{leafopt} with: \code{refl}/\code{tran} (leaf reflectance/transmittance per canopy layer x wavelength, \code{[nl x length(wlP)]}), \code{kChlrel}/\code{kCarrel} (relative contribution of chlorophyll/carotenoids to absorption, same dimensions), \code{Mb}/\code{Mf} (backward/forward fluorescence excitation-emission matrices, replicated per canopy layer), \code{phiI}/\code{phiII} (PSI/PSII relative quantum yield spectra at the fluorescence wavelengths).
#' @export
#'
#' @author 	Christiaan van der Tol  (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' leafopt <- get.fluspect_mSCOPE(mly, spectral, leafbio, soil, optipar, nl,
#'                                 step = 5, get.plots = FALSE)
#' }
get.fluspect_mSCOPE <- function(mly, spectral, leafbio, soil, optipar, nl,step,get.plots=T) {


  if (missing(get.plots)){
    get.plots = FALSE
  }

  leafopt<-list()
  # leaf reflectance, transmittance and the excitation-fluorescence matrices calculation
  # for 60 sublayers
  indStar <- c(1, floor(cumsum(mly[['pLAI']]/sum(mly[['pLAI']])) * nl))  # index of starting for each different layer


  if (missing(step)){
    step_to_model <- 5
    ## not use  this print for get.SCOPE reasons
    #print(paste('fluspect model was running at ',step_to_model,' nm ',sep=''))
    wle <- spectral[['wlE']]  # excitation wavelengths, transpose to column
    wlf <- spectral[['wlF']]  # fluorescence wavelengths, transpose to column
  } else {
    step_to_model <- step

    if (step_to_model >= 8){

      stop('The fluspect model needs to run with a step lower than 8 nm')

    } else if (step_to_model == 1) {

      ## not use  this print for get.SCOPE reasons
      #print('fluspect model was running at 1 nm ')

      wle <- seq(min(spectral[['wlE']]),max(spectral[['wlE']]),step)
      wlf <- seq(min(spectral[['wlF']]),max(spectral[['wlF']]),step)

    } else {
      ## not use  this print for get.SCOPE reasons
      #print(paste('fluspect model was running at ',step_to_model,' nm ',sep=''))

      wle <- seq(min(spectral[['wlE']]),max(spectral[['wlE']]),step)
      wlf <- seq(min(spectral[['wlF']]),max(spectral[['wlF']]),step-1) # fluorescence wavelengths, transpose to column
    }

  }




  leafopt <- list(refl = matrix(nrow=mly$nly, ncol=length(spectral$wlP)),
                  tran = matrix(nrow=mly$nly, ncol=length(spectral$wlP)),
                  Mb = array(dim=c(length(wlf), length(wle), nl)),
                  Mf = array(dim=c(length(wlf), length(wle), nl)),
                  kChlrel = matrix(nrow=mly$nly, ncol=length(spectral$wlP)),
                  kCarrel = matrix(nrow=mly$nly, ncol=length(spectral$wlP)))



  rho_temp <- matrix(NA, nrow = nl, ncol = length(spectral$wlP))
  tau_temp <- matrix(NA, nrow = nl, ncol = length(spectral$wlP))
  kChlrel_temp <- matrix(NA, nrow = nl, ncol = length(spectral$wlP))
  kCarrel_temp <- matrix(NA, nrow = nl, ncol = length(spectral$wlP))

  for (i in 1:mly$nly) { ### before mly$nly , here nly =1
    leafbio$Cab <- mly$pCab[i]
    leafbio$EWT <- mly$pEWT[i]
    leafbio$Car <- mly$pCar[i]
    leafbio$LMA <- mly$pLMA[i]
    leafbio$Cs <- mly$pCs[i]
    leafbio$N <- mly$pN[i]

    leafopt_ml <- SCOPEinR::getFluspect.Cx.SCOPE(inputsLeaf=leafbio,inputsOptipar=SCOPEinR::optipar2021.Pro.CX,
                                                 version = 'SCOPE', step=step_to_model)


    #plot(leafopt_ml$refl,type='l',col='red')
    #plot(leafopt_ml$tran,type='l',col='navyblue')
    #plot(leafopt_ml$kChlrel,type='l',col='forestgreen')
    #plot(leafopt_ml$kCarrel,type='l',col='brown')
    #plot(leafopt_ml$Mb[,10],type='l',col='red')

    leafopt$refl[i,] <- leafopt_ml$refl
    leafopt$tran[i,] <- leafopt_ml$tran

    leafopt$kChlrel[i,] <- leafopt_ml$kChlrel
    leafopt$kCarrel[i,] <- leafopt_ml$kCarrel


    in1 <- indStar[i]
    in2 <- indStar[i+1]

    rho_temp[in1:in2,] <- matrix(rep(leafopt$refl[i,], in2-in1+1), ncol=length(spectral$wlP), byrow=TRUE)
    tau_temp[in1:in2,] <- matrix(rep(leafopt$tran[i,], in2-in1+1), ncol=length(spectral$wlP), byrow=TRUE)

    kChlrel_temp[in1:in2,] <- matrix(rep(leafopt$kChlrel[i,], in2-in1+1), ncol=length(spectral$wlP), byrow=TRUE)
    kCarrel_temp[in1:in2,] <- matrix(rep(leafopt$kCarrel[i,], in2-in1+1), ncol=length(spectral$wlP), byrow=TRUE)

    # replicate this profile layer's Mb/Mf across the canopy sublayers it spans (in1:in2);
    # previously only the last profile layer's Mb/Mf survived for the whole canopy
    n_rep <- in2 - in1 + 1
    leafopt$Mb[,,in1:in2] <- array(rep(leafopt_ml$Mb, n_rep), dim = c(dim(leafopt_ml$Mb)[1], dim(leafopt_ml$Mb)[2], n_rep))
    leafopt$Mf[,,in1:in2] <- array(rep(leafopt_ml$Mf, n_rep), dim = c(dim(leafopt_ml$Mf)[1], dim(leafopt_ml$Mf)[2], n_rep))
  }



  leafopt[['refl']]=rho_temp
  leafopt[['tran']]=tau_temp
  leafopt[['kChlrel']] <- kChlrel_temp
  leafopt[['kCarrel']] <-kCarrel_temp

  wlS <- spectral[['wlS']] # SCOPE wavelengths, make column vectors
  wlF <- spectral[['wlF']]

  iw_coincidents <- match(wlF, wlS)

  leafopt[['phiI']] <- SCOPEinR::optipar2021.Pro.CX$phiI[iw_coincidents]
  leafopt[['phiII']] <- SCOPEinR::optipar2021.Pro.CX$phiII[iw_coincidents]


  if (get.plots == T) {


    ########################################################################
    ### Get leaf reflectance Z
    ########################################################################

    wave.fluspect = c(spectral[['reg1']],spectral[['reg2']],spectral[['reg3']])
    rfl.fuspect = c(leafopt$refl[1,], rep(soil[['rs_thermal']],length(spectral[['IwlT']])))
    rfl.kChlrel = c(leafopt$kChlrel[1,], rep(soil[['rs_thermal']],length(spectral[['IwlT']])))
    rfl.kCarrel = c(leafopt$kCarrel[1,], rep(soil[['rs_thermal']],length(spectral[['IwlT']])))

    df.rfl <- data.frame(wave.fluspect=wave.fluspect, rfl.fuspect = rfl.fuspect,
                         kChlrel=rfl.kChlrel, kCarrel=rfl.kCarrel)


    p.z0 <- ggplot(data = df.rfl, aes(x = wave.fluspect, y = rfl.fuspect)) +
      labs(y= " leaf reflectance Z0 (fluspect-Cx)", x = "")+ xlim(400,2499) +
      geom_line()+ theme_bw()

    print(p.z0)


    ########################################################################
    ### Get Relative portion of pigments contribution
    ########################################################################

    p.k <- ggplot(data = df.rfl, aes(x = wave.fluspect)) +
      labs(y= "Relative portion of pigments contribution", x = "") +
      geom_line(aes(y = kChlrel, color = "chlorophylls"), linewidth = 0.5) +
      geom_line(aes(y = kCarrel, color = "carotenoids"), linewidth = 0.5) +
      scale_color_manual(name = "Irradiance type",
                         values = c("chlorophylls" = "forestgreen", "carotenoids" = "brown")) +
      theme_bw() + xlim(400,800) +
      theme(legend.position = "top") +
      guides(color = guide_legend(title = NULL))
    print(p.k)

    ########################################################################
    ### Get PhiI and PhiII
    ########################################################################


    df.phi<- data.frame(wave=spectral$wlF, phiI = leafopt$phiI,
                        phiII=leafopt$phiII, phiSum=leafopt$phiII+leafopt$phiI)

    p.phi <- ggplot(data = df.phi, aes(x = wave)) +
      labs(y= "phi contribution", x = "") +
      geom_line(aes(y = phiI, color = "phiI"), linewidth = 0.5) +
      geom_line(aes(y = phiII, color = "phiII"), linewidth = 0.5) +
      geom_line(aes(y = phiSum, color = "phiSum"), linewidth = 0.5) +
      scale_color_manual(name = "phi",
                         values = c("phiI" = "forestgreen", "phiII" = "brown",
                                    'phiSum' = 'black')) +
      theme_bw()  +
      theme(legend.position = "top") +
      guides(color = guide_legend(title = NULL))
    print(p.phi)

    ########################################################################
    ### Get MB and Mf
    ########################################################################

    wlf <- seq(640, 850, by = step-1) #
    df.Mb_Mf <- data.frame(wave.f=wlf, Mb = leafopt$Mb[,1,1],
                           Mf=leafopt$Mf[,1,1])

    p.mb <- ggplot(data = df.Mb_Mf, aes(x = wave.f)) +
      labs(y= "Relative portion of pigments contribution", x = "") +
      geom_line(aes(y = Mb, color = "chlorophylls"), linewidth = 0.5) +
      geom_line(aes(y = Mf, color = "carotenoids"), linewidth = 0.5) +
      scale_color_manual(name = "Irradiance type",
                         values = c("chlorophylls" = "forestgreen", "carotenoids" = "brown")) +
      theme_bw()  +
      theme(legend.position = "top") +
      guides(color = guide_legend(title = NULL))
    print(p.mb)

  }


  return(leafopt)

}
