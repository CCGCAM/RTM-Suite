#' Performs fourSAIL2 + PROSPECT simulation based on a set of combinations of input parameters
# ============================================================================= =
#' @param LUT_GB dataframe Includes distribution of biophysical parameters used as green vegetation  (first column).
#' Includes distribution of biophysical parameters used  (second column)
#       - Leafbrown  includes reflectance and 
#' @param inputLUT LUT table with distribution of biophysical parameters used as input parameters in the model
#' @param rsoil numeric. Soil reflectance
#' @param PROSPECTversion Legacy switch, kept for backward compatibility: 'PRO' or anything else ('D'). Ignored if \code{LeafModel} is given.
#' @param LeafModel character. One of 'PROSPECT-PRO', 'PROSPECT-D', 'Liberty', 'Fluspect-B', 'Fluspect-B-Cx' -- same 5 leaf models \code{\link{foursail}} and \code{\link{inform}} support. Defaults to \code{PROSPECTversion}-derived value if not given.
#' @param FieldObserv logical/NULL. If \code{NULL} (default), the "green vegetation" leaf optics are taken from \code{inputLUT} itself (i.e. the same leaf spectrum used for the main canopy), overriding whatever \code{LUT_GB[1,]} would otherwise give. If not \code{NULL}, the green/brown vegetation leaf optics are taken as computed from \code{LUT_GB} instead -- for matching a field-observed green/brown fraction defined independently of \code{inputLUT}.
#'
############################################################################################################
############################################################################################################
#       - TypeLidf  = Type of leaf inclination distribution function
#       - LIDFa = Parameter a.
#         if TypeLidf ==1, controls the average leaf slope
#	        LIDF type 		a 		 b
#	        Planophile 		1		 0
#	        Erectophile    -1	 	 0
#	        Plagiophile 	0		-1
#	        Extremophile 	0		 1
#	        Spherical 	   -0.35 	-0.15
#	        Uniform 0 0
#
#         if TypeLidf ==2, corresponds to average leaf angle
#       - LIDFb = Parameter b
#         if TypeLidf ==1, unused
#         if TypeLidf ==2, controls the distribution's bimodality
#         LIDFa	= average leaf angle (degrees) 0 = planophile	/	90 = erectophile

#       - lai = Leaf Area Index
#       - hot = Hot Spot parameter = ratio of the correlation length of leaf projections in the horizontal plane and the canopy height (doi:10.1016/j.rse.2006.12.013)
#       - tts = Sun zeith angle
#       - tto = Observer zeith angle
#       - psi = Azimuth Sun / Observer
#       - rsoil = Soil reflectance
#       - fraction_brown  = Fraction of brown leaf area
#       - diss  = Layer dissociation factor
#       - Cv = vertical crown cover percentage = % ground area covered with crowns as seen from nadir direction
#       - Zeta = Tree shape factor =  = ratio of crown diameter to crown height
#       - Leafgreen  includes reflectance and transmittance for green vegetation
#       - Leafbrown  includes reflectance and transmittance for brown vegetation
############################################################################################################
############################################################################################################
#'
#' @return list. rdot,rsot,rddt,rsdt
#' rdot: hemispherical-directional reflectance factor in viewing direction
#' rsot: bi-directional reflectance factor
#' rsdt: directional-hemispherical reflectance factor for solar incident flux
#' rddt: bi-hemispherical reflectance factor
#' alfast: canopy absorptance for direct solar incident flux
#' alfadt: canopy absorptance for hemispherical diffuse incident flux
#' @export
#' 
#' @references
#' 
#' Verhoef W & Bach H, 2007. Coupled soil–leaf-canopy and atmosphere radiative transfer modeling to simulate hyperspectral multi-angular surface reflectance and TOA radiance data. Remote Sensing of Environment, 109:166-182. doi:10.1016/j.rse.2006.12.013
#' 
#' Verhoef W, Jia L, Xiao Q & Su Z, 2007. Unified optical-thermal four-stream radiative transfer theory for homogeneous vegetation canopies. IEEE Transactions in Geosciences and Remote Sensing, 45:1808–1822. https://doi.org/10.1109/TGRS.2007.895844
#'
#' Jacquemoud S, Verhoef W, Baret F, Bacour C, Zarco-Tejada PJ, Asner GP, François C & Ustin SL, 2009. PROSPECT+ SAIL models: A review of use for vegetation characterization. Remote Sensing of Environment, 113:S56–S66. https://doi.org/doi:10.1016/j.rse.2008.01.026
#' 
#' Berger K, Atzberger C, Danner M, D’Urso G, Mauser W, Vuolo F & Hank T 2018. Evaluation of the PROSAIL Model Capabilities for Future Hyperspectral Model Environments: A Review Study. Remote Sensing, 10:85. https://doi.org/10.3390/rs10010085
#' 
#' Authors: 
#' 
#' Verhoef W.
#' 
#' Bach H.
#' 
#' Authors of the R version:
#' 
#' Jean-Baptiste Feret
#' 
#' The fourSAIL model is based on a version provided by	Wout Verhoef et al. (2007)
#' 
#' original version downloadable at http://teledetection.ipgp.jussieu.fr/prosail/
#' 
#' Improved and extended version of SAILH model that avoids numerical singularities
#' and works more efficiently if only few parameters change.

foursail2 <- function(LUT_GB=NULL, inputLUT,rsoil, PROSPECTversion='PRO',FieldObserv=NULL,
                      LeafModel=NULL){

#define alll inputs in the models. Retrieved from LUT tables

## fourSAIL
LIDFa=inputLUT[,'LIDFa']; LIDFb=inputLUT[,'LIDFb']; TypeLidf=inputLUT[,'TypeLidf']; lai=inputLUT[,'LAI']
hotspot=inputLUT[,'hspot']; tts=inputLUT[,'tts']; tto=inputLUT[,'tto']; psi=inputLUT[,'psi']
fraction_brown = inputLUT[,'fraction_brown']; diss = inputLUT[,'diss']; Cv = inputLUT[,'Cv'];Zeta = inputLUT[,'Zeta']

# LeafModel supersedes the legacy PROSPECTversion ('PRO'/anything-else) switch,
# so foursail2 can use any of the 5 leaf models fourSAIL/INFORM already support,
# instead of being limited to PROSPECT-PRO/PROSPECT-D.
if (is.null(LeafModel)) {
  LeafModel <- if (identical(PROSPECTversion, 'PRO')) 'PROSPECT-PRO' else 'PROSPECT-D'
}

if (is.null(LUT_GB)){
  message('Please define same spectral domain for GreenVegetation and BrownVegetation and SpecPROSPECT')
  #message('Please define same spectral domain for GreenVegetation and BrownVegetation and SpecPROSPECT')
  # Illustrative green/brown (senescent) defaults, covering the columns needed
  # by any of the 5 leaf models -- override with a real LUT_GB for actual use,
  # exactly as the message above already asks callers to do.
  LUT_GB<-data.frame(N=c(1.5, 2), Cab=c(40,5),Car=c(8,5),Anth=c(0,0),Cbrown=c(0,0),
                     EWT=c(0.01, 0.005), LMA=c(0.009,0.008), alpha=c(40,40),
                     Prot=c(0 , 0),CBC=c(0,0),
                     # Fluspect-B/-Cx extras
                     Cs=c(0, 0.1), fqe=c(0.01, 0.01), Cx=c(0, 0),
                     # Liberty extras (typical values per liberty.R's own docs)
                     cell.d=c(40, 40), inter.c=c(0.045, 0.045), baseline.abs=c(0.0006, 0.0012),
                     leaf.thick=c(1.6, 1.6), albino.abs=c(0, 2), lign.cell=c(2, 6), Nitrogen=c(1, 0.5))
}

## Dispatch on LeafModel, matching foursail()'s pattern exactly so the same
## 5 leaf models work here too (previously only PROSPECT-PRO/-D did).
LRT             <- get.foursail2.leafopt(inputLUT, LeafModel)
GreenVegetation <- get.foursail2.leafopt(LUT_GB[1, , drop = FALSE], LeafModel)
BrownVegetation <- get.foursail2.leafopt(LUT_GB[2, , drop = FALSE], LeafModel)

# NB (fix): Fluspect-B/-Cx's optipar table only spans 400-2400nm (2001 pts),
# so leafgreen/leafbrown's Reflectance/Transmittance below are 2001 long --
# but rsoil (typically the caller's full 400-2500nm, 2101-long spectrum) was
# never truncated to match, unlike foursail()'s equivalent Fluspect branches
# (see foursail.R's "### reduce the wavelength for r soil"). That length
# mismatch propagated into rddsoil/rdosoil/rsdsoil/rsosoil below and from
# there into the whole SAIL2 radiative-transfer core, silently recycling
# rsoil against 2001-long vectors throughout instead of erroring or matching
# domains, corrupting the soil-interaction terms for both leaf models.
if (LeafModel %in% c('Fluspect-B', 'Fluspect-B-Cx')) {
  rsoil <- rsoil[1:2001]
}


###force to use different Green vegetation
if (is.null(FieldObserv)){
  GreenVegetation<-LRT
  BrownVegetation<-BrownVegetation
}
  
  ### Asign in a list the two references spectrum for Green Vegetation and Brown Vegetation
leafgreen<-list()
leafgreen$Reflectance<-GreenVegetation[[2]]
leafgreen$Transmittance<-GreenVegetation[[3]]
leafbrown <- list()
leafbrown$Reflectance<-BrownVegetation[[2]]
leafbrown$Transmittance<-BrownVegetation[[3]]
#rho	 <- 	LRT[[2]] #rho Reflectance
#tau	 <- 	LRT[[3]] #tau Transmittance

########################################
#	1.2 Geometric quAnthities
#########################################

#	This version does not include non-Lambertian soil properties.
#	original codes do, and only need to add the following variables as input
rddsoil <- rdosoil <- rsdsoil <- rsosoil <- rsoil

#	Geometric quAnthities
rd <- pi / 180

#	Generate leaf angle distribution from average leaf angle (ellipsoidal) or (a,b) parameters
if (TypeLidf == 1){
  LeafDistribution <- dladgen(LIDFa,LIDFb)
  lidf <- LeafDistribution$lidf
  litab <- LeafDistribution$litab
  
} else if (TypeLidf == 2){
  LeafDistribution <- campbell(ala = LIDFa)
  lidf <- LeafDistribution$lidf
  litab <- LeafDistribution$litab
}

if (lai < 0){
  message('Please define positive LAI value')
  rddt <- rsdt <- rdot <- rsost <- rsot <- rsoil
  alfast <- alfadt <- 0 * rsoil
} else if (lai == 0){
  tss <- too <- tsstoo <- tdd <- 1.0
  rdd <- rsd <- tsd <- rdo <- tdo <- 0.0
  rso <- rsos <- rsod <- rsodt <- 0.0
  rddt <- rsdt <- rdot <- rsost <- rsot <- rsoil
  alfast <- alfadt <- 0 * rsoil
} else if (lai > 0){
  cts <- cos(rd * tts)
  cto <- cos(rd * tto)
  ctscto <- cts * cto
  tAnths <- tan(rd * tts)
  tAntho <- tan(rd * tto)
  cospsi <- cos(rd * psi)
  dso <- sqrt(tAnths * tAnths + tAntho * tAntho - 2.0 * tAnths * tAntho * cospsi)
  ### ### ### ### ### ### ### ### ### ### ### ### ###
  ## Crown and vegatation clumping effects
  ## similar to flim implementation
  # Clumping effects
  ### ### ### ### ### ### ### ### ### ### ### ### ###
  Cs <- Co <- 1.0
  if (Cv <= 1.0){
    Cs <- 1.0 - (1.0 - Cv)^(1.0 / cts)
    Co <- 1.0 -(1.0 - Cv)^(1.0 / cto)
  }

  Overlap <- 0.0
  if (Zeta > 0.0){
    Overlap <- min(Cs * (1.0 - Co), Co * (1.0 - Cs)) * exp(-dso / Zeta)
  }

  
  Fcd <- Cs * Co + Overlap
  Fcs <- (1.0 - Cs) *Co - Overlap
  Fod <- Cs*(1.0 - Co) - Overlap
  Fos <- (1.0 - Cs)*(1.0 - Co) + Overlap
  Fcdc <- 1.0 - (1.0 - Fcd)^(0.5 / cts + 0.5 / cto)
  
  #	Part depending on diss, fraction_brown, and leaf optical properties
  #	First save the input fraction_brown as the old fraction_brown, as the following change is only artificial
  # Better define an fraction_brown that is actually used: fb, so that the input is not modified!
  
  fb <- fraction_brown
  # if only green leaves
  if (fraction_brown == 0.0){
    fb <- 0.5
    leafbrown$Reflectance <- leafgreen$Reflectance
    leafbrown$Transmittance <- leafgreen$Transmittance
  }
  if (fraction_brown == 1.0){
    fb <- 0.5
    leafgreen$Reflectance <- leafbrown$Reflectance
    leafgreen$Transmittance <- leafbrown$Transmittance
  }
  s <- (1.0 - diss) * fb * (1.0 - fb)
  # rho1 & tau1 : green foliage
  # rho2 & tau2 : brown foliage (bottom layer)
  rho1 <- ((1 - fb - s) * leafgreen$Reflectance + s * leafbrown$Reflectance) / (1 - fb)
  tau1 <- ((1 - fb - s) * leafgreen$Transmittance + s * leafbrown$Transmittance) / (1 - fb)
  rho2 <- (s * leafgreen$Reflectance + (fb - s) * leafbrown$Reflectance) / fb
  tau2 <- (s * leafgreen$Transmittance +(fb - s)* leafbrown$Transmittance) / fb
  
  # angular distance, compensation of shadow length
  #	Calculate geometric factors associated with extinction and scattering
  #	Initialise sums
  ks <- ko <- bf <- sob <- sof <- 0
  
  # Weighted sums over LIDF
  
  for (i in 1:length(litab)){
    ttl <- litab[i]
    ctl <- cos(rd * ttl)
    # SAIL volscatt function gives interception coefficients
    # and two portions of the volume scattering phase function to be
    # multiplied by rho and tau, respectively
    resVolscatt <- volscatt(tts,tto,psi,ttl)

    chi_s <- resVolscatt[[1]] ##chi_s
    chi_o <- resVolscatt[[2]] ##chi_o
    frho <- resVolscatt[[3]] ##frho
    ftau <- resVolscatt[[4]] ##ftau
    # Extinction coefficients
    ksli <- chi_s / cts
    koli <- chi_o / cto
    # Area scattering coefficient fractions
    sobli <- frho * pi / ctscto
    sofli <- ftau *pi / ctscto
    bfli <- ctl * ctl
    ks <- ks + ksli * lidf[i]
    ko <- ko + koli * lidf[i]
    bf <- bf + bfli * lidf[i]
    sob <- sob + sobli * lidf[i]
    sof <- sof + sofli * lidf[i]
  }
  # Geometric factors to be used later in combination with rho and tau
  sdb <- 0.5 * (ks + bf)
  sdf <- 0.5 * (ks - bf)
  dob <- 0.5 * (ko + bf)
  dof <- 0.5 * (ko - bf)
  ddb <- 0.5 * (1. + bf)
  ddf <- 0.5 * (1.- bf)
  
  # LAIs in two layers
  lai1 <- (1 - fb) * lai
  lai2 <- fb * lai
  
  tss <- exp(-ks * lai)
  ck <- exp(-ks * lai1)
  alf <- 1e6
  if (hotspot > 0.0){
    alf <- (dso / hotspot) *2.0 / (ks + ko)
  }
  if (alf > 200.0){
    alf <- 200.0     # inserted H. Bach 1/3/04
  }
  if (alf == 0.0){
    # The pure hotspot
    tsstoo <- tss
    s1 <- (1 - ck) / (ks * lai)
    s2 <- (ck - tss)/(ks * lai)
  } else {
    # Outside the hotspot
    fhot <- lai * sqrt(ko * ks)
    # Integrate 2 layers by exponential simpson method in 20 steps
    # the steps are arranged according to equal partitioning
    # of the derivative of the joint probability function
    x1 <- y1 <- 0.0
    f1 <- 1.0
    ca <- exp( alf * ( fb - 1.0) )
    fint <- (1.0 - ca) * 0.05
    s1 <- 0.0
    for (istep in 1:20){
      if (istep < 20){
        x2 <- -log(1. - istep * fint) / alf
      } else {
        x2 <- 1. - fb
      }
      y2 <- -(ko + ks) * lai * x2 + fhot * (1.0 - exp(-alf * x2)) / alf
      f2 <- exp(y2)
      s1 <- s1 + (f2 - f1) * (x2 - x1) / (y2 - y1)
      x1 <- x2
      y1 <- y2
      f1 <- f2
    }
    fint <- (ca - exp(-alf)) * 0.05
    s2 <- 0.0
    for (istep in 1:20){
      if (istep < 20){
        x2 <- -log(ca - istep * fint) / alf
      } else {
        x2 <- 1.0
      }
      y2 <- -(ko + ks) * lai * x2 + fhot * (1.0 - exp(-alf * x2)) / alf
      f2 <- exp(y2)
      s2 <- s2 + (f2 - f1) * (x2 - x1) / (y2 - y1)
      x1 <- x2
      y1 <- y2
      f1 <- f2
    }
    tsstoo <- f1
  }
  
  # Calculate reflectances and transmittances
  # Bottom layer
  tss <- exp(-ks * lai2)
  too <- exp(-ko * lai2)
  sb <- sdb * rho2 + sdf * tau2
  sf <- sdf * rho2 + sdb * tau2
  
  vb <- dob * rho2 + dof * tau2
  vf <- dof * rho2 + dob * tau2
  
  w2 <- sob * rho2 + sof * tau2
  
  sigb <- ddb * rho2 + ddf * tau2
  sigf <- ddf * rho2 + ddb * tau2
  att <- 1.0 - sigf
  m2 <- (att + sigb) * (att - sigb)
  m2[m2<0] <- 0
  m <- sqrt(m2)
  ### Non Conservative scattering
  f_Non_ConS <- which(m > 0.01)
  ## Conservative scattering
  f_ConS <- which(m <= 0.01)
  
  tdd <- rdd <- tsd <- rsd <- tdo <- rdo <- 0 * m
  rsod <- 0 * m
  if (length(f_Non_ConS) > 0){
    resNCS <- ToolsRTM::NonConservativeScattering(m[f_Non_ConS],lai2,att[f_Non_ConS],sigb[f_Non_ConS],
                                        ks,ko,sf[f_Non_ConS],sb[f_Non_ConS],vf[f_Non_ConS],vb[f_Non_ConS],tss,too)
    tdd[f_Non_ConS] <- resNCS$tdd
    rdd[f_Non_ConS] <- resNCS$rdd
    tsd[f_Non_ConS] <- resNCS$tsd
    rsd[f_Non_ConS] <- resNCS$rsd
    tdo[f_Non_ConS] <- resNCS$tdo
    rdo[f_Non_ConS] <- resNCS$rdo
    rsod[f_Non_ConS] <- resNCS$rsod
  }
  if (length(f_ConS) > 0){
    resCS <- ToolsRTM::ConservativeScattering(m[f_ConS],lai2,att[f_ConS],sigb[f_ConS],
                                    ks,ko,sf[f_ConS],sb[f_ConS],vf[f_ConS],vb[f_ConS],tss,too)
    tdd[f_ConS] <- resCS$tdd
    rdd[f_ConS] <- resCS$rdd
    tsd[f_ConS] <- resCS$tsd
    rsd[f_ConS] <- resCS$rsd
    tdo[f_ConS] <- resCS$tdo
    rdo[f_ConS] <- resCS$rdo
    rsod[f_ConS] <- resCS$rsod
  }
  
  # Set background properties equal to those of the bottom layer on a black soil
  rddb <- rdd
  rsdb <- rsd
  rdob <- rdo
  rsodb <- rsod
  tddb <- tdd
  tsdb <- tsd
  tdob <- tdo
  toob <- too
  tssb <- tss
  # Top layer
  tss <- exp(-ks * lai1)
  too <- exp(-ko * lai1)
  
  sb <- sdb * rho1 + sdf * tau1
  sf <- sdf * rho1 + sdb * tau1
  
  vb <- dob * rho1 + dof * tau1
  vf <- dof * rho1 + dob * tau1
  
  w1 <- sob * rho1 + sof * tau1
  
  sigb <- ddb * rho1 + ddf * tau1
  sigf <- ddf * rho1 + ddb * tau1
  att <- 1.0 - sigf
  
  m2 <- (att + sigb) * (att - sigb)
  m2[m2<0] <- 0
  m <- sqrt(m2)
  f_Non_ConS <- which(m > 0.01)
  f_ConS <- which(m <= 0.01)
  
  tdd <- rdd <- tsd <- rsd <- tdo <- rdo <- 0 * m
  rsod <- 0 * m
  if (length(f_Non_ConS) > 0){
    resNCS <- ToolsRTM::NonConservativeScattering(m[f_Non_ConS],lai1,att[f_Non_ConS],sigb[f_Non_ConS],
                                        ks,ko,sf[f_Non_ConS],sb[f_Non_ConS],vf[f_Non_ConS],vb[f_Non_ConS],tss,too)
    tdd[f_Non_ConS] <- resNCS$tdd
    rdd[f_Non_ConS] <- resNCS$rdd
    tsd[f_Non_ConS] <- resNCS$tsd
    rsd[f_Non_ConS] <- resNCS$rsd
    tdo[f_Non_ConS] <- resNCS$tdo
    rdo[f_Non_ConS] <- resNCS$rdo
    rsod[f_Non_ConS] <- resNCS$rsod
  }
  if (length(f_ConS) > 0){
    resCS <- ToolsRTM::ConservativeScattering(m[f_ConS],lai1,att[f_ConS],sigb[f_ConS],
                                    ks,ko,sf[f_ConS],sb[f_ConS],vf[f_ConS],vb[f_ConS],tss,too)
    tdd[f_ConS] <- resCS$tdd
    rdd[f_ConS] <- resCS$rdd
    tsd[f_ConS] <- resCS$tsd
    rsd[f_ConS] <- resCS$rsd
    tdo[f_ConS] <- resCS$tdo
    rdo[f_ConS] <- resCS$rdo
    rsod[f_ConS] <- resCS$rsod
  }
  
  # Combine with bottom layer reflectances and transmittances (adding method)
  rn <- 1.0 - rdd * rddb
  tup <- (tss * rsdb + tsd * rddb) / rn
  tdn <- (tsd + tss * rsdb * rdd) / rn
  rsdt <- rsd + tup * tdd
  rdot <- rdo + tdd * (rddb * tdo + rdob * too) / rn
  rsodt <- rsod + (tss * rsodb + tdn * rdob) * too + tup * tdo
  
  rsost <- (w1 * s1 + w2 * s2) * lai
  
  rsot <- rsost + rsodt
  
  # Diffuse reflectances at the top and the bottom are now different
  rddt_t <- rdd + tdd * rddb * tdd / rn
  rddt_b <- rddb + tddb * rdd * tddb / rn
  
  # Transmittances of the combined canopy layers
  tsst <- tss * tssb
  toot <- too * toob
  tsdt <- tss * tsdb + tdn * tddb
  tdot <- tdob * too + tddb *(tdo + rdd * rdob * too) / rn
  tddt <- tdd * tddb / rn
  
  # Apply clumping effects to vegetation layer
  rddcb <- Cv * rddt_b
  rddct <- Cv * rddt_t
  tddc <- 1 - Cv + Cv * tddt
  rsdc <- Cs * rsdt
  tsdc <- Cs * tsdt
  rdoc <- Co * rdot
  tdoc <- Co * tdot
  tssc <- 1 - Cs + Cs * tsst
  tooc <- 1 - Co + Co * toot
  
  # New weight function Fcdc for crown contribution (W. Verhoef, 22-05-08)
  rsoc <- Fcdc * rsot
  tssooc <- Fcd * tsstoo + Fcs * toot + Fod * tsst + Fos
  # Canopy absorptance for black background (W. Verhoef, 02-03-04)
  alfas <- 1. - tssc - tsdc - rsdc
  alfad <- 1. - tddc - rddct
  # Add the soil background
  rn <- 1 - rddcb * rddsoil
  tup <- (tssc * rsdsoil + tsdc * rddsoil) / rn
  tdn <- (tsdc + tssc * rsdsoil * rddcb) / rn
  
  rddt <- rddct + tddc * rddsoil * tddc / rn
  rsdt <- rsdc + tup * tddc
  rdot <- rdoc + tddc * (rddsoil * tdoc + rdosoil * tooc) / rn
  rsot <- rsoc + tssooc * rsosoil + tdn * rdosoil * tooc + tup * tdoc
  
  # Effect of soil background on canopy absorptances (W. Verhoef, 02-03-04)
  alfast <- alfas + tup * alfad
  alfadt <- alfad * (1. + tddc * rddsoil / rn)
}

LSTa<- list("rdot" = rdot,"rsot" =rsot,"rddt" =rddt,"rsdt" =rsdt,
            "alfast" = alfast, "alfadt" = alfadt)
return(LSTa)
}

#' Leaf optics dispatch used internally by \code{foursail2}
#'
#' Computes leaf reflectance/transmittance for a single leaf-parameter row
#' using whichever of the 5 supported leaf models is requested, mirroring
#' the dispatch already used in \code{\link{foursail}} and \code{\link{inform}}
#' (so \code{foursail2} is no longer limited to PROSPECT-PRO/PROSPECT-D).
#'
#' @param inputRow data.frame. One row of leaf parameters (column names depend on \code{LeafModel}).
#' @param LeafModel character. One of \code{"PROSPECT-PRO"}, \code{"PROSPECT-D"}, \code{"Liberty"}, \code{"Fluspect-B"}, \code{"Fluspect-B-Cx"}.
#' @return A list whose 2nd element is reflectance and 3rd element is transmittance (same convention as \code{prospect_PRO}/\code{prospect_DB}/\code{liberty}/\code{getFluspect.B}/\code{getFluspect.Cx}).
#' @keywords internal
get.foursail2.leafopt <- function(inputRow, LeafModel) {

  if (LeafModel == 'PROSPECT-PRO') {
    N=inputRow[,'N']; Cab=inputRow[,'Cab']; Car=inputRow[,'Car']; Anth=inputRow[,'Anth']; Cbrown=inputRow[,'Cbrown']
    EWT=inputRow[,'EWT']; LMA=inputRow[,'LMA']; alpha=inputRow[,'alpha']
    Prot=inputRow[,'Prot']; CBC=inputRow[,'CBC']
    LRT <- prospect_PRO(N,Cab,Car,Anth,Cbrown,EWT,LMA,alpha,Prot,CBC)
    message('SAIL2 with PROSPECT-PRO is processing')

  } else if (LeafModel == 'PROSPECT-D') {
    N=inputRow[,'N']; Cab=inputRow[,'Cab']; Car=inputRow[,'Car']; Anth=inputRow[,'Anth']; Brown=inputRow[,'Cbrown']
    EWT=inputRow[,'EWT']; LMA=inputRow[,'LMA']; alpha=inputRow[,'alpha']
    LRT <- prospect_DB(N,Cab,Car,Anth,Brown,EWT,LMA,alpha)
    message('SAIL2 with PROSPECT-D is processing')

  } else if (LeafModel == 'Liberty') {
    LRT <- liberty(inputRow)
    message('SAIL2 with Liberty model is processing')

  } else if (LeafModel == 'Fluspect-B') {
    LRT <- ToolsRTM::getFluspect.B(inputsLeaf = inputRow, inputsOptipar = ToolsRTM::optipar, version = 'D')
    message('SAIL2 with Fluspect-B model is processing')

  } else if (LeafModel == 'Fluspect-B-Cx') {
    LRT <- ToolsRTM::getFluspect.Cx(inputsLeaf = inputRow, inputsOptipar = ToolsRTM::optipar, version = 'Cx')
    message('SAIL2 with Fluspect-B-Cx model is processing')

  } else {
    stop("Unknown LeafModel: '", LeafModel, "'. Use one of 'PROSPECT-PRO', 'PROSPECT-D', ",
         "'Liberty', 'Fluspect-B', 'Fluspect-B-Cx'.")
  }

  LRT
}

