#' fourSAIL model coupled with several leaf models
#' 
#' \code{foursail} fourSAIL simulation based on a set of combinations of input parameters
#' 
#' @param rsoil numeric. Soil reflectance
#' @param inputLUT LUT table with distribution of biophysical parameters used as input parameters in the model
#' @param LeafModel Version of PROSPECT model, Liberty model and  fluspect model. For PROSPECT model: 'PROSPECT-PRO' or 'PROSPECT-D' is accepted. By default 'PROSPECT-PRO' is used.
#'  fluspect model is valid in two options: 'Fluspect-B' and 'Fluspect-B-Cx'. Liberty as 'Liberty'
#' @param spectrum.all a boolean value, False is for SPART and Fluspect Models (400-2400 nm), True for the PROSPECT and Liberty models (400-2500 nm)
#' Liberty model is leaf radiative transfer model designed for conifer needles. Uses 'Liberty'
#' Fluspect-B model is leaf radiative transfer model designed for adding fluorescence emission. Uses 'Fluspect-B' or 'Fluspect-B-Cx'.
#' @return list. rdot,rsot,rddt,rsdt
#' 
#' rdot: hemispherical-directional reflectance factor in viewing direction
#' rsot: bi-directional reflectance factor
#' rsdt: directional-hemispherical reflectance factor for solar incident flux
#' rddt: bi-hemispherical reflectance factor
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
#' @author Wout Verhoef, Bach H. , JJean-Baptiste Feret (Original version in Matlab) 
#' This version is also included in prospect package
#' @author Carlos Camino (Ported version into R wit modification from fourSAIL model in prospect package)
#' 
#' The fourSAIL model is based on a version provided by	Wout Verhoef et al. (2007)
#' 
#' original version downloadable at http://teledetection.ipgp.jussieu.fr/prosail/
#' 
#' Improved and extended version of SAILH model that avoids numerical singularities
#' and works more efficiently if only few parameters change.

#' 
foursail <- function(inputLUT,rsoil, LeafModel='PROSPECT-PRO',spectrum.all = T){

  
  ## for using the range 400-2400 nm (FLuspect and SPART models)
  if (missing(spectrum.all) || spectrum.all) {
    spectrum.all = T ### when is true is using (400-2500 nm) for PROSPECT and Liberty
  }
    
  ## parameters for fourSAIL
  LIDFa=inputLUT[,'LIDFa']; LIDFb=inputLUT[,'LIDFb']; TypeLidf=inputLUT[,'TypeLidf']; lai=inputLUT[,'LAI']
  hotspot=inputLUT[,'hspot']; tts=inputLUT[,'tts']; tto=inputLUT[,'tto']; psi=inputLUT[,'psi']
  
########################################
#	1.1 Leaf optical properties
#########################################
if (LeafModel == 'PROSPECT-PRO') {
  
  #define alll inputs in the models. retreived from LUT tables
  N=inputLUT[,'N']; Cab=inputLUT[,'Cab']; Car=inputLUT[,'Car']; Anth=inputLUT[,'Anth']; Cbrown=inputLUT[,'Cbrown']
  EWT=inputLUT[,'EWT']; LMA=inputLUT[,'LMA']; alpha=inputLUT[,'alpha']
  Prot=inputLUT[,'Prot'];CBC=inputLUT[,'CBC']
  # run LeafModel ='PROSPECT-PRO'
  LRT <- prospect_PRO(N,Cab,Car,Anth,Cbrown,EWT,LMA,alpha,Prot,CBC)
  message('SAIL with PROSPECT-PRO is processing')
  if (spectrum.all == F){
    ### reduce the wavelength for r soil
    rsoil =  rsoil[1:2001]
  }
  
  } else if (LeafModel == 'Liberty'){
    
  #define alll inputs in the models. retreived from LUT tables
  # run LeafModel ='PRO'
  LRT <- liberty(inputLUT)
  message('SAIL with Liberty model is processing')
  if (spectrum.all == F){
    ### reduce the wavelength for r soil
    rsoil =  rsoil[1:2001]
  }
  
  } else if (LeafModel == 'Fluspect-B'){

  # run LeafModel ='fluspect-D'
  LRT <- ToolsRTM::getFluspect.B(inputsLeaf= inputLUT, inputsOptipar =ToolsRTM::optipar, version = 'D')
  message('SAIL with Fluspect-D model is processing')
  ### reduce the wavelength for r soil
  rsoil =  rsoil[1:2001]
    
  
  
  } else if (LeafModel == 'Fluspect-B-Cx'){
 
  # run LeafModel ='fluspect-Cx'
  LRT <- ToolsRTM::getFluspect.Cx(inputsLeaf= inputLUT, inputsOptipar =ToolsRTM::optipar, version = 'Cx')
  message('SAIL with Fluspect-PRO model is processing')
  ### reduce the wavelength for r soil
  rsoil =  rsoil[1:2001]
  
  } else if (LeafModel == 'PROSPECT-D') {
  
  #define alll inputs in the models. retreived from LUT tables
  N=inputLUT[,'N']; Cab=inputLUT[,'Cab']; Car=inputLUT[,'Car']; Anth=inputLUT[,'Anth']; Brown=inputLUT[,'Cbrown']
  EWT=inputLUT[,'EWT']; LMA=inputLUT[,'LMA']; alpha=inputLUT[,'alpha']
  # run LeafModel ='PROSPECT-D'
  LRT <- prospect_DB(N,Cab,Car,Anth,Brown,EWT,LMA,alpha)
  message('SAIL with PROSPECT-D is processing')
  if (spectrum.all == F){
    ### reduce the wavelength for r soil
    rsoil =  rsoil[1:2001]
  }

  } else {
    
  stop('a leaf model is needed')
    
  }
  
if (spectrum.all == T) {
  rho	 <- 	LRT[[2]]
  tau	 <- 	LRT[[3]]
} else {
  rho	 <- 	LRT[[2]][1:2001]
  tau	 <- 	LRT[[3]][1:2001]
}


########################################
#	1.2 Geometric quAnthities
#########################################

rd <- pi / 180
cts		 <-  cos(rd * tts)
cto		 <-  cos(rd * tto)
ctscto	 <-  cts * cto
tAnths	 <-  tan(rd * tts)
tAntho	 <-  tan(rd * tto)
cospsi	 <-  cos(rd * psi)
dso		 <-  sqrt(tAnths * tAnths + tAntho * tAntho - 2 * tAnths * tAntho * cospsi)



###########################################################################################################################
#	1.3 Generate leaf angle distribution from average leaf angle (ellipsoidal) or (a,b) parameters
###########################################################################################################################

if (TypeLidf==1){
LeafDistribution <- dladgen(LIDFa,LIDFb)
lidf <- LeafDistribution$lidf
litab <- LeafDistribution$litab

} else if (TypeLidf==2){
  LeafDistribution <- campbell(LIDFa)
  lidf <- LeafDistribution$lidf
  litab <- LeafDistribution$litab
}


  # angular distance, compensation of shadow length
	#	Calculate geometric factors associated with extinction and scattering
	#	Initialise sums
	ks	 <-  0
	ko	 <-  0
	bf	 <-  0
	sob	 <-  0
	sof	 <-  0

	#	Weighted sums over LIDF
    na <- length(litab)
    #ksli<-rep(NA,13)
	for (i in 1:na){
	  
		ttl <- litab[i]# leaf inclination discrete values
		ctl <- cos(rd * ttl)
		#	SAIL volume scattering phase function gives interception and portions to be
		#	multiplied by rho and tau

		chi_s_chi_o_frho_ftau <- ToolsRTM::volscatt(tts,tto,psi,ttl)
		chi_s<-chi_s_chi_o_frho_ftau[[1]]
		chi_o<-chi_s_chi_o_frho_ftau[[2]]
		frho<-chi_s_chi_o_frho_ftau[[3]]
		ftau<-chi_s_chi_o_frho_ftau[[4]] 

		#********************************************************************************
		#*                   SUITS SYSTEM COEFFICIENTS
		#*
		#*	ks  : Extinction coefficient for direct solar flux
		#*	ko  : Extinction coefficient for direct observed flux
		#*	att : Attenuation coefficient for diffuse flux
		#*	sigb : Backscattering coefficient of the diffuse downward flux
		#*	sigf : Forwardscattering coefficient of the diffuse upward flux
		#*	sf  : Scattering coefficient of the direct solar flux for downward diffuse flux
		#*	sb  : Scattering coefficient of the direct solar flux for upward diffuse flux
		#*	vf   : Scattering coefficient of upward diffuse flux in the observed direction
		#*	vb   : Scattering coefficient of downward diffuse flux in the observed direction
		#*	w   : Bidirectional scattering coefficient
		#********************************************************************************

		#	Extinction coefficients ksli
		ksli <- chi_s / cts
		koli <- chi_o / cto

		#	Area scattering coefficient fractions
		sobli	 <-  frho * pi / ctscto
		sofli	 <-  ftau * pi / ctscto
		bfli	 <-  ctl * ctl
		ks	 <-  ks + ksli * lidf[i]
	
		ko	 <-  ko + koli * lidf[i]
		bf	 <-  bf + bfli * lidf[i]
		sob	 <-  sob + sobli * lidf[i]
		sof	 <-  sof + sofli * lidf[i]
		######
		
	}
    
  ######################################################
	#	Geometric factors to be used later with rho and tau
  #####################################################
    
	sdb	 <-  0.5 * (ks + bf) #
	sdf	 <-  0.5 * (ks - bf) # weight of specular2diffuse     foward  scatter coefficient

	ddb	 <-  0.5 * (1.+ bf) #
	ddf	 <-  0.5 * (1.- bf) # weight of diffuse2diffuse back scatter coefficient
	
	dob	 <-  0.5 * (ko + bf) # weight of diffuse2directional  back    scatter coefficient
	dof	 <-  0.5 * (ko - bf) # weight of diffuse2directional  forward scatter coefficient

	#	Here rho and tau come in
	sigb <-  ddb * rho + ddf * tau
	sigf <-  ddf * rho + ddb * tau
	att	 <-  1 - sigf
	m2  <- (att + sigb) * (att - sigb)
	m2[which(m2 <= 0)] <- 0
	
	m    <- sqrt(m2)

	sb  <- sdb * rho + sdf *tau
	sf	 <-  sdf * rho + sdb * tau
	vb	 <-  dob * rho + dof * tau
	vf	 <-  dof * rho + dob * tau
	w	 <-  sob * rho + sof * tau

	
	######################################################
	#	Here the LAI comes in
	#   Outputs for the case LAI = 0
	######################################################
	
	if (lai<0){
		tss		 <-  1
		too		 <-  1
		tsstoo	 <-  1
		rdd		 <-  0
		tdd		 <-  1
		rsd		 <-  0
		tsd		 <-  0
		rdo		 <-  0
		tdo		 <-  0
		rso		 <-  0
		rsos	 <-  0
		rsod	 <-  0

		rddt	 <-  rsoil
		rsdt	 <-  rsoil
		rdot	 <-  rsoil
		rsodt	 <-  0 * rsoil
		rsost	 <-  rsoil
		rsot	 <-  rsoil

	} else {
    

	######################################################
	#	Other cases (LAI > 0)
	######################################################
  ###########
  tss	 <-  exp( -ks * lai)
  too	 <-  exp( -ko * lai)
  
  NonConScatt <- ToolsRTM::NonConservativeScattering(m,lai,att,sigb,ks,ko,sf,sb,vf,vb,tss,too)
  tdd=unlist(NonConScatt[[1]]) #tdd
  rdd=unlist(NonConScatt[[2]]) #rdd
  tsd=unlist(NonConScatt[[3]]) #tsd
  rsd=unlist(NonConScatt[[4]]) #rsd
  tdo=unlist(NonConScatt[[5]]) #tdo
  rdo=unlist(NonConScatt[[6]]) #rdo
  rsod=unlist(NonConScatt[[7]]) #rsod

	#############################################################################
	#	Treatment of the hotspot-effect
	#############################################################################
	
	alf <- 1e6

	######################################################################
	#	Apply correction 2/(K+k) suggested by F-M Br?on
	if (hotspot > 0 ){
		alf <- (dso / hotspot) * 2 / (ks + ko)
	}
	######################################################################
	
	
	######################################################################
	if (alf > 200){# inserted H Bach 1/3/04
		alf <- 200
	}
	######################################################################
	
	######################################################################
	if (alf == 0) {
  		#	The pure hotspot - no shadow
  		tsstoo <- tss
  		sumint <- (1 - tss) / (ks * lai)
  		
      } else {
      
      ######################################################################
  		#	Outside the hotspot
      ######################################################################  
  		fhot <- lai * sqrt(ko * ks)
  		#	Integrate by exponential Simpson method in 20 steps
  		#	the steps are arranged according to equal partitioning
  		#	of the slope of the joint probability function
  		x1 <- 0
  		y1 <- 0
      f1 <- 1
      fint <- (1. - exp(-alf)) * 0.05
  		sumint <- 0

    		for (i in 1:20){
    		  if (i < 20){ 
    		  x2 <- -log(1 - i * fint) / alf 
    		  } else {
    		  x2 <- 1 }
    			
    		  y2 <- -(ko + ks) * lai * x2 + fhot *(1 - exp(-alf * x2)) / alf
    		  #print(y2)
    			f2 <- exp(y2)
    			#print(sumint)
    			
    	    sumint <- sumint + (f2 - f1) * (x2 - x1) / (y2 - y1)
    		
    			x1 <- x2
    			y1 <- y2
    			f1 <- f2
         }
		  tsstoo <- f1
      }
  
	######################################################################
	
	
	
#	Bidirectional reflectance
#	Single scattering contribution
	rsos <- w * lai * sumint
#	Total canopy contribution
	rso <- rsos + rsod

  #	Interaction with the soil
  dn <- 1- rsoil * rdd
  # rddt: bi-hemispherical reflectance factor
  rddt <- rdd + tdd * rsoil * tdd / dn
  # rsdt: directional-hemispherical reflectance factor for solar incident flux
  rsdt <- rsd + (tsd + tss) * rsoil * tdd / dn
  # rdot: hemispherical-directional reflectance factor in viewing direction
  rdot <- rdo + tdd * rsoil * (tdo + too) / dn
  # rsot: bi-directional reflectance factor
  rsodt <- rsod + ((tss + tsd) * tdo + (tsd + tss * rsoil * rdd) * too) * rsoil / dn
  rsost <- rsos + tsstoo * rsoil
  rsot <- rsost + rsodt
}
LSTa<- list(rdot=rdot,rsot=rsot,rddt=rddt,rsdt=rsdt)
return(LSTa)

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
############################################################################################################
############################################################################################################

}
