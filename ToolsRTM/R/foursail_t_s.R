#' Performs PROSAIL simulation based on a set of combinations of input parameters
#' 
#' for estimating Crown transmittance in sun direction (t_s)
#' 
#' @param inputLUT LUT table with distribution of biophysical parameters used as input parameters in the model. Must include a 'skyl' column (fraction of diffuse sky illumination, 0-1) in addition to the standard LAI/LIDFa/LIDFb/TypeLidf/tts/tto/psi columns — this is not part of the standard PROSAIL LUT structure and must be added explicitly.
#' @param rsoil numeric. Soil reflectance (here the reflectance of understorey)
#' @param rleaf leaf reflectance  
#' @param tleaf leaf transmittance 
#' 
#' @return Crown transmittance in sun direction (t_s)
#' 
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

#' 
foursail_t_s<- function(inputLUT,rsoil=r_understorey, rleaf=r_leaf,tleaf=t_leaf){

# Crown transmittance in observation direction (t_o)
  
lai=inputLUT[,'LAI']
hotspot=0
LIDFa=inputLUT[,'LIDFa']; LIDFb=inputLUT[,'LIDFb']; TypeLidf=inputLUT[,'TypeLidf']; 
tts=inputLUT[,'tts']; tto=inputLUT[,'tto']; psi=inputLUT[,'psi']
skyl=inputLUT[,'skyl']

# reflectance and transmittance

rho	 <- rleaf
tau	 <- 	tleaf

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

if (TypeLidf == 1){
#LeafDistribution <- dladgen(LIDFa,LIDFb)
#lidf <- LeafDistribution$lidf
#litab <- LeafDistribution$litab

lidf = c(0.015192247821401716, 0.04511513125626976, 0.07366721933874043, 0.09998095795183792, 0.12325683701156054, 0.14278760558390557, 0.1579798610972386, 0.16837196070089022, 0.034475081609994906, 0.034644632694447175, 0.03477199446646273, 0.03485697201080851, 0.03489949845644191,
                      0.00028378419874902573, 0.0020294058823425495, 0.005751638473401871, 0.01200202743766678, 0.02190796087454504, 0.03787410106306353, 0.06589705963563201, 0.12690963854753423, 0.04115845607276447, 0.05181790293362443, 0.06942677499148037, 0.106277087297172, 0.4586641625920237,
                      0.0011565951416267486, 0.008876829634519247, 0.029891034176080182, 0.09640335028409708, 0.7273443849356243, 0.09640334777606041, 0.029891033573473225, 0.008876829401686992, 0.000569757676969429, 0.0003409834748974161, 0.000173365433308037, 6.345407373831158e-05, 9.034417918662996e-06,
                      0.42712701495904715, 0.051885576869933114, 0.016954996965587388, 0.003890519113906976, 0.00028378419874897087, 0.0038905192522350474, 0.016954997336226962, 0.05188557794898474, 0.019952583631293264, 0.02655232968560639, 0.0375291397600469, 0.0606223355617268, 0.2824706247166563,
                      0.03789183438420043, 0.043216944609998746, 0.05490877250501834, 0.07525964096053808, 0.10741859382700167, 0.14983497911429472, 0.1813652598851876, 0.18082510675991648, 0.03458082011568342, 0.03411370192811436, 0.033742780197829614, 0.0334862403361027, 0.033355325376113854,
                      0.11111111416724702, 0.11111110780104928, 0.11111111416724703, 0.11111110780104927, 0.111111114167247, 0.11111110780104927, 0.11111111416724706, 0.11111110780104927, 0.022222225379928462, 0.022222219013730782, 0.022222225379928573, 0.02222221901373067, 0.022222223339496305)
tx1 <- c(10,20,30,40,50,60,70,80,82,84,86,88,90)
tx2 <- c(0,10,20,30,40,50,60,70,80,82,84,86,88)
litab<- (tx1 + tx2) / 2 

} else if (TypeLidf == 2){
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

		chi_s_chi_o_frho_ftau <- volscatt(tts,tto,psi,ttl)
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

		#	Extinction coefficientsksli
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
    
	sdb	 <-  0.5 * (ks + bf)
	sdf	 <-  0.5 * (ks - bf)
	dob	 <-  0.5 * (ko + bf)
	dof	 <-  0.5 * (ko - bf)
	ddb	 <-  0.5 * (1 + bf)
	ddf	 <-  0.5 * (1 - bf)

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
		rsodt	 <-  0*rsoil
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
      fint <- (1 - exp(-alf)) * 0.05
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
LSTa<- list(rdot,rsot,rddt,rsdt)

PARdifo = skyl;
PARdiro = 1 - skyl;

trans_hemi = (rdot * PARdiro + rsot * PARdifo) / (PARdiro + PARdifo)


return(trans_hemi)

}
