#' The INvertible FOrest Reflectance Model coupled with leaf models
#'
#' \code{inform} INFORM  simulation based on a set of combinations of inpuparameters
#' @param inputLUT LUT table with distribution of biophysical parameters used as input parameters in the model
#' @param rsoil numeric. Soil reflectance
#' @param LeafModel Version of PROSPECT model; 'PRO' or 'D' is accepted. By default 'PROSPECT-PRO' is used. Liberty model
#' and Fluspect-B and Fluspect-B-Cx is also provided. Options are : 'PROSPECT-PRO','PROSPECT-D', 'Liberty', 'Fluspect-B' and 'Fluspect-B-Cx'.
#'
#' @return rfl 
#' @export
#'
#' 
#' @references
#' 
#' INFORM (Atzberger, 2000; Schlerf&Atzberger, 2006) simulates the bi-directional reflectance
#' of forest stands between 400 and 2500 nm. 
#' 
#' INFORM is essentially an innovative combination of FLIM (Rosema et al., 1992), SAIL (Verhoef, 1984), and PROSPECT (Jacquemoud et al. 1996)
#' 
#' Atzberger, C. 2000: Development of an invertible forest reflectance model: The INFOR-Model.
#' In: Buchroithner (Ed.): A decade of trans-european remote sensing cooperation. Proceedings
#' of the 20th EARSeL Symposium Dresden, Germany, 14.-16. June 2000: 39-44.
#' 
#' Schlerf, M. & Atzberger, C. (2006): Inversion of a forest reflectance model to estimate biophysical
#' canopy variables from hyperspectral remote sensing data. Remote Sensing of Environment, 100: 281-294
#' 
#' Rosema, A., Verhoef, W., Noorbergen, H. 1992: A new forest light interaction model in support of forest
#' monitoring. Remote Sensing of Environment, 42: 23-41.
#' 
#' Jacquemoud S., Ustin S.L., Verdebout J., Schmuck G., Andreoli G., Hosgood B. (1996): Estimating leaf
#' biochemistry using the PROSPECT leaf optical properties model, Remote Sens. Environ., 56:194-202.
#' 
#' Verhoef, W. 1984: Light scattering by leaf layers with application to canopy reflectance modeling: The
#' SAIL model. Remote Sensing of Environment, 16: 125-141.
#' 
#' Basic version of INFORM: Clement Atzberger, 1999
#' INFORM modifications and validation: Martin Schlerf, 2004-2007
#' 
#' @author Clement Atzberger and Martin Schlerf (Original version in Matlab and python)
#' @author Carlos Camino (Ported version into R)
#' 

inform<- function(inputLUT=NULL,rsoil=rsoil, LeafModel='PROSPECT-PRO'){
  
  # _____________________________________________________________________________________________________________
  
  # SUBROUTINES
  #
  # M-File                    Function
  # inform_prospect.R         Main INFORM code
  # msail_background.R         SAIL-PROSPECT-SOIL to compute background reflectance (adapted with PRO-D)
  # msail_inf.R                SAIL-PROSPECT-SOIL to compute infinite crown reflectance (adapted with PRO-D)
  # sail_t_s.R                SAIL-PROSPECT-SOIL to compute crown transmittance for sun direction (adapted with PRO-D)
  # sail_t_o.R                SAIL-PROSPECT-SOIL to compute crown transmittance for observation direction (adapted with PRO-D)
  # prospect_PRO.R             Leaf reflecance model PROSPECT-PRO
  # prospect_D.R             Leaf reflecance model PROSPECT-D
  # s13aaf.M                  Integral for PROSPECT no uses in this code
  # calctav                     Refraction index for PROSPECT
  # data                  Absorption coefficients for PROSPECT
  
  # _____________________________________________________________________________________________________________
  
  # INPUT VARIABLES
  
  # PARA_PROSPECT: Leaf Input Parameters
  # Variable	                         Designation	    Unit	    Default value
  # leaf structure parameter              N               -           2
  # chlorophyll a+b content in            Cab             ug/cm-2     60
  # equivalent water thickness            Cw              g/cm-2      0.025
  # dry matter content                    Cm              g/cm-2       0.025
  
  # PARA_INFORM: Canopy Input Parameters
  # Variable	                          Designation       Unit	    Default value
  # Scale factor for soil reflectance     scale                       1
  # Single tree leaf area index           lai             m2 m-2      7
  # Leaf area index of understorey	    laiu	        m2 m-2	    0.1
  # Stem density                          sd              ha-1        650
  # Tree height                           h               m           20
  # Crown diameter                        cd              m           4.5
  # Average leaf angle of tree canopy	    LIDFa	            deg	        55
  
  
  # External Input Parameters
  # Variable	                          Designation       Unit	    Default value
  # Sun zenith angle 	                    tts	        deg	        30
  # Observation zenith angle 	            tto	        deg	        0
  # Azimuth angle	                        phi	            deg	        0
  # Fraction of diffuse radiation	        skyl	        fraction	0.1
  
  # Other Input Data
  # Variable	                          Designation
  # rsoil                                Soil spectrum
  # _____________________________________________________________________________________________________________
  
  # OUTPUT VARIABLES
  # Variable	                          Designation
  # Forest reflectance                    r_forest
  # Soil reflectance                      rsoil
  # Understorey reflectance               r_understorey
  # Infinite canopy reflectance           r_sail_inf
  # Crown closure                         co
  # Crown factor                          C
  # Ground factor                         G
  # Leaf reflectance                      r_leaf
  # Leaf transmittance                    t_leaf
  # Crown transmittance for teta_s        t_s
  # Crown transmittance for teta_o        t_o
  # _____________________________________________________________________________________________________________
  
if ((is.null(inputLUT)) | (is.null(rsoil)) ){
  message('Please add the missing parameters')
  stop()
}

## fourSAIL
LIDFa=inputLUT[,'LIDFa']; LIDFb=inputLUT[,'LIDFb']; TypeLidf=inputLUT[,'TypeLidf']; lai=inputLUT[,'LAI']
hotspot=inputLUT[,'hspot']; tts=inputLUT[,'tts']; tto=inputLUT[,'tto']; psi=inputLUT[,'psi']

## INform model

laiu=inputLUT[,'LAIu']
sd=inputLUT[,'sd']; cd=inputLUT[,'cd']; h=inputLUT[,'h']; psi=inputLUT[,'psi']
skyl=inputLUT[,'skyl']

# _____________________________________________________________________________________________________________

# Computing spectral components (SAIL, PROSPECT, LIBERTY models)

# Scaling of soil spectrum (to account for effects due to shadow and soil moisture)
#refl_soil <- scale*rsoil ## individual 

refl_soil = rsoil  ### we give refl_soil based on scale factor and rsoil

# Computing of understorey reflectance for dicotyledoneae (ala=45,hotspot=0, N=2, Cab=30, Cw=0.05 CM=0,05)
#r_understorey <- ToolsRTM::msail_background(inputLUT=inputLUT,rsoil=refl_soil,LeafModel= LeafModel, typeLAI = 'understorey') 
if (LeafModel == 'Fluspect-B' | LeafModel == 'Fluspect-B-Cx' ) {
  

  r_understorey <- ToolsRTM::foursail.inform(inputLUT=inputLUT,rsoil=refl_soil, typeLAI = 'understorey') 
  r_understorey = r_understorey[1:2001]  ### we give refl_soil based on scale factor and rsoil
} else {

  r_understorey <- ToolsRTM::foursail.inform(inputLUT=inputLUT,rsoil=refl_soil, typeLAI = 'understorey') 
  
}


if (LeafModel == 'PROSPECT-PRO') {
  ## Prospect inputs
  N=inputLUT[,'N']; Cab=inputLUT[,'Cab']; Car=inputLUT[,'Car']; Anth=inputLUT[,'Anth']; Cbrown=inputLUT[,'Cbrown']
  EWT=inputLUT[,'EWT']; LMA=inputLUT[,'LMA'];alpha=inputLUT[,'alpha']
  Prot=inputLUT[,'Prot'];CBC=inputLUT[,'CBC']
  
  #LeafModel = 'PRO'
  LRT<- prospect_PRO(N,Cab,Car,Anth,Cbrown,EWT,LMA,alpha,Prot,CBC)
  # Computing of leaf reflecance and transmittance
  r_leaf <- LRT[[2]] #rho Reflectance
  t_leaf <- LRT[[3]] #tau Transmittance
} else if (LeafModel == 'Liberty'){
  #define alll inputs in the models. retreived from LUT tables
  
  # run LeafModel ='Liberty'
  LRT <- ToolsRTM::liberty(inputLUT)
  r_leaf <- LRT[[2]] #rho Reflectance
  t_leaf <- LRT[[3]] #tau Transmittance
} else if (LeafModel == 'Fluspect-B'){
  #define alll inputs in the models. retreived from LUT tables
  
  # run LeafModel ='Fluspect-B'
  LRT <- getFluspect.B(inputsLeaf= inputLUT, inputsOptipar =ToolsRTM::optipar,version='D')
  message('SAIL with fluspect model is processing')
  r_leaf <- LRT[[2]] #rho Reflectance
  t_leaf <- LRT[[3]] #tau Transmittance
  
} else if (LeafModel == 'Fluspect-B-Cx'){
  #define alll inputs in the models. retreived from LUT tables
  
  # run LeafModel ='Fluspect-Cx'
  LRT <- getFluspect.Cx(inputsLeaf= inputLUT, inputsOptipar =ToolsRTM::optipar,version='Cx')
  message('SAIL with fluspect model is processing')
  r_leaf <- LRT[[2]] #rho Reflectance
  t_leaf <- LRT[[3]] #tau Transmittance
  
} else if (LeafModel == 'PROSPECT-D'){
  
  #LeafModel ='PROSPECT-D'
  N=inputLUT[,'N']; Cab=inputLUT[,'Cab']; Car=inputLUT[,'Car']; Anth=inputLUT[,'Anth']; Brown=inputLUT[,'Cbrown']
  EWT=inputLUT[,'EWT']; LMA=inputLUT[,'LMA'];alpha=inputLUT[,'alpha']
  
  LRT <- LRT <- prospect_DB(N,Cab,Car,Anth,Brown,EWT,LMA,alpha)
  # Computing of leaf reflecance and transmittance
  r_leaf <- LRT[[2]] #rho Reflectance
  t_leaf <- LRT[[3]] #tau Transmittance
}

# Computing of infinitive crown reflectance for a very dense forest canopy (LAI=15, hot=0.04, N=1.5)
if (LeafModel == 'Fluspect-B' | LeafModel == 'Fluspect-B-Cx' ) {
  r_sail_inf <- foursail.inf(inputLUT=inputLUT,rsoil=r_understorey[1:2001],rleaf=r_leaf[1:2001],tleaf=t_leaf[1:2001], short.waves=T)
  r_sail_inf<-r_sail_inf[1:2001]
} else {
  r_sail_inf <- ToolsRTM::foursail.inf(inputLUT=inputLUT,rsoil=r_understorey,rleaf=r_leaf,tleaf=t_leaf)
}


# _____________________________________________________________________________________________________________

# Ground coverage (FLIM model)

# Computation of the average horizontal area of a single tree crown in hectare (k) corrected
# by the factor 'adapt'. 'adapt' was adjusted in such a way, that modelled values of 'co'
# agreed with measured values of 'co' for a given stand structure (as defined by 'cd' and
# 'sd'). To compare measured with modelled stands for certain values of canopy lai (lai_c),
# an empirical relation between 'lai_c' and 'co' was derived from field data at Idarwald site
# (co=0.078*lai_c+0.23, n=40).

# adapt=1; # no correction
# adapt=0.6;

adapt <- 1
k <- adapt * (pi * (cd / 2)^2) / 10000

# angles (degree) to angles (radian)
tto_ <- tto * pi / 180 #teta_o
tts_ <- tts * pi / 180 #teta_s
psi_ <- psi * pi / 180



## Coverage and Shadowing 
# Observed ground coverage  by crowns (co) under observation zenith angle teta_o

co <- 1 - exp(-k * sd / cos(tto_) ) # eq 1 from Rosema et al 1992
# Ground coverage by shadow (cs) under a solar zenith angle teta_s
cs <- 1 - exp(-k * sd / cos(tts_) )
# Geometrical factor (g) depending on the illumination and viewing geometry
g <- ((tan(tto_)^2 + tan(tts_))^2 - 2 * tan(tto_) * tan(tts_) * cos(psi_))^(0.5)

# Correlation coefficient (p)

p <- exp(-g * h / cd)



# _____________________________________________________________________________________________________________

# Ground surface fractions (FLIM model)

# Tree crowns with shadowed background (Fcd)
Fcd <- co * cs + p* (co * (1 - co) * cs * (1 - cs))^0.5

# Tree crowns with sunlit background (Fcs)
Fcs <- co * (1 - cs) - p * (co * (1 - co)* cs * (1 - cs) )^(0.5)

# Shadowed open space (Fod)
Fod <- (1 - co ) * cs - p * (co * (1 - co) * cs * (1 - cs) )^(0.5)

# Sunlit open space (Fos)
Fos <- (1 - co) * (1 - cs) + p * (co * (1 - co) * cs * (1 - cs) )^(0.5)


# _____________________________________________________________________________________________________________

# Crown transmittance (SAIL model)

# Angles (radian) to angles (degree)
#tts <- tts * 180/pi
#tto <- tto * 180/pi
#psi <- psi * 180/pi
# Crown transmittance in sun direction (t_s)
t_s <- ToolsRTM::foursail_t_s(inputLUT=inputLUT,rsoil=r_understorey,rleaf=r_leaf,tleaf=t_leaf)


# Crown transmittance in observation direction (t_o)
t_o <- ToolsRTM::foursail_t_o(inputLUT=inputLUT,rsoil=r_understorey,rleaf=r_leaf,tleaf=t_leaf)


# _____________________________________________________________________________________________________________

# Forest reflectance (FLIM model)

# Ground factor (G), that is ground contribution to scene reflectance
G <- Fcd * t_s * t_o + Fcs * t_o + Fod * t_s + Fos

# Crown factor (C), that is crown contribution to scene reflectance
# C=(1-t_s.*t_o)*cs*co;  # Original formula

C <- Fcd*(1 - t_s * t_o) # Similar to original formula

# C=Fcd*(1-t_s.*t_o)+Fcs*(1-t_o); # Modification, rejected by W. Verhoef
# Der Faktor 'Fcs*(1-t_o)' wurde erg??????nzt. Er beschreibt den Beitrag, der von
# beleuchteten Kronen (sunlit crowns, Fcs) erfolgt abz??????glich des Beitrags
# des darunterliegenden Untergrundes (Fcs*t_o), der schon in G enthalten
# ist.

# Forest reflectance
r_forest = (r_sail_inf * C) + (r_understorey * G) #*10

return(r_forest)
}


