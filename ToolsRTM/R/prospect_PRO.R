# _______________________________________________________________________
#
# prospect_PRO.m
# version 7.0 (January, 7th 2020)
# subroutines required: calctav.m, dataSpec_PRO.m
# _______________________________________________________________________
#
# Plant leaf reflectance and transmittance are calculated from 400 nm to
# 2500 nm (1 nm step) with the following parameters:
#
#       - N   = leaf structure parameter
#       - Cab = chlorophyll a+b content in ?g/cm2
#       - Car = carotenoids content in ?g/cm?
#       - Anth = Anthocyanin content in nmol/cm2
#       - Cbrown= brown pigments content in arbitrary units
#       - EWT  = equivalent water thickness in g/cm? or cm
#       - LMA  = dry matter content in g/cm2
#       - Prot = protein content g/cm2
#       - CBC= non protein dry matter content in g/cm2
#
# Here are some examples observed during the LOPEX'93 experiment on
# fresh (F) and dry (D) leaves :
#
# ---------------------------------------------
#                N     Cab     EWT        LMA    
# ---------------------------------------------
# min          1.000    0.0  0.004000  0.001900
# max          3.000  100.0  0.040000  0.016500
# corn (F)     1.518   58.0  0.013100  0.003662
# rice (F)     2.275   23.7  0.007500  0.005811
# clover (F)   1.875   46.7  0.010000  0.003014
# laurel (F)   2.660   74.1  0.019900  0.013520
# ---------------------------------------------
# min          1.500    0.0  0.000063  0.0019
# max          3.600  100.0  0.000900  0.0165
# bamboo (D)   2.698   70.8  0.000117  0.009327
# lettuce (D)  2.107   35.2  0.000244  0.002250
# walnut (D)   2.656   62.8  0.000263  0.006573
# chestnut (D) 1.826   47.7  0.000307  0.004305
# ---------------------------------------------
# _______________________________________________________________________

#'  PROSPECT-PRO
#'  
#' @param N numeric. Leaf structure parameter
#' @param Cab numeric. Chlorophyll content (microg.cm-2)
#' @param Car numeric. Carotenoid content (microg.cm-2)
#' @param Anth numeric. Anthocyain content (microg.cm-2)
#' @param Cbrown numeric. Brown pigment content (Arbitrary units)
#' @param EWT numeric. Equivalent Water Thickness (g.cm-2)
#' @param LMA numeric. Leaf Mass per Area (g.cm-2)
#' @param alpha numeric. Solid angle for incident light at surface of leaf
#' @param Prot numeric. protein content  (g.cm-2)
#' @param CBC numeric. NonProtCarbon-based constituent content (g.cm-2)

#'
#' @return List of lambda with leaf directional-hemisphrical reflectance and transmittance 
#' @importFrom expint expint
#' @export
#' 
#' @references
#' 
#' Féret J-B, Gitelson AA, Noble SD & Jacquemoud S, 2017. PROSPECT-D: Towards modeling leaf optical properties through a complete lifecycle. Remote Sensing of Environment, 193, 204–215. https://doi.org/10.1016/j.rse.2017.03.004
#'
#' Féret, J.B., Berger, K., de Boissieu, F., Malenovský, Z., 2021. PROSPECT-PRO for estimating content of nitrogen-containing leaf proteins and other carbon-based constituents. Remote Sens. Environ. 252. https://doi.org/10.1016/j.rse.2020.112173
#'
#' Jacquemoud S, Baret F, Hanocq J-F, 1992. Modeling spectral and bidirectional soil reflectance. Remote Sensing of Environment, 41, 123–132. https://doi.org/10.1016/0034-4257(92)90072-R
#'
#' Jacquemoud, S., Baret, F., 1990. PROSPECT: a model of leaf optical properties spectra. Remote Sens. Environ. 34, 75–91. https://doi.org/10.1016/0034-4257 (90)90100-Z.
#' 
#'Authors:
#'
#'Jean-Baptiste Feret (jb.feret@teledetection.fr)
#'
#'Frédéric Baret (baret@avignon.inra.fr)
#'
#'Stephane Jacquemoud  (jacquemoud@ipgp.fr)
#'
#'This function includes numerical optimizations proposed in the FLUSPECT code
#'
#'Authors: 
#'
#'Wout Verhoef
#'
#'Christiaan van der Tol (c.vandertol@utwente.nl)
#'
#'Joris Timmermans
#'
#' 
prospect_PRO<-function(N,Cab,Car,Anth,Cbrown,EWT,LMA,alpha,Prot,CBC){  

data <- ToolsRTM::dataSpec_PRO # read.table('parameters/dataSpec_PRO.csv',header = T, sep=',')
lambda  <- data[,1] ##wavelenght
nr      <- data[,2] ## refractive index of leaf material
Kab     <- data[,3] ## specific absorption coefficient of chlorophyll (a+b) (cm2.microg-1) 
Kcar    <- data[,4] ## specific absorption coefficient of carotenoids (cm2.microg-1)   
Kant    <- data[,5] ## specific absorption coefficient of Anthocyanins (cm2.nmol-1) 
KBrown  <- data[,6] ## specific absorption coefficient of brown pigments (arbitrary units) 
Kw      <- data[,7] ## specific absorption coefficient of water (cm-1)   
Km      <- data[,8] ## specific absorption coefficient of dry matter (cm2.g-1)  
Kprot   <- data[,9] ## specific absorption coefficient of proteins (cm2.g-1)   
Knonprot<- data[,10] ## specific absorption coefficient of non proteic dry matter (cm2.g-1)  

if (LMA > 0 & (Prot > 0 | CBC > 0)) {
  message('With LMA=0 proteins and CBC should be set to 0')
  message('LMA = PROT + CBC')
  message("The  prospect-pro model is processing with LMA=0")
  LMA <- 0
}


Kall    <- (Cab * Kab + Car * Kcar + Anth * Kant + Cbrown * KBrown + EWT * Kw + LMA * Km + Prot * Kprot + CBC * Knonprot) / N

############################################################################################ 
## This is for INFORM model
############################################################################################ 
# ## When N = 0 and is used for estimating reflectance understory
# Kall_t<-c()
# if (Kall[1] == 'Inf') { ## when dividing by cero
#   Kall[which(Kall == 'Inf')]<-1
# } else {
#     for (i in c(1:length(Kall))){
#     Kall_t[i]=(1-Kall[i])*exp(-Kall[i])+Kall[i]^(2)*ToolsRTM::s13aaf(Kall[i])
#     }
#   Kall = Kall_t
# }
############################################################################################ 
############################################################################################ 
j       <- which(Kall > 0)# Non-conservative scattering (normal case)
t1      <- (1 - Kall)*exp(-Kall)
t2      <- Kall^2*expint::expint(Kall)
#t2      <- Kall^2*expint(Kall)
tau     <- rep(1, length(t1))
tau[j]  <- t1[j] + t2[j]

# ***********************************************************************
# reflectance and transmittance of one layer
# ***********************************************************************
# Allen W.A., Gausman H.W., Richardson A.J., Thomas J.R. (1969),
# Interaction of isotropic ligth with a compact plant leaf, J. Opt.
# Soc. Am., 59(10):1376-1379.
# ***********************************************************************
# reflectivity and transmissivity at the interface
#-------------------------------------------------
#talf    <- calctav(40,nr) ##default alpha=40
talf    <- calctav(alpha,nr)
ralf    <- 1 - talf
t12     <- calctav(90,nr)
r12     <- 1 - t12
t21     <- t12 / (nr^2)
r21     <- 1 - t21


# top surface side
denom   <- 1 - (r21 * r21 * (tau^2))
Ta      <- (talf * tau * t21) / denom
Ra      <- ralf + (r21 * tau * Ta)

# bottom surface side
t       <- t12 * tau * t21 / denom
r       <- r12 + (r21 * tau * t)


# ***********************************************************************
# reflectance and transmittance of N layers
# Stokes equations to compute properties of next N-1 layers (N real)
# Normal case
# ***********************************************************************
# Stokes G.G. (1862), On the intensity of the light reflected from
# or transmitted through a pile of plates, Proc. Roy. Soc. Lond.,
# 11:545-556.
# ***********************************************************************
D       <- sqrt((1 + r + t) * (1 + r - t) * (1 - r + t) * (1 - r - t))
rq      <- r^2
tq      <- t^2
a       <- (1 + rq - tq + D)/(2 * r)
b       <- (1-rq + tq + D)/(2 * t)


bNm1    <- b^(N - 1)#
bN2     <- bNm1^2
a2      <- a^2
denom   <- a2 * bN2 - 1
Rsub    <- a * (bN2 - 1) / denom
Tsub    <- bNm1 * (a2 - 1) / denom

# Case of zero absorption
j       <- which(r + t >= 1)
Tsub[j] <- t[j] / (t[j] + (1 - t[j]) * (N - 1))
Rsub[j]	 <-  1 - Tsub[j]

# Reflectance and transmittance of the leaf: combine top layer with next N-1 layers
denom   <- 1 - Rsub * r
tran    <- Ta * Tsub / denom
refl    <- Ra + (Ta * Rsub * t) / denom

LRT <- list(lambda = lambda, refl = refl, tran = tran)
return(LRT)
}
