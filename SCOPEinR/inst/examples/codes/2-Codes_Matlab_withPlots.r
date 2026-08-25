
rm(list=ls())

##############################################################################################################################
#	0. Libraries   -----    
##############################################################################################################################
if (!require("hsdar")) { install.packages("hsdar"); require("hsdar") }  ### hsdar for PROSAIL
if (!require("RColorBrewer")) { install.packages("RColorBrewer"); require("RColorBrewer") }  ### colors
if (!require("pls")) { install.packages("pls"); require("pls") }  ### PLSR
if (!require("signal")) { install.packages("signal"); require("signal") }  ### interpolations
if (!require("prospectr")) { install.packages("prospectr"); require("prospectr") }  ## for resample2
if (!require("MASS")) { install.packages("MASS"); require("MASS") }  ## for smoth t the rfl from leaves
if (!require("caret")) { install.packages("caret"); require("caret") }  ##  models (random forest)

if (!require("nnet")) { install.packages("nnet"); require("nnet") }  ## nNEt models
if (!require("NeuralNetTools")) { install.packages("NeuralNetTools"); require("NeuralNetTools") }  ## nNEt models
if (!require("neuralnet")) { install.packages("neuralnet"); require("neuralnet") }  ## nNEt models
if (!require("e1071")) { install.packages("e1071"); require("e1071") }  ### SVM model
if (!require("parallel")) { install.packages("parallel"); require("parallel") }  ### Paralell
if (!require("doParallel")) { install.packages("doParallel"); require("doParallel") }  ### Paralell foreach and caret


#######################################################################################################################################
######################## 0.   Get spectra from HyperVNIR ----     
##############################################################################################################################
dataset<-read.table('Tables/DataBase/7-Database_Mallorca_SVMachine.csv', head=T,sep="," ,na.strings = "NA")

dim(dataset)
head(dataset)
dataset[1:10,1:10]

##### RFL
rfl.bands<-names(dataset[,grep(colnames(dataset),pattern="RFL.",fixed = TRUE)])
length(rfl.bands)

dataset.matrix<-as.matrix(dataset[,rfl.bands])
bandVar=names(dataset)[grep("^RFL.",names(dataset))]
wave.sensor<-as.numeric(gsub('RFL.', '', bandVar))

IDs<-dataset$ID
SpecDataset <- speclib(dataset.matrix, wave.sensor)
idSpeclib(SpecDataset) <- as.character(IDs)
SI(SpecDataset) <- dataset[,c(2:98)]
mask(SpecDataset)<-c(801,990,1098, 1169,1311,1498)
plot(SpecDataset)

##### RAD
rad.bands<-names(dataset[,grep(colnames(dataset),pattern="L.7",fixed = TRUE)])
rad.bands<-rad.bands[11:38]
length(rad.bands)
dataset.matrix.rad<-as.matrix(dataset[,rad.bands])
bandVar.rad=names(dataset)[grep("^L.",names(dataset))]
bandVar.rad<-bandVar.rad[4:31]
wave.sensor.rad<-as.numeric(gsub('L.', '', bandVar.rad))


IDs<-dataset$ID
SpecDataset.rad <- speclib(dataset.matrix.rad*1000, wave.sensor.rad)
idSpeclib(SpecDataset.rad) <- as.character(IDs)
SI(SpecDataset.rad) <- dataset[,c(2:98)]
plot(SpecDataset.rad, ylab='Radiance sensor')

##########################################################################################
### 1. Read tables -----  
##########################################################################################
 
### Read 20  Tables/SCOPE_outputs//verification_run_2020-11-03-1932
List<-list.files("Tables/SCOPE_outputs/",full.names =T)
print(List)
output_code<-3 ### 20k with RFL soil corrected (good)
Data_input<-read.table(paste(List[output_code],'/pars_and_input.dat',sep=""),header=T)
dim(Data_input)
######################################################################################
## units in nm for RFL
wave=c(seq(from=400, to=2400,by=1),seq(from=2500, to=15000,by=100),seq(from=16000, to=50000,by=1000))

######################################################################################
###Solar(Esun) and sky (Esky) irradiance above the canopy
#irradiance in (W m-2 um-1)
#Rin*(fEsun+fEsky)
######################################################################################
Esky_<-as.matrix(read.table(paste('Tables/SCOPE_outputs/data/input/radiationdata/Esky_.dat',sep=""),header=F,skip=0,na.strings = "0000"))
Esun_<-as.matrix(read.table(paste('Tables/SCOPE_outputs/data/input/radiationdata/Esun_.dat',sep=""),header=F,skip=0,na.strings = "0000"))
wave=c(seq(from=400, to=2400,by=1),seq(from=2500, to=15000,by=100),seq(from=16000, to=50000,by=1000))

Data_Irrad<-as.matrix(read.table(paste(List[output_code],'/irradiance_spectra.dat',sep=""),header=F,skip=2,na.strings = "0000"))
SpecIrrad <- speclib(Data_Irrad[,c(1:2001)], wave[c(1:2001)])
plot(SpecIrrad,ylab='Irradiance')

##  Plot Irradiance

axis_x<-expression(bold('Wavelength (nm)'))
axis_irrad<-expression(bold('E'~ '('~W~m^-2~mu~m^-1~sr^-1~ ')'))

par(mfrow=c(1,1),  mar = c(5,5,1.1,1),bg= "white", 
    font.main=1.0, cex.main=1.2,font.axis=2, cex.axis=1.0, las=1,
    font.lab=2, cex.lab=1.0)

plot(NA,NA, type='l',lwd=2,lty=1, xlim=c(400,2000),ylim=c(0,1500),xlab=axis_x,ylab=axis_irrad)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
par(new=T)
plot(wave,Esun_, type='l',lwd=2,lty=1, xlim=c(400,2000),ylim=c(0,1500),xlab=axis_x,ylab=axis_irrad, col='tomato')
par(new=T)
plot(wave,Esky_, type='l',lwd=2,lty=1, xlim=c(400,2000),ylim=c(0,1500),xlab=axis_x,ylab=axis_irrad,col='turquoise3')
par(new=T)
plot(SpecIrrad,lwd=2,lty=2,type='l',col='navyblue',xlim=c(400,2000),ylim=c(0,1500),xlab=axis_x,ylab=axis_irrad)

legend(x="topright",ncol=1,cex=1.2, 
       c('Eo','Esun', 'Esky'),
       fill =c('navyblue','tomato','turquoise3'),text.col="black",border = "black", bty = "y")



#######################################################################################
### Radiance
##radiance W m-2 sr-1 um-1  (rad.Lo) spectrum of reflected radiance Wm-2 microm-1 sr-1
#######################################################################################

rad.Lo<-as.matrix(read.table(paste(List[output_code],'/spectrum_obsdir_optical.dat',sep=""),header=F,skip=2,na.strings = "0000"))
SpecLo <- speclib(rad.Lo[,c(1:2001)], wave[c(1:2001)])
plot(SpecLo,ylab='radiance')

# rad.LoF_.('fluorescence.dat) 
#fluorescence per simulation for wavelengths of 640 to 850 nm, with 1 nm resolution 
#Rad.LoF_, this is the spectrum of fluorescence radiance (Wm-2 microm-1 sr-1), from 640-850 nm
rad.LoF<-as.matrix(read.table(paste(List[output_code],'/fluorescence.dat',sep=""),header=F,skip=2,na.strings = "0000"))
wave_rad.LoF=seq(from=640, to=850,by=1)
SpecradLoF <- speclib(rad.LoF, wave_rad.LoF)
plot(SpecradLoF, ylab='fluorescence')

##fluorescence coming directly from sunlit leaves
#Units in W m-2um-1 sr-1
rad.LoF_sunlit_<-as.matrix(read.table(paste(List[output_code],'/fluorescence_sunlit.dat',sep=""),header=F,skip=2,na.strings = "0000"))
wave_rad.LoF=seq(from=640, to=850,by=1)
SpecradLoF_sunlit <- speclib(rad.LoF_sunlit_, wave_rad.LoF)
plot(SpecradLoF_sunlit, ylab='fluorescence sunlit leaves')

##fluorescence coming directly from shaded leaves
#Units in W m-2um-1 sr-1
rad.LoF_shaded_<-as.matrix(read.table(paste(List[output_code],'/fluorescence_shaded.dat',sep=""),header=F,skip=2,na.strings = "0000"))
wave_rad.LoF=seq(from=640, to=850,by=1)
SpecradLoF_shaded <- speclib(rad.LoF_shaded_, wave_rad.LoF)
plot(SpecradLoF_shaded, ylab='fluorescence shaded leaves')

##fluorescence coming directly from PS
#Units in W m-2um-1 sr-1
rad.LoF_PS_<-as.matrix(read.table(paste(List[output_code],'/fluorescence_emitted_by_all_photosystems.dat',sep=""),header=F,skip=2,na.strings = "0000"))
wave_rad.LoF=seq(from=640, to=850,by=1)
SpecradLoF_PS <- speclib(rad.LoF_PS_, wave_rad.LoF)
plot(SpecradLoF_PS, ylab='fluorescence PSs')

## soil and leaves fluorescence_scattered
rad.LoF_sc_<-as.matrix(read.table(paste(List[output_code],'/fluorescence_scattered.dat',sep=""),header=F,skip=2,na.strings = "0000"))
wave_rad.LoF=seq(from=640, to=850,by=1)
SpecradLoF_sc <- speclib(rad.LoF_sc_, wave_rad.LoF)
plot(SpecradLoF_sc, ylab='fluorescence soil and leaves scattered')

## Emitted fluorescence irradiance by all leaves together
#Units in W m-2um-1
Fem_<-as.matrix(read.table(paste(List[output_code],'/fluorescence_emitted_by_all_leaves.dat',sep=""),header=F,skip=2,na.strings = "0000"))
wave_rad.LoF=seq(from=640, to=850,by=1)
SpecradFem_ <- speclib(Fem_, wave_rad.LoF)
plot(SpecradFem_, ylab='Emitted fluorescence irradiance ')

## Total fluorescence by all leaves together
SpecradLoF_all<-as.matrix(spectra(SpecradLoF_sc))+as.matrix(spectra(SpecradLoF_shaded))+as.matrix(spectra(SpecradLoF_sunlit)) 
SpecradLoF_all <- speclib(SpecradLoF_all, wave_rad.LoF)


######################################################################################
#reflectance  #fraction of radiation in observation direction *pi / irradiance 
# Lo_*pi/ Rin*(fEsun+fEsky)
######################################################################################

### Reflectance 
Data_rfl<-as.matrix(read.table(paste(List[output_code],'/reflectance.dat',sep=""),header=F,skip=2,na.strings = "0000"))
SpecRfl <- speclib(Data_rfl[,c(1:2001)], wave[c(1:2001)])
plot(SpecRfl,ylab='reflectance')

#######################################################################################
### Reflectance pure
##########################################################################################

Esun_sp <- speclib(Esun_[c(1:2001),1], wave[c(1:2001)],)
Esky_sp <- speclib(Esky_[c(1:2001),], wave[c(1:2001)],)
Esun_m<-Esun_[,]
Esun_m<-matrix(Esun_m, nrow=20000, ncol=length(Esun_m), byrow=TRUE)
Esky_m<-Esky_[,]
Esky_m<-matrix(Esky_m, nrow=20000, ncol=length(Esky_m), byrow=TRUE)

SpecRfl_pure<-pi*(rad.Lo)/(Esun_m+Esky_m) ##
SpecRfl_pure <- speclib(SpecRfl_pure[,c(1:2001)], wave[c(1:2001)])
mask(SpecRfl_pure)<-c(1311,1498,1800,1940)
plot(SpecRfl_pure,ylab='reflectance pure')

#######################################################################################
### Reflectance pure
##########################################################################################
v1<-c(rep(0,240))
v1_m<-matrix(v1, nrow=20000, ncol=length(v1), byrow=TRUE)
v2<-c(rep(0,1711))
v2_m<-matrix(v2, nrow=20000, ncol=length(v2), byrow=TRUE)
rad.LoF_all.range<-cbind(v1_m,as.matrix(rad.LoF),v2_m)
dim(rad.LoF_all.range)

Total_rad<-rad.Lo+rad.LoF_all.range
Spec_Total_rad<-speclib(Total_rad[,c(1:2001)], wave[c(1:2001)])
plot(rad.LoF_all.range[1,],ylab='Radiance (total)')
plot(Spec_Total_rad,ylab='Radiance (total)')

axis_x<-expression(bold('Wavelength (nm)'))
axis_y<-expression(bold('Reflectance'))
plot(wave,rad.Lo[1,],type='l',lwd=2,col='navyblue',lty=1, xlim=c(740,800),ylim=c(40,150),xlab=axis_x,ylab=axis_y)
par(new=T)
plot(wave,Total_rad[1,],type='l',lwd=2,col='red',lty=1, xlim=c(740,800),ylim=c(40,150),xlab=axis_x,ylab=axis_y)


rfl_apparent<-pi*(Total_rad)/(Esun_m+Esky_m)
SpecRfl_apparant <- speclib(rfl_apparent[,c(1:2001)], wave[c(1:2001)])
mask(SpecRfl_apparant)<-c(1311,1498,1800,1940)
plot(SpecRfl_apparant,ylab='reflectance apparent')

axis_x<-expression(bold('Wavelength (nm)'))
axis_y<-expression(bold('Reflectance'))

##Plot RFL, RFL pure (free FS) and RFL Apparent(polluted with Fs)
plot(SpecRfl[2,],type='l',lwd=2,col='navyblue',lty=1, xlim=c(400,800),ylim=c(0,0.8),
     xlab=axis_x,ylab=axis_y)
par(new=T)
plot(SpecRfl_pure[2,],type='l',lwd=2,col='tomato',lty=1, xlim=c(400,800),ylim=c(0,0.8),
     xlab=axis_x,ylab=axis_y)
par(new=T)
plot(SpecRfl_apparant[2,],type='l',lwd=2,col='turquoise3',lty=1, xlim=c(400,800),ylim=c(0,0.8),
     xlab=axis_x,ylab=axis_y)

legend(x="topleft",ncol=1,cex=0.8, 
       c('rfl','rfl pure', 'rfl apparent'),
       fill =c('navyblue','tomato','turquoise3'),text.col="black",border = "black", bty = "y")





##########################################################################################
### 2. Plots Fluorescence ----- 
##########################################################################################

library(RColorBrewer)
color.d = brewer.pal(7, "Blues")
axis_x<-expression(bold('wavelength (nm)'))
axis_y<-expression(bold('Fluorescence emission'))
par(mfrow=c(1,1),  mar = c(5,5,1.1,1),bg= "white", 
    font.main=1.0, cex.main=1.0,font.axis=2, cex.axis=1.0, las=1,
    font.lab=2, cex.lab=1.0)

plot(NA,NA, lwd=2,lty=2,type='l',col='forestgreen',ylim=c(0,8),xlim=c(640,850),xlab=axis_x,ylab=axis_y)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
par(new=T)
plot(SpecradLoF_all, lwd=2,lty=2,type='l',col='black',ylim=c(0,8),xlim=c(640,850),xlab=axis_x,ylab=axis_y)
par(new=T)
plot(SpecradLoF_shaded, lwd=2,lty=2,type='l',col=color.d[1],ylim=c(0,8),xlim=c(640,850),xlab=axis_x,ylab=axis_y)
par(new=T)
plot(SpecradLoF_sunlit, lwd=2,lty=2,type='l',col=color.d[7],ylim=c(0,8),xlim=c(640,850),xlab=axis_x,ylab=axis_y)
par(new=T)
plot(SpecradLoF_sc, lwd=2,lty=2,type='l',col=color.d[3],ylim=c(0,8),xlim=c(640,850),xlab=axis_x,ylab=axis_y)

legend("topleft", legend = c(expression(bold('total')),expression(bold('shaded leaves')),
                              expression(bold('sunlit leaves')),expression(bold('scattered leaves'))),
       fill=c('black',color.d[c(1,7,3)]),cex=0.6)

##########################################################################################
### 3. Plots Reflectance ----- 
##########################################################################################

color.d = brewer.pal(7, "Blues")
axis_x<-expression(bold('wavelength (nm)'))
axis_y<-expression(bold('Reflectance'))
par(mfrow=c(1,1),  mar = c(5,5,1.1,1),bg= "white", 
    font.main=1.0, cex.main=1.2,font.axis=2, cex.axis=1.0, las=1,
    font.lab=2, cex.lab=1.0)

plot(NA,NA, lwd=2,lty=2,type='l',col='forestgreen',ylim=c(0,0.6),xlim=c(400,1700),xlab=axis_x,ylab=axis_y)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
par(new=T)
plot(SpecDataset, lwd=2,lty=2,type='l',col='black',ylim=c(0,0.6),xlim=c(400,1700),xlab=axis_x,ylab=axis_y)
par(new=T)
plot(SpecRfl, fun='mean',lwd=2,lty=2,type='l',col='forestgreen',ylim=c(0,0.6),xlim=c(400,1700),xlab=axis_x,ylab=axis_y)

legend("topright", legend = c(expression(bold('Hyperspectral')),expression(bold('SCOPE 20k'))),
       fill=c('black','forestgreen'),cex=0.8)
##########################################################################################
##########################################################################################
color.d = brewer.pal(7, "Blues")
axis_x<-expression(bold('wavelength (nm)'))
axis_y<-expression(bold('Reflectance'))
plot(NA,NA, lwd=2,lty=2,type='l',col='forestgreen',ylim=c(0,0.8),xlim=c(400,1700),xlab=axis_x,ylab=axis_y)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
par(new=T)
for (i in c(1:100)){
  par(new=T)
  plot(SpecRfl[i], fun='max',lwd=2,lty=2,type='l',col=color.d[3],ylim=c(0,0.8),xlim=c(400,1700),xlab=axis_x,ylab=axis_y)

  
}

par(new=T)
plot(SpecDataset, fun='max',lwd=2,lty=2,type='l',col='black',ylim=c(0,0.8),xlim=c(400,1700),xlab=axis_x,ylab=axis_y)
legend("topleft", legend = c(expression(bold('Hyperspectral')),expression(bold('SCOPE 20k'))),
       fill=c('black',color.d[3]),cex=0.8)


#save.image(paste('Tables/SCOPE_outputs/Z-Inverion_save.RData',sep=''))

##########################################################################################
### 4. RAD Resampled from Model  ----- 
##########################################################################################

## Savitzky-Golay
Total_rad.sgolay <- noiseFiltering(Spec_Total_rad, method="sgolay", n=25)
## 
## Resampling to fwhm VNIR sensor
resample.fun<-function(center, wl, fwhm)
{
  a <- dnorm(wl, mean = center, sd = fwhm/2)
  a <- (a-min(a))/(max(a) - min(a))
  return(a)
}

##### VNIR (RAD)
center <- wave.sensor.rad
fwhm   <- 6.4
wl     <- c(740:850)

## Create spectral response with gaussian density function
response.vnir <- speclib(t(sapply(center, resample.fun, wl, fwhm)), wl)
## Perform resampling
Total_rad.resample.vnir <- spectralResampling(Spec_Total_rad, response_function = response.vnir)
Total_rad.sgolayvnir <- spectralResampling(Total_rad.sgolay, response_function = response.vnir)

rfl.scope.sgolay<- as.matrix(spectra(Total_rad.sgolayvnir))
wave.hyper<-wave.sensor
wave.sb.sgolay<-Total_rad.sgolayvnir@wavelength
rfl.sc.sgolay<-list()
for (m in c(1:dim(rfl.scope.sgolay)[1])) {
  rfl.sc.sgolay[[m]]<-interp1(wave.sb.sgolay, rfl.scope.sgolay[m,], wave.hyper, extrap = T) ##radiance
}
rfl.sc.sgolay<-do.call(rbind, rfl.sc.sgolay)
dim(rfl.sc.sgolay)


##########################################################################################
##########################################################################################

color.d = brewer.pal(7, "Reds")
axis_x<-expression(bold('wavelength (nm)'))
axis_y<-expression(bold('Radiance'))

par(mfrow=c(1,1),  mar = c(5,5,0.4,0.2),bg= "white", 
    font.main=1.0, cex.main=1.2,font.axis=2, cex.axis=1.0, las=1,
    font.lab=2, cex.lab=1.0)

plot(NA,NA, lwd=2,lty=2,type='l',col='forestgreen',ylim=c(0,200),xlim=c(740,800),xlab=axis_x,ylab=axis_y)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")

for (i in c(1:25)){
par(new=T)
plot(Total_rad.resample.vnir[i],lwd=1,lty=2,type='b',col=color.d[3], xlim=c(740,800),ylim=c(0,200),xlab=axis_x,ylab=axis_y)
}
for (i in c(1:25)){
  par(new=T)
  plot(wave.hyper,rfl.sc.sgolay[i,],lwd=1,lty=2,type='b',col='blue', xlim=c(740,800),ylim=c(0,200),xlab=axis_x,ylab=axis_y)
}
par(new=T)
plot(SpecDataset.rad, fun='mean',col='black',ylim=c(0,200),xlim=c(740,800),xlab=axis_x,ylab=axis_y)

legend("topleft", legend = c(expression(bold('Hyperspectral')),expression(bold('SCOPE 20k'))),
       fill=c('black',color.d[3]),cex=0.8)

##########################################################################################
### 4. Get SIF from Model  ----- 
##########################################################################################

source("Codes_R/Functions/getFLD_Scope.R")

Vcmax<-c()
SIF_sim<-c()
SIF_sim.resample<-c()
SIF_sim.sgolay<-c()
SIF_hyper<-c()

n_cases<-20000
progress_bar = txtProgressBar(min=0, max=n_cases, style = 3, char="=")
for (m in c(1:n_cases)) {
  #print(m)
  Vcmax[m]<-Data_input[m,'Vcmo']
  rad.total.i<-Spec_Total_rad@spectra[m][341:451]
  wave.i<-SpecIrrad@wavelength[341:451]##wave from 740:850
  irrad.i<-data.frame( wave=wave.i,Eo=SpecIrrad@spectra[m][341:451])
  SIF_sim[m]<-getFLD2(rad.total.i,wave.i,irrad.i)
  setTxtProgressBar(progress_bar, value = m)
}

close(progress_bar)

#################################################################################
#################################################################################

color.d = brewer.pal(7, "Blues")
axis_x<-expression(bold('Vcmax '~ '('~mu~mol~m^-2~s^-1~')'))
axis_y<-expression(bold('SIF '~ '('~W~m^-2~nm^-1~sr^-1~ ')'))
axis_A<-expression(bold('A '~ '('~mu~mol~m^-2~s^-1~')'))

r2<-round(cor(Vcmax,SIF_sim,use='pairwise.complete.obs')^2,2)
mylabel.r = bquote(bold(r)^2 == .(format(r2, digits = 3)))
plot(NA,NA,col='black',ylim=c(0.5,8),xlim=c(0.5,150),xlab=axis_x,ylab=axis_y)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
par(new=T)
plot(Vcmax,SIF_sim,pch=18,col=color.d[3], ylim=c(0.5,8),xlim=c(0,150),xlab=axis_x,ylab=axis_y)
abline(lm(SIF_sim~Vcmax),lwd=2,lty=2,col='black')
text(x = 25, y = 6, labels = mylabel.r, col='black')

##########################################################################################
##########################################################################################

###### SIF using resampled radiance
#### Here the range is 740:850
source("Codes_R/Functions/getFLD_Scope_resample.R")
progress_bar = txtProgressBar(min=0, max=n_cases, style = 3, char="=")
for (m in c(1:n_cases)) {
  #print(m)
  Vcmax[m]<-Data_input[m,'Vcmo']
  rad.total.i<-Total_rad.resample.vnir@spectra[m]
  wave.i<-SpecIrrad@wavelength[341:451]
  irrad.i<-data.frame(wave=wave.i,Eo=SpecIrrad@spectra[m][341:451])
  SIF_sim.resample[m]<-getFLD2_resample(rad.total.i, wave.sensor.rad,irrad.i)
  setTxtProgressBar(progress_bar, value = m)
}

close(progress_bar)

#################################################################################
#################################################################################

r2<-round(cor(Vcmax,SIF_sim.resample,use='pairwise.complete.obs')^2,2)
mylabel.r = bquote(bold(r)^2 == .(format(r2, digits = 3)))
plot(NA,NA,col='black',ylim=c(0,20),xlim=c(0.5,150),xlab=axis_x,ylab=axis_y)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
par(new=T)
plot(Vcmax,SIF_sim.resample,pch=18,col=color.d[3], ylim=c(0,20),xlim=c(0,150),xlab=axis_x,ylab=axis_y)
abline(lm(SIF_sim.resample~Vcmax),lwd=2,lty=2,col='black')
text(x = 25, y = 6, labels = mylabel.r, col='black')


##########################################################################################
##########################################################################################

##########################################################################################
### 5. Inversions with RAD Hyperspectral  ----- 
##########################################################################################

source("Codes_R/Functions/InversionOpt_v2.R")

### Method 1 using resmapling size of SCOPE
rfl.scope<- as.matrix(spectra(Total_rad.resample.vnir))
rfl.sensor<- as.matrix(spectra(SpecDataset.rad))
wave.sb<-Total_rad.resample.vnir@wavelength
rfl.sensor_m <- rfl.sensor

for (m in c(1:dim(rfl.sensor)[1])) {
  rfl.sensor_m[m,]<-interp1(wave.sensor.rad, rfl.sensor[m,], wave.sb, extrap = T) ##radiance
}
LUT<-Data_input

inv.RMSE_v1<-InversionOpt_nOpt(rfl.sensor=rfl.sensor_m[1:100,],rfl.prosail=rfl.scope,wave=wave.sb, 
                            n=n_cases,nOpt=20,method='merit-RMSE')
## for the best spectra from PROSAIL ( with min RMSE)
Table_lut_v1<-inv.RMSE_v1[[1]]
filename=paste('Tables/SCOPE_outputs/outs_inversions/1-Inversion_vcmax_v1.csv')
write.table(Table_lut_v1, file = filename, sep=",", row.names = F, col.names = T,append = F)
data.inv.v1<-cbind(dataset[1:100,],Table_lut_v1)
plot(data.inv.v1$Vcmo,data.inv.v1$SIF2)
##########################################################################################
##########################################################################################
#### best agrrement
### Method 2 using resmapling from VNIR
rfl.scope<- rfl.sc.sgolay
rfl.sensor<- as.matrix(spectra(SpecDataset.rad))
wave.sb<-wave.hyper
rfl.sensor_m <- rfl.sensor

LUT<-Data_input

inv.RMSE_v2<-InversionOpt_nOpt(rfl.sensor=rfl.sensor_m,rfl.prosail=rfl.scope,wave=wave.sb, 
                            n=n_cases,nOpt=20,method='merit-RMSE')
## for the best spectra from PROSAIL ( with min RMSE)
Table_lut_v2<-inv.RMSE_v2[[1]]
filename=paste('Tables/SCOPE_outputs/outs_inversions/1-Inversion_vcmax_v2_Rad_Sgolayn20k.csv')
write.table(Table_lut_v2, file = filename, sep=",", row.names = F, col.names = T,append = F)

data.inv.v2<-cbind(dataset,Table_lut_v2)
cols.remove<-'DNCabxc'
data.inv.v2 <- data.inv.v2[, ! names(data.inv.v2) %in% cols.remove, drop = F]
#data.inv.v2<-na.omit(data.inv.v2)
dim(data.inv.v2)
plot(data.inv.v2$Vcmo,data.inv.v2$SIF2)
##########################################################################################
### 5. Check the inversions ----- 
##########################################################################################

L.bands<-names(data.inv.v2[,grep(colnames(data.inv.v2),pattern="R.",fixed = TRUE)])
L.bands<-L.bands[3:30]
wave.sim<-as.numeric(gsub('R.', '', L.bands))
L.bands.sensor<-names(data.inv.v2[,grep(colnames(data.inv.v2),pattern="L.7",fixed = TRUE)])
L.bands.sensor<-L.bands.sensor[11:38]
wave.sensor<-as.numeric(gsub('L.', '', L.bands.sensor))
data.sim<-as.matrix(data.inv.v2[,L.bands])
data.hyper<-as.matrix(data.inv.v2[,L.bands.sensor])

##########################################################################################
##########################################################################################

axis_x<-expression(bold('Wavelength'))
axis_y<-expression(bold('Radiance '~ '('~W~m^-2~nm^-1~sr^-1~ ')'))
plot(NA,NA,col='black',ylim=c(50,150),xlim=c(740,800),xlab=axis_x,ylab=axis_y)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")

for (j in c(1)){
  par(new=T)
  plot(wave.sim,c(data.sim[j,]),type='l',col=color.d[3], ylim=c(50,150),xlim=c(740,800),xlab=axis_x,ylab=axis_y)
  par(new=T)
  plot(wave.sensor,c(data.hyper[j,]*1000),type='l',col='black', ylim=c(50,150),xlim=c(740,800),xlab=axis_x,ylab=axis_y)
}
##########################################################################################
##########################################################################################

######### Adding inversion to main dataset


 #data.inv.v2<-subset(data.inv.v2, Method == 'SVM')
data.inv.v2$Incidence[(data.inv.v2$SEV_July== 0)] <- 0
data.inv.v2$Incidence[!(data.inv.v2$SEV_July == 0)] <- 1
### Parcels 4, 11-14,16 (irrigated)
### Parcels 8-10,17 (rainfed)

data.inv.v2$Water[(data.inv.v2$Parcel ==8 | data.inv.v2$Parcel == 9 | data.inv.v2$Parcel == 10 | data.inv.v2$Parcel ==17)] <- 'Rainfed'
data.inv.v2$Water[!(data.inv.v2$Parcel ==8 | data.inv.v2$Parcel == 9 | data.inv.v2$Parcel == 10 | data.inv.v2$Parcel ==17)] <- 'Irrigated'

data.inv.v2$SEV[data.inv.v2$SEV_July == 0] <- 0
data.inv.v2$SEV[(data.inv.v2$SEV_July > 0   & data.inv.v2$SEV_July <= 1.5)] <- 1
data.inv.v2$SEV[(data.inv.v2$SEV_July > 1.5 & data.inv.v2$SEV_July <= 2.5)] <- 2
data.inv.v2$SEV[(data.inv.v2$SEV_July > 2.5 & data.inv.v2$SEV_July <= 3.5)] <- 3
data.inv.v2$SEV[data.inv.v2$SEV_July > 3.5] <- 4
##########################################################################################
##### Method 1 for estimating SIF using a single irradiance file1( 1 ASD file)
##########################################################################################

SIF_inverted_v1<-c()
Vcmax_inverted<-c()
data.ASD<-read.table(paste('Tables/ASD/smarts.Simulado/SMART_Correcion.csv',sep=""),sep=',',header=T)
data.ASD<-data.ASD[,c(1,3)]
progress_bar = txtProgressBar(min=0, max=n_cases, style = 3, char="=")
for (m in c(1:dim(data.inv.v2)[1])) {
  
  Vcmax_inverted[m]<-data.inv.v2[m,'Vcmo']
  rad.total.i<-c(data.sim[m,L.bands])
  wave.i<-wave.sim ##wave from 740:850
  SIF_inverted_v1[m]<-getFLD2_resample(rad.total.i,wave.i,data.ASD)
  setTxtProgressBar(progress_bar, value = m)
}

close(progress_bar)
data.inv.v2$SIF_inverted_m1<-SIF_inverted_v1/10
##########################################################################################
##########################################################################################

######### Adding SIF v1 to main dataset
axis_vcmax<-expression(bold('Vcmax '~ '('~mu~mol~m^-2~s^-1~')'))
axis_sif<-expression(bold('SIF '~ '('~W~m^-2~nm^-1~sr^-1~ ')'))


plot(NA,NA,col='black',ylim=c(0,8),xlim=c(0.5,150),xlab=axis_vcmax,ylab=axis_sif)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
par(new=T)
plot(data.inv.v2$Vcmo,data.inv.v2$SIF_inverted_m1,pch=18,col=color.d[3], ylim=c(0,8),xlim=c(0,150),xlab=axis_vcmax,ylab=axis_sif)
abline(lm(SIF_inverted_m1~Vcmo, data=data.inv.i),lwd=2,lty=2,col='black')
###############################################################
###### plots by Parcel Vcmo vs SIF_inverted_m1
###############################################################
for (k in unique(data.inv.v2$Parcel)){
  ii<-k
  data.inv.i<-subset(data.inv.v2, Parcel == ii)
  r2<-round(cor(data.inv.i$Vcmo,data.inv.i$SIF_inverted_m1,use='pairwise.complete.obs')^2,2)
  mylabel.r = bquote(bold(r)^2 == .(format(r2, digits = 3)))
  plot(NA,NA,col='black',ylim=c(0,8),xlim=c(0.5,150),xlab=axis_vcmax,ylab=axis_sif)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
  par(new=T)
  plot(data.inv.i$Vcmo,data.inv.i$SIF_inverted_m1,pch=18,col=color.d[3], ylim=c(0,8),xlim=c(0,150),xlab=axis_vcmax,ylab=axis_sif)
  abline(lm(SIF_inverted_m1~Vcmo, data=data.inv.i),lwd=2,lty=2,col='black')
  text(x = 25, y = 6, labels = mylabel.r, col='black')
}

##########################################################################
#######################   Boxplots Vcmax   #########################
##########################################################################
library(dplyr)
require(ggplot2)

color.d = brewer.pal(7, "Reds")

axis_vcmax<-expression(bold('Vcmax '~ '('~mu~mol~m^-2~s^-1~')'))
axis_sif<-expression(bold('SIF '~ '('~W~m^-2~nm^-1~sr^-1~ ')'))


data.inv.v2 %>% 
  ggplot(aes(x = SEV, y = Vcmo)) +
  theme_bw()+
  #geom_point(aes(color = factor(SEV)))    +
  geom_boxplot(aes(fill=factor(SEV)), width=0.3) +
  #scale_color_manual(values = color.d)+ 
  scale_fill_manual(values = color.d) +
  ylab (axis_vcmax) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14, face = "bold", colour = 'Black'),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.y = element_text(size = 14, face = "bold",colour = 'Black'),
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.box.background = element_rect(colour = "black"))+
  labs(fill = "SEV") +
  scale_y_continuous(limits = c(0,150))

##########################################################################
#######################   Boxplots SIF   #########################
##########################################################################
library(dplyr)
color.d = brewer.pal(7, "Reds")
axis_y<-expression(bold('V'['cmo']))
axis_x<-expression(bold('SEV'))
require(ggplot2)

data.inv.v2 %>% 
  ggplot(aes(x = SEV, y = SIF_inverted_m1)) +
  theme_bw()+
  #geom_point(aes(color = factor(SEV)))    +
  geom_boxplot(aes(fill=factor(SEV)), width=0.3) +
  #scale_color_manual(values = color.d)+ 
  scale_fill_manual(values = color.d) +
  ylab (axis_sif) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14, face = "bold", colour = 'Black'),
        axis.title.y = element_text(size=14, face="bold"),
        axis.text.y = element_text(size = 14, face = "bold",colour = 'Black'),
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.box.background = element_rect(colour = "black"))+
  labs(fill = "SEV") +
  scale_y_continuous(limits = c(0,7))

##########################################################################################
##########################################################################################
##########################################################################################
##### Method 2 for estimating SIF using irradiance files ( ASD files in each parcel))
##########################################################################################

path_irrad<-'Tables/ASD'
Files_ASD <- list.files(paste(path_irrad), pattern = "IRRAD", full.names = T, recursive = F)
#print(Files_ASD)
Parcels<-c(1,10,15,6,7,8,9,4,16,11,12,13,14,17)
#Parcels<-c(1,10,11)
files_irradiance<-c(5,5,5,2,2,1,1,4,4,6,6,6,6,4)
#files_irradiance<-c(5,5,6)

data_inverted<-cbind(data.inv.v2,data.sim[,L.bands])

Vcmax.inv.v2<-list()
SIF.inv.v2<-list()
ID.inv.v2<-list()
SEV.inv.v2<-list()
Parcel.inv.v2<-list()
Water.inv.v2<-list()

progress_bar = txtProgressBar(min=0, max=n_cases, style = 3, char="=")
for (p in c(1:length(Parcels))){
 # print(p)
  pp<-Parcels[p]
  data.ASD<-read.table(Files_ASD[files_irradiance[p]],sep=",",header = T)
  colnames(data.ASD)<-c("wave","E")
  data.ASD$E<-data.ASD$E 
  data.inv.p<-subset(data_inverted, Parcel == pp)
  L.bands.sb<-names(data.inv.p[,grep(colnames(data.inv.p),pattern="R.7",fixed = TRUE)])
  wave.sim.sb<-wave.sim
  data.sim.p<-as.matrix(data.inv.p[,L.bands])

  SIF_v2<-c()
  Vcmax_v2<-c()
    for (m in c(1:dim(data.inv.p)[1])) {
      
      Vcmax_v2[m]<-data.inv.p[m,'Vcmo']
    
      rad.total.m<-data.sim.p[m,]
      wave.m<-wave.sim ##wave from 740:850
      SIF_v2[m]<-getFLD2_resample(rad.total.m,wave.m,data.ASD)/10
    
    }
  ID.inv.v2[[p]]<-data.inv.p$ID
  SEV.inv.v2[[p]]<-data.inv.p$SEV
  Parcel.inv.v2[[p]]<-data.inv.p$Parcel
  Water.inv.v2[[p]]<-data.inv.p$Water
  Vcmax.inv.v2[[p]]<-Vcmax_v2
  SIF.inv.v2[[p]]<-SIF_v2
  data.parcel<-data.frame(Vcmo=Vcmax.inv.v2[[p]],
                          SIF_inverted=SIF.inv.v2[[p]],
                          SEV=data.inv.p$SEV )
  plot_p1<-ggplot(data.parcel, aes(x=Vcmo, y=SIF_inverted,color=as.factor(SEV))) + 
    scale_color_brewer(type = 'qual') +
    geom_point(shape=18)+    theme_bw()
  #print(plot_p1)
  setTxtProgressBar(progress_bar, value = p)
}

Vcmax_inverted_v2 <- unlist(Vcmax.inv.v2, recursive = F)
SIF_inverted_v2 <- unlist(SIF.inv.v2, recursive = F)
ID_inverted_v2 <- unlist(ID.inv.v2, recursive = F)
SEV_inverted_v2 <- unlist(SEV.inv.v2, recursive = F)
Parcel_inverted_v2 <- unlist(Parcel.inv.v2, recursive = F)
Water_inverted_v2 <- unlist(Water.inv.v2, recursive = F)
data.inverted<-data.frame(ID=ID_inverted_v2,
                          Vcmo=Vcmax_inverted_v2,
                          SIF_inverted=SIF_inverted_v2,
                          SEV=SEV_inverted_v2,
                          Parcel=Parcel_inverted_v2,
                          Water=Water_inverted_v2)

data.inv.v2$SIF_inverted_m2<-SIF_inverted_v2

dataArea<-read.table('Tables/DataBase/ML-quartilesTc/Database_Npredictions_SVM_CARS_Inversions_complete_withquartileTc_with_Area.csv', head=T,sep="," ,na.strings = "NA")
dataArea<-dataArea[,c('ID','Area','X0.','X25.','X50.','X75.','X100.','SD_Celsius','N_predict','Cab','Car','Anth','Cw','Cm','LAI','LIDFa','Proteins','CBC')]
colnames(dataArea)<-c('ID','Area','X0.','X25.','X50.','X75.','X100.','SD_Celsius','N_predict','Cab_ML','Car_ML','Anth_ML','Cw_ML','Cm_ML','LAI_ML','LIDFa_ML','Proteins_ML','CBC_ML')

data.inv.v3<-merge(data.inv.v2,dataArea, by.x='ID',by.y='ID',all.x=T,all.y=F)
dim(data.inv.v3)

data.inv.RF<-subset(data.inv.v3, Water == 'Rainfed')
data.inv.IR<-subset(data.inv.v3, Water == 'Irrigated')
filename=paste('Tables/SCOPE_outputs/outs_inversions/2-Table_withInversions_SCOPE_all.csv')

write.table(data.inv.v3, file = filename, sep=",", row.names = F, col.names = T,append = F)


#########################################################################
### Vcmax vs SIF2 (method 1 single ASD files) by Parcels
#######################################################################
axis_vcmax<-expression(bold('Vcmax '~ '('~mu~mol~m^-2~s^-1~')'))
axis_sif<-expression(bold('SIF '~ '('~W~m^-2~nm^-1~sr^-1~ ')'))

data.inv.v3%>%
  ggplot(aes(x=Vcmo, y=SIF_inverted_m1,color=Parcel)) + 
 # scale_color_gradient(low="blue", high="red")  +
  scale_color_gradientn(colours = rainbow(12)) +  
  ylab (axis_sif) +  xlab (axis_vcmax)+
  geom_point(shape=18)+    theme_bw() +
  geom_smooth(method=lm,  linetype="dashed",
              color="black") +
  theme( legend.position = "top",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.box.background = element_rect(colour = "black")) +
  scale_x_continuous(limits = c(0,150)) +
  scale_y_continuous(limits = c(2,5))


#########################################################################
### Vcmax vs SIF2 (method 1 single ASD files) by SEV
#######################################################################

axis_vcmax<-expression(bold('Vcmax '~ '('~mu~mol~m^-2~s^-1~')'))
axis_sif<-expression(bold('SIF '~ '('~W~m^-2~nm^-1~sr^-1~ ')'))

data.inv.v3 %>%
  ggplot(aes(x=Vcmo, y=SIF_inverted_m1,color=SEV)) + 
  scale_color_gradientn(colours = rainbow(5)) +  
  ylab (axis_sif)+  xlab (axis_vcmax)+
  geom_point(shape=18)+      theme_bw()+
  geom_smooth(method=lm,  linetype="dashed",
              color="black")+
  theme(axis.text.x = element_text(size = 10, face = "bold", colour = 'Black'),
        axis.title.y = element_text(size=10, face="bold"),
        axis.text.y = element_text(size = 10, face = "bold",colour = 'Black'),
        legend.position = "top") +
  scale_y_continuous(limits = c(2,5)) +
  scale_x_continuous(limits = c(0,150))



#########################################################################
### Vcmax vs Nitrogen (method 1 single ASD files) by SEV
#######################################################################

axis_vcmax<-expression(bold('Vcmax '~ '('~mu~mol~m^-2~s^-1~')'))
axis_n<-expression(bold('N '~ '(%)'))

data.inv.v3 %>%
  ggplot(aes(x=Vcmo, y=N_predict,color=SEV)) + 
  scale_color_gradientn(colours = rainbow(5)) +  
  ylab (axis_n)+  xlab (axis_vcmax)+
  geom_point(shape=18)+      theme_bw()+
  geom_smooth(method=lm,  linetype="dashed",
              color="black")+
  theme(axis.text.x = element_text(size = 10, face = "bold", colour = 'Black'),
        axis.title.y = element_text(size=10, face="bold"),
        axis.text.y = element_text(size = 10, face = "bold",colour = 'Black'),
        legend.position = "top") +
  scale_y_continuous(limits = c(1,3)) +
  scale_x_continuous(limits = c(0,150))
data.inv.RF%>%
  ggplot(aes(x=Vcmo, y=N_predict,color=SEV)) + 
  scale_color_gradientn(colours = rainbow(5)) +  
  ylab (axis_n)+  xlab (axis_vcmax)+
  geom_point(shape=18)+      theme_bw()+
  geom_smooth(method=lm,  linetype="dashed",
              color="black")+
  theme(axis.text.x = element_text(size = 10, face = "bold", colour = 'Black'),
        axis.title.y = element_text(size=10, face="bold"),
        axis.text.y = element_text(size = 10, face = "bold",colour = 'Black'),
        legend.position = "top") +
  scale_y_continuous(limits = c(1,3)) +
  scale_x_continuous(limits = c(0,150))


#########################################################################
### Vcmax vs SIF2 (method 1 single ASD files) by SEV
#######################################################################

axis_vcmax<-expression(bold('Vcmax '~ '('~mu~mol~m^-2~s^-1~')'))
axis_sif<-expression(bold('SIF '~ '('~W~m^-2~nm^-1~sr^-1~ ')'))
axis_Cab<-expression(bold('Cab '~ '('~u~gr~cm^-2~')'))

data.inv.v3 %>%
  ggplot(aes(x=Cab, y=SIF_inverted_m1,color=SEV)) + 
  scale_color_gradientn(colours = rainbow(5)) +  
  ylab (axis_sif)+  xlab (axis_Cab)+
  geom_point(shape=18)+      theme_bw()+
  geom_smooth(method=lm,  linetype="dashed",
              color="black")+
  theme(axis.text.x = element_text(size = 10, face = "bold", colour = 'Black'),
        axis.title.y = element_text(size=10, face="bold"),
        axis.text.y = element_text(size = 10, face = "bold",colour = 'Black'),
        legend.position = "top") +
  scale_y_continuous(limits = c(2,5)) +
  scale_x_continuous(limits = c(0,60))




###########################################################################
###########################################################################
data.inv.outliers<-subset(data.inv.v3, (Parcel <= 14 ))
data.inv.outliers<-subset(data.inv.outliers, (Area <=19 ))
unique(data.inv.outliers$Parcel)

axis_tc<-expression(bold('Tree-Crown Tc'))

ggplot(data.inv.outliers, aes(x=Vcmo, y=X50.,color=Parcel)) + 
  scale_color_gradientn(colours = rainbow(5)) +  
  ylab (axis_tc)+  xlab (axis_vcmax)+
  geom_point(shape=18)+    theme_bw()+
  geom_smooth(method=lm,  linetype="dashed",
              color="black")+
  theme(axis.text.x = element_text(size = 10, face = "bold", colour = 'Black'),
        axis.title.y = element_text(size=10, face="bold"),
        axis.text.y = element_text(size = 10, face = "bold",colour = 'Black'),
        legend.position = "top",
        legend.title = element_text(colour="black", size=10, 
                                    face="bold"),
        legend.text = element_text(colour="black", size=10, 
                                   face="bold"),
        strip.text.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white")) +
  scale_y_continuous(limits = c(30,45)) +
  scale_x_continuous(limits = c(0,150))


##################################################################################### 
### Plot Vcmax /SIF By Means Parcels and SEV
##################################################################################### 


axis_vcmax<-expression(bold('Vcmax '~ '('~mu~mol~m^-2~s^-1~')'))
axis_sif<-expression(bold('SIF '~ '('~W~m^-2~nm^-1~sr^-1~ ')'))

dataset.mean <-aggregate(data.inv.v3[,c('Vcmo','SIF2','SIF_inverted_m1','SIF_inverted_m2')], by=list(Parcel=data.inv.v3$Parcel),
                         FUN=mean, na.rm=TRUE)
dataset.sd <-aggregate(data.inv.v3[,c('Vcmo','SIF2','SIF_inverted_m1','SIF_inverted_m2')], by=list(Parcel=data.inv.v3$Parcel),
                         FUN=sd, na.rm=TRUE)
na.omit(dataset.mean)
na.omit(dataset.mean)
ii<-'Vcmo'
jj<-'SIF_inverted_m1'
# Default line plot
xsd<-dataset.sd[,ii]/2

ysd<-dataset.sd[,jj]/2


X<-dataset.mean[,ii]
Y<-dataset.mean[,jj]

### stats 
r2<-round(cor(X, Y,use='pairwise.complete.obs')^2,2)
mylabel.r2 = bquote(bold(r)^2 == .(format(r2, digits = 3)))
fit<-lm(Y~X)
inter=as.vector(fit$coefficients[1])
slope<-as.vector(fit$coefficients[2])
p<-ggplot(dataset.mean, aes_string(x=ii, y=jj, group='Parcel')) + 
  geom_point()+ theme_bw()+  ylab (axis_sif) +   xlab (axis_vcmax) + 
  scale_color_brewer(type = 'qual') +
  geom_errorbar(aes(ymin=Y-ysd, ymax=Y+ysd), width=.2,
                position=position_dodge(0.05)) +
  geom_errorbarh(aes(xmin = X-xsd, xmax = X+xsd)) +
  theme(legend.position = "top") +
  scale_y_continuous(limits = c(0,7)) +
  scale_x_continuous(limits = c(75,150))  +
  geom_abline(intercept = inter, slope = slope, color ='black',linetype = "dashed") + 
  annotate(geom="text", size=3.0,x=38, y=49, label=mylabel.r2,color="black")
print(p)













###################################################################################
### 6. Inversions with RFL Hyperspectral  ----- 
##########################################################################################
source("Codes_R/Functions/InversionOpt_v2.R")

dataset.matrix<-as.matrix(dataset[,rfl.bands])
bandVar=names(dataset)[grep("^RFL.",names(dataset))]
wave.sensor<-as.numeric(gsub('RFL.', '', bandVar))
##### VNIR (RFL)
center <- wave.sensor[1:44]
fwhm   <- 6.4
wl     <- c(400:801)
SpecRfl_apparant.vnir<-SpecRfl_apparant
mask(SpecRfl_apparant.vnir)<-c(810,2501)
Rfl_apparant.vnir<-as.matrix(spectra(SpecRfl_apparant.vnir))
SpecRfl_apparant.vnir <- speclib(Rfl_apparant.vnir, SpecRfl_apparant.vnir@wavelength)
plot(SpecRfl_apparant.vnir,ylab='reflectance apparent')

## Create spectral response with gaussian density function
response.vnir <- speclib(t(sapply(center, resample.fun, wl, fwhm)), wl)
## Perform resampling
SpecRfl_apparant.vnir.resampling <- spectralResampling(SpecRfl_apparant.vnir, response_function = response.vnir)

rfl.scope.resampling<- as.matrix(spectra(SpecRfl_apparant.vnir.resampling))
rfl.scope<- as.matrix(spectra(SpecRfl_apparant.vnir))

plot(rfl.scope.resampling[1,])
plot(rfl.scope[1,])

rfl.sensor<- as.matrix(spectra(SpecDataset))

wave.sb.resampling<-SpecRfl_apparant.vnir.resampling@wavelength
wave.sb<-SpecRfl_apparant.vnir@wavelength

rfl.sensor_m <- rfl.sensor
### selected 505-565 [13:18]
#for resampling -->wave.sb.resampling
#wave:505-565 (green area)
#wave: 640_850
wave_inversions<-c(400:800)
wave_inversions<-wave.sensor[1:44]
wave_inversions<-c(505:565)
rfl.sensor.inv<-list()
for (m in c(1:dim(rfl.sensor)[1])) {
  rfl.sensor.inv[[m]]<-interp1(wave.sensor[1:44], rfl.sensor[m,1:44], wave_inversions, extrap = T) ##radiance
}
rfl.sensor_m<-do.call(rbind, rfl.sensor.inv)


rfl.scope.inv<-list()
for (m in c(1:dim(rfl.scope.resampling)[1])) {
  rfl.scope.inv[[m]]<-interp1(wave.sb.resampling, rfl.scope.resampling[m,], wave_inversions, extrap = T) ##radiance
}
rfl.scope_m<-do.call(rbind, rfl.scope.inv)

##########################################################################################
##########################################################################################

LUT<-Data_input

inv.RMSE.v3<-InversionOpt_nOpt(rfl.sensor=rfl.sensor_m[c(1:500),],rfl.prosail=rfl.scope_m,wave=wave_inversions, 
                            n=n_cases,nOpt=20,method='merit-RMSE')

## for the best spectra from PROSAIL ( with min RMSE)
Table_lut_v3<-inv.RMSE.v3[[1]]
filename=paste('Tables/SCOPE_outputs/outs_inversions/1-Inversion_vcmax_rfl_app_Hyper_Nopt20_greenArea.csv')
write.table(inv.RMSE.v3, file = filename, sep=",", row.names = F, col.names = T,append = F)

##########################################################################################
##########################################################################################
data.inv.v3<-cbind(dataset[1:500,],Table_lut_v3)
data.inv.v3<-subset(data.inv.v3, Method == 'SVM')

data.inv.v3$Incidence[(data.inv.v3$SEV_July== 0)] <- 0
data.inv.v3$Incidence[!(data.inv.v3$SEV_July == 0)] <- 1
### Parcels 4, 11-14,16 (irrigated)
### Parcels 8-10,17 (rainfed)
data.inv.v3$Water[(data.inv.v3$Parcel ==8 | data.inv.v3$Parcel == 9 | data.inv.v3$Parcel == 10 | data.inv_v3$Parcel ==17)] <- 'Rainfed'
data.inv.v3$Water[!(data.inv.v3$Parcel ==8 | data.inv.v3$Parcel == 9 | data.inv.v3$Parcel == 10 | data.inv_v3$Parcel ==17)] <- 'Irrigated'

data.inv.v3$SEV[data.inv.v3$SEV_July == 0] <- 0
data.inv.v3$SEV[(data.inv.v3$SEV_July > 0   & data.inv.v3$SEV_July <= 1.5)] <- 1
data.inv.v3$SEV[(data.inv.v3$SEV_July > 1.5 & data.inv.v3$SEV_July <= 2.5)] <- 2
data.inv.v3$SEV[(data.inv.v3$SEV_July > 2.5 & data.inv.v3$SEV_July <= 3.5)] <- 3
data.inv.v3$SEV[data.inv.v3$SEV_July > 3.5] <- 4

##########################################################################################
##########################################################################################

axis_x<-expression(bold('Vcmax '~ '('~mu~mol~m^-2~s^-1~')'))
axis_y<-expression(bold('SIF '~ '('~W~m^-2~nm^-1~sr^-1~ ')'))


ggplot(data.inv.v3, aes(x=Vcmo, y=SIF2,color=as.factor(SEV))) + 
  scale_color_brewer(type = 'qual') +
  geom_point(shape=18)+    theme_bw()


plot(data.inv.v3$Vcmo,data.inv.v3$SIF2)

##########################################################################################
### 5. Check the inversions ----- 
##########################################################################################

L.bands<-names(data.inv.v3[,grep(colnames(data.inv.v3),pattern="R.",fixed = TRUE)])
L.bands<-L.bands[3:46]
wave.sim<-as.numeric(gsub('R.', '', L.bands))
L.bands.sensor<-names(data.inv.v3[,grep(colnames(data.inv.v3),pattern="RFL.",fixed = TRUE)])
L.bands.sensor<-L.bands.sensor[1:44]
wave.sensor<-as.numeric(gsub('RFL.', '', L.bands.sensor))
data.sim<-as.matrix(data.inv.v3[,L.bands])
data.hyper<-as.matrix(data.inv.v3[,L.bands.sensor])

##########################################################################################
##########################################################################################

axis_x<-expression(bold('Wavelength'))
axis_y<-expression(bold('RFL'))
plot(NA,NA,col='black',ylim=c(0,0.6),xlim=c(400,800),xlab=axis_x,ylab=axis_y)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")

for (j in c(1)){
  par(new=T)
  plot(wave.sim,c(data.sim[j,]),type='l',col=color.d[3],lwd=3,ylim=c(0,0.6),xlim=c(400,800),xlab=axis_x,ylab=axis_y)
  par(new=T)
  plot(wave.sensor,c(data.hyper[j,]),type='l',col='black', ylim=c(0,0.6),xlim=c(400,800),xlab=axis_x,ylab=axis_y)
}
