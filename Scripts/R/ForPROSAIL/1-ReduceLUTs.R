
rm(list=ls())

##############################################################################################################################
#	0. load main Libraries   -----
##############################################################################################################################

if (!require("RColorBrewer")) { install.packages("RColorBrewer"); require("RColorBrewer") }  ### colors
if (!require("parallel")) { install.packages("parallel"); require("parallel") }  ### Paralell
if (!require("doParallel")) { install.packages("doParallel"); require("doParallel") }  ### Paralell foreach and caret
# My packages in R
if (!require("ToolsRTM")) { install.packages("ToolsRTM"); require("ToolsRTM") }

if (!require("dplyr")) { install.packages("dplyr"); require("dplyr") }  ###
if (!require("data.table")) { install.packages("data.table"); require("data.table") }  ###
if (!require("zoo")) { install.packages("zoo"); require("zoo") }  ###
# NOTE: signal / caretEnsemble / neuralnet / sigmoid were required here but
# never actually used anywhere in this script -- dropped to cut install weight.
# hsdar (speclib/spectralResampling) was also dropped: it was archived from
# CRAN on 2023-07-06 and can no longer be installed via install.packages() on
# current R versions. The two operations this script actually needed from it
# (Gaussian-response spectral resampling, wavelength-range masking) are
# reimplemented below as plain-matrix helpers -- see 'resample_spectra()' and
# 'mask_wavelengths()'.

use.dev.source <- TRUE  # TRUE = load ToolsRTM from local source (ToolsRTM/R); FALSE = use the installed package
if (use.dev.source) devtools::load_all("../../../ToolsRTM/R")

out_root <- "../../../outs/ForPROSAIL"  # outputs always go under <repo root>/outs/, never inside Scripts/
dir.create(file.path(out_root, "LUTs"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_root, "Sims"), showWarnings = FALSE, recursive = TRUE)

## --- Lightweight replacement for hsdar::speclib()/spectralResampling() -----
## Resample a reflectance matrix (rows = spectra, columns = `from_wl`) onto a
## new wavelength grid `to_wl`, using a Gaussian response function of the
## given fwhm centred at each `to_wl` value -- the same principle hsdar's
## spectralResampling() used internally, without the speclib object wrapper.
resample.fun <- function(center, wl, fwhm) {
  a <- dnorm(wl, mean = center, sd = fwhm / 2)
  a <- (a - min(a)) / (max(a) - min(a))
  a
}
resample_spectra <- function(refl_matrix, from_wl, to_wl, fwhm) {
  response <- t(sapply(to_wl, resample.fun, wl = from_wl, fwhm = fwhm))
  response <- response / rowSums(response)
  resampled <- as.matrix(refl_matrix) %*% t(response)
  colnames(resampled) <- to_wl
  resampled
}
## Replacement for hsdar::mask(x) <- c(a1,b1,a2,b2,...): TRUE/FALSE keep-mask
## over `wl`, excluding each paired [a,b] range (e.g. water-absorption bands).
mask_wavelengths <- function(wl, ranges) {
  keep <- rep(TRUE, length(wl))
  for (k in seq(1, length(ranges), by = 2)) keep[wl >= ranges[k] & wl <= ranges[k + 1]] <- FALSE
  keep
}

##############################################################################################################################
#	1. Open  LUT with correlations  -----
##############################################################################################################################

#LUT = read.table(file.path(out_root, 'LUTs/INFORM_LUT50k_v23.csv'), sep=',', header = T)
LUT = read.table(file.path(out_root, 'LUTs/LUT_PROSAIL_with_n20k.csv'), sep=',', header = T)
nSamples =dim(LUT)[1]
print(nSamples)
head(LUT)


LUT.sb <-LUT[!duplicated(LUT), ]
dim(LUT.sb)


df.plot <-LUT
plot.rfl <-cbind(df.plot[c(1,3:13)], stack(df.plot[22:188]))
colnames(plot.rfl)<-c(names(plot.rfl)[1:12],'RFL','wave')
plot.rfl$wave =as.numeric(sub("R.", "", plot.rfl$wave, fixed = TRUE))
head(plot.rfl)

plot.rfl$GroupID <- cut(plot.rfl$Cab, breaks=c(0,20,30,40,60,80,Inf),labels = c('0-20', '20-30', '30-40','40-60','60-80','>80'))
GroupID<-plot.rfl$GroupID

plot.rfl.mean <- aggregate(plot.rfl[,c(1:13)], by = list(GroupID=plot.rfl[,'GroupID'], wave=plot.rfl[,'wave']), FUN = "mean")
head(plot.rfl.mean)

##############################################################################################################################
# 2.   Reflectance plots   ----
##############################################################################################################################

# output folder
paths.outs <- paste0(out_root, "/Sims/")
ifelse(!dir.exists(paths.outs), dir.create(paths.outs, recursive = TRUE), FALSE)


spectral_plot<-ggplot(plot.rfl.mean, aes(x=wave, y=RFL, group=GroupID)) +theme_bw() +
  geom_line(aes(color=as.factor(GroupID)),linetype = "dashed", size=0.7) +
  labs(color='Cab values',x='Wavelength',y='Reflectance')+  ggtitle("Degradation by chlorophyll")+
  geom_point(aes(color=as.factor(GroupID)), size=0.4) +
  theme(plot.title = element_text(hjust = 0.5, size=12)) #legend.title=element_blank()

print(spectral_plot)
ggsave(paste(paths.outs,'1-Reflectance_spectra_byCab.png',sep=''),
      width = 12, height = 8,  dpi = 300,units = "cm")



##############################################################################################################################
# 3.   Filtering LUT based on the observations  ----
##############################################################################################################################

field_data_path <- file.path(out_root, "FieldData/2019_2020RawData_Vcmax_withRFL_resampled.csv")
if (!file.exists(field_data_path)) {
  stop("Field data not found at '", field_data_path, "'. This is produced by running ",
       "1-GetSimulationsLUTs.R through step 9, which itself needs Carlos' private 2019/2020 ",
       "Vcmax field campaign data (see Scripts/ForPROSAIL/data/README.md).")
}
data.sensor = read.table(field_data_path, sep=',', header = T)
dim(data.sensor)

bands.rfl<- names(data.sensor[,grep(colnames(data.sensor),pattern="R.",fixed = TRUE)])
wave_=as.numeric(substr(bands.rfl,3,6))
dim(data.sensor)

data.simulations<-LUT
dim(data.simulations)


min_=as.data.frame(apply(data.simulations[,bands.rfl],2,min))
max_=as.data.frame(apply(data.simulations[,bands.rfl],2,max))
mean_=as.data.frame(apply(data.simulations[,bands.rfl],2,mean))
rfl.plot.sim=cbind(min_, max_,mean_)
names(rfl.plot.sim)=c('rfl.min', 'rfl.max','rfl.mean')
rfl.plot.sim[,'wave']=wave_
head(rfl.plot.sim)

min_=as.data.frame(apply(data.sensor[,bands.rfl],2,min))
max_=as.data.frame(apply(data.sensor[,bands.rfl],2,max))
mean_=as.data.frame(apply(data.sensor[,bands.rfl],2,mean))
rfl.plot.sensor=cbind(min_, max_,mean_)
names(rfl.plot.sensor)=c('rfl.min', 'rfl.max','rfl.mean')
rfl.plot.sensor[,'wave']=wave_


axis_x<-bquote(bold('reflectance')) # axis x
axis_y<-bquote(bold('wavelegth (nm)')) # axis x
statsLabel = paste0('')

#Plot
spectral.plot<-ggplot(rfl.plot.sim , aes(x = wave, y = rfl.min))+
  theme_bw()+
  geom_ribbon(data=rfl.plot.sim,aes(x=wave, ymin=rfl.min, ymax=rfl.max), fill = "grey80",alpha=0.6) +
  geom_ribbon(data=rfl.plot.sensor, aes(x=wave, ymin=rfl.min, ymax=rfl.max),  fill = "black",alpha=0.2) +

  geom_line(data=rfl.plot.sim, aes(x=wave, y=rfl.min),color = "navyblue",alpha=0.5,lty=2) +
  geom_line(data=rfl.plot.sim, aes(x=wave, y=rfl.mean), color = "navyblue",alpha=0.5,lty=2) +
  geom_line(data=rfl.plot.sim, aes(x=wave, y=rfl.max), color = "navyblue",alpha=0.5,lty=2) +

  geom_line(data=rfl.plot.sensor, aes(x=wave, y=rfl.max), color = "black") +
  geom_line(data=rfl.plot.sensor, aes(x=wave, y=rfl.mean), color = "black") +
  geom_line(data=rfl.plot.sensor, aes(x=wave, y=rfl.min), color = "black") +
  labs(color = '', x='Wavelength',y='Reflectance')+  ggtitle("spectral domain") +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size=10,face="bold"),
        axis.title = element_text(face="bold", size=10),
        axis.text.y=element_text(hjust = 0.5, size=10,face="bold"),
        axis.text.x=element_text(hjust = 0.5, size=10,face="bold"),
        legend.title=element_blank())
print(spectral.plot)

ggsave(paste0(out_root, '/Sims/1-Comparison_PROSAIL_Observations_withMinMax.png'),
       width = 10, height = 10,  dpi = 300,units = "cm")

##############################################################################################################################
# 4.   Filtering LUT based on the observations  ----
##############################################################################################################################

data.simulations<-LUT
dim(data.simulations)

for (i in 1:length(bands.rfl)) {
  print(bands.rfl[i])
  bands_i=bands.rfl[i]
  wave=as.numeric(substr(bands.rfl[i],3,6))

  if (wave <= 520){

    print('filter in the blue region ...')
    min=0.5*min(data.sensor[bands_i])

    max=2.2*max(data.sensor[bands_i])

  } else if (between(wave,520+1, 670)){

    print('filter in the green-red regions...')
    min=0.7*min(data.sensor[bands_i])
    max=1.4*max(data.sensor[bands_i])
  } else if (between(wave,670+1, 750)){

    print('filter in the red-edge region ...')
    min=0.005*min(data.sensor[bands_i])
    max=1.9*max(data.sensor[bands_i])
  } else {
    print('filter in the NIR ...')
    min=0.75*min(data.sensor[bands_i])
    max=1.2*max(data.sensor[bands_i])

  }
  data.simulations = data.simulations[data.simulations[bands_i] > min & data.simulations[bands_i] < max,]

  print(paste('dimensions is: ',dim(data.simulations)[1]))
}

min_=as.data.frame(apply(data.simulations[,bands.rfl],2,min))
max_=as.data.frame(apply(data.simulations[,bands.rfl],2,max))
mean_=as.data.frame(apply(data.simulations[,bands.rfl],2,mean))
rfl.plot.sim=cbind(min_, max_,mean_)
names(rfl.plot.sim)=c('rfl.min', 'rfl.max','rfl.mean')
rfl.plot.sim[,'wave']=wave_



#Plot
axis_x<-bquote(bold('reflectance')) # axis x
axis_y<-bquote(bold('wavelegth (nm)')) # axis x
statsLabel = paste0('')

spectral.plot<-ggplot(rfl.plot.sim , aes(x = wave, y = rfl.min))+
  theme_bw()+
  geom_ribbon(data=rfl.plot.sim,aes(x=wave, ymin=rfl.min, ymax=rfl.max), fill = "grey80",alpha=0.6) +
  geom_ribbon(data=rfl.plot.sensor, aes(x=wave, ymin=rfl.min, ymax=rfl.max),  fill = "black",alpha=0.2) +

  geom_line(data=rfl.plot.sim, aes(x=wave, y=rfl.min),color = "navyblue",alpha=0.5,lty=2) +
  geom_line(data=rfl.plot.sim, aes(x=wave, y=rfl.mean), color = "navyblue",alpha=0.5,lty=2) +
  geom_line(data=rfl.plot.sim, aes(x=wave, y=rfl.max), color = "navyblue",alpha=0.5,lty=2) +

  geom_line(data=rfl.plot.sensor, aes(x=wave, y=rfl.max), color = "black") +
  geom_line(data=rfl.plot.sensor, aes(x=wave, y=rfl.mean), color = "black") +
  geom_line(data=rfl.plot.sensor, aes(x=wave, y=rfl.min), color = "black") +
  labs(color = '', x='Wavelength',y='Reflectance')+  ggtitle("spectral domain") +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size=10,face="bold"),
        axis.title = element_text(face="bold", size=10),
        axis.text.y=element_text(hjust = 0.5, size=10,face="bold"),
        axis.text.x=element_text(hjust = 0.5, size=10,face="bold"),
        legend.title=element_blank())
print(spectral.plot)

ggsave(paste0(out_root, '/Sims/1-Comparison_PROSAIL_Observations_withMinMax_filtered_n5-10k.png'),
       width = 10, height = 10,  dpi = 300,units = "cm")



##############################################################################################################################
# 5.   Generated a LUT based on the observations  ----
##############################################################################################################################


LUT.new=ToolsRTM::getLUTs(inputs=data.simulations[,c(1,3:21)], dependencies='Car',nLUT=20000,setseed = F)
head(LUT.new)
dim(LUT.new)
names(LUT.new)

LUT.new['Cab']<- scales::rescale(LUT.new[,'Cab'], to = c(5, 60))
LUT.new['Car']<- scales::rescale(LUT.new[,'Car'], to = c(0, 10))
LUT.new['LIDFa']<- scales::rescale(LUT.new[,'LIDFa'], to = c(60, 70))

ggplot(LUT.new , aes(x = Cab, y = Car)) +
  theme_bw() + geom_point(aes(), size=0.4)

data <- ToolsRTM::dataSpec_PDB
Rsoil1  <- data[,11]  # rsoil1 = dry soil
Rsoil2 <- data[,12]
rsoil.filter<-list()
for (k in c(1:dim(LUT.new)[1])){
  rsoil<- c(LUT.new[k,'psoil']*Rsoil1+(1-LUT.new[k,'psoil'])*Rsoil2)
  rsoil.filter[[k]]<-rsoil#

}
##############################################################################################################################
# 3.Get Simulations  from a single day ----
##############################################################################################################################

## choose number of processors/cores
use.parallel <- TRUE  # FALSE = run sequentially (safer on machines where spawning worker processes is unstable)
no_cores <- max(1, detectCores() - 2)
if (use.parallel) {
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
} else {
  foreach::registerDoSEQ()
}

start_time <- Sys.time()
sim.rfl<-list()
sims<-foreach(i=1:nSamples) %dopar% {
  ###

  data.foursail_pro<-ToolsRTM::foursail(inputLUT=LUT.new[i,],rsoil=rsoil.filter[[i]],LeafModel = 'PROSPECT-PRO')

  rdot<-data.foursail_pro[[1]]
  rsot<-data.foursail_pro[[2]]
  data.foursail_pro<-ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT.new[i,'tts'],data.light =ToolsRTM::dataSpec_PDB)

  sim.rfl[[i]]<-data.foursail_pro


} ##end paralle

if (use.parallel) stopCluster(cl)
end_time <- Sys.time()
print(end_time - start_time)


sim.canopy<-do.call(rbind,sims)
data <- ToolsRTM::dataSpec_PDB
wave<-data[,1]

# hsdar::speclib()/idSpeclib()/SI()/mask() replaced with plain matrices + the
# resample_spectra()/mask_wavelengths() helpers defined near the top of this
# script (hsdar was archived from CRAN in 2023, see note above). SI(Spec.simula)
# <- LUT in the original was dead metadata -- LUT.new (not LUT) is what's
# actually joined onto the resampled reflectance below.
IDs<-c(1:nSamples)
keep.simula <- mask_wavelengths(wave, c(300,470,841,2600))

##############################################################################################################################
# 4.   Resample simulatons to specfic sensors  ----
############################################################################################################################### Usage of Speclib with spectral response values

##### VNIR
data.spectra<-as.matrix(data.sensor[,bands.rfl])

center <- wave_
fwhm   <- mean(diff(sort(unique(wave_))))  # approximate FWHM as mean band spacing (no explicit sensor FWHM metadata available)
wl     <- seq(468,842,2)

## Resample the simulated canopy spectra (masked to the sensor's VNIR range) onto the sensor's band grid
Spec.simula_resample.vnir <- resample_spectra(sim.canopy[, keep.simula, drop = FALSE], wave[keep.simula], wl, fwhm)

##############################################################################################################################
# 5.   Save Simulations and LUTs  ----
##############################################################################################################################

rfl.hypertoExport<-as.data.frame(Spec.simula_resample.vnir)
colnames(rfl.hypertoExport)<-paste('R.',wl,sep='')
rfl.bands = paste('R.',wl,sep='')

LUT.rfl.added<-cbind(ID=IDs, LUT.new,rfl.hypertoExport)
dim(LUT.rfl.added)

LUT.rfl.added.Indices<-getIndices(data=LUT.rfl.added,  pattern.rfl = "R.")
n_indices=dim(LUT.rfl.added)[2]+1

indices.names<-names(LUT.rfl.added.Indices)[n_indices:dim(LUT.rfl.added.Indices)[2]]

write.table(LUT.rfl.added.Indices,file=paste0(out_root,'/LUTs/LUT_PROSAIL_with_n',nSamples/1000,'k_filtered.csv'),sep=',',row.names = F)

variable_inputs <- c('ID','Cab','Car','Anth','CBC','EWT','Prot','N','LIDFa','LAI')
LUT.rfl.added.Indices.sb <- LUT.rfl.added.Indices[,c(variable_inputs,rfl.bands,indices.names)]
write.table(LUT.rfl.added.Indices.sb,file=paste0(out_root,'/LUTs/LUT_PROSAIL_with_n',nSamples/1000,'k_sb_MainTraits_filtered.csv'),sep=',',row.names = F)


##############################################################################################################################
# 6.   Reflectance plot ----
##############################################################################################################################



LUT.to_plot<-LUT.rfl.added.Indices #subset(LUT.rfl.added.Indices, Car  > 10)# & Cab >=15)
hist(LUT.to_plot[,'Cab'])
dim(LUT.to_plot)
min_=as.data.frame(apply(LUT.to_plot[,bands.rfl],2,min))
max_=as.data.frame(apply(LUT.to_plot[,bands.rfl],2,max))
mean_=as.data.frame(apply(LUT.to_plot[,bands.rfl],2,mean))
rfl.plot.sim=cbind(min_, max_,mean_)
names(rfl.plot.sim)=c('rfl.min', 'rfl.max','rfl.mean')
rfl.plot.sim[,'wave']=wave_



#Plot
axis_x<-bquote(bold('reflectance')) # axis x
axis_y<-bquote(bold('wavelegth (nm)')) # axis x
statsLabel = paste0('')

spectral.plot<-ggplot(rfl.plot.sim , aes(x = wave, y = rfl.min))+
  theme_bw()+
  geom_ribbon(data=rfl.plot.sim,aes(x=wave, ymin=rfl.min, ymax=rfl.max), fill = "grey80",alpha=0.6) +
  geom_ribbon(data=rfl.plot.sensor, aes(x=wave, ymin=rfl.min, ymax=rfl.max),  fill = "black",alpha=0.2) +

  geom_line(data=rfl.plot.sim, aes(x=wave, y=rfl.min),color = "navyblue",alpha=0.5,lty=2) +
  geom_line(data=rfl.plot.sim, aes(x=wave, y=rfl.mean), color = "navyblue",alpha=0.5,lty=2) +
  geom_line(data=rfl.plot.sim, aes(x=wave, y=rfl.max), color = "navyblue",alpha=0.5,lty=2) +

  geom_line(data=rfl.plot.sensor, aes(x=wave, y=rfl.max), color = "black") +
  geom_line(data=rfl.plot.sensor, aes(x=wave, y=rfl.mean), color = "black") +
  geom_line(data=rfl.plot.sensor, aes(x=wave, y=rfl.min), color = "black") +
  labs(color = '', x='Wavelength',y='Reflectance')+  ggtitle("spectral domain") +
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size=10,face="bold"),
        axis.title = element_text(face="bold", size=10),
        axis.text.y=element_text(hjust = 0.5, size=10,face="bold"),
        axis.text.x=element_text(hjust = 0.5, size=10,face="bold"),
        legend.title=element_blank())
print(spectral.plot)

ggsave(paste0(out_root, '/Sims/2-Comparison_PROSAIL_Observations_withMinMax_fromLUT_filtered_20k.png'),
       width = 10, height = 10,  dpi = 300,units = "cm")




df.plot <-LUT.rfl.added.Indices
plot.rfl <-cbind(df.plot[c(1,3:13)], stack(df.plot[22:188]))

colnames(plot.rfl)<-c(names(plot.rfl)[1:12],'RFL','wave')
plot.rfl$wave =as.numeric(sub("R.", "", plot.rfl$wave, fixed = TRUE))
head(plot.rfl)

plot.rfl$GroupID <- cut(plot.rfl$LIDFa, breaks=c(30,40,50,60,Inf),labels = c('30-40', '40-50','50-60','>60'))
GroupID<-plot.rfl$GroupID

plot.rfl.mean <- aggregate(plot.rfl[,c(1:13)], by = list(GroupID=plot.rfl[,'GroupID'], wave=plot.rfl[,'wave']), FUN = "mean")
head(plot.rfl.mean)


spectral_plot<-ggplot(plot.rfl.mean, aes(x=wave, y=RFL, group=GroupID)) +theme_bw() +
  geom_line(aes(color=as.factor(GroupID)),linetype = "dashed", size=0.7) +
  labs(color='LIDFa values',x='Wavelength',y='Reflectance')+  ggtitle("Degradation by LIDFa")+
  geom_point(aes(color=as.factor(GroupID)), size=0.4) +
  theme(plot.title = element_text(hjust = 0.5, size=12)) #legend.title=element_blank()

print(spectral_plot)
ggsave(paste(paths.outs,'2-Reflectance_spectra_byLIDFa_usingLUT_filtered.png',sep=''),
       width = 12, height = 8,  dpi = 300,units = "cm")

