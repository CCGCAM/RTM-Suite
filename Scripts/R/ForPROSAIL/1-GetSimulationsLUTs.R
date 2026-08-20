#adjust view and ilumination angles, chl, LAI and Prot values, date for soil subset and date for sentinel-2 subset

rm(list=ls())

##############################################################################################################################
#	0. load main Libraries   -----
##############################################################################################################################

if (!require("RColorBrewer")) { install.packages("RColorBrewer"); require("RColorBrewer") }  ### colors
if (!require("signal")) { install.packages("signal"); require("signal") }  ### Savitzky-Golay smoothing
if (!require("parallel")) { install.packages("parallel"); require("parallel") }  ### Paralell
if (!require("doParallel")) { install.packages("doParallel"); require("doParallel") }  ### Paralell foreach and caret
# NOTE: rstanarm / ADGofTest / pspline were required here but never actually used
# anywhere in this script -- dropped to cut install weight. hsdar (speclib/
# spectralResampling/mask/noiseFiltering) was also dropped: it was archived
# from CRAN on 2023-07-06 and can no longer be installed via install.packages()
# on current R versions. What this script needs from it is covered directly
# by ToolsRTM's own functions below (get.spectral.convolution.gaussian() for
# resampling) or by base R / signal:: calls inline -- no local reimplementation.

library(tidyverse)
library(ggplot2)

# ToolsRTM isn't on CRAN -- install from GitLab if missing.
if (!requireNamespace("ToolsRTM", quietly = TRUE)) {
  remotes::install_gitlab("caminoccg/toolsrtm")
}
library(ToolsRTM)

##############################################################################################################################
#	1. Generate LUTs with correlations  -----
##############################################################################################################################

# Bundled with the ToolsRTM package (was 'Tables/LUTs/inputs/inputs_PROSAIL.csv', which
# doesn't exist in this repo -- the real file ships inside the package instead).
inputsPRO = read.table(system.file("LUts", "inputs_PROSAIL.csv", package = "ToolsRTM"),
                       sep=',', header = T, fileEncoding="UTF-8-BOM")
nSamples =2000


LUT<-as.data.frame(getLUT(inputs = inputsPRO, nLUT=nSamples, setseed = 1234))
ID<-c(1:nSamples)
LUT<-cbind(ID, LUT)
head(LUT)
dim(LUT)

#variable_inputs <- c('ID','Cab','Car','Anth','LMA','EWT','N','LAI','LIDFa','tts','tto')
#LUT.sb <- LUT[,c(variable_inputs)]

LUT.pigments<-ToolsRTM::getCor(n_inputs = 4,setseed = 1234,distribution = 'Uniform',nLUT = nSamples,rho = 0.95,
                               Varnames = c('Cab','Car','Anth','LAI'),
                               MinRange = c(0.5,0.1,0,0.0), MaxRange = c(90,25,7,7))

summary(LUT.pigments$LUT)
LUT.pigments$LUT[,'Cab']<- scales::rescale(LUT.pigments$LUT[,'Cab'], to = c(0.5, 80))
LUT.pigments$LUT[,'Car']<- scales::rescale(LUT.pigments$LUT[,'Car'], to = c(0.5, 20))
LUT.pigments$LUT[,3]<- scales::rescale(LUT.pigments$LUT[,3], to = c(0, 7)) #Anth
LUT.pigments$LUT[,'LAI']<- scales::rescale(LUT.pigments$LUT[,'LAI'], to = c(0.5, 7))

LUT$Cab<-LUT.pigments$LUT$Cab
LUT$Car<-LUT.pigments$LUT$Car
LUT$Anth<-LUT.pigments$LUT$Anth
LUT$LAI<-LUT.pigments$LUT$LAI


summary(LUT)

#set.seed(12345)
LUT$N<-runif(nSamples, min = 1.5, max = 3)
LUT$hspot<-runif(nSamples, min = 0, max = 1)
LUT$LMA<-0
LUT$EWT<- 0.0090
LUT$alpha <- 40 # default
LUT$Cbrown<-runif(nSamples, min = 0, max = 0.15)

#set.seed(12345)
LUT$tts<-runif(nSamples, min = 0, max = 0)
LUT$tto<-runif(nSamples, min = 0, max = 10)
LUT$psi<-runif(nSamples, min = 0, max = 45)


LUT$CBC<-runif(nSamples, min = 0.001,max = 0.03)
LUT$Prot<-0.0015 #runif(nSamples, min = 0.001,max = 0.03)
LUT$LIDFa<-runif(nSamples, min = 30,max = 70)

summary(LUT)
dim(LUT)


##############################################################################################################################
# 2.Get soil reflectance ----
##############################################################################################################################

data <- ToolsRTM::dataSpec_PDB
Rsoil1  <- data[,11]  # rsoil1 = dry soil
Rsoil2 <- data[,12]  # rsoil2 = wet soil

###### ### ### ### ### ### ### ### ###
### This spte is for varying psoil
###### ### ### ### ### ### ### ### ###
#psoil	 <-  1    # soil factor (psoil=0: wet soil / psoil=1: dry soil)

psoil	 <-  runif(nSamples, 0.5, 1)
rsoil0<-list()
for (k in c(1:nSamples)){
 rsoil<- c(psoil[k]*Rsoil1+(1-psoil[k])*Rsoil2)
 rsoil0[[k]]<-rsoil#
 #print(plot(rsoil0[[k]]))
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
  ### using SEntinel-2 spectra from soil

  data.foursail_pro<-ToolsRTM::foursail(inputLUT=LUT[i,],rsoil=rsoil0[[i]],LeafModel = 'PROSPECT-PRO') ##rsoil=rsoil0 [[i]] Si hay diferente RFL suelo en la LUT. SI fijamos la misma rsoil para toda la LUT quitar [[i]]

  rdot<-data.foursail_pro[[1]]
  rsot<-data.foursail_pro[[2]]
  data.foursail_pro<-ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT[i,'tts'],data.light =ToolsRTM::dataSpec_PDB)

  #sim.rfl[[i]]<-data.inform[[1]]
  sim.rfl[[i]]<-data.foursail_pro
  #sim.rfl[[i]]<-data.foursail2_pro

} ##end paralle

if (use.parallel) stopCluster(cl)
end_time <- Sys.time()
print(end_time - start_time)

##############################################################################################################################
# 4. Reflectance matrix + wavelength mask  ----
##############################################################################################################################
sim.canopy<-do.call(rbind,sims)
data <- ToolsRTM::dataSpec_PDB
wave<-data[,1]
IDs<-c(1:nSamples)
# keep only the 470-841nm VNIR range (excludes [300,470] and [841,2600])
keep.simula <- !(wave >= 300 & wave <= 470) & !(wave >= 841 & wave <= 2600)
#plot at 1 nm
#matplot(wave[keep.simula], t(sim.canopy[, keep.simula, drop = FALSE]), type = "l", lty = 1, xlim=c(440,850),ylim=c(0,0.6))

##############################################################################################################################
# 4b.  Save the native-resolution simulated LUT + indices  ----
##############################################################################################################################
# Doesn't need field data at all -- full-resolution (masked to ~470-841nm
# VNIR) simulated PROSAIL reflectance + vegetation indices. Saved regardless
# of whether Carlos' private field dataset (steps 5-10 below) is present, so
# outs/ForPROSAIL/ always has real output from a run of this script.

# Outputs always go under <repo root>/outs/, never inside Scripts/.
out_root <- "outs/ForPROSAIL"
dir.create(file.path(out_root, "LUTs"),     showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_root, "Sensor"),   showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_root, "Sims"),     showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_root, "FieldData"),showWarnings = FALSE, recursive = TRUE)

rfl.native <- as.data.frame(sim.canopy[, keep.simula, drop = FALSE])
colnames(rfl.native) <- paste0('R.', wave[keep.simula])

LUT_rfl.native <- cbind(ID = IDs, LUT, psoil = psoil, rfl.native)
LUT_rfl.native.withIndices <- getIndices(data = LUT_rfl.native, pattern.rfl = "R.")

native_lut_path <- paste0(out_root, "/LUTs/LUT_PROSAIL_native_n", nSamples / 1000, "k.csv")
write.table(LUT_rfl.native.withIndices, file = native_lut_path, sep = ',', row.names = F)
cat("Saved native-resolution simulated LUT (no field data needed) to '", native_lut_path, "'\n", sep = '')

##############################################################################################################################
# 4c.  Diagnostic plot: native reflectance spectra by Cab  ----
##############################################################################################################################
# Quick visual sanity check that the simulation looks physically reasonable
# (reflectance should drop in the visible as Cab increases, rise sharply at
# the red edge, etc.) -- print()ed so it displays if you're running this
# interactively (RStudio/console), and always saved to outs/ regardless.

rfl.cols.native <- grep("^R\\.", names(LUT_rfl.native.withIndices), value = TRUE)
rfl.cols.plot   <- rfl.cols.native[seq(1, length(rfl.cols.native), by = 5)]  # subsample wavelengths, this is just a quick-look plot
plot.rfl.native <- cbind(LUT_rfl.native.withIndices[, c("ID", "Cab")],
                          stack(LUT_rfl.native.withIndices[, rfl.cols.plot]))
colnames(plot.rfl.native) <- c("ID", "Cab", "RFL", "wave")
plot.rfl.native$wave <- as.numeric(sub("R.", "", plot.rfl.native$wave, fixed = TRUE))
plot.rfl.native$GroupID <- cut(plot.rfl.native$Cab, breaks = c(0, 20, 30, 40, 60, 80, Inf),
                                labels = c('0-20', '20-30', '30-40', '40-60', '60-80', '>80'))
plot.rfl.native.mean <- aggregate(RFL ~ GroupID + wave, data = plot.rfl.native, FUN = "mean")

spectral_plot_native <- ggplot(plot.rfl.native.mean, aes(x = wave, y = RFL, group = GroupID)) +
  theme_bw() +
  geom_line(aes(color = GroupID), linetype = "dashed", size = 0.7) +
  labs(color = 'Cab (ug/cm2)', x = 'Wavelength (nm)', y = 'Reflectance') +
  ggtitle(sprintf("Native simulated reflectance by Cab (n=%d)", nSamples)) +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

print(spectral_plot_native)
native_plot_path <- paste0(out_root, "/Sims/0-Native_reflectance_spectra_byCab.png")
ggsave(native_plot_path, plot = spectral_plot_native, width = 12, height = 8, dpi = 300, units = "cm")
cat("Saved diagnostic plot to '", native_plot_path, "'\n", sep = '')

##############################################################################################################################
# 5.   load dataset   ----
##############################################################################################################################

# This is Carlos' own 2019/2020 Vcmax field campaign data -- proprietary, not
# bundled in this repo. Drop your own copy at the path below (same columns:
# 'X<wavelength>' reflectance columns + a 'Vcmax' column) to run steps 5-10
# (resample the simulated spectra onto the field spectrometer's own bands,
# and compare against the measurements). Without it, only the native-resolution
# LUT saved above (step 4b) is produced -- everything below is field-data-only.
field_data_path <- "data/2019_2020RawData_Vcmax.csv"
if (!file.exists(field_data_path)) {
  message("Field data not found at '", field_data_path, "' -- skipping steps 5-10 (field ",
          "comparison, sensor-band resampling). The native-resolution LUT from step 4b above ",
          "was still saved. See Scripts/ForPROSAIL/data/README.md to enable steps 5-10.")
} else {

dataset = read.table(field_data_path, sep=',', header = T)
head(dataset)

wave.bands<- round(as.numeric(substr(names(dataset[,grep(colnames(dataset),pattern="X",fixed = TRUE)]),2,8)),4)
rfl.bands<- names(dataset[,grep(colnames(dataset),pattern="X",fixed = TRUE)])

data.spectra<-as.matrix(dataset[,rfl.bands])
wave_ =c(wave.bands)


jpeg(file.path(out_root, "Sensor/1-Check_RFL_Sensor.jpg"), units="in", width=5, height=5, res=300)

matplot(wave_, t(data.spectra), type = "l", lty = 1, xlab = "wave (nm)", ylab = "Reflectance",
        col='black', lwd=2, xlim=c(400,1000))
rect(470, 0.6, 840, 0, col = "grey", density = 60)
dev.off()

# keep only the 470-841nm VNIR range (excludes [300,470] and [841,2600])
keep.data <- !(wave_ >= 300 & wave_ <= 470) & !(wave_ >= 841 & wave_ <= 2600)
## predefined filter functions
## Savitzky-Golay (row-wise, via signal::sgolayfilt directly)
data.spectra.sgolay <- t(apply(as.matrix(data.spectra), 1, signal::sgolayfilt, p = 3, n = 25))
# same keep.data mask applies -- sgolay smoothing doesn't change the wavelength grid

jpeg(file.path(out_root, "Sensor/1-Check_RFL_Sensor_withSgolay.jpg"), units="in", width=5, height=5, res=300)
matplot(wave_[keep.data], t(data.spectra.sgolay[, keep.data, drop = FALSE]), type = "l", lty = 1,
        xlab = "wave (nm)", ylab = "Reflectance", col='black', lwd=2, xlim=c(440,850))
dev.off()
##############################################################################################################################
# 6.   Resample simulatons to specfic sensors  ----
############################################################################################################################### uses ToolsRTM::get.spectral.convolution.gaussian()'s bulk (matrix) mode

##### VNIR
center <- wave_[keep.data] #wave.bands[5:190] ## 490-815 nm
fwhm   <- mean(diff(sort(unique(center))))  # approximate FWHM as mean band spacing (no explicit sensor FWHM metadata available)
wl     <- seq(468,842,2)

## Resample both the simulated canopy spectra and the field sensor spectra onto the same wl grid, so they're directly comparable
Spec.simula_resample.vnir <- get.spectral.convolution.gaussian(
  centers = wl, fwhm = fwhm, rfl = sim.canopy[, keep.simula, drop = FALSE], wave = wave[keep.simula])
Spec.sensor.vnir <- get.spectral.convolution.gaussian(
  centers = wl, fwhm = fwhm, rfl = data.spectra.sgolay[, keep.data, drop = FALSE], wave = center)


##############################################################################################################################
# 7.   Comparison Plots  ----
##############################################################################################################################

color.d = brewer.pal(7, "Blues")
axis_x<-expression(bold('wave (nm)'))
axis_y<-expression(bold('Reflectance'))


jpeg(file.path(out_root, "Sims/1-Comparison_PROSAIL_Observations.jpg"), units="in", width=5, height=5, res=300)
par(mfrow=c(1,1),  mar = c(5,5,1.1,1),bg= "white",
    font.main=1, cex.main=1,font.axis=2, cex.axis=1, las=1,
    font.lab=2, cex.lab=1.0)
plot(NA,NA, lwd=0,lty=2,type='l',col='forestgreen',ylim=c(0,0.6),xlim=c(450,850),xlab=axis_x,ylab=axis_y)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
par(new=T)

matplot(wl, t(Spec.simula_resample.vnir), type="l", lwd=1,lty=2,col='black',ylim=c(0,0.6),xlim=c(450,850),xlab=axis_x,ylab=axis_y)
par(new=T)
matplot(wl, t(Spec.sensor.vnir), type="l", lwd=1,lty=2,col='forestgreen',ylim=c(0,0.6),xlim=c(450,850),xlab=axis_x,ylab=axis_y)

legend("topleft", legend = c(expression(bold('PROSAIL')),expression(bold('Obser.'))),
       fill=c('black','forestgreen'),cex=0.8)
dev.off()

##############################################################################################################################
# 8.   Save Simulations and LUTs  ----
##############################################################################################################################

rfl.hypertoExport<-as.data.frame(Spec.simula_resample.vnir)
colnames(rfl.hypertoExport)<-paste('R.',wl,sep='')
rfl.bands = paste('R.',wl,sep='')

LUT_rfl.hyper<-cbind(ID=IDs, LUT,psoil=psoil,rfl.hypertoExport)
head(LUT_rfl.hyper)

LUT_rfl.hyper.withIndices<-getIndices(data=LUT_rfl.hyper,  pattern.rfl = "R.")
n_indices=dim(LUT_rfl.hyper)[2]+1
indices.names<-names(LUT_rfl.hyper.withIndices)[n_indices:dim(LUT_rfl.hyper.withIndices)[2]]


write.table(LUT_rfl.hyper.withIndices,file=paste(out_root,'/LUTs/LUT_PROSAIL_with_n',nSamples/1000,'k.csv',sep=''),sep=',',row.names = F)


variable_inputs <- c('ID','Cab','Car','Anth','CBC','EWT','Prot','N','LIDFa','LAI')
LUT_rfl.hyper.sb <- LUT_rfl.hyper.withIndices[,c(variable_inputs,rfl.bands,indices.names)]
write.table(LUT_rfl.hyper.sb,file=paste(out_root,'/LUTs/LUT_PROSAIL_with_n',nSamples/1000,'k_sb_MainTraits.csv',sep=''),sep=',',row.names = F)


##############################################################################################################################
# 9.   Save Dataset  ----
##############################################################################################################################


rfl.hyper.resampled<-as.data.frame(Spec.sensor.vnir)
colnames(rfl.hyper.resampled)<-paste('R.',wl,sep='')
dataset.resampled<-cbind(ID=c(1:dim(dataset)[1]), dataset[,c(1:5)],rfl.hyper.resampled)
head(dataset.resampled)

dataset.resampled.withIndices<-getIndices(data=dataset.resampled,  pattern.rfl = "R.")

write.table(dataset.resampled.withIndices,file=file.path(out_root, "FieldData/2019_2020RawData_Vcmax_withRFL_resampled.csv"),sep=',',row.names = F)


##############################################################################################################################
# 10.   Reflectance plots and more  ----
##############################################################################################################################

# output folder
paths.outs <- paste0(out_root, "/FieldData/")


##############################################################################################################################
####################### Reflectance plot

df.plot <-dataset.resampled.withIndices
plot.rfl <-cbind(df.plot[1:7], stack(df.plot[10:173]))
colnames(plot.rfl)<-c(names(plot.rfl)[1:7],'RFL','wave')
plot.rfl$wave =as.numeric(sub("R.", "", plot.rfl$wave, fixed = TRUE))

plot.rfl$GroupID <- cut(plot.rfl$Plot, breaks=c(0,100,200,300,400,500,Inf),labels = c('0-100', '100-200', '200-300','300-400','400-500','>500'))
GroupID<-plot.rfl$GroupID
unique(plot.rfl$Group)
spectral_plot<-ggplot(plot.rfl, aes(x=wave, y=RFL, group=GroupID)) +theme_bw() +
  geom_line(aes(color=as.factor(GroupID)),linetype = "dashed", size=0.7) +
  labs(color = GroupID, x='Wavelength',y='Reflectance')+  ggtitle("spectral signatures")+
  geom_point(aes(color=as.factor(GroupID)), size=0.4) +
  theme(plot.title = element_text(hjust = 0.5, size=12), legend.title=element_blank())
print(spectral_plot)
ggsave(paste(paths.outs,'1-Reflectance_spectra.png',sep=''),
       width = 12, height = 8,  dpi = 300,units = "cm")


##############################################################################################################################
####################### Reflectance plot


axis_x<-bquote(bold('measured Vcmax')) # axis x
axis_y<-bquote(bold('NDVI')) # axis x
statsLabel = paste0('')

scatter.vcmax <-  ggplot(dataset.resampled.withIndices, aes_string(x='Vcmax', y='NDVI')) +
  geom_point(alpha=0.4,shape = 16,aes(), size=1.5) +  #geom_smooth(method=lm, formula = 'y ~ x', se=F,lty=2, lwd=1) +
  theme_bw() + ylim(0.5,0.8) + xlim(25,90)  + ggtitle('ANN model')  +

  scale_color_gradient(low = "#0091ff", high = "#f0650e") +
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, size=10,face="bold"),
        axis.title = element_text(face="bold", size=10),
        axis.text.y=element_text(hjust = 0.5, size=10,face="bold"),
        axis.text.x=element_text(hjust = 0.5, size=10,face="bold"),
        legend.title=element_blank()) +
  labs(title = statsLabel, x=axis_x,y=axis_y, size=10,face="bold") +
  stat_smooth(method = "lm",formula = y ~ x,geom = "smooth",col='red',lty=2,se=T)
print(scatter.vcmax)

ggsave(paste(paths.outs,'1-Relationship_Vcmax_vs_NDVI.png',sep=''),
       width = 8, height = 8,  dpi = 300,units = "cm")

}  # end if (file.exists(field_data_path)) -- steps 5-10



