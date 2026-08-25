
rm(list= ls())

##############################################################################################################################
#	0. load main Libraries   -----    
##############################################################################################################################

if (!require("RColorBrewer")) { install.packages("RColorBrewer"); require("RColorBrewer") }  ### colors
if (!require("ggplot2")) { install.packages("ggplot2"); require("ggplot2") }  ### Paralell foreach and caret

if (!require("signal")) { install.packages("signal"); require("signal") }  ### interpolations
if (!require("parallel")) { install.packages("parallel"); require("parallel") }  ### Paralell
if (!require("doParallel")) { install.packages("doParallel"); require("doParallel") }  ### Paralell foreach and caret
# My packages in R
if (!require("ToolsRTM")) { install.packages("ToolsRTM"); require("ToolsRTM") }  ### Paralell foreach and caret
if (!require("SCOPEinR")) { install.packages("SCOPEinR"); require("SCOPEinR") }  ### Paralell foreach and caret

if (!require("expint")) { install.packages("expint"); require("expint") }  ### Needed for  INFORM model

if (!require("parallel")) { install.packages("parallel"); require("parallel") }  ### Paralell
if (!require("doParallel")) { install.packages("doParallel"); require("doParallel") }  ### Paralell foreach and caret
if (!require("doFuture")) { install.packages("doFuture"); require("doFuture") }  ### Paralell foreach and caret

##############################################################################################################################
# 1.1. Get Soil reflectance from PROSAIL model -----    
##############################################################################################################################
nSamples<-200
ID=(1:nSamples)    


data <- ToolsRTM::dataSpec_PDB
Rsoil1  <- data[,11]  # rsoil1 = dry soil
Rsoil2 <- data[,12]  # rsoil2 = wet soil 
j=2
set.seed(j*1256)
#psoil	 <-  1    # soil factor (psoil=0: wet soil / psoil=1: dry soil)
psoil	 <-  runif(nSamples, 0, 1) 
rsoil0<-list()
for (k in c(1:nSamples)){
  
  rsoil<- c(psoil[k]*Rsoil1+(1-psoil[k])*Rsoil2)
  rsoil0[[k]]<-rsoil#
}
set.seed(Sys.time())

##############################################################################################################################
# 1.2. Simulated SPART  -----    
##############################################################################################################################
inputs.SPART = ToolsRTM::inputsSPART

LUT<-as.data.frame(ToolsRTM::getLUT(inputs = inputs.SPART, nLUT=1, setseed = 1234))




LUT$Cab=40
LUT$Car= 10
LUT$Anth= 10

LUT$N = 1.5
LUT$EWT =0.02
LUT$LMA =0.01
LUT$TypeLidf=1
LUT$LIDFa= -0.35
LUT$LIDFb=-0.15
LUT$LAI=3
LUT$hspot= 0.05

LUT$tto=0
LUT$tts = 40
LUT$psi= 0

LUT$aot550=0.78#0.3246
LUT$uo3 = 0.348
LUT$uh2o = 1.4116
LUT$alt_m =0 
LUT$Pa0 = 1.00133e-03


sim.spart <-  SPART(inputLUT=LUT[1,], optipar=SCOPEinR::optipar2021.Pro.CX, 
                    CanopyModel = 'fourSAIL',
                    LeafModel='PROSPECT-D',
                    df.irradiance = NULL,
                    sensor.i = ToolsRTM::TerraAqua.MODIS,get.plots=F)

#head(sim.spart)

plot.rfl <- ggplot(data = sim.spart$output, aes(x = wave)) +
  labs(y = "reflectance", x = "") +
  geom_point(aes(y = rfl.toa, color = "TOA rfl."), size = 2) +
  geom_line(aes(y = rfl.toa, color = "TOA rfl.")) +

  geom_point(aes(y = rfl.toc, color = "TOC rfl. (SMAC)"), size = 2) +
  geom_line(aes(y = rfl.toc, color = "TOC rfl. (SMAC)")) +

#  geom_point(aes(y = rfl.toc.BRDF, color = "TOC rfl. (BRDF)"), size = 1) +
 # geom_line(aes(y = rfl.toc.BRDF, color = "TOC rfl. (BRDF)")) +
  theme_bw() +
  guides(color = guide_legend(title = "Reflectance:"), linetype = guide_legend(title = "Reflectance:"), shape = guide_legend(title = "Reflectance:")) +
  theme(legend.position="top",
        strip.text.x = element_text(size = 12, color = "black", face = "bold"),
        plot.title = element_text(hjust = 0.5, size=12,face="bold"),
        panel.background = element_rect(fill="grey90"),
        panel.grid.minor = element_line(colour = "grey90"),
        panel.grid.major = element_line(colour = "grey90"),
        axis.title = element_text(face="bold", size=14),
        legend.text=element_text(size=12,face="bold"),
        panel.spacing.x = unit(4, "mm"),
        axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
        axis.text.x=element_text(angle=0,hjust = 0.5, size=12,face="bold"),
        legend.title = element_blank(),
        legend.key = element_rect(fill = "transparent", color = "transparent"),
        legend.background = element_rect(fill = "transparent", color = "transparent"))

print(plot.rfl)

##############################################################################################################################
# 2. Run RTM   ----   
##############################################################################################################################
nSamples = 100
inputs.SPART = ToolsRTM::inputsSPART

LUT<-as.data.frame(ToolsRTM::getLUT(inputs = inputs.SPART, nLUT=nSamples, setseed = 1234))

## choose number of processors/cores
no_cores <- detectCores() - 2 
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# Explicitly export the LUT variable
clusterExport(cl, "LUT")
start_time <- Sys.time()
sim.rfl<-list()
sims<-foreach(i=1:nSamples) %dopar% {
  
  data.spart<- ToolsRTM::SPART.simN(inputLUT = LUT[i,],optipar=SCOPEinR::optipar2021.Pro.CX, 
     CanopyModel = 'fourSAIL',
     LeafModel='PROSPECT-D',
     df.irradiance = NULL,
     sensor.i = ToolsRTM::Sentinel2A.MSI)
  rfl.toa <- data.spart$output$rfl.toa
  rfl.toc <- data.spart$output$rfl.toc
  rfl.toc.brdf <- data.spart$output$rfl.toc.BRDF
  sim.rfl[[i]]<-rfl.toc
  

}
##end paralle
stopCluster(cl)
end_time <- Sys.time()
print(end_time - start_time)

sim.canopy<-do.call(rbind,sims)
wave<-c(ToolsRTM::Sentinel2A.MSI$center_wvl)
wave



# Convert matrix to data frame
df <- data.frame(sim.canopy)
df$row <- 1:nrow(df)  # Add a row identifier
head(df)

# Reshape the data to long format
df_long <- tidyr::gather(df, key = "band", value = "value", -row)
# Make 'band' an ordered factor with desired order
df_long$band <- factor(df_long$band, levels = paste0("X", 1:ncol(df)))


# Plot using ggplot2
ggplot(df_long, aes(x = band, y = value, group = row, color = as.factor(row))) +
  geom_line() +
  geom_point() +
  scale_x_discrete(labels = wave, breaks = unique(df_long$band)) +
  guides(color = guide_legend(title = "Reflectance:"), linetype = guide_legend(title = "Reflectance:"), shape = guide_legend(title = "Reflectance:")) +
  
  labs(title = "SPART Simulations", x = "Bands", y = "reflectance") + 
  theme_bw() +
  theme(legend.position = "none")

##############################################################################################################################
# 1.3. Simulated liberty  -----    
##############################################################################################################################
i=1
inputs.liberty = ToolsRTM::inputsLiberty
LUT<-as.data.frame(ToolsRTM::getLUT_liberty(inputs = inputs.liberty, nLUT=nSamples, setseed = 1234))


sim.liberty<-ToolsRTM::liberty(inputLUT = LUT[i,])


sim.prosail.liberty<-ToolsRTM::foursail(inputLUT=LUT[i,],rsoil=rsoil0[[1]],LeafModel = 'Liberty')
rdot<-sim.prosail.liberty[[1]]
rsot<-sim.prosail.liberty[[2]]
rfl.prosail.liberty<-ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT[2,'tts'],data.light=ToolsRTM::dataSpec_PDB)
plot(rfl.prosail.liberty)

df.rfl <- data.frame(wave =c(400:2500), rfl.toc = rfl.prosail.liberty,rdot = rdot)

plot.Liberty <- ggplot(data = df.rfl, aes(x = wave)) +
  labs(y= " reflectance", x = "") +
  geom_line(aes(y = rfl.toc, color = "rfl toc")) +
  geom_line(aes(y = rdot, color = "rdot")) + theme_bw() +
  guides(color = guide_legend(title = "reflectance"), linetype = guide_legend(title = "reflectance"), shape = guide_legend(title = "reflectance")) 

print(plot.Liberty)
##############################################################################################################################
# 1.4. Simulated PROSAIL-PRO -----    
##############################################################################################################################
inputs.sail = ToolsRTM::inputsPROSAIL
LUT<-as.data.frame(getLUT(inputs = inputs.sail, nLUT=nSamples, setseed = 1234))

sim.prosail.pro<-foursail(inputLUT=LUT[i,],rsoil=rsoil0[[1]],LeafModel = 'PROSPECT-PRO')

rdot<-sim.prosail.pro[[1]]
rsot<-sim.prosail.pro[[2]]
rfl.prosail.pro<-ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT[2,'tts'],data.light=ToolsRTM::dataSpec_PDB)


df.rfl <- data.frame(wave =c(400:2500), rfl.toc = rfl.prosail.pro,rdot = rdot)

plot.prosail <- ggplot(data = df.rfl, aes(x = wave)) +
  labs(y= " reflectance", x = "") +
  geom_line(aes(y = rfl.toc, color = "rfl toc")) +
  geom_line(aes(y = rdot, color = "rdot")) + theme_bw() +
  guides(color = guide_legend(title = "reflectance"), linetype = guide_legend(title = "reflectance"), shape = guide_legend(title = "reflectance")) 

print(plot.prosail)

##############################################################################################################################
# 1.5. Simulated PROSAIL-D  -----    
##############################################################################################################################

sim.prosail.d<-foursail(inputLUT=LUT[i,],rsoil=rsoil0[[1]],LeafModel = 'PROSPECT-D')
rdot<-sim.prosail.d[[1]]
rsot<-sim.prosail.d[[2]]
rfl.prosail.d<-ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT[2,'tts'],data.light=ToolsRTM::dataSpec_PDB)

df.rfl <- data.frame(wave =c(400:2500), rfl.toc = rfl.prosail.d,rdot = rdot)

plot.prosail <- ggplot(data = df.rfl, aes(x = wave)) +
  labs(y= " reflectance", x = "") +
  geom_line(aes(y = rfl.toc, color = "rfl toc")) +
  geom_line(aes(y = rdot, color = "rdot")) + theme_bw() +
  guides(color = guide_legend(title = "reflectance"), linetype = guide_legend(title = "reflectance"), shape = guide_legend(title = "reflectance")) 

print(plot.prosail)

##############################################################################################################################
# 1.3. Simulated SAILH-Fluspect  -----    
##############################################################################################################################
inputs.fluspect = ToolsRTM::inputsFlUSPECT
LUT<-as.data.frame(getLUT(inputs = inputs.fluspect, nLUT=nSamples, setseed = 1234))

sim.prosail.fluspect<-foursail(inputLUT=LUT[i,],rsoil=rsoil0[[1]],LeafModel = 'Fluspect-B')
rdot<-sim.prosail.fluspect[[1]]
rsot<-sim.prosail.fluspect[[2]]
rfl.prosail.fluspect<-ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT[2,'tts'],data.light=ToolsRTM::dataSpec_PDB,short.waves = T)

df.rfl <- data.frame(wave =c(400:2400), rfl.toc = rfl.prosail.fluspect,rdot = rdot)

plot.sail.fluspect <- ggplot(data = df.rfl, aes(x = wave)) +
  labs(y= " reflectance", x = "") +
  geom_line(aes(y = rfl.toc, color = "rfl toc")) +
  geom_line(aes(y = rdot, color = "rdot")) + theme_bw() +
  guides(color = guide_legend(title = "reflectance"), linetype = guide_legend(title = "reflectance"), shape = guide_legend(title = "reflectance")) 

print(plot.sail.fluspect)


##############################################################################################################################
# 1.5. Simulated SAILH-Fluspect-B-Cx  -----    
##############################################################################################################################
inputs.fluspect = ToolsRTM::inputsRTMs
LUT<-as.data.frame(getLUT(inputs = inputs.fluspect, nLUT=nSamples, setseed = 1234))

sim.prosail.fluspect<-foursail(inputLUT=LUT[i,],rsoil=rsoil0[[1]],LeafModel = 'Fluspect-B-Cx')
rdot<-sim.prosail.fluspect[[1]]
rsot<-sim.prosail.fluspect[[2]]
rfl.prosail.fluspect<-ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT[2,'tts'],data.light=ToolsRTM::dataSpec_PDB,short.waves = T)

df.rfl <- data.frame(wave =c(400:2400), rfl.toc = rfl.prosail.fluspect,rdot = rdot)

plot.sail.fluspect <- ggplot(data = df.rfl, aes(x = wave)) +
  labs(y= " reflectance", x = "") +
  geom_line(aes(y = rfl.toc, color = "rfl toc")) +
  geom_line(aes(y = rdot, color = "rdot")) + theme_bw() +
  guides(color = guide_legend(title = "reflectance"), linetype = guide_legend(title = "reflectance"), shape = guide_legend(title = "reflectance")) 

print(plot.sail.fluspect)



##############################################################################################################################
# 1.2. Simulated different models  -----    
##############################################################################################################################

inputs = ToolsRTM::inputsRTMs
LUT<-as.data.frame(getLUT(inputs = inputs, nLUT=nSamples, setseed = 1234))


sim.sail.fluspect.Cx<-foursail(inputLUT=LUT[i,],rsoil=rsoil0[[1]],LeafModel = 'Fluspect-B-Cx')
rdot<-sim.sail.fluspect.Cx[[1]]
rsot<-sim.sail.fluspect.Cx[[2]]
rfl.sail.fluspect.Cx<-ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT[2,'tts'],data.light=ToolsRTM::dataSpec_PDB,short.waves = T)


sim.sail.fluspect.B<-foursail(inputLUT=LUT[i,],rsoil=rsoil0[[1]],LeafModel = 'Fluspect-B')
rdot<-sim.sail.fluspect.B[[1]]
rsot<-sim.sail.fluspect.B[[2]]
rfl.sail.fluspect.B<-ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT[2,'tts'],data.light=ToolsRTM::dataSpec_PDB,short.waves = T)


sim.prosail.pro<-foursail(inputLUT=LUT[i,],rsoil=rsoil0[[1]],LeafModel = 'PROSPECT-PRO')
rdot<-sim.prosail.pro[[1]]
rsot<-sim.prosail.pro[[2]]
rfl.prosail.pro<-ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT[2,'tts'],data.light=ToolsRTM::dataSpec_PDB,short.waves = T)


sim.prosail.d<-foursail(inputLUT=LUT[i,],rsoil=rsoil0[[1]],LeafModel = 'PROSPECT-D')
rdot<-sim.prosail.d[[1]]
rsot<-sim.prosail.d[[2]]
rfl.prosail.d<-ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT[2,'tts'],data.light=ToolsRTM::dataSpec_PDB,short.waves = F)

sim.sail.liberty<-foursail(inputLUT=LUT[i,],rsoil=rsoil0[[1]],LeafModel = 'Liberty')
rdot<-sim.sail.liberty[[1]]
rsot<-sim.sail.liberty[[2]]
rfl.sail.liberty<-ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT[2,'tts'],data.light=ToolsRTM::dataSpec_PDB,short.waves = F)


sim.inform.pro<-inform(inputLUT = LUT[i,],rsoil=rsoil0[[1]],LeafModel = 'PROSPECT-PRO')
sim.inform.d<-inform(inputLUT = LUT[i,],rsoil=rsoil0[[1]],LeafModel = 'PROSPECT-D')
sim.inform.liberty<-inform(inputLUT = LUT[i,],rsoil=rsoil0[[1]],LeafModel = 'Liberty')
sim.inform.fluspect<-inform(inputLUT = LUT[i,],rsoil=rsoil0[[1]],LeafModel = 'Fluspect-B-Cx')


sims<-as.data.frame(cbind(c(400:2400),
                          rfl.sail.liberty[1:2001],rfl.prosail.pro[1:2001],
                          rfl.prosail.d[1:2001],rfl.sail.fluspect.Cx,
                          
                          sim.inform.pro[1:2001],sim.inform.d[1:2001],
                          sim.inform.liberty[1:2001],sim.inform.fluspect))
colnames(sims)<-c('wave','rfl.prosail.liberty','rfl.prosail.pro','rfl.prosail.d','rfl.prosail.f',
                  'rfl.inform.pro','rfl.inform.d','rfl.inform.liberty','rfl.inform.f')
head(sims)

##############################################################################################################################
# 1.3.  Plots different models  -----    
##############################################################################################################################

# Create plot
ggplot(sims, aes(x = wave)) + 
  
  geom_line(aes(y = rfl.inform.liberty, color = "INFORM + Liberty")) +
  geom_line(aes(y = rfl.inform.f, color = "INFORM + FLUSPECT-Cx-B")) +
  geom_line(aes(y = rfl.inform.d, color = "INFORM + PROSPECT-D")) +
  geom_line(aes(y = rfl.inform.pro, color = "INFORM + PROSPECT-PRO")) +
  
  geom_line(aes(y = rfl.prosail.liberty, color = "SAILH + Liberty")) +
  geom_line(aes(y = rfl.prosail.f, color = "SAILH + FLUSPECT-Cx-B")) +
  geom_line(aes(y = rfl.prosail.d, color = "SAILH + PROSPECT-D")) +
  geom_line(aes(y = rfl.prosail.pro, color = "SAILH + PROSPECT-PRO ")) +

  # Set axis labels and title
  labs(x = "Wavelength (nm)", y = "Reflectance", 
       color = "RTM models",
       title = "")  + theme_bw() + ylim(0,0.7) +
  theme(legend.position = c(0.75, 0.8),
        legend.box.background = element_rect(color = "black",linewidth=1),
        plot.title = element_text(hjust = 0.5, size=14,face="bold"),
        axis.title = element_text(face="bold", size=14),
        legend.text = element_text(face="bold", size=10),
        axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
        axis.text.x=element_text(hjust = 0.5, size=12,face="bold"),
       legend.title=element_blank()) +
       # legend.title = element_text(face = "bold", size = 14)) +
guides(color = guide_legend(title.position = "top", title.hjust = 0.5))

