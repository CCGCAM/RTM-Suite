
rm(list= ls())

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#	0. load main Libraries   -----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (!require("RColorBrewer")) { install.packages("RColorBrewer"); require("RColorBrewer") }  ### colors
if (!require("signal")) { install.packages("signal"); require("signal") }  ### interpolations
if (!require("parallel")) { install.packages("parallel"); require("parallel") }  ### Paralell
if (!require("doParallel")) { install.packages("doParallel"); require("doParallel") }  ### Paralell foreach and caret
if (!require("expint")) { install.packages("expint"); require("expint") }
if (!require("ggplot2")) { install.packages("ggplot2"); require("ggplot2") }

# My packages in R
if (!require("ToolsRTM")) { install.packages("ToolsRTM"); require("ToolsRTM") }  ### Paralell foreach and caret
if (!require("dplyr")) { install.packages("dplyr"); require("dplyr") }  ### dataframes

dir.create('../../outs/Sensibility/', showWarnings = FALSE, recursive = TRUE)
pdf('../../outs/Sensibility/Rplots.pdf')  # catches any stray plot()/print() call, keeps it out of Scripts/

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 1. Get spectra from GetLUT ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Define number of simulations
nSamples = 15 ###  length(seq(from = 0.5, to = 7.5, by = 0.5))

inputsPRO = ToolsRTM::inputsPROSAIL
LUT<-as.data.frame(getLUT(inputs = inputsPRO, nLUT=nSamples, setseed = 1234))

plant.traits <- c('LAI','Cab','Cbrown', 'EWT')

for (depVar in plant.traits){
  print(depVar)


  #Inputs by default for Leaf model (PROSPECT)
  if (depVar == 'LAI') {
    LUT$N = 2.5;  LUT$Cab = 40; LUT$Car = 5; LUT$Anth = 0; LUT$Cbrown = 0.5
    LUT$EWT = 0.09; LUT$LMA = 0.012; LUT$alpha = 40
    LUT$Prot = 0; LUT$CBC = 0
    LUT$hspot = 0; LUT$LIDFa = 45;
    LUT$tts = 0; LUT$tto = 0; LUT$psi = 0
    LUT$LAI = seq(from = 0.5, to = 7.5, by = 0.5)
  } else if (depVar == 'Cab') {

    LUT$N = 2.5;  LUT$Cab = seq(5,90,6); LUT$Car = 5; LUT$Anth = 0; LUT$Cbrown = 0.5
    LUT$EWT = 0.09; LUT$LMA = 0.012; LUT$alpha = 40
    LUT$Prot = 0; LUT$CBC = 0
    LUT$hspot = 0; LUT$LIDFa = 45;
    LUT$tts = 0; LUT$tto = 0; LUT$psi = 0
    LUT$LAI = 3
  }  else if (depVar == 'Cbrown') {

    LUT$N = 2.5;  LUT$Cab = 40; LUT$Car = 5; LUT$Anth = 0; LUT$Cbrown =seq(0,1,0.07)
    LUT$EWT = 0.09; LUT$LMA = 0.012; LUT$alpha = 40
    LUT$Prot = 0; LUT$CBC = 0
    LUT$hspot = 0; LUT$LIDFa = 45;
    LUT$tts = 0; LUT$tto = 0; LUT$psi = 0
    LUT$LAI = 3
  } else if (depVar == 'EWT') {

    LUT$N = 2.5;  LUT$Cab = 40; LUT$Car = 5; LUT$Anth = 0; LUT$Cbrown = 0.5
    LUT$EWT = seq(0.001,0.5,0.035); LUT$LMA = 0.012; LUT$alpha = 40
    LUT$Prot = 0; LUT$CBC = 0
    LUT$hspot = 0; LUT$LIDFa = 45;
    LUT$tts = 0; LUT$tto = 0; LUT$psi = 0
    LUT$LAI = 3
  }
  head(LUT)




  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # 2. Get soil reflectance
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  dataSpec <- ToolsRTM::dataSpec_PDB # Spectral library
  Rsoil1  <- dataSpec[,11]  # rsoil1 = dry soil
  Rsoil2 <- dataSpec[,12]  # rsoil2 = wet soil
  psoil	 <-  0.5    # soil factor (psoil=0: wet soil / psoil=1: dry soil)
  rsoil <- c(psoil*Rsoil1+(1-psoil)*Rsoil2)


  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # 3. Get Simulations  ----
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


  ## choose number of processors/cores
  no_cores <- detectCores() - 2
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)

  start_time <- Sys.time()
  sim.rfl <- list()
  sims <- foreach(i=1:nSamples) %dopar% {

    data.foursail_pro <- ToolsRTM::foursail(inputLUT=LUT[i,],rsoil=rsoil,LeafModel = 'PROSPECT-PRO')
    rdot <- data.foursail_pro[[1]]
    rsot <- data.foursail_pro[[2]]
    data.foursail_pro <- ToolsRTM::Compute_BRF(rdot=rdot,rsot=rsot,tts=LUT[i,'tts'],data.light =ToolsRTM::dataSpec_PDB)

    sim.rfl[[i]]<-data.foursail_pro


  } ##end parallel

  stopCluster(cl)
  end_time <- Sys.time()
  print(end_time - start_time)

  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # 4.   Generate Plot object  ----
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  # Resample to SE2

  sim.canopy<-as.matrix(do.call(rbind,sims))
  wave<-c(400:2500)


  waves.se2 <- subset(ToolsRTM::sensor.characteristics, Sensor == 'Sentinel2a')
  center_wvl <- waves.se2$average

  # Initialize an empty list to store interpolated reflectance data
  SE2_df <- list()
  refl.se2 <- list()
  # Perform linear interpolation for each row in sim.canopy
  for (i in 1:nrow(sim.canopy)) {
    # Perform linear interpolation to resample to Sentinel-2 bands
    reflectance.i <- signal::interp1(wave, sim.canopy[i, ], center_wvl, method = "spline")
    # Store the interpolated reflectance in the list
    SE2_df[[i]] <- data.frame(center_wvl = center_wvl, reflectance = reflectance.i, n.sim=i,trait=LUT[,depVar][i])
    refl.se2[[i]] <- reflectance.i
  }

  # Convert the list to a matrix
  df.SE2 <- do.call(rbind, SE2_df)
  df.rfl.SE2 <- do.call(rbind, refl.se2)

  # Calculate average, 25th percentile, and 50th percentile for each band
  stats.SE2 <- df.SE2 %>%
    group_by(center_wvl) %>%
    summarise(
      average = mean(reflectance,na.rm=T),
      median = median(reflectance,na.rm=T),
      percentile_05 = quantile(reflectance, 0.05,na.rm=T),
      percentile_50 = quantile(reflectance, 0.50,na.rm=T),
      percentile_95 = quantile(reflectance, 0.95,na.rm=T)
    )
  # Define the colors for the greens
  gradient_colors <- colorRampPalette(c("cornsilk4",'gold3', "darkolivegreen4"))(15)


  #::::::::::::::::::::::::::::::::::::::::::
  # 4.1   Generate Plot   ----
  #::::::::::::::::::::::::::::::::::::::::::

  # Create a ggplot scatter plot with modifications
  plot.ind <- ggplot(df.SE2, aes(x = center_wvl, y = reflectance,color=as.factor(trait))) + #ylim(0,1) +
    #geom_tile(aes(color = trait), width = 1) +
    geom_line(size = 1,alpha=0.5) +
    scale_color_manual(values = gradient_colors) +
    theme_bw() +
    theme(legend.position="right",
          plot.title = element_text(hjust = 0.5, size=10,face="bold"),
          panel.background = element_rect(fill="grey87"),
          axis.title = element_text(face="bold", size=14),
          axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
          axis.text.x=element_text(hjust = 0.5, size=12,face="bold"),
          legend.title=element_blank())
  #print(plot.ind)


  #::::::::::::::::::::::::::::::::::::::::::
  # 4.2   Generate Plot   ----
  #::::::::::::::::::::::::::::::::::::::::::

  plot_ <- ggplot(stats.SE2, aes(x = center_wvl)) + ylim(0,0.5) +
    geom_line(aes(y = average), size = 0.8) +
    geom_line(aes(y = median),  linetype = "dashed", linewidth = 0.8) +
    geom_ribbon(aes(ymin = percentile_05, ymax = percentile_95),fill = "black", alpha = 0.3) +
    labs(
      title = paste('',sep=''),
      x = "wavelength (nm)",
      y = "Reflectance"
    ) +
    theme_bw() +  theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12,face="bold"),
                        panel.background = element_rect(fill="grey90"),
                        axis.title = element_text(face="bold", size=14),
                        axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
                        axis.text.x=element_text(hjust = 0.5, size=12,face="bold"),
                        legend.title=element_blank())

  #plot(plot_)

  ### for plots (comparison with Field data)
  paths.plots = '../../outs/Sensibility/'
  ifelse(!dir.exists(paths.plots), dir.create(paths.plots, recursive = TRUE), FALSE)

  ggsave(paste(paths.plots,'1-Figure-',depVar,'.png',sep=''),
         plot=plot_,
         width = 12, height = 12,  dpi = 300,units = "cm")

  #::::::::::::::::::::::::::::::::::::::::::
  # 4.3   Generate Plot   ----
  #::::::::::::::::::::::::::::::::::::::::::

  sim.canopy.ggplot <-data.frame(Trait= LUT[,depVar],ID=1:nrow(sim.canopy), sim.canopy)
  colnames(sim.canopy.ggplot)<-c('Trait','ID',paste('R',c(400:2499),sep='.'))
  sim.canopy.ggplot <-cbind(sim.canopy.ggplot[1], stack(sim.canopy.ggplot[3:dim(sim.canopy.ggplot)[2]]))
  sim.canopy.ggplot$ind =as.numeric(sub("R.", "", sim.canopy.ggplot$ind, fixed = TRUE))

  spectral_plot<-ggplot(sim.canopy.ggplot, aes(x=ind, y=values, group=Trait)) +
    geom_line(aes( color=Trait), size=0.8) +
    labs(x='Wavelength',y='Reflectance') +  ggtitle("")+
    viridis::scale_color_viridis(direction = -1) +
    theme_bw() +
    theme(legend.position="right",
          strip.text.x = element_text(size = 12, color = "black", face = "bold"),
          plot.title = element_text(hjust = 0.5, size=10,face="bold"),
          panel.background = element_rect(fill="grey75"),
          axis.title = element_text(face="bold", size=14),
          axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
          legend.title =element_text(hjust = 0.0,size=12,face="bold"),
          axis.text.x=element_text(hjust = 0.5, size=12,face="bold"))
  print(spectral_plot)
  ggsave(paste(paths.plots,'2-Figure-spectral-',depVar,'.png',sep=''),
         plot=spectral_plot,
         width = 12, height = 12,  dpi = 300,units = "cm")

  #::::::::::::::::::::::::::::::::::::::::::
  # 4.4   Generate Plot   ----
  #::::::::::::::::::::::::::::::::::::::::::

  sim.canopy.ggplot.se2 <-data.frame(Trait= LUT[,depVar],ID=1:nrow(df.rfl.SE2), df.rfl.SE2)
  colnames(sim.canopy.ggplot.se2)<-c('Trait','ID',paste('R',center_wvl,sep='.'))
  sim.canopy.ggplot.se2 <-cbind(sim.canopy.ggplot.se2[1], stack(sim.canopy.ggplot.se2[3:dim(sim.canopy.ggplot.se2)[2]]))
  sim.canopy.ggplot.se2$ind =as.numeric(sub("R.", "", sim.canopy.ggplot.se2$ind, fixed = TRUE))


  spectral_plot.se<-ggplot(sim.canopy.ggplot.se2, aes(x=ind, y=values, group=Trait)) +
    geom_line(aes( color=Trait), size=0.8) +
    labs(x='Wavelength',y='Reflectance') +  ggtitle("")+
    viridis::scale_color_viridis(direction = -1) +
    theme_bw() +
    theme(legend.position="right",
          strip.text.x = element_text(size = 12, color = "black", face = "bold"),
          plot.title = element_text(hjust = 0.5, size=10,face="bold"),
          panel.background = element_rect(fill="grey75"),
          axis.title = element_text(face="bold", size=14),
          axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
          legend.title =element_text(hjust = 0.0,size=12,face="bold"),
          axis.text.x=element_text(hjust = 0.5, size=12,face="bold"))
  print(spectral_plot.se)
  ggsave(paste(paths.plots,'2-Figure-spectral-',depVar,'_withSE2.png',sep=''),
         plot=spectral_plot.se,
         width = 12, height = 12,  dpi = 300,units = "cm")


  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # 4.   Get LUT w/ Indices ----
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


  # Convert resampled list to wide format df

  SE2_wide <- tidyr::pivot_wider(df.SE2, names_from = c("center_wvl"), values_from = c("reflectance"))
  SE2_wide <- tidyr::unnest(SE2_wide,names(SE2_wide))

  # Rename bands
  names(SE2_wide) <- c('n.sim','trait','B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12')
  # Calculate indices
  SE2.bands <- c('B1','B2','B3','B4','B5','B6','B7','B8','B8A','B9','B10','B11','B12')
  SE2.ind <- getIndicesSE2(SE2_wide[,SE2.bands], sensor = "Sentinel-2a", df.data = NULL, fast.process =T)

  # Combine LUT, bands and indices
  LUT_SE2 <- cbind(LUT,SE2_wide)
  LUT_SE2_Indices <- cbind(LUT,SE2_wide,SE2.ind)

  sigma <- 0.15
  knr <- exp(-(LUT_SE2_Indices['B8'] - LUT_SE2_Indices['B4'])^2/(2*sigma^2))
  LUT_SE2_Indices['kNDVI'] <- (1-knr) / (1+knr)
  # Calculate EVI using bands B2 and B8A from the dataframe
  LUT_SE2_Indices <- LUT_SE2_Indices %>%
    mutate(EVI = 2.5 * ((B8A - B4) / (B8A + 6.0 * B4 - 7.5 * B2+ 1.0))) %>%
    mutate(NDVIv = ((B8 - B4) / (B8 + B4)) * B8)  %>%
    rename(CR.re.nir = CR.red.nir.6)


  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # 4.   Get Plots ----
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # Define the colors for the greens
  gradient_colors <- colorRampPalette(c("cornsilk4",'gold3', "darkolivegreen4"))(15)

  # Create a ggplot scatter plot with modifications
  plot.ind <- ggplot(LUT_SE2_Indices, aes_string(x = depVar, y = 'kNDVI',color=depVar)) + ylim(0,1) +
    geom_point(size = 5,alpha=0.8) +
    geom_line(size = 1,alpha=0.2, color='darkolivegreen4') +
    viridis::scale_color_viridis() +

    # scale_color_gradientn(colors = gradient_colors)  +
    theme_bw() +
    theme(legend.position="right",
          plot.title = element_text(hjust = 0.5, size=10,face="bold"),
          panel.background = element_rect(fill="grey87"),
          axis.title = element_text(face="bold", size=14),
          axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
          axis.text.x=element_text(hjust = 0.5, size=12,face="bold"),
          legend.title=element_blank())
  #print(plot.ind)


  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # 5.   Get Plots ----
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



  axis_y<- expression(bold("VIs"))

  if (depVar == 'Cab') {
    axis_x<- expression(bold("Cab (" ~ μg ~ cm^-2 ~")"))
    # Select the columns Trait, LAI, NDVI, NDVIv, kNDVI, EVI, CR.red.nir.6
    selected_cols <- c("trait", "LAI",'Cab', 'Cbrown','EWT',"NDRE", "NDVIv", "kNDVI", "EVI", 'CIre',"CR.re.nir")
  } else if (depVar == 'LAI')  {
    axis_x<- expression(bold("LAI (" ~ m^2 ~ m^-2 ~")"))
    # Select the columns Trait, LAI, NDVI, NDVIv, kNDVI, EVI, CR.red.nir.6
    selected_cols <- c("trait", "LAI",'Cab', 'Cbrown','EWT',"NDVI", "NDVIv", "kNDVI", "EVI", 'CIre',"CR.re.nir")
  } else if (depVar == 'Cbrown') {
    axis_x<- expression(bold("Cbrown pigments"))
    # Select the columns Trait, LAI, NDVI, NDVIv, kNDVI, EVI, CR.red.nir.6
    selected_cols <- c("trait", "LAI",'Cab', 'Cbrown','EWT',"ARI", "NDVIv", "kNDVI", "EVI", 'CIre',"CR.re.nir")
  }  else if (depVar == 'EWT')  {
    axis_x<- expression(bold("EWT (" ~ g ~ cm^-2 ~")"))
    # Select the columns Trait, LAI, NDVI, NDVIv, kNDVI, EVI, CR.red.nir.6
    selected_cols <- c("trait", "LAI",'Cab', 'Cbrown','EWT',"NDWI", "WDRVI", "NBR", "CR.SWIR", 'WET',"SBI")
  }


  # Convert the selected columns to long format
  lut.long <- LUT_SE2_Indices %>%
    dplyr::select(all_of(selected_cols)) %>%
    tidyr::pivot_longer(cols = c(-trait, -LAI, -Cab,-Cbrown, -EWT),
                        names_to = "Variable",
                        values_to = "Value")


  head(lut.long)

  # Create a ggplot scatter plot with modifications
  plot.ind <- ggplot(lut.long, aes_string(x = depVar, y = 'Value',color=depVar)) +
    geom_point(size = 4,alpha=0.8) +
    viridis::scale_color_viridis(direction=-1) +
    facet_wrap(~Variable, ncol = 2,  scales = "free") +
    labs(title = "", x = axis_x, y = axis_y, size = 12, face = "bold") + theme_bw()+
    # scale_color_gradientn(colors = gradient_colors)  +
    theme_bw() +
    theme(legend.position="right",
          strip.text.x = element_text(size = 12, color = "black", face = "bold"),
          plot.title = element_text(hjust = 0.5, size=10,face="bold"),
          panel.background = element_rect(fill="grey87"),
          axis.title = element_text(face="bold", size=14),
          axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
          legend.title =element_text(hjust = 0.0, size=12,face="bold"),
          axis.text.x=element_text(hjust = 0.5, size=12,face="bold"))
  print(plot.ind)

  ggsave(paste(paths.plots,'3-Figure-',depVar,'-sensibility.png',sep=''),
         plot=plot.ind,
         width = 18, height = 26,  dpi = 300,units = "cm")
}






dev.off()  # close the Rplots.pdf capture device opened near the top
