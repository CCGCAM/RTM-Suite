
rm(list=ls())

##############################################################################################################################
#	0. load main Libraries   -----
##############################################################################################################################

if (!require("sensitivity")) { install.packages("sensitivity"); require("sensitivity") }  ### for corrplots
if (!require("sensobol")) { install.packages("sensobol"); require("sensobol") }  ### for corrplots

# My packages in R
if (!require("ToolsRTM")) { install.packages("ToolsRTM"); require("ToolsRTM") }  ### Paralell foreach and caret
if (!require("RColorBrewer")) { install.packages("RColorBrewer"); require("RColorBrewer") }  ### colors
if (!require("ggplot2")) { install.packages("ggplot2"); require("ggplot2") }  ### ggplot2
if (!require("dplyr")) { install.packages("dplyr"); require("dplyr") }  ###


if (!require("ToolsRTM")) { install.packages("ToolsRTM"); require("ToolsRTM") }  ### Paralell foreach and caret
if (!require("dplyr")) { install.packages("dplyr"); require("dplyr") }  ### dataframes
if (!require("parallel")) { install.packages("parallel"); require("parallel") }  ### Paralell
if (!require("doParallel")) { install.packages("doParallel"); require("doParallel") }  ### Paralell foreach and caret
if (!require("signal")) { install.packages("signal"); require("signal") }  ### interpolations

source('CodesR/functions/get.sobol.indices.R')

version <-'PROSAIL'

paths.models=paste('Tables/Models/',version,'/',sep='')
ifelse(!dir.exists(paths.models), dir.create(paths.models), FALSE)

##############################################################################################################################
#	1. Open  LUT with correlations  -----
##############################################################################################################################


LUT.1<- data.table:: fread(paste('Tables/LUTs/',version,'/LUT5k_PROSAIL_SE2_Autumn_wIndices.csv',sep='')) ##
LUT.2<- data.table:: fread(paste('Tables/LUTs/',version,'/LUT5k_PROSAIL_SE2_Winter_wIndices.csv',sep='')) ##
LUT.3<- data.table:: fread(paste('Tables/LUTs/',version,'/LUT5k_PROSAIL_SE2_Summer_wIndices.csv',sep='')) ##
LUT.4<- data.table:: fread(paste('Tables/LUTs/',version,'/LUT5k_PROSAIL_SE2_Spring_wIndices.csv',sep='')) ##
LUT <- as.data.frame(rbind(LUT.4,LUT.3,LUT.2,LUT.1))
dim(LUT)
names(LUT)[1:22]

LUT['BF.Anth'] <- (LUT['B3'] - LUT['B2']) / (LUT['B3'] + LUT['B2'])
LUT['CR.red.nir']  <- LUT['B7']  / (LUT['B8'] + ( (LUT['B6'] - LUT['B8']) / (740.5- 833) ) * (782.8 - 832.8))



##############################################################################################################################
#	2. Get Sobol indices  -----
##############################################################################################################################
if (version == 'PROSAIL') {
  traits <- c('Cab', 'Car', 'Anth', 'N', 'Cbrown', 'CBC', 'Prot', 'EWT', 'LIDFa', 'LAI', 'hspot','tto', 'tts')
} else {
  traits <- c('Cab', 'Car', 'Anth', 'N', 'Cbrown', 'CBC', 'Prot', 'EWT', 'LIDFa', 'LAIu', 'sd', 'cd', 'hspot', 'tto', 'tts')

}

spectral.indices<-c('ARI','Greeness','TCARI','TCARI_OSAVI','GNDVI',
                    'CCCI','REP','Datt6','CIre','CR.red.nir','NDRE','CR.red.nir.6',
                    'NDVI','MTVI2','Redness','MTVI1',
                    'CR.SWIR','CR.SWIR.2','NDWI','WET','NBR.2')
bands.SE2<-c('B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A','B11', 'B12')


# Find the starting index of "NDVI" columns
start_index <- grep("^NDVI", names(LUT))
indices.names <- start_index:length(names(LUT))
spectra.all <- names(LUT)[indices.names]
#print(spectra.all)
# Get the lower triangle


cor.amtrix <- round(cor(LUT[,c(traits,spectral.indices)]), 1)
plot.cor <- ggcorrplot::ggcorrplot(cor.amtrix,
           hc.order = TRUE,
           type = "lower",
           outline.color = "white") +theme_bw()
# Seleect the inputs for make the Analysis of Senisibility
indices<-spectral.indices


num_samples <- dim(LUT)[1] / 2 #length(input_factors) * 10

# Iterate over each spectral traits and compute Sobol indices
# Initialize a list to store Sobol indices for each spectral traits
df.sobol<-list()
# Iterate over each spectral i.inputs and compute Sobol indices
for (i.index in indices) {
  df.sobol[[i.index]] <- get.sobol.indices(data = LUT[, c(traits, i.index)],
                                           output = i.index, N = num_samples, normalize = T)
}

# Convert the list to a dataframe
df.sobol <- do.call(rbind, df.sobol)
df.sobol$Band <- factor(df.sobol$Band, levels = indices)

# Convert negative values to positive in the dataframe
df.sobol[, c("Si", "STi")] <- abs(df.sobol[, c("Si", "STi")])
rownames(df.sobol) <-NULL
head(df.sobol)




#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## :::::: Define colors
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Define colors for the parameters
my_colors <- rainbow(length(unique(df.sobol$Parameter)))

# Define custom colors
my_colors <- c(Cab = "seagreen4", Car = "#A1D99B", Anth = "red",N='magenta',
               Cbrown = "orange4",  CBC = "lightgoldenrod2", Prot = "gold1",
               LIDFa = "seashell2", LAI = "seashell4",
               EWT = "steelblue",
               hspot='grey11',tto='lightcyan2',tts='cyan4',
               LAIu = "seashell4", cd='darkslategray2', sd='darksalmon')

df.sobol$Parameter <- factor(df.sobol$Parameter, levels = traits)


##############################################################################################################################
#	3. Get plots  -----
##############################################################################################################################

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::: 3.1. Barplot df.sobol (100)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

plot.barplot<-ggplot(df.sobol, aes(x = Band, y = Si_norm, fill = Parameter)) +
  geom_bar(stat = "identity", position = "stack",colour = "black") +
  scale_fill_manual(values = my_colors) +
  labs(title = "",
       x = "",
       y = "Si Sobol",
       fill = "Parameter") +
  theme_bw() + theme(legend.position="right",
                     plot.title = element_text(hjust = 0.5, size=12,face="bold"),
                     panel.background = element_rect(fill="white"),
                     plot.background = element_rect(fill = 'white', color = 'white'),
                     legend.key = element_rect(fill = "white", color = "white"),
                     axis.title = element_text(face="bold", size=12),
                     axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
                     axis.text.x=element_text(angle = 90, vjust = 0.5,hjust = 1, size=12,face="bold"),
                     legend.title=element_blank())

##print(plot.barplot)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::: 3.2. Rose-Barplot Si_norm (100)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# Create the rose plot
plot.sobol<-ggplot(df.sobol, aes(x = Band, y = Si_norm, fill = Parameter)) +
  geom_bar(stat = "identity", position = "stack",colour = "black") +
  coord_polar() +  # Convert to polar coordinates
  scale_fill_manual(values = my_colors) +
 # scale_fill_brewer(palette="Greens")+xlab("")+ylab("") +
  labs(title = "Si Sobol",
       x = "",
       y = "Si Sobol",
       fill = "Parameter") +
  theme_bw() +  xlab("")+ylab("") +
  theme(legend.position="right",
        plot.title = element_text(hjust = 0.5, size=12,face="bold"),
        panel.background = element_rect(fill="white"),
        plot.background = element_rect(fill = 'white', color = 'white'),
        legend.key = element_rect(fill = "white", color = "white"),
        axis.title = element_text(face="bold", size=12),
        axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
        axis.text.x=element_text(angle = 0, size=10,color = "black",
                                 face="bold"),
        legend.title=element_blank()) +  # Replace with your desired labels
  scale_x_discrete(labels = c('ARI','Greeness','TCARI','       TCARI-OSAVI','GNDVI',
                              'CCCI','REP','Datt6','CIre','CR.red.nir(8)','NDRE','CR.red.nir(8A)',
                              'NDVI','MTVI2','Redness','MTVI1',
                              'CR-SWIR','CR-SWIR (8)','NDWI','WET','NBR-2'))

print(plot.sobol)


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#:::::: 3.2. Rose-Barplot I.Johnson_norm (100)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


# Create the rose plot
plot.Johnson <-ggplot(df.sobol, aes(x = Band, y =  I.Johnson_norm, fill = Parameter)) +
  geom_bar(stat = "identity", position = "stack",colour = "black") +
  coord_polar() +  # Convert to polar coordinates
  scale_fill_manual(values = my_colors) +
  # scale_fill_brewer(palette="Greens")+xlab("")+ylab("") +
  labs(title = "Johnson",
       x = "",
       y = "",
       fill = "Parameter") +
  theme_bw() +  xlab("")+ylab("") +
  theme(legend.position="right",
        plot.title = element_text(hjust = 0.5, size=12,face="bold"),
        panel.background = element_rect(fill="white"),
        plot.background = element_rect(fill = 'white', color = 'white'),
        legend.key = element_rect(fill = "white", color = "white"),
        axis.title = element_text(face="bold", size=12),
        axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
        axis.text.x=element_text(angle = 0,  color = "black", size=10,face="bold"),
        legend.title=element_blank()) +
  # Replace with your desired labels
  scale_x_discrete(labels = c('ARI','Greeness','TCARI','TCARI-OSAVI','GNDVI',
                              'CCCI','REP','Datt6','CIre','CR.red.nir(8)','NDRE','CR.red.nir(8A)',
                              'NDVI','MTVI2','Redness','MTVI1',
                              'CR-SWIR','CR-SWIR (8)','NDWI','WET','NBR-2'))



print(plot.Johnson)

## for plots (comparison with Field data)
paths.plots = 'Plots/Sensibility/'
ifelse(!dir.exists(paths.plots), dir.create(paths.plots), FALSE)

ggsave(paste(paths.plots,'1-Johnson-indices-',version,'.png',sep=''),
       plot=plot.Johnson,
       width = 20, height = 20,  dpi = 300,units = "cm")

