
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

##############################################################################################################################
#	2. Bands and variables  -----
##############################################################################################################################

inputs <- c('Cab', 'Car', 'Anth', 'CBC','Prot','Cbrown', 'N', 'LIDFa', 'LAI', 'EWT', 'hspot','psi')
bands.SE2<-c('B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A','B11', 'B12')
# Datos de salida
bands <- LUT[, bands.SE2]

#Datos de entrada
X <- as.data.frame(LUT[,inputs])
# Define the parameter ranges
parameters <- lapply(X, function(x) c(min(x), max(x)))
ranges <- apply(LUT[,inputs], 2, range)
# Now, create the q.arg list using these ranges
q_arg <- apply(ranges, 2, function(x) list(min = x[1], max = x[2]))

summary(LUT[,1:18])
# Define the names of the input factors
input_factors <- names(X)
# This should be a multiple of the number of input variables
num_samples <- dim(X)[1] / 2 #length(input_factors) * 10

##############################################################################################################################
#	3. Sensitivity analysis -----
##############################################################################################################################


df.johnson<- data.frame(Band = character(), Parameter=character(), Index = numeric())
#df.SmthSpl <- data.frame(Band = character(), Parameter=character(), Index = numeric())
df.src <- data.frame(Band = character(), Parameter=character(), Index = numeric())
df.fast <- data.frame(Band = character(), Parameter=character(), Index = numeric())

for (band in colnames(bands)) {
  print(band)
  y <- bands[, band]

  # Perform the  the Johnson indices
  ind.johnson <- sensitivity:: johnson(X, y)
  df.johnson <- rbind(df.johnson, data.frame(Band = band, Parameter =rownames(ind.johnson$johnson), Index = ind.johnson$johnson$original ))

  # Perform the Sobol' First Order Indices with B-spline Smoothing
  #sobol.sm <-sensitivity::sobolSmthSpl(y,X)
  #sobol.sm <- as.data.frame(sobol.sm$S)
  #df.SmthSpl <- rbind(df.SmthSpl,data.frame(Band = band, Parameter = rownames(sobol.sm), Index = sobol.sm$Si))

  src_ <-sensitivity::src(X,y, rank = F, logistic = FALSE, nboot = 0, conf = 0.95)
  df.src <- rbind(df.src,data.frame(Band = band, Parameter = rownames(src_$SRC), Index = src_$SRC$original))

  # Perform extended-FAST" method (Saltelli et al. 1999)
  #fast99_ <- sensitivity::fast99(model = NULL, X = X, y = y,factors = input_factors,
   #                              q = rep('qunif', length(input_factors)), q.arg = q_arg,
    #                             n = num_samples,nboot = 1000)

}



# Adjust variable to Levels

df.johnson$Band <- factor(df.johnson$Band, levels = bands.SE2)
#df.SmthSpl$Band <- factor(df.SmthSpl$Band, levels = bands.SE2)
df.src$Band <- factor(df.src$Band, levels = bands.SE2)


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## :::::: Define colors
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Define custom colors
my_colors <- c(Cab = "seagreen3", Car = "darkseagreen4", Anth = "tomato1",N='red',
               Cbrown = "orange3", EWT = "steelblue", hspot='grey',psi='mediumpurple2',
               LIDFa = "seashell2", LAI = "lightsteelblue3", CBC = "lightgoldenrod2", Prot = "brown4")




#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## :::::: plot df.src
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

ggplot(df.src, aes(x = Band, y = Index, fill = Parameter)) +
  geom_bar(stat = "identity", position = "stack",colour = "black") +
  scale_fill_manual(values = my_colors) +
  labs(title = "",
       x = "",
       y = "Standardized Regression Coefficients",
       fill = "Parameter") +
  theme_bw() + theme(legend.position="right",
                     plot.title = element_text(hjust = 0.5, size=12,face="bold"),
                     panel.background = element_rect(fill="white"),
                     plot.background = element_rect(fill = 'white', color = 'white'),
                     legend.key = element_rect(fill = "white", color = "white"),
                     axis.title = element_text(face="bold", size=12),
                     axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
                     axis.text.x=element_text(angle = 0, vjust = 0.5,hjust = 1, size=12,face="bold"),
                     legend.title=element_blank())
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## :::::: plot df.johnson
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

ggplot(df.johnson, aes(x = Band, y = Index, fill = Parameter)) +
  geom_bar(stat = "identity", position = "stack",colour = "black") +
  scale_fill_manual(values = my_colors) +
  labs(title = "",
       x = "",
       y = "Johnson",
       fill = "Parameter") +
  theme_bw() + theme(legend.position="right",
                     plot.title = element_text(hjust = 0.5, size=12,face="bold"),
                     panel.background = element_rect(fill="white"),
                     plot.background = element_rect(fill = 'white', color = 'white'),
                     legend.key = element_rect(fill = "white", color = "white"),
                     axis.title = element_text(face="bold", size=12),
                     axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
                     axis.text.x=element_text(angle = 0, vjust = 0.5,hjust = 1, size=12,face="bold"),
                     legend.title=element_blank())



#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## :::::: plot df.SmthSpl
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


ggplot(df.SmthSpl, aes(x = Band, y = Index, fill = Parameter)) +
  geom_bar(stat = "identity", position = "stack",colour = "black") +
  scale_fill_manual(values = my_colors) +
  labs(title = "",
       x = "",
       y = " Sobol' First Order Indices with B-spline Smoothing",
       fill = "Parameter") +
  theme_bw() + theme(legend.position="right",
                     plot.title = element_text(hjust = 0.5, size=12,face="bold"),
                     panel.background = element_rect(fill="white"),
                     plot.background = element_rect(fill = 'white', color = 'white'),
                     legend.key = element_rect(fill = "white", color = "white"),
                     axis.title = element_text(face="bold", size=12),
                     axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
                     axis.text.x=element_text(angle = 0, vjust = 0.5,hjust = 1, size=12,face="bold"),
                     legend.title=element_blank())






##############################################################################################################################
#	4. plots by wavelength  -----
##############################################################################################################################
center_wvl <-c(ToolsRTM::Sentinel2A.MSI$center_wvl)[c(2:9,12:13)]
center_wvl<-c(520,560,654,701,743,779,840,865,1625,2260)
# Function to convert Sentinel-2 bands to wavelengths
band_to_wavelength <- function(band) {
  wavelengths <- c('B2' = 520, 'B3' = 560, 'B4' = 654, 'B5' = 701, 'B6' = 743,
                   'B7' = 779, 'B8' = 842, 'B8A' = 865, 'B11' = 1625, 'B12' = 2260)
  return(wavelengths[band])
}

# Add wavelength column
df.SmthSpl <- df.SmthSpl %>%
  mutate(Wavelength = sapply(Band, band_to_wavelength))

# Normalize Index to 100% by band
df.SmthSpl <- df.SmthSpl %>%
  group_by(Band) %>%
  mutate(Normalized_Index = (Index - min(Index)) / (max(Index) - min(Index)) * 100) %>%
  mutate(Normalized_Index.v2 = Index / sum(Index) * 100)  %>%
  ungroup() %>%
  arrange(Band, Parameter) %>%
  group_by(Band) %>%
  mutate(Cumulative_Index = cumsum(Normalized_Index.v2))



df.SmthSpl <- as.data.frame(df.SmthSpl)

# Plot
plot.sobol<- ggplot(df.SmthSpl, aes(x = Wavelength, y = Normalized_Index.v2, fill = Parameter)) +
  geom_area(color = "black",alpha = 0.8) +   # Add area under the line
  scale_fill_manual(values = my_colors) +  # Fill area with specified colors
  labs(title = "",
       x = "Wavelength (nm)",
       y = "Relative Contribution (%)",
       color = "") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        panel.background = element_rect(fill = "grey87"),
        plot.background = element_rect(fill = 'white', color = 'white'),
        legend.key = element_rect(fill = "white", color = "white"),
        axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle=90, vjust = 0.5),
        axis.text = element_text(size = 10,face = "bold"),
        legend.title = element_text(face = "bold", size = 12)) +
  scale_x_continuous(breaks = center_wvl,label=bands.SE2) +  # Adjust x-axis breaks
  scale_y_continuous(breaks = seq(0, 100, by = 10))  # Adjust y-axis breaks

print(plot.sobol)
######################################################################################################################
#	5. plots by wavelength  -----
##############################################################################################################################


# Add wavelength column
df.johnson <- df.johnson %>%
  mutate(Wavelength = sapply(Band, band_to_wavelength))

# Normalize Index to 100% by band
df.johnson <- df.johnson %>%
  group_by(Band) %>%
  mutate(Normalized_Index = (Index - min(Index)) / (max(Index) - min(Index)) * 100) %>%
  mutate(Normalized_Index.v2 = Index / sum(Index) * 100)  %>%
  ungroup() %>%
  arrange(Band, Parameter) %>%
  group_by(Band) %>%
  mutate(Cumulative_Index = cumsum(Normalized_Index.v2))



df.johnson <- as.data.frame(df.SmthSpl)

# Plot
plot.johnson <-ggplot(df.johnson, aes(x = Wavelength, y = Normalized_Index.v2, fill = Parameter)) +
  geom_area(color = "black",alpha = 0.8) +   # Add area under the line
  scale_fill_manual(values = my_colors) +  # Fill area with specified colors
  labs(title = "",
       x = "Wavelength (nm)",
       y = "Relative Contribution (%)",
       color = "") +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
        panel.background = element_rect(fill = "grey87"),
        plot.background = element_rect(fill = 'white', color = 'white'),
        legend.key = element_rect(fill = "white", color = "white"),
        axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle=90, vjust = 0.5),
        axis.text = element_text(size = 10,face = "bold"),
        legend.title = element_text(face = "bold", size = 12)) +

  scale_x_continuous(breaks = center_wvl,label=bands.SE2) +  # Adjust x-axis breaks
  scale_y_continuous(breaks = seq(0, 100, by = 10))  # Adjust y-axis breaks

print(plot.johnson)

