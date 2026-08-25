
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
LUT <- read.csv("Tables/Tables_Erika/LUTVM-20k-SAIL-PROSPECT-PRO_withIndices.csv", sep=',', header = TRUE, nrows=2000)
#LUT <- read.csv("Tables/Tables_Erika/LUTVS-20k-SAIL-PROSPECT-PRO_withIndices.csv", sep=',', header = TRUE, nrows=2000)

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


df.johnson <- as.data.frame(df.johnson)

# Plot
plot.johnson <-ggplot(df.johnson, aes(x = Wavelength, y = Normalized_Index.v2, fill = Parameter)) +
  geom_area(color = "black",alpha = 0.8) +   # Add area under the line
  scale_fill_manual(values = my_colors) +  # Fill area with specified colors
  labs(title = "",
       x = "Longitud de onda (nm)",
       y = "Contribución relativa (%)",
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

ggsave("Plots/Sensibility/plot.johnson_VM.svg", plot.johnson, width = 10, height = 6, units = "in", dpi = 300)

##############################################################################################################################
#	5. Rose-Barplot I.Johnson_norm (100)
##############################################################################################################################

# Create the rose plot
plot.Johnson <-ggplot(df.johnson, aes(x = Band, y =  Normalized_Index.v2, fill = Parameter)) +
  geom_bar(stat = "identity", position = "stack",colour = "black") +
  coord_polar() +  # Convert to polar coordinates
  scale_fill_manual(values = my_colors) +
  # scale_fill_brewer(palette="Greens")+xlab("")+ylab("") +
  labs(title = "",
       x = "",
       y = "",
       fill = "Parámetro") +
  theme_bw() +  xlab("")+ylab("") +
  theme(legend.position="right",
        plot.title = element_text(hjust = 0.5, size=12,face="bold"),
        panel.background = element_rect(fill="white"),
        plot.background = element_rect(fill = 'white', color = 'white'),
        legend.key = element_rect(fill = "white", color = "white"),
        axis.title = element_text(face="bold", size=12),
        axis.text.y=element_text(hjust = 0.5, size=12,face="bold"),
        axis.text.x=element_text(angle = 0, vjust = -1,hjust = -1, size=10,
                                 color = "black",face="bold"),
        legend.title=element_blank())


print(plot.Johnson)
ggsave("Plots/Sensibility/plot.johnson2_VM.svg", plot.Johnson, width = 10, height = 6, units = "in", dpi = 300)

##############################################################################################################################
#	6. Plots  1 nm -----
##############################################################################################################################

# Agrupar los datos por el parámetro
grouped_data <- df.johnson %>%
  group_by(Parameter)

# Secuencia de longitudes de onda
longitudes_onda_interp <- seq(520, 2256, by = 1)

# Crear una lista para almacenar los resultados
interpolated_data_list <- list()

# Iterar sobre cada parámetro único
for (parametro in unique(df.johnson$Parameter)) {
  # Filtrar los datos para el parámetro actual
  current_data <- subset(df.johnson, Parameter == parametro)
  
  # Interpolar los valores de porcentaje utilizando spline
  indice_interp <- signal::interp1(
    x = current_data$Wavelength,
    y = current_data$Index,
    xi = longitudes_onda_interp,
    method = "spline"
  )
  
  # Crear un tibble con los resultados interpolados utilizando spline
  interpolated_data <- tibble(
    Long_onda = longitudes_onda_interp,
    Indice = indice_interp
  )
  
  # Almacenar los resultados en la lista
  interpolated_data_list[[parametro]] <- interpolated_data
}


# Combinar todos los dataframes de la lista en uno solo
interpolated_data_combined <- bind_rows(interpolated_data_list, .id = "Parameter")

# Visualizar los primeros registros del nuevo dataframe
head(interpolated_data_combined)

#Pasar los índices a valores de porcentaje

interpolated_data_combined <- interpolated_data_combined %>%
  group_by(Long_onda) %>%
  mutate(porcentaje = Indice / sum(Indice) * 100) %>%
  mutate(porcentaje2= (Indice - min(Indice))/(max(Indice)-min(Indice))*100) %>%
  ungroup() %>%
  arrange(Long_onda, Parameter) %>%
  group_by(Long_onda) %>%
  mutate(Cumulative_Index = cumsum(porcentaje2))

head(interpolated_data_combined)

# Crear un gráfico de líneas con áreas sombreadas
plot.Johnson1nm<- ggplot(interpolated_data_combined, aes(x = Long_onda, y = porcentaje, fill = Parameter)) +
  geom_line() +
  geom_ribbon(aes(ymin = 0, ymax = porcentaje), alpha = 0.3) +
  scale_fill_manual(values = my_colors) +
  labs(x = "Longitud de onda", y = "Índices de Johnson (%)") +
  scale_color_discrete(name = "Parameter") +
  scale_fill_discrete(name = "Parameter") +
  theme_minimal() + scale_y_continuous(limits = c(0, 105))


print(plot.Johnson1nm)
ggsave("Plots/Sensibility/plot.johnson1nm_VM.svg", plot.Johnson1nm, width = 10, height = 6, units = "in", dpi = 300)
