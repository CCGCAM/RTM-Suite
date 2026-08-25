################################################################################
############################## ÍNDICES DE JOHNSON
################################################################################
rm(list= ls())
library(sensitivity)
library(ggplot2)
library(dplyr)
library(rlang)
library(signal)

datos <- read.csv("Tables/Tables_Erika/LUTVM-20k-SAIL-PROSPECT-PRO_withIndices.csv", sep=',', header = TRUE, nrows=2000)
summary(datos)
head(datos)

# Datos de salida
bands <- datos[, c("B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A","B11", "B12")]

#Datos de entrada
X <- as.data.frame(datos[,c('Cab', 'Car', 'Anth', 'CBC','Prot','Cbrown', 'N', 'LIDFa', 'LAI', 'EWT', 'hspot','psi')])


indices_df <- data.frame(Banda = character(), Parametro=character(), Indice = numeric())

for (banda in colnames(bands)) {
  y <- bands[, banda]
  x <- johnson(X, y)
  
  # Agregar filas al data frame
  indices_df <- rbind(indices_df, data.frame(Banda = banda, Parametro = x$johnson, Indice = x$johnson$original )) 
}


indices_df$Parametro <- rownames(indices_df)  # Extrae los nombres de los parámetros
rownames(indices_df) <- NULL

# Reemplazando los valores 'Car1', 'Car2', 'Car3' por 'Car'

indices_df <- indices_df %>%
  mutate(Parametro = case_when(
    grepl("^Cab\\d+$", Parametro) ~ "Cab",
    grepl("^Car\\d+$", Parametro) ~ "Car",
    grepl("^Anth\\d+$", Parametro) ~ "Anth",
    grepl("^CBC\\d+$", Parametro) ~ "CBC",
    grepl("^Prot\\d+$", Parametro) ~ "Prot",
    grepl("^Cbrown\\d+$", Parametro) ~ "Cbrown",
    grepl("^N\\d+$", Parametro) ~ "N",
    grepl("^LIDFa\\d+$", Parametro) ~ "LIDFa",
    grepl("^LAI\\d+$", Parametro) ~ "LAI",
    grepl("^EWT\\d+$", Parametro) ~ "EWT",
    grepl("^hspot\\d+$", Parametro) ~ "hspot",
    grepl("^psi\\d+$", Parametro) ~ "psi",
    TRUE ~ Parametro
  ))

# Convertir la variable Banda a factor con el orden deseado
indices_df$Banda <- factor(indices_df$Banda, levels = c("B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A","B11", "B12")) 
head(indices_df)

ggplot(indices_df, aes(x = Banda, y = Indice, fill = Parametro)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "",
       x = "Banda",
       y = "Índices de Johnson",
       fill = "Parámetro") +
  theme_bw()

###########################################################################################

ind_johnson <- indices_df %>%
  group_by(Banda) %>%
  mutate(
    porcentaje = original / sum(original) * 100,
    Long_onda = case_when(
      Banda == "B2" ~ 520,
      Banda == "B3" ~ 560,
      Banda == "B4" ~ 654,
      Banda == "B5" ~ 701,
      Banda == "B6" ~ 743,
      Banda == "B7" ~ 779,
      Banda == "B8" ~ 789,
      Banda == "B8A" ~ 871,
      Banda == "B11" ~ 1639,
      Banda == "B12" ~ 2256,
      TRUE ~ NA_real_
    )
  )
head(ind_johnson)

ggplot(ind_johnson, aes(x = Banda, y = porcentaje, fill = Parametro)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "",
       x = "Banda",
       y = "Índices de Johnson",
       fill = "Parámetro") +
  theme_bw()


################################################################################

# Agrupar los datos por el parámetro
grouped_data <- ind_johnson %>%
  group_by(Parametro)

# Secuencia de longitudes de onda
longitudes_onda_interp <- seq(520, 2256, by = 1)

# Crear una lista para almacenar los resultados
interpolated_data_list <- list()

# Iterar sobre cada parámetro único
for (parametro in unique(ind_johnson$Parametro)) {
  # Filtrar los datos para el parámetro actual
  current_data <- subset(ind_johnson, Parametro == parametro)
  
  # Interpolar los valores de porcentaje utilizando spline
  indice_interp <- signal::interp1(
    x = current_data$Long_onda,
    y = current_data$Indice,
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
interpolated_data_combined <- bind_rows(interpolated_data_list, .id = "Parametro")

# Visualizar los primeros registros del nuevo dataframe
head(interpolated_data_combined)

#Pasar los índices a valores de porcentaje

interpolated_data_combined <- interpolated_data_combined %>%
  group_by(Long_onda) %>%
  mutate(
    porcentaje = Indice / sum(Indice) * 100,
  )

head(interpolated_data_combined)


# Crear un gráfico de líneas con áreas sombreadas
ggplot(interpolated_data_combined, aes(x = Long_onda, y = porcentaje, color = Parametro, fill = Parametro)) +
  geom_line() +
  geom_ribbon(aes(ymin = 0, ymax = porcentaje), alpha = 0.3) +
  labs(x = "Longitud de onda", y = "Índice de Johnson (%)") +
  scale_color_discrete(name = "Parámetro") +
  scale_fill_discrete(name = "Parámetro") +
  theme_minimal() + scale_y_continuous(limits = c(0, 100))



# Crear un gráfico de barras apiladas
ggplot(interpolated_data_combined, aes(x = Long_onda, y = porcentaje, fill = Parametro)) +
  geom_bar(stat = "identity") +
  labs(x = "Longitud de onda", y = "Índice") +
  scale_fill_discrete(name = "Parámetro") +
  theme_minimal()+ scale_y_continuous(limits = c(0, 110))

######################################################################
#----Grafico de un solo parámetro


interpolated_data_combinedCab <- subset(interpolated_data_combined, Parametro == "LAI")

ggplot(interpolated_data_combinedCab, aes(x = Long_onda, y = porcentaje )) +
  geom_line() +
  geom_area(fill = "Blue", alpha = 0.5) +  # Agregar área bajo la curva
  labs(x = "Longitud de Onda (nm)", y = "Porcentaje") +
  ggtitle("") +
  theme_minimal() + scale_y_continuous(limits = c(0, 100))


#############################################################################################
###Aqui está solo de un parámetro

ind_johnson_N<- subset(ind_johnson,Parametro == "N")

longitudes_onda_interp <- seq(520, 2256, by = 1)
# Interpolar los valores de porcentaje utilizando spline
porcentajes_interp_splineN <- interp1(
  x = ind_johnson_N$Long_onda,
  y = ind_johnson_N$porcentaje,
  xi = longitudes_onda_interp,
  method = "spline"
)

# Crear un tibble con los resultados interpolados utilizando spline
interpolated_data_splineN <- tibble(
  Long_onda = longitudes_onda_interp,
  Porcentaje = porcentajes_interp_splineN
)

# Graficar los datos interpolados utilizando spline
ggplot(interpolated_data_splineN, aes(x = Long_onda, y = Porcentaje)) +
  geom_line() +
  geom_area(fill = "lightblue", alpha = 0.5) +  # Agregar área bajo la curva
  labs(x = "Longitud de Onda (nm)", y = "Porcentaje") +
  ggtitle("") +
  theme_minimal()



