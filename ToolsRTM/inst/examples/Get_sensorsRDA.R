
require(dplyr)
## Sentinel 2A
Sensors <-c('ALI', 'Hyperion', 'Landsat4', 'Landsat5', 'Landsat7', 'Landsat8',
            'MODIS', 'Quickbird', 'RapidEye', 'Sentinel2a', 'Sentinel2b', 'WorldView2-4', 'WorldView2-8')
ls.sensors<-list()
for (i in c(1:length(Sensors))){
  print(i)
  data_s <-hsdar::get.sensor.characteristics(Sensors[i], response_function=FALSE)
  if(Sensors[i] == 'EnMAP'){
    ls.sensors[[i]]  <- data_s %>%
      mutate(Sensor = Sensors[i]) %>%
      dplyr::select(Sensor, channel, center, fwhm)
  } else {
    # Calculate spectral average for each channel using dplyr
    data_s <- data_s %>%
      mutate(average = (lb + ub) / 2)  %>%
      mutate(Sensor = Sensors[i])
    # Reorder the columns to have the sensor column before the others
    ls.sensors[[i]]  <- data_s %>%
      dplyr::select(Sensor, channel, lb, ub, average)
  }

}

# Combine all data frames into one
sensor.characteristics <- do.call(rbind, ls.sensors)
save(sensor.characteristics, file = "data/sensor.characteristics.rda")


# EnMap
data_s <-hsdar::get.sensor.characteristics('EnMAP', response_function=FALSE)
EnMap.characteristics  <- data_s %>%
  mutate(Sensor = 'EnMAP') %>%
  dplyr::select(Sensor, channel, center, fwhm)
head(EnMap.characteristics)
save(EnMap.characteristics, file = "data/EnMap.characteristics.rda")




