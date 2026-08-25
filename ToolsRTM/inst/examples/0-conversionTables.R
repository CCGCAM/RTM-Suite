
srf.sentinel2a <- read.csv("Tables/srf_sentinel2a.csv")
save(srf.sentinel2a, file = "data/srf.sentinel2a.rda")

srf.sentinel2b <- read.csv("Tables/srf_sentinel2b.csv")
save(srf.sentinel2b, file = "data/srf.sentinel2b.rda")

srf.prisma <- read.csv("Tables/srf_prisma.csv")
save(srf.prisma, file = "data/srf.prisma.rda")

fwhm.prisma <- read.csv("Tables/prisma_fwhm.csv")
save(fwhm.prisma, file = "data/fwhm.prisma.rda")



# library to read matlab data formats into R
library(R.matlab)
require(reticulate)
source_python("examples/SPART_v2/python/Read_pkl.py")

Sentinel2A.MSI <- read_pickle_file("examples/SPART_v2/Coef_SMAC/Sentinel2A-MSI.pkl")
save(Sentinel2A.MSI, file = "data/Sentinel2A.MSI.rda")

Sentinel2B.MSI <- read_pickle_file("examples/SPART_v2/Coef_SMAC/Sentinel2B-MSI.pkl")
save(Sentinel2B.MSI, file = "data/Sentinel2B.MSI.rda")

LANDSAT4.TM <- read_pickle_file("examples/SPART_v2/Coef_SMAC/LANDSAT4-TM.pkl")
save(LANDSAT4.TM, file = "data/LANDSAT4.TM.rda")

LANDSAT5.TM <- read_pickle_file("examples/SPART_v2/Coef_SMAC/LANDSAT5-TM.pkl")
save(LANDSAT5.TM, file = "data/LANDSAT5.TM.rda")

LANDSAT7.ETM <- read_pickle_file("examples/SPART_v2/Coef_SMAC/LANDSAT7-ETM.pkl")
save(LANDSAT7.ETM, file = "data/LANDSAT7.ETM.rda")

LANDSAT8.OLI <- read_pickle_file("examples/SPART_v2/Coef_SMAC/LANDSAT8-OLI.pkl")
save(LANDSAT8.OLI, file = "data/LANDSAT8.OLI.rda")

Sentinel3B.OLCI <- read_pickle_file("examples/SPART_v2/Coef_SMAC/Sentinel3B-OLCI.pkl")
save(Sentinel3B.OLCI, file = "data/Sentinel3B.OLCI.rda")

Sentinel3A.OLCI <- read_pickle_file("examples/SPART_v2/Coef_SMAC/Sentinel3A-OLCI.pkl")
save(Sentinel3A.OLCI, file = "data/Sentinel3A.OLCI.rda")

TerraAqua.MODIS <- read_pickle_file("examples/SPART_v2/Coef_SMAC/TerraAqua-MODIS.pkl")
save(TerraAqua.MODIS, file = "data/TerraAqua.MODIS.rda")


optical.parameters <- read_pickle_file("examples/SPART_v2/optical_params_prospect-d.pkl")
save(TerraAqua.MODIS, file = "data/TerraAqua.MODIS.rda")




inputsRTMs <- read.csv("examples/LUTs/inputs_Leaf_canopy_v2.csv")
save(inputsRTMs, file = "data/inputsRTMs.rda")

inputsPROSAIL <- read.csv("examples/LUTs/inputs_PROSAIL.csv")
save(inputsPROSAIL, file = "data/inputsPROSAIL.rda")

inputsLiberty <- read.csv("examples/LUTs/inputs_Liberty.csv")
save(inputsLiberty, file = "data/inputsLiberty.rda")

inputsINFORM <- read.csv("examples/LUTs/inputs_INFORM.csv")
save(inputsINFORM, file = "data/inputsINFORM.rda")

inputsFlUSPECT <- read.csv("examples/LUTs/inputs_FLUSPECT.csv")
save(inputsFlUSPECT, file = "data/inputsFlUSPECT.rda")

inputsSPART <- read.csv("examples/LUTs/inputs_SPART.csv")
save(inputsSPART, file = "data/inputsSPART.rda")



# read in our data
inputs.SPART <- read.csv("examples/SPART_v2/inputs_SPART.csv")

save(inputs.SPART, file = "data/inputs.SPART.rda")

# read in our data
EIrrad <- readMat("examples/SPART_v2/TOC2TOA_inputs/Extraterrestrial_irradiance.mat")
Extraterrestrial_irradiance <- data.frame(wave = EIrrad$wl.Ea, EIrrad = EIrrad$Ea)
save(Extraterrestrial_irradiance, file = "data/Extraterrestrial_irradiance.rda")

ggplot(data = Extraterrestrial_irradiance, aes(x = wave, y = EIrrad)) +
  labs(y= "soil reflectance", x = "")+ xlim(400,2500) +
  geom_line() + theme_bw()


# read in our data
Optipar2020.ProspectD.BSM2019 <- readMat("examples/SPART_v2/PROSPECT_BSM_inputs/Optipar2020_ProspectD_BSM2019.mat")
Optipar2020_ProspectD_BSM2019


# read in our data
Landsat4 <- readMat("examples/SPART_v2/TOC2TOA_inputs/sensors_config_SMAC/spart_sensor_info_LANDSAT4-TM.mat")
sensorLandsat4 <-list()
sensorLandsat4[['mission']] <- c(Landsat4$sensor[[1]])
sensorLandsat4[['name']] <- "TM"  
sensorLandsat4[['band.id.all']] <- c(1:7)
sensorLandsat4[['res.spatials']] <- as.data.frame(Landsat4$sensor[[4]])
sensorLandsat4[['rang.wvls']] <- as.data.frame(Landsat4$sensor[[5]])
sensorLandsat4[['swath.widths']] <- as.data.frame(Landsat4$sensor[[6]])
sensorLandsat4[['revisit.days']] <- as.data.frame(Landsat4$sensor[[7]])
sensorLandsat4[['band.width']] <- as.data.frame(Landsat4$sensor[[8]])
sensorLandsat4[['center.wvl']] <- as.data.frame(Landsat4$sensor[[9]])
sensorLandsat4[['band.id.smac']] <- as.data.frame(Landsat4$sensor[[10]])
sensorLandsat4[['SMAC.coef']] <- as.data.frame(Landsat4$sensor[[11]])
sensorLandsat4[['wl.smac']] <- as.data.frame(Landsat4$sensor[[12]])
sensorLandsat4[['p.srf']] <- as.data.frame(Landsat4$sensor[[13]])
sensorLandsat4[['wl.srf']] <- as.data.frame(Landsat4$sensor[[14]])
sensorLandsat4[['id.smac.in.all']] <- as.data.frame(Landsat4$sensor[[15]])
sensorLandsat4[['p.srf.smac ']] <- as.data.frame(Landsat4$sensor[[16]])
sensorLandsat4[['wl.srf.smac ']] <- as.data.frame(Landsat4$sensor[[17]])
save(sensorLandsat4, file = "data/sensorLandsat4.rda")

# read in our data
Landsat5 <- readMat("examples/SPART_v2/TOC2TOA_inputs/sensors_config_SMAC/spart_sensor_info_LANDSAT5-TM.mat")
sensorLandsat5 <-list()
sensorLandsat5[['mission']] <- c(Landsat5$sensor[[1]])
sensorLandsat5[['name']] <- "TM"  
sensorLandsat5[['band.id.all']] <- c(1:7)
sensorLandsat5[['res.spatials']] <- as.data.frame(Landsat5$sensor[[4]])
sensorLandsat5[['rang.wvls']] <- as.data.frame(Landsat5$sensor[[5]])
sensorLandsat5[['swath.widths']] <- as.data.frame(Landsat5$sensor[[6]])
sensorLandsat5[['revisit.days']] <- as.data.frame(Landsat5$sensor[[7]])
sensorLandsat5[['band.width']] <- as.data.frame(Landsat5$sensor[[8]])
sensorLandsat5[['center.wvl']] <- as.data.frame(Landsat5$sensor[[9]])
sensorLandsat5[['band.id.smac']] <- as.data.frame(Landsat5$sensor[[10]])
sensorLandsat5[['SMAC.coef']] <- as.data.frame(Landsat5$sensor[[11]])
sensorLandsat5[['wl.smac']] <- as.data.frame(Landsat5$sensor[[12]])
sensorLandsat5[['p.srf']] <- as.data.frame(Landsat5$sensor[[13]])
sensorLandsat5[['wl.srf']] <- as.data.frame(Landsat5$sensor[[14]])
sensorLandsat5[['id.smac.in.all']] <- as.data.frame(Landsat5$sensor[[15]])
sensorLandsat5[['p.srf.smac ']] <- as.data.frame(Landsat5$sensor[[16]])
sensorLandsat5[['wl.srf.smac ']] <- as.data.frame(Landsat5$sensor[[17]])
save(sensorLandsat5, file = "data/sensorLandsat5.rda")



# read in our data
Landsat7 <- readMat("examples/SPART_v2/TOC2TOA_inputs/sensors_config_SMAC/spart_sensor_info_LANDSAT7-ETM.mat")
sensorLandsat7 <-list()
sensorLandsat7[['mission']] <- c(Landsat7$sensor[[1]])
sensorLandsat7[['name']] <- "ETM"  
sensorLandsat7[['band.id.all']] <- c(1:7)
sensorLandsat7[['res.spatials']] <- as.data.frame(Landsat7$sensor[[4]])
sensorLandsat7[['rang.wvls']] <- as.data.frame(Landsat7$sensor[[5]])
sensorLandsat7[['swath.widths']] <- as.data.frame(Landsat7$sensor[[6]])
sensorLandsat7[['revisit.days']] <- as.data.frame(Landsat7$sensor[[7]])
sensorLandsat7[['band.width']] <- as.data.frame(Landsat7$sensor[[8]])
sensorLandsat7[['center.wvl']] <- as.data.frame(Landsat7$sensor[[9]])
sensorLandsat7[['band.id.smac']] <- as.data.frame(Landsat7$sensor[[10]])
sensorLandsat7[['SMAC.coef']] <- as.data.frame(Landsat7$sensor[[11]])
sensorLandsat7[['wl.smac']] <- as.data.frame(Landsat7$sensor[[12]])
sensorLandsat7[['p.srf']] <- as.data.frame(Landsat7$sensor[[13]])
sensorLandsat7[['wl.srf']] <- as.data.frame(Landsat7$sensor[[14]])
sensorLandsat7[['id.smac.in.all']] <- as.data.frame(Landsat7$sensor[[15]])
sensorLandsat7[['p.srf.smac ']] <- as.data.frame(Landsat7$sensor[[16]])
sensorLandsat7[['wl.srf.smac ']] <- as.data.frame(Landsat7$sensor[[17]])
save(sensorLandsat7, file = "data/sensorLandsat7.rda")


# read in our data
Landsat8 <- readMat("examples/SPART_v2/TOC2TOA_inputs/sensors_config_SMAC/spart_sensor_info_LANDSAT8-OLI.mat")
sensorLandsat8 <-list()
sensorLandsat8[['mission']] <- c(Landsat8$sensor[[1]])
sensorLandsat8[['name']] <- "OLI"  
sensorLandsat8[['band.id.all']] <- c(1:7)
sensorLandsat8[['res.spatials']] <- as.data.frame(Landsat8$sensor[[4]])
sensorLandsat8[['rang.wvls']] <- as.data.frame(Landsat8$sensor[[5]])
sensorLandsat8[['swath.widths']] <- as.data.frame(Landsat8$sensor[[6]])
sensorLandsat8[['revisit.days']] <- as.data.frame(Landsat8$sensor[[7]])
sensorLandsat8[['band.width']] <- as.data.frame(Landsat8$sensor[[8]])
sensorLandsat8[['center.wvl']] <- as.data.frame(Landsat8$sensor[[9]])
sensorLandsat8[['band.id.smac']] <- as.data.frame(Landsat8$sensor[[10]])
sensorLandsat8[['SMAC.coef']] <- as.data.frame(Landsat8$sensor[[11]])
sensorLandsat8[['wl.smac']] <- as.data.frame(Landsat8$sensor[[12]])
sensorLandsat8[['p.srf']] <- as.data.frame(Landsat8$sensor[[13]])
sensorLandsat8[['wl.srf']] <- as.data.frame(Landsat8$sensor[[14]])
sensorLandsat8[['id.smac.in.all']] <- as.data.frame(Landsat8$sensor[[15]])
sensorLandsat8[['p.srf.smac ']] <- as.data.frame(Landsat8$sensor[[16]])
sensorLandsat8[['wl.srf.smac ']] <- as.data.frame(Landsat8$sensor[[17]])
save(sensorLandsat8, file = "data/sensorLandsat8.rda")


# read in our data
Sentinel3A <- readMat("examples/SPART_v2/TOC2TOA_inputs/sensors_config_SMAC/spart_sensor_info_Sentinel3A-OLCI.mat")
sensorSentinel3A <-list()
sensorSentinel3A[['mission']] <- c(Sentinel3A$sensor[[1]])
sensorSentinel3A[['name']] <- "OLCI"  
sensorSentinel3A[['band.id.all']] <- c(1:7)
sensorSentinel3A[['res.spatials']] <- as.data.frame(Sentinel3A$sensor[[4]])
sensorSentinel3A[['rang.wvls']] <- as.data.frame(Sentinel3A$sensor[[5]])
sensorSentinel3A[['swath.widths']] <- as.data.frame(Sentinel3A$sensor[[6]])
sensorSentinel3A[['revisit.days']] <- as.data.frame(Sentinel3A$sensor[[7]])
sensorSentinel3A[['band.width']] <- as.data.frame(Sentinel3A$sensor[[8]])
sensorSentinel3A[['center.wvl']] <- as.data.frame(Sentinel3A$sensor[[9]])
sensorSentinel3A[['band.id.smac']] <- as.data.frame(Sentinel3A$sensor[[10]])
sensorSentinel3A[['SMAC.coef']] <- as.data.frame(Sentinel3A$sensor[[11]])
sensorSentinel3A[['wl.smac']] <- as.data.frame(Sentinel3A$sensor[[12]])
sensorSentinel3A[['p.srf']] <- as.data.frame(Sentinel3A$sensor[[13]])
sensorSentinel3A[['wl.srf']] <- as.data.frame(Sentinel3A$sensor[[14]])
sensorSentinel3A[['id.smac.in.all']] <- as.data.frame(Sentinel3A$sensor[[15]])
sensorSentinel3A[['p.srf.smac ']] <- as.data.frame(Sentinel3A$sensor[[16]])
sensorSentinel3A[['wl.srf.smac ']] <- as.data.frame(Sentinel3A$sensor[[17]])
save(sensorSentinel3A, file = "data/sensorSentinel3A.rda")





# read in our data
Sentinel3B <- readMat("examples/SPART_v2/TOC2TOA_inputs/sensors_config_SMAC/spart_sensor_info_Sentinel3B-OLCI.mat")
sensorSentinel3B<-list()
sensorSentinel3B[['mission']] <- c(Sentinel3B$sensor[[1]])
sensorSentinel3B[['name']] <- "OLCI"  
sensorSentinel3B[['band_id_all']] <- c(1:7)
sensorSentinel3B[['res_spatials']] <- as.data.frame(Sentinel3B$sensor[[4]])
sensorSentinel3B[['rang_wvls']] <- as.data.frame(Sentinel3B$sensor[[5]])
sensorSentinel3B[['swath_widths']] <- as.data.frame(Sentinel3B$sensor[[6]])
sensorSentinel3B[['revisit_days']] <- as.data.frame(Sentinel3B$sensor[[7]])
sensorSentinel3B[['band_width']] <- as.data.frame(Sentinel3B$sensor[[8]])
sensorSentinel3B[['center_wvl']] <- as.data.frame(Sentinel3B$sensor[[9]])
sensorSentinel3B[['band_id_smac']] <- as.data.frame(Sentinel3B$sensor[[10]])
sensorSentinel3B[['SMAC_coef']] <- as.data.frame(Sentinel3B$sensor[[11]])
sensorSentinel3B[['wl_smac']] <- as.data.frame(Sentinel3B$sensor[[12]])
sensorSentinel3B[['p_srf']] <- as.data.frame(Sentinel3B$sensor[[13]])
sensorSentinel3B[['wl_srf']] <- as.data.frame(Sentinel3B$sensor[[14]])
sensorSentinel3B[['id_smac_in_all']] <- as.data.frame(Sentinel3B$sensor[[15]])
sensorSentinel3B[['p_srf_smac ']] <- as.data.frame(Sentinel3B$sensor[[16]])
sensorSentinel3B[['wl_srf_smac ']] <- as.data.frame(Sentinel3B$sensor[[17]])
save(sensorSentinel3B, file = "data/sensorSentinel3B.rda")



# read in our data
MODIS <- readMat("examples/SPART_v2/TOC2TOA_inputs/sensors_config_SMAC/spart_sensor_info_TerraAqua-MODIS.mat")
sensorMODIS<-list()
sensorMODIS[['mission']] <- c(MODIS$sensor[[1]])
sensorMODIS[['name']] <- "TerraAqua"  
sensorMODIS[['band.id.all']] <- c(1:7)
sensorMODIS[['res.spatials']] <- as.data.frame(MODIS$sensor[[4]])
sensorMODIS[['rang.wvls']] <- as.data.frame(MODIS$sensor[[5]])
sensorMODIS[['swath.widths']] <- as.data.frame(MODIS$sensor[[6]])
sensorMODIS[['revisit.days']] <- as.data.frame(MODIS$sensor[[7]])
sensorMODIS[['band.width']] <- as.data.frame(MODIS$sensor[[8]])
sensorMODIS[['center.wvl']] <- as.data.frame(MODIS$sensor[[9]])
sensorMODIS[['band.id.smac']] <- as.data.frame(MODIS$sensor[[10]])
sensorMODIS[['SMAC.coef']] <- as.data.frame(MODIS$sensor[[11]])
sensorMODIS[['wl.smac']] <- as.data.frame(MODIS$sensor[[12]])
sensorMODIS[['p.srf']] <- as.data.frame(MODIS$sensor[[13]])
sensorMODIS[['wl.srf']] <- as.data.frame(MODIS$sensor[[14]])
sensorMODIS[['id.smac.in.all']] <- as.data.frame(MODIS$sensor[[15]])
sensorMODIS[['p.srf.smac ']] <- as.data.frame(MODIS$sensor[[16]])
sensorMODIS[['wl.srf.smac ']] <- as.data.frame(MODIS$sensor[[17]])
save(sensorMODIS, file = "data/sensorMODIS.rda")
