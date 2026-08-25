
soil.rfl<-read.table('examples/SCOPE_v2.1/input/soil_spectra/soilnew.txt',header =F)
colnames(soil.rfl) <- c('wave','rfl.soil1','rfl.soil2','rfl.soil3')
save(soil.rfl, file=paste('data/','soil.rfl','.rda',sep=''))

#Get options for running SCOPE mode
data.opts<-read.table('input/setoptions.csv',header=T, sep=',')
save(data.opts, file=paste('data/','data.opts','.rda',sep=''))

#Get LUT default for running SCOPE mode
SCOPE.LUT.default=read.table('input/LUT_input.csv',header=T,sep=',')
save(SCOPE.LUT.default, file=paste('data/','SCOPE.LUT.default','.rda',sep=''))



constants= load(file='data/constants.rda')
#constants <- data.constant[1:14,]
save(constants, file=paste('data/','constants','.rda',sep=''))

###### directional angles
files<-list.files(path='examples/SCOPE_v2.1/input/directional',pattern = '.dat',full.names = T)
n_files<-list.files(path='examples/SCOPE_v2.1/input/directional',pattern = '.dat',full.names = F)
for (i in 1:3){

  anglesfile<-read.table(files[i],header=F)
  print(head(anglesfile))
  name_table<-sub('\\.dat$', '', n_files[i])
  if (i == 1){
    brdf_angles_no_oversampling = anglesfile
    save(brdf_angles_no_oversampling, file=paste('data/',name_table,'.rda',sep=''))
  } else if( i == 2 ){
    brdf_angles = anglesfile
    save(brdf_angles, file=paste('data/',name_table,'.rda',sep=''))
  } else if( i == 3 ){
    brdf_angles2 = anglesfile
    save(brdf_angles2, file=paste('data/',name_table,'.rda',sep=''))
  }

}


###### dataset data
files<-list.files(path='examples/SCOPE_v2.1/input/dataset dat',pattern = '.dat',full.names = T)
n_files<-list.files(path='examples/SCOPE_v2.1/input/dataset dat',pattern = '.dat',full.names = F)
for (i in 1:10){

  anglesfile<-read.table(files[i],header=F)
  print(head(anglesfile))
  name_table<-sub('\\.dat$', '', n_files[i])
  if (i == 1){
    ea_ = anglesfile
    save(ea_, file=paste('data/',name_table,'.rda',sep=''))
  } else if( i == 2 ){
    p_ = anglesfile
    save(p_, file=paste('data/',name_table,'.rda',sep=''))
  } else if( i == 3 ){
    Rin_ = anglesfile
    save(Rin_, file=paste('data/',name_table,'.rda',sep=''))
  } else if( i == 4 ){
    Rli_ = anglesfile
    save(Rli_, file=paste('data/',name_table,'.rda',sep=''))
  } else if( i == 5 ){
    t_ = anglesfile
    save(t_, file=paste('data/',name_table,'.rda',sep=''))
  } else if( i == 6 ){
    Ta_ = anglesfile
    save(Ta_, file=paste('data/',name_table,'.rda',sep=''))
  } else if( i == 7 ){
    table_Jmax_ = anglesfile
    save(table_Jmax_, file=paste('data/',name_table,'.rda',sep=''))
  } else if( i == 8 ){
    table_Vcmax_ = anglesfile
    save(table_Vcmax_, file=paste('data/',name_table,'.rda',sep=''))
  } else if( i == 9 ){
    u_ = anglesfile
    save(u_, file=paste('data/',name_table,'.rda',sep=''))
  } else if( i == 10 ){
    year_ = anglesfile
    save(year_, file=paste('data/',name_table,'.rda',sep=''))
  }

}

###### Radiation files
files<-list.files(path='examples/SCOPE_v2.1/input/radiationdata',pattern = '.dat',full.names = T)
n_files<-list.files(path='examples/SCOPE_v2.1/input/radiationdata',pattern = '.dat',full.names = F)
for (i in 1:2){
  atmfile<-read.table(files[i],header=F)
  print(head(atmfile))
  name_table<-sub('\\.dat$', '', n_files[i])
  if (i == 1){
    Esky_ = atmfile
    save(Esky_, file=paste('data/',name_table,'.rda',sep=''))
  } else if( i == 2 ){
    Esun_ = atmfile
    save(Esun_, file=paste('data/',name_table,'.rda',sep=''))
  }


}


constants= load(file='data/constants.rda')

library(R.matlab)
### this model is from SPART
Optipar2020_prospectD_BSM2019 <- readMat("examples/SPART_v1.21/Optipar2020_ProspectD_BSM2019.mat")
optipar2020.prospectD.BSM2019<-list()
optipar2020.prospectD.BSM2019$wl <- c(Optipar2020_prospectD_BSM2019$optipar[[1]])
optipar2020.prospectD.BSM2019$nr <- c(Optipar2020_prospectD_BSM2019$optipar[[2]])
optipar2020.prospectD.BSM2019$Kab <- c(Optipar2020_prospectD_BSM2019$optipar[[3]])
optipar2020.prospectD.BSM2019$Kca <- c(Optipar2020_prospectD_BSM2019$optipar[[4]])
optipar2020.prospectD.BSM2019$Ks <- c(Optipar2020_prospectD_BSM2019$optipar[[5]])
optipar2020.prospectD.BSM2019$Kw <- c(Optipar2020_prospectD_BSM2019$optipar[[6]])
optipar2020.prospectD.BSM2019$Kdm <- c(Optipar2020_prospectD_BSM2019$optipar[[7]])
optipar2020.prospectD.BSM2019$phiI <- c(Optipar2020_prospectD_BSM2019$optipar[[8]])
optipar2020.prospectD.BSM2019$phiII <- c(Optipar2020_prospectD_BSM2019$optipar[[9]])
optipar2020.prospectD.BSM2019$KcaV <- c(Optipar2020_prospectD_BSM2019$optipar[[10]])
optipar2020.prospectD.BSM2019$KcaZ <- c(Optipar2020_prospectD_BSM2019$optipar[[11]])
optipar2020.prospectD.BSM2019$Kant <- c(Optipar2020_prospectD_BSM2019$optipar[[12]])
optipar2020.prospectD.BSM2019$KcaV2 <- c(Optipar2020_prospectD_BSM2019$optipar[[13]])
optipar2020.prospectD.BSM2019$phi <- c(Optipar2020_prospectD_BSM2019$optipar[[14]])
optipar2020.prospectD.BSM2019$GSV <- Optipar2020_prospectD_BSM2019$optipar[[15]]
optipar2020.prospectD.BSM2019$nw <- c(Optipar2020_prospectD_BSM2019$optipar[[16]])


save(optipar2020.prospectD.BSM2019, file = "data/optipar2020.prospectD.BSM2019.rda")

# library to read matlab data formats into R
library(R.matlab)
# read in our data
optipar2017_ProspectD <- readMat("examples/SCOPE_v2.1/input/fluspect_parameters/Optipar2017_ProspectD.mat")
head(optipar2017_ProspectD)
optipar2017.ProspectD<-list()
optipar2017.ProspectD$wl <- c(optipar2017_ProspectD$optipar[[1]])
optipar2017.ProspectD$nr <- c(optipar2017_ProspectD$optipar[[2]])
optipar2017.ProspectD$Kab <- c(optipar2017_ProspectD$optipar[[3]])
optipar2017.ProspectD$Kca <- c(optipar2017_ProspectD$optipar[[4]])
optipar2017.ProspectD$Ks <- c(optipar2017_ProspectD$optipar[[5]])
optipar2017.ProspectD$Kw <- c(optipar2017_ProspectD$optipar[[6]])
optipar2017.ProspectD$Kdm <- c(optipar2017_ProspectD$optipar[[7]])
optipar2017.ProspectD$phiI <- c(optipar2017_ProspectD$optipar[[8]])
optipar2017.ProspectD$phiII <- c(optipar2017_ProspectD$optipar[[9]])
optipar2017.ProspectD$KcaV <- c(optipar2017_ProspectD$optipar[[10]])
optipar2017.ProspectD$KcaZ <- c(optipar2017_ProspectD$optipar[[11]])
optipar2017.ProspectD$Kant <- c(optipar2017_ProspectD$optipar[[12]])
optipar2017.ProspectD$phi <- c(optipar2017_ProspectD$optipar[[13]])
optipar2017.ProspectD$GSV <- optipar2017_ProspectD$optipar[[14]]
optipar2017.ProspectD$nw <- c(optipar2017_ProspectD$optipar[[15]])
save(optipar2017.ProspectD, file = "data/optipar2017_ProspectD.rda")

optipar2021_Pro <- readMat("examples/SCOPE_v2.1/input/fluspect_parameters/Optipar2021_ProspectPRO_CX.mat")
head(optipar2021_Pro)
optipar2021.Pro.CX<-list()
optipar2021.Pro.CX$wl <- c(optipar2021_Pro$optipar[[1]])
optipar2021.Pro.CX$nr <- c(optipar2021_Pro$optipar[[2]])
optipar2021.Pro.CX$Kab <- c(optipar2021_Pro$optipar[[3]])
optipar2021.Pro.CX$Kca <- c(optipar2021_Pro$optipar[[4]])
optipar2021.Pro.CX$Ks <- c(optipar2021_Pro$optipar[[5]])
optipar2021.Pro.CX$Kw <- c(optipar2021_Pro$optipar[[6]])
optipar2021.Pro.CX$Kdm <- c(optipar2021_Pro$optipar[[7]])
optipar2021.Pro.CX$phiI <- c(optipar2021_Pro$optipar[[8]])
optipar2021.Pro.CX$phiII <- c(optipar2021_Pro$optipar[[9]])
optipar2021.Pro.CX$KcaV <- c(optipar2021_Pro$optipar[[10]])
optipar2021.Pro.CX$KcaZ <- c(optipar2021_Pro$optipar[[11]])
optipar2021.Pro.CX$Kant <- c(optipar2021_Pro$optipar[[12]])
optipar2021.Pro.CX$phi <- c(optipar2021_Pro$optipar[[13]])
optipar2021.Pro.CX$GSV <- optipar2021_Pro$optipar[[14]]
optipar2021.Pro.CX$nw <- c(optipar2021_Pro$optipar[[15]])
optipar2021.Pro.CX$Kp <- c(optipar2021_Pro$optipar[[16]])
optipar2021.Pro.CX$Kcbc <- c(optipar2021_Pro$optipar[[17]])
optipar2021.Pro.CX$phiE <- c(optipar2021_Pro$optipar[[18]])
save(optipar2021.Pro.CX, file = "data/optipar2021_ProspectPRO_CX.rda")

