
#' Get CWSI for almonds
#'
#' @param df  dataframe
#' @param Ta  air temperature in degrees
#' @param HR  relative humidity in %
#' @param df.data datafrma ewith indicators
#'
#' @return CWSI
#' @export
#'
#' 
getCWSI <- function(df, Ta, HR,df.data) {

  cwsi.list = list()
  for (i in c(1:dim(df)[1])){
  values<-df[i,]
  ##correccion Apoggie
  values<-0.8463*values+7.7403
  ##la presion de vapor en saturacion 
  esTa<-0.61078*exp((17.269*Ta)/(Ta+237.3))
  ea=HR/100*esTa
 
  #El DPV es la diferencia entre la presi?n de saturaci?n y la presi?n actual de vapor
  #DPV=es-ea
  #DPV<- ((100-HR)/100)*esTa
  DPV=esTa-ea 
  
  ##### Linea Base para el Almendro
  #BaseLine<--1.23*DPV+3.65 #(almedros v1)
  ### (Tc-TaLL) --> +-1.0078*DPV+3.765
  VPD_0<-3.765
  BaseLine<-+-1.0078*DPV+ VPD_0  ###(almendros v2) 
  #Tc-Ta (VPD=0) es 3.765
  
  Tc_Ta<-values-Ta
  
  ###la presi?n de vapor en saturaci?n 
  esTc<-0.61078*exp((17.269*(Ta+VPD_0))/((Ta+VPD_0)+237.3))
  AEs<-+esTc-esTa 
  Tc_Ta_UL<-1.0078*AEs+VPD_0
  CWSI<-(Tc_Ta-BaseLine)/(Tc_Ta_UL-BaseLine)
  CWSI<-as.vector(CWSI)
  cwsi.list[[i]] <- CWSI
  }
  
  df.cwsi <- data.frame(matrix(unlist(cwsi.list), nrow=length(cwsi.list), byrow=T))
  colnames(df.cwsi)<-c('CWSI')
  df.final<-cbind(df.data,Tc=df,df.cwsi, Ta=Ta)
  return(df.final)
  
}



