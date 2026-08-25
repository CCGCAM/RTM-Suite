
#' Get LUT table for SCOPE model 
#'
#' @param inputs a dataframe with the ranges and the distribution
#' @param dataICOS Flux network with radiation, air temperature ...
#' @param timeStart the initial time step in this format '2018-01-01'
#' @param timeEnd the initial time step in this format '2018-01-30'
#' @param freq.hour the diurnal time or hourly time options: 'diurnal', 'hourly' NULL at 12.00 GTM
#' @param nLUT number of rows for each time step
#' @param SetFixedSeed  Fixed seed
#'
#' @return a dataframe with the LUT for Time series
#' @export
#'
#' @examples here adding examples ....
#' 
getLUT_time<-function(inputs=NULL, dataICOS= NULL, timeStart='2018-01-01', timeEnd='2018-01-30', freq.hour = 'diurnal', nLUT=100,SetFixedSeed=T){
  
  
  if (!require("lubridate")) { install.packages("lubridate"); require("lubridate") }  ### for time spteps
  #inputs =read.table('examples/SCOPE/input_bordersM.csv',sep=',',header=T)
  hours= seq(from=as.POSIXct(timeStart, tz="GMT"),to=as.POSIXct(as.Date(timeStart)+1, tz="GMT"), by="hour" ) 
  diurnal= format(hours, format="%H:%M:%S")[9:18]
  all.times<-hours[1:24]
  no_hour = hours[13]
  if (freq.hour == 'diurnal'){
    time.step = diurnal
  } else if (is.null(freq.hour)){
    time.step = no_hour
  } else if (freq.hour == 'hourly'){
    time.step = all.times
  }
  
  
  if (dim(inputs)[2] != 10) {
    message(' Please provide a LUT table with 10 columns')
    print('with theses columns and names: variable, lower,upper,units,Distribution,Mean_D,Std_D,Dependecies,include,default')
    stop()
  } 
  
  if (is.null(nLUT)){
    message('number of varitions for each input is fixed to 100')
    nLUT = 100
  }

  if (is.null(timeStart) | is.null(timeEnd)){
    time=seq(1:nLUT)
    
  } else if(is.null(timeEnd)){
    message('only 30 days from intial date will be computed')
    time_nLUT=sort(rep(seq(as.Date(timeStart),as.Date(timeStart) + 100, "days"),nLUT))
    Year<-substr(as.character(time_nLUT),1,4)
    time = lubridate::ymd_hms(paste(time_nLUT,time.step,sep=''),tz='GMT')
    
  } else{
    
    time_nLUT=sort(rep(seq(as.Date(timeStart),as.Date(timeEnd), "days"),nLUT))
    Year<-substr(as.character(time_nLUT),1,4)
    time = lubridate::ymd_hms(paste(time_nLUT,time.step,sep=''),tz='GMT')
  }
  
  message('step 1: processing LUT table for each time step')
  pb <- progress::progress_bar$new(total = length(time),  format = "  Processing [:bar] :percent eta: :eta",
                                  clear = FALSE, width= 60)
  list.time<-list()
  var_to_mean <- c('windS','SW_in','LW_in','pressure','Tair', 'RH')
  
  data.ICOS$time <- round(data.ICOS$time)
  data.ICOS.t<-aggregate(data.ICOS[,var_to_mean],by=list(time=data.ICOS$time.UT,Day=data.ICOS$Day,Month=data.ICOS$Month, Year = data.ICOS$Year  ),mean)
  data.ICOS.t$Date = as.POSIXct(paste(data.ICOS.t$Year, data.ICOS.t$Month, data.ICOS.t$Day,sep='-'), tz="GMT")
  data.ICOS.t$DateS <- data.ICOS.t$Date + data.ICOS.t$time * 3600 
  minute(data.ICOS.t$DateS) <-0
  data.ICOS.t<-aggregate(data.ICOS.t[,c('Year', 'Month', 'Day',var_to_mean)],by=list(Date=data.ICOS.t$DateS),mean)
  #head(data.ICOS.t)
  data.ICOS.t$Date <-as.POSIXct(data.ICOS.t$Date, tz='GMT')
  
  ### Get variation for all parameters
  for (i in c(1:length(time))){
    pb$tick()
      
      ## Get LUT for all inputs
      if (SetFixedSeed == T | is.null(SetFixedSeed)){
           list.time[[i]]<-data.frame(t=time[i], getLUT(inputs = inputs, nLUT = nLUT, setseed = 123))
         }  else{
          list.time[[i]]<-data.frame(t=time[i], getLUT(inputs = inputs, nLUT = nLUT, setseed = i*123))
         } 
    
      ## For those time[i] where ICOS site has value of Rin, Rli, Ta and p and Wind and RH
       df.icos <- subset(data.ICOS.t, Date == time[i])
       df.inputs <-as.data.frame(list.time[[i]])
       if (rlang::is_empty(df.icos)) {
         list.time[[i]] <- list.time[[i]]
       } else{
         #broadband incoming shortwave radiation (0.4-2.5 um)
         if (rlang::is_empty(df.icos$SW_in)){
           df.inputs$Rin <- df.inputs$Rin
         } else {
           df.inputs$Rin <- rep(df.icos$SW_in,nLUT)
         }
         #broadband incoming longwave radiation (2.5-50 um)
         if (rlang::is_empty(df.icos$LW_in)){
           df.inputs$Rli <- df.inputs$Rli
         } else {
           df.inputs$Rli <- rep(df.icos$LW_in,nLUT)
         }
         # Wind velocity
         if (rlang::is_empty(df.icos$u)){
           df.inputs$u <- df.inputs$u
         } else {
           df.inputs$u <- rep(df.icos$windS,nLUT)
         }
         # Pressure
         if (rlang::is_empty(df.icos$pressure)){
           df.inputs$p <- df.inputs$p
         } else {
           df.inputs$p <- rep(df.icos$pressure*10,nLUT)
         }
         ## Air temperature
         if (rlang::is_empty(df.icos$Tair)){
           df.inputs$Ta <- df.inputs$Ta
         } else {
           df.inputs$Ta <- rep(df.icos$Tair,nLUT)
         }
         ## Relative humidity
         if (rlang::is_empty(df.icos$RH)){
           df.inputs$RH <- df.inputs$RH
         } else {
           df.inputs$RH <- rep(df.icos$RH/ 100,nLUT)
         }
         
         list.time[[i]] <-df.inputs
       }
       
       
      
  } # end for


 
  df <- do.call("rbind",list.time)
  dim(df)
  Year<-substr(as.character(df[,'t']),1,4)
  Month<-substr(as.character(df[,'t']),6,7)
  Day<-substr(as.character(df[,'t']),9,10)
  Hour<-substr(as.character(df[,'t']),12,13)
  time_<-as.numeric(paste(Year,Month,Day,Hour,'00',sep=''))
  df<- cbind(time_,df)
  
  return(df)
}
