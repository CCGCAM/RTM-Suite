

#' Get FLD2
#'
#' @param values (vector with radiance) 
#' @param wavelengths  vector with wavelengths
#' @param irradiance  irradiance
#'
#' @return FLD2
#' @export
#'
#' 
getFLD2_SE <- function(values, wavelengths, irradiance) {

  Eout <- irradiance[which(irradiance[,1] == 750),2] #a750
  Ein <- irradiance[which(irradiance[,1] == 761),2] #b 761 
  v740_782 = signal::interp1(wavelengths, values, c(740, 782.8), extrap = F) ##radiance
  Lout <- v740_782[1] #c
  Lin <- v740_782[2] #d
  #R<-((c-d)/(a-b)) 
  R<-(Lout-Lin)/(Eout-Ein)
  #f=d-Rb 
  FLD2 <- (Lin-R*Ein)
 #print(FLD2)
  return(FLD2)
}

