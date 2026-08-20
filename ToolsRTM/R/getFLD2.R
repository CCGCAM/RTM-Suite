

#' Get FLD2
#'
#' @param values (vector with radiance) 
#' @param wavelengths  vector with wavelengths
#' @param irradiance  irradiance
#'
#' @return get FLD
#' @export
#'
#' 
getFLD2 <- function(values, wavelengths, irradiance) {

  Eout <- irradiance[which(irradiance[,1] == 750),2] #a750
  Ein <- irradiance[which(irradiance[,1] == 761),2] #b 761 
  v750_762 = signal::interp1(wavelengths, values, c(750, 761), extrap = F) ##radiance
  Lout <- v750_762[1] #c
  Lin <- v750_762[2] #d
  #R<-((c-d)/(a-b)) 
  R<-(Lout-Lin)/(Eout-Ein)
  #f=d-Rb 
  FLD2 <- (Lin-R*Ein)
 #print(FLD2)
  return(FLD2)
}

