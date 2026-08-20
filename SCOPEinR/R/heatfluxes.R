#' @title get.heatfluxes
#' \code{get.heatfluxes} Calculates  latent and sensible heat flux
#'
#' @param ra   aerodynamic resistance for heat         s m-1
#' @param rs   stomatal resistance                     s m-1
#' @param Tc   leaf temperature                        oC
#' @param ea   vapour pressure above canopy            hPa
#' @param Ta   air temperature above canopy            oC
#' @param e_to_q  conv. from vapour pressure to abs hum   hPa-1
#' @param Ca  ambient CO2 concentration               umol m-3
#' @param Ci  intercellular CO2 concentration         umol m-3
#'
#' @return a list with
#'  lEc         latent heat flux of a leaf              W m-2
#   Hc          sensible heat flux of a leaf            W m-2
#   ec          vapour pressure at the leaf surface     hPa
#   Cc          CO2 concentration at the leaf surface   umol m-3
#' @export
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
# last updates
# First: 7 Dec 2007
# updated: 15 Apr 2009 CvdT     changed layout
# updated: 14 Sep 2012 CvdT     added ec and Cc to output
# updated: 09 Dec 2019 CvdT     modified for computational efficiency
#' @examples
#' \dontrun{
#' out <- get.heatfluxes(ra, rs, Tc, ea, Ta, e_to_q, Ca, Ci)
#' }


get.heatfluxes<-function(ra,rs,Tc,ea,Ta,e_to_q,Ca,Ci){

  # dependencies functions: ebal.R

  # functions for saturated vapour pressure
  #function for saturated pressure function es(hPa)=f(T(C))
  es_fun <- function(Temp) (6.107 * 10 ^(7.5 * Temp /(237.3 + Temp)))

  # function for slope of the saturated pressure function (s(hPa/C) = f(T(C), es(hPa))
  s_fun <- function(ei,Temp) (ei * 2.3026 * 7.5 * 237.3 / (237.3 + Temp)^2)

  # output:
  #   lEc         latent heat flux of a leaf              W m-2
  #   Hc          sensible heat flux of a leaf            W m-2
  #   ec          vapour pressure at the leaf surface     hPa
  #   Cc          CO2 concentration at the leaf surface   umol m-3

  heatfluxes<-list()


  rhoa = subset(SCOPEinR::constants,constant == 'rhoa')[[2]]
  cp = subset(SCOPEinR::constants,constant == 'cp')[[2]]
  #MH2O =  subset(SCOPEinR::constants,constant == 'MH2O')[[2]]
  #R  =  subset(SCOPEinR::constants,constant == 'R')[[2]]

  # Evapor. heat in J kg-1
  lambda = (2.501 - 0.002361 * Tc) * 1e6

  heatfluxes[['lambda']] <- lambda

  ei <- es_fun(Tc)

  s <- s_fun(ei, Tc)

  heatfluxes[['s']] <- s

  ##################################################################
  ## Old function use PSI as main argument in the get.heatfluxes
  #ei <- es * exp(1E-3 * PSI * MH2O /R / (Tc + 273.15))
  ##################################################################

  qi <- ei * e_to_q

  qa <-  ea * e_to_q

  # Latent heat flux in W m-2
  heatfluxes[['lE']] <- rhoa / (ra + rs) * lambda * (qi - qa)

  # Sensible heat flux in W m-2
  heatfluxes[['H']] <- (rhoa * cp) /ra *(Tc - Ta)

  # [W m-2] vapour pressure at the leaf surface
  heatfluxes[['ec']] <- ea + (ei - ea) * ra / (ra + rs)

  # CO2 concentration at the leaf surface in umol m-2 s-1
  heatfluxes[['Cc']] <- Ca - (Ca - Ci) * ra / (ra + rs)


  return(heatfluxes)

}

