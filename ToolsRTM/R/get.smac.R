
#' Apply an atmospheric correction using SMAC method (Rahman and Dedieu, 1994)
#'
#' This function applies the Simplified Method for Atmospheric Correction (SMAC) to satellite
#' imagery based on the provided LUT (Look-Up Table) and sensor characteristics.
#'
#' @param inputLUT A data frame containing the LUT used for atmospheric correction.
#' @param sensor A dataset with sensor information for which the calculations are performed. Options include "LANDSAT4.TM", "LANDSAT5.TM", 
#'               "LANDSAT5.TM, "LANDSAT7.ETM", "LANDSAT8.OLI","Sentinel2A.MSI", "Sentinel2B.MSI", "Sentinel3A.OLCI",
#'               "Sentinel3B.OLCI"and "TerraAqua.MODIS"
#'
#' @return A data frame with atmospheric corrected values.
#' @export
#'
#' @examples
#' \dontrun{
#' inputLUT = ToolsRTM::inputsSPART
#' LUT <- as.data.frame(ToolsRTM::getLUT(inputs = ToolsRTM::inputsSPART, nLUT = 1, setseed = 1234))
#' corrected_data <- get.smac(inputLUT, sensor = ToolsRTM::TerraAqua.MODIS)
#' 
#' }
get.smac <- function(inputLUT = LUT,sensor)   {
  
  # Translate by Carlos Camino from Matlab to R cran 
  # Improved by Peiqi Yang, ITC, University of Twente, from the original SMAC
  
  # Main Author: Peiqi Yang  (p.yang@utwente.nl)
  # Last update: 26/Aug/2019
  #
  
  ## get spectral ranges for SMAC and SCOPE
  #data.spectral<-ToolsRTM::get.spectra.spart(getSpectral = T)
  data.spectral<-get.spectra.spart(getSpectral = T)
  ## get angles: 
  
  tts <- inputLUT[1,'tts']
  tto <- inputLUT[1,'tto']
  psi <- inputLUT[1,'psi']
  
  ## get atm properties: 
  
  Pa <- inputLUT[1,'Pa']

  taup550 <- inputLUT[1,'aot550']
  uo3 <- inputLUT[1,'uo3']
  uh2o <- inputLUT[1,'uh2o']
  alt_m <- inputLUT[1,'alt_m']
  Pa0 <- inputLUT[1,'Pa0']
  DOY <- inputLUT[1,'DOY']
  FWHM <- inputLUT[1,'FWHM']
  
  
  if ((is.null(Pa[1]) | Pa == -999.000)) {
    #Pa = ToolsRTM::get.Altitude2Pa(altitude= alt_m,Pa0=Pa0)
    Pa = get.Altitude2Pa(altitude= alt_m,Pa0=Pa0)
  }
  
  #sensors.properties = ToolsRTM::get.coef.SMAC(sensor = sensor)
  sensors.properties = get.coef.SMAC(sensor = sensor)
  
  wlSensor =  sensors.properties[['wl.smac']]
  coefs.SMAC  = sensors.properties[['coefs.SMAC']]
  Sensor.name = sensors.properties[['Sensor.name']]
 
  ah2o <- c(unlist(coefs.SMAC['ah2o']))
  nh2o <-  c(unlist(coefs.SMAC['nh2o']))
  ao3 <-  c(unlist(coefs.SMAC['ao3']))
  no3 <-  c(unlist(coefs.SMAC['no3']))
  ao2 <-  c(unlist(coefs.SMAC['ao2']))
  no2 <-  c(unlist(coefs.SMAC['no2']))
  po2 <-  c(unlist(coefs.SMAC['po2']))
  aco2 <-  c(unlist(coefs.SMAC['aco2']))
  nco2 <-  c(unlist(coefs.SMAC['nco2']))
  pco2 <-  c(unlist(coefs.SMAC['pco2']))
  ach4 <-  c(unlist(coefs.SMAC['ach4']))
  nch4 <-  c(unlist(coefs.SMAC['nch4']))
  pch4 <-  c(unlist(coefs.SMAC['pch4']))
  ano2 <-  c(unlist(coefs.SMAC['ano2']))
  nno2 <-  c(unlist(coefs.SMAC['nno2']))
  pno2 <-  c(unlist(coefs.SMAC['pno2']))
  aco <-  c(unlist(coefs.SMAC['aco']))
  nco <-  c(unlist(coefs.SMAC['nco']))
  pco <-  c(unlist(coefs.SMAC['pco']))
  a0s <-  c(unlist(coefs.SMAC['a0s']))
  a1s <-  c(unlist(coefs.SMAC['a1s']))
  a2s <-  c(unlist(coefs.SMAC['a2s']))
  a3s <-  c(unlist(coefs.SMAC['a3s']))
  a0T <-  c(unlist(coefs.SMAC['a0T']))
  a1T <-  c(unlist(coefs.SMAC['a1T']))
  a2T <-  c(unlist(coefs.SMAC['a2T']))
  a3T <-  c(unlist(coefs.SMAC['a3T']))
  taur <-  c(unlist(coefs.SMAC['taur']))
  sr <-  c(unlist(coefs.SMAC['sr']))
  a0taup <-  c(unlist(coefs.SMAC['a0taup']))
  a1taup <-  c(unlist(coefs.SMAC['a1taup']))
  wo <-  c(unlist(coefs.SMAC['wo']))
  gc <-  c(unlist(coefs.SMAC['gc']))
  a0P <-  c(unlist(coefs.SMAC['a0P']))
  a1P <-  c(unlist(coefs.SMAC['a1P']))
  a2P <-  c(unlist(coefs.SMAC['a2P']))
  a3P <-  c(unlist(coefs.SMAC['a3P']))
  a4P <-  c(unlist(coefs.SMAC['a4P']))
  Rest1 <-  c(unlist(coefs.SMAC['Rest1']))
  Rest2 <-  c(unlist(coefs.SMAC['Rest2']))
  Rest3 <-  c(unlist(coefs.SMAC['Rest3']))
  Rest4 <-  c(unlist(coefs.SMAC['Rest4']))
  Resr1 <-  c(unlist(coefs.SMAC['Resr1']))
  Resr2 <-  c(unlist(coefs.SMAC['Resr2']))
  Resr3 <-  c(unlist(coefs.SMAC['Resr3']))
  Resa1 <-  c(unlist(coefs.SMAC['Resa1']))
  Resa2 <-  c(unlist(coefs.SMAC['Resa2']))
  Resa3 <-  c(unlist(coefs.SMAC['Resa3']))
  Resa4 <-  c(unlist(coefs.SMAC['Resa4']))
  

  # --------------calculate the reflectance at TOA------------------------
  
  cdr <- pi / 180
  crd <- 180 / pi
  
  us <- cos(tts * cdr)         # cos(tts)
  uv <- cos(tto * cdr)         # cos(tto)
  Peq <- Pa / 1013.25          # air press/ std air pressure
  
  # ------------ 1) air mass Eq(5) --------------------------------------
  
  m <- 1 / us + 1 / uv        #
  
  # ----- 2) aerosol optical depth in the spectral band, taup Eq. (16)----
  
  taup <- a0taup + a1taup * taup550 # a0 + a1 * tau550
  
  # --3) gaseous transmissions (downward and upward paths) (Eq. 5 and 6)----
  
  # to3 = 1; th2o= 1; to2 = 1; tco2= 1; tch4= 1;
  # UH2O is the water vapour integrated content in g/cm2
  
  uo2 <- (Peq^po2)    # Vertically integrated absorber amount
  uco2 <- (Peq^pco2)
  uch4 <- (Peq^pch4)
  uno2 <- (Peq^pno2)
  uco <- (Peq^pco)
  
  # ao2[11:15] <- ao2[11:15] * 2
  
  to3 <- exp(ao3 * (uo3 * m)^no3)
  th2o <- exp(ah2o * (uh2o * m)^nh2o)
  to2 <- exp(ao2 * (uo2 * m)^no2)
  tco2 <- exp(aco2 * (uco2 * m)^nco2)
  tch4 <- exp(ach4 * (uch4 * m)^nch4)
  tno2 <- exp(ano2 * (uno2 * m)^nno2)
  tco <- exp(aco * (uco * m)^nco)
  
  tg <- th2o * to3 * to2 * tco2 * tch4 * tco * tno2  # Eq. 6
  
  # ------- 5) spherical albedo of the atmosphere Eq. 7------------------------
  
  s <- a0s * Peq + a3s + a1s * taup550 + a2s * (taup550^2)  # modification of Eq.8
  
  # ------- 6) Total scattering transmission Eq. 9 ----------------------------
  
  ttetas <- a0T + a1T * taup550/us + ((a2T) * Peq + a3T) / (1 + us)  # downward
  ttetav <- a0T + a1T * taup550/uv + ((a2T) * Peq + a3T) / (1 + uv)  # upward
  
  # ------- 7) scattering angle cosine  Eq.14 -------------------------------
  
  cksi <- - (us * uv + sqrt(1 - us^2) * sqrt(1 - uv^2) * cos(psi * cdr))
  
  if (cksi < -1) {  # it seems cksi is always >= -1
    cksi <- -1.0
  }
  
  # -------------- 8) scattering angle in degree ---------------------------
  
  ksiD <- crd * acos(cksi)
  
  # -------------- 9) rayleigh atmospheric reflectance ----------------------
  
  ray_phase <- 0.7190443 * (1 + cksi^2) + 0.0412742  # Eq. 13
  ray_ref <- (taur * ray_phase) / (4 * us * uv)  # Eq. 11
  ray_ref <- ray_ref * Pa / 1013.25  # correction for pressure variation (not in the paper)
  taurz <- taur * Peq  # Eq. 12
  
  # -------------- 10) aerosol atmospheric reflectance (3.4.2) ----------------
  
  aer_phase <- a0P + a1P * ksiD + a2P * ksiD^2 + a3P * (ksiD^3) + a4P * (ksiD^4)  # extension of Eq. 17 aerosol phase function
  ak2 <- (1 - wo) * (3 - wo * 3 * gc)  # Eq. 15 k^2
  ak <- sqrt(ak2)  # Eq. 15 k
  
  # -------------- X Y Z Appendix --------------------------------------------
  
  e <- -3 * us^2 * wo / (4 * (1 - ak2 * us^2))  # E = -3 * XMUS
  f <- -(1 - wo) * 3 * gc * us^2 * wo / (4 * (1 - ak2 * us^2))
  dp <- e / (3 * us) + us * f
  d <- e + f
  b <- 2 * ak / (3 - wo * 3 * gc)
  delta <- exp(ak * taup) * (1 + b)^2 - exp(-ak * taup) * (1 - b)^2
  ww <- wo / 4
  ss <- us / (1 - ak2 * us^2)
  q1 <- 2 + 3 * us + (1 - wo) * 3 * gc * us * (1 + 2 * us)
  q2 <- 2 - 3 * us - (1 - wo) * 3 * gc * us * (1 - 2 * us)
  q3 <- q2 * exp(-taup / us)
  c1 <- ((ww * ss) / delta) * (q1 * exp(ak * taup) * (1 + b) + q3 * (1 - b))
  c2 <- -((ww * ss) / delta) * (q1 * exp(-ak * taup) * (1 - b) + q3 * (1 + b))
  cp1 <- c1 * ak / (3 - wo * 3 * gc)
  cp2 <- -c2 * ak / (3 - wo * 3 * gc)
  z <- d - wo * 3 * gc * uv * dp + wo * aer_phase / 4
  x <- c1 - wo * 3 * gc * uv * cp1
  y <- c2 - wo * 3 * gc * uv * cp2
  aa1 <- uv / (1 + ak * uv)
  aa2 <- uv / (1 - ak * uv)
  aa3 <- us * uv / (us + uv)
  
  aer_ref1 <- x * aa1 * (1 - exp(-taup / aa1))  # Eq.13(1)
  aer_ref2 <- y * aa2 * (1 - exp(-taup / aa2))  # Eq.13(2)
  aer_ref3 <- z * aa3 * (1 - exp(-taup / aa3))  # Eq.13(3)
  aer_ref <- (aer_ref1 + aer_ref2 + aer_ref3) / (us * uv)
  
  # -------------- 11) Residue Rayleigh (not in the paper) --------------------
  
  Res_ray <- Resr1 + Resr2 * taur * ray_phase / (us * uv) + Resr3 * ((taur * ray_phase / (us * uv))^2)
  
  # -------------- 12) Residue Aerosol ----------------------------------------
  
  Res_aer <- (Resa1 + Resa2 * (taup * m * cksi) + Resa3 * ((taup * m * cksi)^2)) + Resa4 * ((taup * m * cksi)^3)
  
  # -------------- 13) Term coupling molecule/aerosol -------------------------
  
  tautot <- taup + taurz
  
  Res_6s <- (Rest1 + Rest2 * (tautot * m * cksi) + Rest3 * ((tautot * m * cksi)^2)) + Rest4 * ((tautot * m * cksi)^3)
  # ----------------- 14) total atmospheric reflectance ----------------------
  
  atm_ref <- ray_ref - Res_ray + aer_ref - Res_aer + Res_6s
  
  # -------------------- 15) TOA reflectance -----------------------------------
  
  # r_toa = r_surf * tg * ttetas * ttetav / (1 - r_surf * s) + (atm_ref * tg)
  
  ## added by Peiqi Yang for non-lambertian surface
  tdir_tts <- exp(-tautot/us)  # downward
  tdir_ttv <- exp(-tautot/uv)  # upward
  tdif_tts <- ttetas - tdir_tts  # downward
  tdif_ttv <- ttetav - tdir_ttv  # upward
  
  # Create a list for atmospheric optics parameters
  data.atm <- list()
  
  data.atm[['Ta_ss']] <- tdir_tts  # directional transmittance for direct incidence
  data.atm[['Ta_sd']] <- tdif_tts  # hemispherical transmittance for direct incidence
  data.atm[['Ta_oo']] <- tdir_ttv  # directional transmittance for direct incidence (in the viewing direction)
  data.atm[['Ta_do']] <- tdif_ttv  # hemispherical transmittance for direct incidence (in the viewing direction)
  data.atm[['Ta_s']] <- ttetas    # directional transmittance for diffuse light
  data.atm[['Ta_o']] <- ttetav    # hemispherical transmittance for diffuse light
  data.atm[['Tg']] <- tg        # total scattering transmission
  
  data.atm[['Ra_dd']] <- s          # hemispherical atmospheric reflectance for diffuse light
  data.atm[['Ra_so']] <- atm_ref    # directional atmospheric reflectance for direct incidence
 
  return(data.atm) 
}

  

