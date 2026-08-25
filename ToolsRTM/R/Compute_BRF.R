
#' Computes bidirectional reflectance factor based on outputs from PRO-4SAIL and sun position
#'
#' Authors of the version:Jean-Baptiste FERET (jb.feret@teledetection.fr)
#' Florian de BOISSIEU (fdeboiss@gmail.com)
#' Copyright 2019/11 Jean-Baptiste FERET
#' 
#' The direct and diffuse light are taken into account as proposed by:
#' Francois et al. (2002) Conversion of 400-1100 nm vegetation albedo
#' measurements into total shortwave broadband albedo using a canopy
#' radiative transfer model, Agronomie
#' 
#' Es = direct
#' Ed = diffuse
#'
#' @param rdot numeric. Hemispherical-directional reflectance factor in viewing direction
#' @param rsot numeric. Bi-directional reflectance factor
#' @param tts numeric. Solar zenith angle
#' @param data.light list. direct and diffuse radiation for clear conditions, is NULL use default values
#' @param short.waves boolean . Is true the outputs is shorted to wave for Fluspect-model.
#' @return BRF numeric. Bidirectional reflectance factor
#' @export
#' 
Compute_BRF  <- function(rdot=NULL,rsot=NULL,tts=NULL,data.light=NULL,short.waves=T){
  
  ############################## #
  ##	direct / diffuse light	##
  ############################## #

  if (is.null(data.light)){
    Es <- ToolsRTM::dataSpec_PDB[,11]
    Ed <- ToolsRTM::dataSpec_PDB[,12]
    rd <- pi / 180
    
  } else{
    Es <- data.light$direct_light ##
    Ed <- data.light$diffuse_light ##
    rd <- pi / 180
   
  }
  ############################## #
  ##	Get only waves for direct / diffuse light
  ## for Fluspect-D and Cx-B
  ############################## #

  # NB (fix): was `if (missing(short.waves)) short.waves = FALSE` -- this
  # ignored whatever value the caller actually passed (or the signature's own
  # declared default, short.waves=T) and forced FALSE whenever the argument
  # was omitted. Because rdot/rsot are 2101 long for PROSPECT-D/-PRO/Liberty
  # but only 2001 long for Fluspect-B/-Cx (its native 400-2400nm domain),
  # skipping truncation for an omitted short.waves silently left Es/Ed at
  # 2101 -- for Fluspect-based rdot/rsot (2001) that length mismatch caused R
  # to silently recycle Es/Ed against rdot/rsot below, corrupting the last
  # ~100nm (roughly 2401-2500nm) of the output with wrapped-around values
  # from the start of the spectrum (400-500nm) instead of an error or the
  # correct truncation. Match Es/Ed to rdot's actual length directly instead
  # of trusting a caller-supplied flag -- this also makes the flag redundant,
  # kept only for backward API compatibility with existing callers.
  n <- length(rdot)
  Es <- Es[seq_len(n)]
  Ed <- Ed[seq_len(n)]
  #
  # if (skyl == 0.1){
  #   #  by default skyl = 0.1
  #   skyl = 0.847 - 1.61 * sin_90tts + 1.04 * sin_90tts ** 2 #this equation return sky=0.1 
  # } else {
  #   skyl=inputLUT[,'skyl']
  # }
  
  skyl <- 0.847- 1.61 * sin((90 - tts) * rd)+ 1.04 * sin((90 - tts) * rd)*sin((90 - tts) * rd) # diffuse radiation (Francois et al., 2002)
  PARdiro <- (1 - skyl) * Es
  PARdifo <- skyl * Ed
  BRF <- (rdot * PARdifo + rsot * PARdiro)/(PARdiro + PARdifo)

  return(BRF)
}