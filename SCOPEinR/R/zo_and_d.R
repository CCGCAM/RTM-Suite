#' @title get.zo_and_d model
#' \code{get.zo_and_d} Calculates roughness length for momentum and zero
#'  plane displacement from vegetation height and LAI
#' @param inputLUT list. Canopy structure inputs (vegetation height, LAI, etc.) for the roughness-length calculation.
#' @param constants list. Physical constants used in the aerodynamic roughness formulation.
#' @param calc.heat logical. Whether to include the heat-flux-related roughness adjustment terms.
#' @param calc.rss_rbs logical. Whether to also compute soil/boundary-layer resistance terms alongside roughness length and displacement height.
#'
#' @return zo_and_d: a list with roughness length for momentum and zero-plane displacement height.
#' @export
#' @author 	A. Verhoef (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' last updates:
#'   - 17 November 2008
#    - 17 April 2013 (structures)

#' @references Verhoef, McNaughton & Jacobs (1997), HESS 1, 81-91
#' @examples
#' \dontrun{
#' get.zo_and_d(inputLUT, constants, calc.heat = TRUE, calc.rss_rbs = TRUE)
#' }
get.zo_and_d<-function(inputLUT,constants,calc.heat,calc.rss_rbs){

  zo_and_d = list()
  ## constants
  # constants used (as global)
  #   kappa       Von Karman's constant

  kappa   =  subset(constants,constant == 'kappa')[[2]]

  canopy<- SCOPEinR::getinputLUT(inputLUT, dataset='canopy')
  soil<- SCOPEinR::getinputLUT(inputLUT, dataset='soil',
                               calc.heat = calc.heat,
                               calc.rss_rbs =  calc.rss_rbs)


  # soil fields used:
  #   Cd          Averaged drag coefficient for the vegetation
  #   CR          Drag coefficient for isolated tree
  #   CSSOIL      Drag coefficient for soil
  #   CD1         Fitting parameter
  #   Psicor      Roughness layer correction
  ## parameters
  CR      = canopy[['CR']];
  CSSOIL  = canopy[['CSSOIL']];
  CD1     = canopy[['CD1']];
  Psicor  = canopy[['Psicor']];

  # canopy fields used as inpuyt:
  #   LAI         one sided leaf area index
  #   hc           vegetation height (m)
  LAI     = canopy[['LAI']];
  hc       = canopy[['hc']];

  ## calculations
  sq      = sqrt(CD1 * LAI / 2);
  G1      = max(3.3, (CSSOIL + CR * LAI / 2)^(-0.5));
  if((LAI > 1e-7) & (hc > 1e-7)){
    d  = hc * (1 - (1 - exp( -sq))/sq);
    zo_and_d$d<-d
  } else {
    ## NB (fix): was `zo_and_d$d <- d`, referencing `d` which is never
    ## assigned on this branch (only the if-branch above defines it) --
    ## errors with "object 'd' not found" on a fresh call, or silently
    ## picks up a stale `d` left over from a previous call in the same R
    ## session/environment (the same stray-variable bug pattern found
    ## elsewhere in this codebase, e.g. RTMo.R's `Rndir <- fs_ * Asun[j]`).
    d <- 0
    zo_and_d$d <- d
  }
  # output:
  #   zom  roughness lenght for momentum (m)
  #   d zero plane displacement (m)

  # Eq 12 in Verhoef et al (1997)
  zom = (hc - d) * exp(-kappa * G1 + Psicor);
  zo_and_d$zom <-zom
  return(zo_and_d)
}

