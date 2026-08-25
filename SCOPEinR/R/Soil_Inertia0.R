#' Calculate the soil thermal inertia from known soil thermal properties
#'
#' \code{Soil_Inertia0} computes the soil thermal inertia \code{GAM} directly
#' from the soil specific heat, bulk density and thermal conductivity, used
#' when \code{options.soil_heat_method == 0} in the SCOPE energy balance.
#'
#' @param cs numeric. Soil specific heat capacity (J kg^-1 K^-1).
#' @param rhos numeric. Soil bulk density (kg m^-3).
#' @param lambdas numeric. Soil thermal conductivity (W m^-1 K^-1).
#'
#' @return numeric. Soil thermal inertia \code{GAM} (J m^-2 K^-1 s^-1/2), used to compute the soil heat flux.
#' @export
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' GAM <- Soil_Inertia0(cs, rhos, lambdas)
#' }
Soil_Inertia0<-function(cs,rhos,lambdas){
  # soil thermal inertia
  GAM = sqrt(cs * rhos * lambdas);  # soil thermal intertia
  return(GAM)
}

