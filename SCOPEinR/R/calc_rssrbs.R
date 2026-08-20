#' Function for calculating the rss and rbs values based on the input parameters (SMC, LAI and rbs)
#'
#' \code{calc_rssrbs} calculates the soil surface resistance to evaporation
#' (\code{rss}) as an empirical function of soil moisture content, and scales
#' the boundary-layer resistance of the soil (\code{rbs}) with leaf area index.
#'
#' @param SMC numeric. Soil moisture content (m^3 m^-3).
#' @param LAI numeric. Leaf area index of the canopy (m^2 m^-2).
#' @param rbs numeric. Boundary-layer resistance of the soil at reference LAI = 3.3 (s m^-1); overwritten in the returned list with the LAI-scaled value.
#'
#' @return A list with:
#' \describe{
#'   \item{rss}{numeric. Soil surface resistance to evaporation (s m^-1), computed as \code{11.2 * exp(42 * (0.22 - SMC))}.}
#'   \item{rbs}{numeric. Soil boundary-layer resistance scaled by LAI (s m^-1), computed as \code{rbs * LAI / 3.3}.}
#' }
#' @export
#'
#' @author Carlos Camino
#'
#' @examples
#' \dontrun{
#' rss_rbs <- calc_rssrbs(SMC, LAI, rbs)
#' }
calc_rssrbs<-function(SMC,LAI,rbs){

rss_rbs<-list()
rss = 11.2 * exp(42 * (0.22 - SMC))
rss_rbs$rss <- rss

rbs = rbs * LAI / 3.3
rss_rbs$rbs <-rbs

return(rss_rbs)

}

