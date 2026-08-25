
#' @title get.aggregator.ebal
#'
#' @description
#' \code{get.aggregator.ebal} aggregates a leaf-scale flux, computed
#' separately for sunlit and shaded leaves, into a single LAI-scaled
#' canopy-integrated total flux, by internally calling
#' \code{\link{meanleaf.v2}} for each contribution and combining them
#' weighted by the fraction of sunlit leaf area \code{Fs}.
#'
#' @param LAI numeric. Total (one-sided) leaf area index of the canopy (m2 m-2).
#' @param sunlit_flux array with dimensions (nlincl, nlazi, nlayers) (or as required by \code{canopy.choice}). Leaf-scale flux for sunlit leaves.
#' @param shaded_flux numeric vector of length nlayers. Leaf-scale flux for shaded leaves (averaged over leaf angle, one value per layer).
#' @param Fs numeric vector of length nlayers. Fraction of sunlit leaf area at each canopy layer (\code{Ps}); the shaded contribution is weighted by \code{1 - Fs}.
#' @param data.canopy list. Canopy structural properties: \code{nlayers} (number of canopy layers), \code{nlincl} (number of leaf inclination classes), \code{nlazi} (number of leaf azimuth classes), \code{lidf} (leaf inclination distribution function weights).
#' @param canopy.choice character. Aggregation mode passed to \code{\link{meanleaf.v2}}: one of \code{'angles'}, \code{'layers'}, or \code{'angles_and_layers'}.
#'
#' @return numeric. LAI-scaled canopy-integrated total flux \code{flux_tot}, combining the sunlit and shaded contributions.
#' @export
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#'
#' @examples
#' \dontrun{
#' flux_tot <- get.aggregator.ebal(LAI, sunlit_flux, shaded_flux, Fs,
#'                                  data.canopy, canopy.choice = 'angles_and_layers')
#' }
get.aggregator.ebal <- function(LAI, sunlit_flux, shaded_flux, Fs, data.canopy, canopy.choice) {

  nl <- data.canopy[['nlayers']]
  nli <- data.canopy[['nlincl']]
  nlazi <- data.canopy[['nlazi']]
  lidf <- data.canopy[['lidf']]

  # if (is.vector(sunlit_flux)){
  #
  #   sunlit_flux_ <- array(0, dim = c(nli, nlazi,nl))
  #   sunlit_flux_[,,1:nl] <- sunlit_flux
  #
  # }
  #
  #
  # if (is.vector(shaded_flux)){
  #
  #   shaded_flux_ <- array(0, dim = c(nli, nlazi,nl))
  #   shaded_flux_[,,1:nl] <- shaded_flux
  # }
  #

  flux_tot <- LAI * (meanleaf.v2(data.canopy, sunlit_flux, canopy.choice, Fs) + meanleaf.v2(data.canopy, shaded_flux, canopy.choice = 'layers', 1 - Fs))
  return(flux_tot)
}


#' @title meanleaf.v2
#'
#' @description
#' \code{meanleaf.v2} averages a leaf-scale flux or property over leaf angle
#' classes and/or canopy layers, weighted by the leaf inclination
#' distribution function (\code{lidf}) and/or a layer weight \code{Ps}
#' (typically the sunlit/shaded fraction). This is a vectorized variant of
#' \code{\link{meanleaf}} used internally by \code{\link{get.aggregator.ebal}}.
#'
#' @param data.canopy list. Canopy structural properties: \code{nlayers} (number of canopy layers), \code{nlincl} (number of leaf inclination classes), \code{nlazi} (number of leaf azimuth classes), \code{lidf} (leaf inclination distribution function weights).
#' @param F_ array or numeric vector. Leaf-scale flux/property to average; expected as an array with dimensions (nlincl, nlazi, nlayers) for \code{canopy.choice = 'angles'} or \code{'angles_and_layers'}, or as a numeric vector of length \code{nlayers} for \code{canopy.choice = 'layers'}.
#' @param canopy.choice character. Averaging mode: \code{'angles'} (average over leaf inclination/azimuth only), \code{'layers'} (weighted average over layers using \code{Ps}), or \code{'angles_and_layers'} (average over both).
#' @param Ps numeric vector of length nlayers, or scalar. Layer weight (e.g. fraction of sunlit leaf area) used when \code{canopy.choice} includes \code{'layers'}.
#'
#' @return numeric. Single averaged value \code{Fout_vertical}.
#' @export
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' Fout <- meanleaf.v2(data.canopy, F_, canopy.choice = 'layers', Ps)
#' }
meanleaf.v2 <- function(data.canopy, F_, canopy.choice, Ps) {
  nl <- data.canopy[['nlayers']]
  nli <- data.canopy[['nlincl']]
  nlazi <- data.canopy[['nlazi']]
  lidf <- data.canopy[['lidf']]

  Fout <- array(0, dim = c(nli, nlazi, nl))

  switch(canopy.choice,

         "angles" = {
           for (j in 1:nli) {
             Fout[j,,] <- F_[j,,] * lidf[j]
           }

           Fout_vertical <- sum(Fout) / nlazi

         },

         "layers" = {

           #this is only for a single vector (nlayers)

           #F_nlayer <- apply(F_, c(3), sum)
           #Fout_vertical <- Ps[1:nl] *  F_nlayer / nl
           #return a single value
           Fout_vertical = sum(Ps * (F_) ) / nl
         },

         "angles_and_layers" = {
           for (j in 1:nli) {
             Fout[j,,] <- F_[j,,] * lidf[j]
           }
           for (j in 1:nl) {
             Fout[,,j] <- Fout[,,j] * Ps[j]
           }
           Fout_vertical <- sum(sum(sum(Fout))) / (nlazi * nl)
         }
  )

  return(Fout_vertical)
}

