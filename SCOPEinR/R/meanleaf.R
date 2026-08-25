


#' Calculates the layer average and the canopy average of leaf properties per layer, per leaf angle and per leaf azimuth (36)
#'
#' \code{meanleaf} averages a leaf-scale property, given as a 3D array over
#' leaf inclination, leaf azimuth and canopy layer, weighted by the leaf
#' inclination distribution function (\code{lidf}) and/or a per-layer weight
#' \code{Ps} (typically the fraction of sunlit leaf area).
#'
#' @param canopy list. Canopy structural properties: \code{nlayers} (number of canopy layers), \code{nlincl} (number of leaf inclination classes), \code{nlazi} (number of leaf azimuth classes), \code{lidf} (leaf inclination distribution function weights).
#' @param F_ numeric array with dimensions (nli, nlazi, nl). Leaf-scale property to average, indexed by leaf inclination class, leaf azimuth class, and canopy layer.
#' @param canopy.choice character. Integration method: \code{'angles'} (integration over leaf angles only, one value per layer), \code{'layers'} (integration over layers only, weighted by \code{Ps}, single value), or \code{'angles_and_layers'} (integration over both leaf angles and layers, weighted by \code{lidf} and \code{Ps}).
#' @param Ps numeric vector of length nl. Fraction of sunlit leaf area per layer, used as the layer weight for \code{canopy.choice \%in\% c('layers', 'angles_and_layers')}.
#'
#' @return numeric vector. \code{Fout_vertical}: for \code{canopy.choice = 'angles'}, one averaged value per layer (length \code{nl}); for \code{'layers'} or \code{'angles_and_layers'}, a single averaged value.
#' @export
#'
#' @details Last update: 7 December 2007; update 11 February 2008 made modular (Joris Timmermans); update 25 Feb 2013 Wout Verhoef: proposed name change, removed globals and used canopy-structure for input.
#'
#' @author 	Christiaan van der Tol (Original version in Matlab)
#' @author Carlos Camino (Ported version into R)
#'
#' @examples
#' \dontrun{
#' Fout <- meanleaf(canopy, F_, canopy.choice = 'angles_and_layers', Ps)
#' }
meanleaf <- function(canopy, F_, canopy.choice, Ps) {

  nl <- canopy[['nlayers']]
  nli <- canopy[['nlincl']]
  nlazi <- canopy[['nlazi']]
  lidf <- canopy[['lidf']]


  # create an empty list
  F_list <- list()
  # loop over the third dimension of F_
  for (i in 1:dim(F_)[3]) {
    # extract a 13 x 36 matrix from F_ and add it to the list
    F_list[[i]] <- F_[,,i]
  }

  # Output:
  #   Fout    in case of choice = 'angles': [nl]
  #           in case of choice = 'angles_and_layers': [1]

  Fout <- array(0, dim = c(nli, nlazi,nl))

  switch(canopy.choice,

         # Integration over leaf angles
         'angles' = {
           for (j in 1:nli) {
             Fout[j,,] <- F_[j,,] * lidf[j]
           }

           Fout_mean <- list()
           for (j in 1:nl) {
             Fout_mean[[j]] <- sum(Fout[, , j]) / nlazi
           }
           # One value per layer: sum over leaf inclination/azimuth,
           # weighted by lidf, divided by nlazi.
           Fout_vertical <- unlist(Fout_mean)
         },

         # Integration over layers only
         'layers' = {


           #this is only for a single vector (nlayers)

           #F_nlayer <- apply(F_, c(3), sum)
           #Fout_vertical <- Ps[1:nl] *  F_nlayer / nl
           #return a single value
           Fout_vertical = sum(Ps * (F_) ) / nl

         },

         # Integration over both leaf angles and layers
         'angles_and_layers' = {
           for (j in 1:nli) {
             Fout[j,,] <- F_[j,,] * lidf[j]
           }
           for (j in 1:nl) {
             Fout[,,j] <- Fout[,,j] * Ps[j]
           }
           # Sums over leaf angle, azimuth, and layer to return a single
           # canopy-wide scalar, matching how callers in ebal.R use this
           # value (e.g. Htot <- Hstot + Hctot, where Hstot is already
           # a scalar).
           Fout_vertical <- sum(Fout) / (nlazi * nl)
         }
  )

  return(Fout_vertical)
}
